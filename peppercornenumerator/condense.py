#
#  peppercornenumerator/condense.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

import operator
import collections
import numpy as np
import itertools as it
from functools import reduce

from .utils import tarjans
from .objects import (SingletonError, 
                      PepperReaction,
                      PepperMacrostate)

class CondensationError(Exception):
    pass

class PepperCondensation:
    def __init__(self, enumerator):
        self.enumerator = enumerator

        # Helper dictionaries
        self.reactions_consuming = get_reactions_consuming(self.complexes, 
                                                           list(self.detailed_reactions))

        # Map complexes to their products from 1-1 fast reactions.
        fast11_products = {c: set().union(*[r.products for r in rxns \
                            if r.arity == (1, 1) and self.is_fast(r)]) \
                                for (c, rxns) in self.reactions_consuming.items()}

        self.SCCs = tarjans(list(self.complexes), fast11_products)
        self.scc_containing = {c: scc for scc in self.SCCs for c in scc}

        # Condensed graph components
        self.set_to_fate = dict()
        self._complex_fates = dict()
        self._condensed_reactions = None

        # Condensed reaction rate calculations
        self.stationary_dist = dict()
        self.cplx_decay_prob = collections.defaultdict(float)
        self.exit_prob = dict()
        self.reaction_decay_prob = collections.defaultdict(float)

    @property
    def complexes(self):
        return self.enumerator.complexes

    @property
    def detailed_reactions(self):
        return self.enumerator.reactions

    @property
    def resting_macrostates(self):
        if self._condensed_reactions is None:
            self.condense()
        return set(self.set_to_fate.values())

    @property
    def condensed_reactions(self):
        if self._condensed_reactions is None:
            self.condense()
        return iter(self._condensed_reactions)

    def is_fast(self, rxn):
        return rxn.arity[0] == 1 and rxn.rate_constant[0] >= max(self.enumerator.k_slow, self.enumerator.k_fast)

    def complex_fates(self, cplx):
        if cplx not in self._complex_fates:
            self.compute_fates(self.scc_containing[cplx])
        return self._complex_fates[cplx]

    def compute_fates(self, scc):
        """ Calculates all fates for a given SCC. """
        cfates = self._complex_fates
        rxncon = self.reactions_consuming

        if scc[0] in cfates:
            return # This SCC has been processed before

        scc_set = frozenset(scc)
        out_rxns = [r for c in scc for r in rxncon[c] \
                        if self.is_fast(r) and is_outgoing(r, scc_set)]

        if len(out_rxns) == 0: # A resting macrostate
            try:
                resting_ms = PepperMacrostate(scc_set)
            except SingletonError as err:
                resting_ms = err.existing
            self.set_to_fate[scc_set] = resting_ms
            # find all reactions between complexes in this SCC (non-outgoing 1,1 reactions)
            internal_reactions = [r for c in scc for r in rxncon[c] 
                    if (r.arity == (1, 1) and not is_outgoing(r, scc_set))]
            self.stationary_dist[resting_ms] = get_stationary_distribution(scc, internal_reactions)
            fate = (resting_ms,)
            for c in scc:
                assert c not in cfates
                cfates[c] = SetOfFates([fate])
                # All complexes in this SCC decay to this SCC with probability 1.
                self.cplx_decay_prob[(c, fate)] = 1.0
        else:
            # All combinations of fates of each product of outgoing reactions.
            reaction_fate_combinations = [list(cartesian_product(
                map(self.complex_fates, rxn.products))) for rxn in out_rxns]

            # Compute the fates of each of the outgoing reactions by summing
            # each element above. Returns a list, with each element
            # corresponding to the frozenset of fates for one reaction.
            reaction_fates = [sorted(map(tuple_sum, combination))
                              for combination in reaction_fate_combinations]
            
            # note that these two are equivalent; the intermediate
            # reaction_fate_combinations is only calculated for rates
            assert reaction_fates == [sorted(cartesian_sum(
                list(map(self.complex_fates, rxn.products)))) for rxn in out_rxns]

            # find all fast reactions involving complexes in this SCC
            fast_reactions = [r for c in scc for r in rxncon[c] if self.is_fast(r)]
            self.exit_prob[scc_set] = get_exit_probabilities(scc, fast_reactions)

            # calculate the decay probability for each reaction fate
            for i, rxn in enumerate(out_rxns):
                # for each possible fate of rxn:
                #  calculate the probability that the products decay to fate.
                for fate in reaction_fates[i]:
                    # get all combinations of product fates that sum to the
                    # same reaction fate.
                    combinations = reaction_fate_combinations[i]
                    self.get_reaction_decay_probabilities(rxn, fate, combinations)
                    for c in scc:
                        # P(c decays to F) = 
                        # = P(SCC exits via rxn | SCC was entered via c) * P(rxn decays to F)
                        self.cplx_decay_prob[(c, fate)] += self.exit_prob[scc_set][
                                (c, rxn)] * self.reaction_decay_prob[(rxn, fate)]
            set_of_fates = frozenset().union(*reaction_fates)
            for c in scc:
                assert c not in cfates
                cfates[c] = SetOfFates(set_of_fates)
        return

    def get_reaction_decay_probabilities(self, rxn, fate, combinations = None):
        """ Calculate the decay probability of a reaction to a particular fate.

        Each combination (`fates`) that sums to `fate` constitutes a possible
        way this reaction can decay to `fate`, and therefore contributes to 
        P(r -> fate). This contribution is the joint probability that each
        product `d` of r decays to the corresponding fate `f` in `fates`.
        """
        if combinations is None:
            combinations = list(cartesian_product(list(map(self.complex_fates, rxn.products))))

        # NOTE: introduced as a consequence of a bug discovered in the
        # cooperative binding case. Before you calculate (add to) reaction
        # decay probabilities, make sure that the dictionary is empty.
        self.reaction_decay_prob = collections.defaultdict(float)

        for fates in combinations:
            if tuple_sum(fates) != fate:
                continue
            self.reaction_decay_prob[(rxn, fate)] += times(
                    self.cplx_decay_prob[(d, f)] for (d, f) in zip(rxn.products, fates))
        return

    def condense(self):
        """ Do the reaction network condensation. """
        if self._condensed_reactions is not None:
            raise CondensationError('condensation called twice')
        else:
            self._condensed_reactions = set()

        for scc in self.SCCs:
            self.compute_fates(scc)

        for rxn in self.detailed_reactions:
            # Fast reactions have been handled by compute_fates
            if self.is_fast(rxn):
                continue
            if rxn.arity == (0, 0):
                # Should be safe to continue here ... but why would you?
                raise CondensationError(f'I refuse to condense {rxn}!')

            # Get the corresponding fates (resting macrostate)
            reactant_fates = [self.complex_fates(x) for x in rxn.reactants]
            product_fates = [self.complex_fates(x) for x in rxn.products]

            # Make sure each reactant is a resting macrostate: 
            assert all(Fc.is_singleton() for Fc in reactant_fates)

            # Generate new reactions with each non-trivial combination
            for (reactants, products) in it.product(cartesian_sum(reactant_fates),
                                                    cartesian_sum(product_fates)):
                if sorted(reactants) == sorted(products):
                    # skip the trivial case
                    continue
                prxn = PepperReaction(reactants, products, rtype = 'condensed')
                const = (self.get_condensed_rate(prxn), '/M' * (len(reactants)-1) + '/s')
                if prxn.rate_constant[0] is None:
                    prxn.rate_constant = const
                self._condensed_reactions.add(prxn)
        return

    def get_condensed_rate(self, rxn):
        """
        Condensed rate from P + Q -> Fate
        Sum_{i,j,l} (p_i * q_j * K_{ij}^l * d^l(Fate)
        """
        products = tuple(rxn.products)

        # calculate reaction rate by summing over all representative detailed reactions
        detailed_reactions_consuming = set(
                [r for rms in rxn.reactants for c in rms.complexes 
                    for r in self.reactions_consuming[c]])

        rate = 0.0
        for r in detailed_reactions_consuming:
            # Calculate the probability that this detailed reaction decays to
            # the products P(r -> F), where F = products.
            self.get_reaction_decay_probabilities(r, products)

            # probability that the products of this detailed reaction decay
            # into fates that yield the condensed reaction
            product_probability = self.reaction_decay_prob[(r, products)]
            assert 0 <= product_probability < 1.000001

            # probability that the resting macrostates comprising the reactants
            # of the condensed reaction will be in the right configuration to
            # perform the detailed reaction
            reactant_probabilities = times(self.stationary_dist[self.set_to_fate[
                frozenset(self.scc_containing[a])]][a] for a in r.reactants)
            assert 0 <= reactant_probabilities < 1.000001

            # rate of the detailed reaction
            k = r.rate_constant[0]
            assert 0 < k

            # overall contribution of detailed reaction r to rate of the condensed reaction 
            rate += reactant_probabilities * k * product_probability

            if isinstance(rate, complex) and rate.imag > 0:
                log.warning(f"Detailed reaction {r} contributes a complex rate: " + \
                            f"{rate.real} + {rate.imag} to condensed reaciton {rxn}.")
        return rate

class SetOfFates:
    """ The set of multisets that can be reached from a given complex. """

    def __init__(self, states):
        self.states = frozenset([tuple(sorted(s)) for s in states])

    def is_singleton(self):
        return len(self.states) == 1 and len(list(self.states)[0]) == 1

    def __iter__(self):
        return iter(self.states)

    def __eq__(self, other):
        return self.states == other.states

    def __len__(self):
        return len(self.states)

    def __repr__(self):
        return f"{self.__class__}({self.states})"

def is_outgoing(reaction, scc):
    """ bool: the reaction leads out of the passed scc. """
    return not (set(reaction.products) <= scc)

def times(iterable):
    """ float: multiply each element of the given `iterable`. """
    return reduce(operator.mul, iterable)

def tuple_sum(iterable):
    """ tuple: combines an iterable of tuples into one sorted tuple. """
    return tuple(sorted(reduce(operator.concat, tuple(iterable), tuple())))

def cartesian_product(iters):
    """ iter: the cartesian product of the input iterables. """
    return (tuple(x) for x in it.product(*iters))

def cartesian_sum(iters):
    """ iter: the cartesian sum of the input iterables. """
    return [tuple_sum(x) for x in cartesian_product(iters)]

def get_reactions_consuming(complexes, reactions):
    """ dict: maps complexes to lists of reactions where they appear as a reactant. """
    return {c: [r for r in reactions if (c in r.reactants)] for c in complexes}

def get_reactions_producing(complexes, reactions):
    """ dict: maps complexes to lists of reactions where they appear as a product. """
    return {c: [r for r in reactions if (c in r.products)] for c in complexes}

def get_stationary_distribution(scc, reactions, epsilon = 1e-5):
    """ Calculate the stationary distribution of complexes in a SCC.

    Returns:
        dict: Mapping of complexes to stationary distribution (flt).
    """
    log.debug(f"Calculating stationary distribution for scc: {[x.name for x in scc]}")

    scc_set = frozenset(scc)
    scc_list = sorted(scc)
    L = len(scc)

    # assign a numerical index to each complex
    complex_indices = {c: i for (i, c) in enumerate(scc_list)}

    # T is the transition rate matrix, defined as follows:
    # T_{ij} = { rate(j -> i)       if  i != j
    #          { - sum_{j'} T_{j'i} if  i == j
    T = np.zeros((L, L))

    # add transition rates for each reaction to T
    for r in reactions:
        assert r.arity[0] == 1 and next(r.reactants) in scc
        assert r.arity[1] == 1 and next(r.products) in scc
        # r : a -> b
        # T_{b,a} = rate(r : a -> b)
        a = next(r.reactants)
        b = next(r.products)
        T[complex_indices[b]][complex_indices[a]] = r.rate_constant[0]

    T0 = np.copy(T)

    # compute diagonal elements of T
    T_diag = np.sum(T, axis=0)  # sum over columns
    for i in range(L):
        T[i][i] = -T_diag[i]

    # calculate eigenvalues
    (w, v) = np.linalg.eig(T)
    # w is array of eigenvalues
    # v is array of eigenvectors, where column v[:,i] is eigenvector 
    # corresponding to the eigenvalue w[i].

    # find eigenvector corresponding to eigenvalue zero (or nearly 0)
    i = np.argmin(np.abs(w))
    if abs(w[i]) > epsilon:
        log.warning("Bad stationary distribution for resting set transition matrix. " + \
                    f"Eigenvalue {w[i]} has magnitude greater than epsilon = {epsilon}. " + \
                    f"Markov chain may be periodic, or epsilon may be too high. Eigenvalues: {w}")

    # Check that the stationary distribution is good.
    s = v[:, i]
    if not ((s >= 0).all() or (s <= 0).all()): 
        log.error('Condensed reaction rates may be incorrect. ' + \
           f'Stationary distribution of resting complexes is an eigenvector of mixed sign: {s}')
    s = s / np.sum(s)
    if not (abs(np.sum(s) - 1) < epsilon):
        log.error('Condensed reaction rates may be incorrect. ' + \
           f'Stationary distribution of resting complexes does not sum to 1: {abs(np.sum(s))}.')
    return {c: s[i] for (c, i) in complex_indices.items()}
    
def get_exit_probabilities(scc, reactions):
    """ Calculate the exit probabilities of complexes via reactions in a SCC.

    Returns:
        dict: Mapping of (incoming complex, outgoing rxn) to exit probability (flt).
    """
    log.debug(f"Calculating exit probabilities for scc: {[x.name for x in scc]}")
    # build set and list of elements in SCC; assign a numerical index to each complex
    scc_set = frozenset(scc)
    scc_list = sorted(scc)
    complex_indices = {c: i for (i, c) in enumerate(scc_list)}

    # sort reactions into internal and outgoing; assign a numerical index to each reaction
    r_outgoing = sorted([r for r in reactions if is_outgoing(r, scc_set)])
    r_internal = sorted([r for r in reactions if not is_outgoing(r, scc_set)])
    exit_indices = {r: i for (i, r) in enumerate(r_outgoing)}

    L = len(scc) # num of complexes in SCC
    e = len(r_outgoing) # num of exit pathways

    # add transition rates for each internal reaction
    T = np.zeros((L, L))
    for r in r_internal:
        assert r.arity == (1, 1)
        a = next(r.reactants)
        b = next(r.products)
        T[complex_indices[a]][complex_indices[b]] = r.rate_constant[0]

    # add transition rates for each outgoing reaction
    Te = np.zeros((L, e))
    for r in r_outgoing:
        a = next(r.reactants)
        Te[complex_indices[a]][exit_indices[r]] = r.rate_constant[0]

    # the full transition matrix P_{L+e x L+e} would be
    #
    # P = ( Q_{LxL}  R_{Lxe}   )
    #     ( 0_{exL}  I_{exe}   )
    #
    # but we really only care about Q to calculate the fundamental matrix,
    # so we construct
    #
    # P = (T_{LxL} Te_{Lxe})
    P = np.hstack((T, Te))

    # then normalize P along each row, to get the overall transition
    # probabilities, e.g. P_ij = P(i -> j), where i,j in 0...L+e
    P = P / np.sum(P, 1)[:, np.newaxis]
    #log.debug(f"{P=}")

    # extract the interior transition probabilities (Q_{LxL})
    Q = P[:, 0:L]
    #log.debug(f"{Q=}")

    # extract the exit probabilities (R_{Lxe})
    R = P[:, L:L + e]
    #log.debug(f"{R=}")

    # calculate the fundamental matrix (N = (I_L - Q)^-1)
    N = np.linalg.inv(np.eye(L) - Q)
    #log.debug(f"{N=}")

    # make sure all elements of fundamental matrix are >= 0
    if not (N >= 0).all():
        log.error('Condensed reaction rates may be incorrect. ' + \
            f'Negative elements in fundamental matrix.')

    # calculate the absorption matrix (B = NR)
    B = np.dot(N, R)

    # map tuples of (incoming complex, outgoing reaction) to exit probabilities.
    return {(c, r): B[i, j] for (c, i) in complex_indices.items()
            for (r, j) in exit_indices.items()}

