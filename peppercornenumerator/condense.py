#
#  peppercornenumerator/condense.py
#  EnumeratorProject
#
from __future__ import absolute_import, print_function, division

import logging
log = logging.getLogger(__name__)

import operator
import collections
import itertools as it
import numpy as np
from functools import reduce

from peppercornenumerator.objects import PepperComplex
from peppercornenumerator.objects import PepperReaction
from peppercornenumerator.objects import PepperMacrostate
from peppercornenumerator.objects import DSDDuplicationError

class CondensationError(Exception):
    pass

class PepperCondensation(object):
    def __init__(self, enumerator):
        self._enumerator = enumerator

        # Condensed graph components
        self._SCCs = None
        self._set_to_fate = dict()
        self._complex_fates  = dict()
        self._condensed_reactions = None

        # Helper dictionaries
        self._scc_containing = None
        self._reactions_consuming = None

        # Condensed reaction rate calculations
        self._exit_probabilities = dict()
        self._stationary_distributions = dict()
        self._decay_probabilities = collections.defaultdict(float)
        self._reaction_decay_probabilities = collections.defaultdict(float)

    @property
    def enumerator(self):
        return self._enumerator

    @property
    def complexes(self):
        return self._enumerator.complexes

    @property
    def resting_sets(self):
        print("# Deprecated function: PepperCondensation.resting_sets. Replace with PepperCondensation.resting_macrostates")
        return self.resting_macrostates

    @property
    def resting_macrostates(self):
        if not self._set_to_fate:
            raise CondensationError('need to call condense first')
        return list(self._set_to_fate.values())

    @property
    def resting_set_representatives(self):
        if not self._set_to_fate:
            raise CondensationError('need to call condense first')
        return [x.canonical for x in list(self._set_to_fate.values())]

    @property
    def set_to_fate(self):
        if not self._set_to_fate:
            raise CondensationError('need to call condense first')
        return self._set_to_fate

    @property
    def cplx_to_fate(self):
        if not self._complex_fates:
            raise CondensationError('need to call condense first')
        return self._complex_fates

    @property
    def reactions(self):
        return self.condensed_reactions + self.detailed_reactions

    @property
    def detailed_reactions(self):
        return self._enumerator.reactions

    @property
    def condensed_reactions(self):
        if self._condensed_reactions is None:
            self.condense()
        return list(self._condensed_reactions)

    @property
    def k_fast(self):
        return self._enumerator.k_fast

    def is_fast(self, rxn):
        if rxn.arity[0] != 1:
            return False
        if self._enumerator.local_elevation:
            from peppercornenumerator.enumerator import local_elevation_rate
            k = local_elevation_rate(rxn)
        else :
            k = rxn.const
        return k >= self.k_fast

    def reactions_consuming(self, cplx):
        if self._reactions_consuming is None:
            self._reactions_consuming = get_reactions_consuming(self.complexes, self.reactions)
        return self._reactions_consuming[cplx]

    @property
    def SCCs(self):
        """
        Return (compute) a list of SCCs
        """
        if self._SCCs is None:
            if not self._reactions_consuming:
                self._reactions_consuming = get_reactions_consuming(
                        self.complexes, self.detailed_reactions)
            self._SCCs = tarjans(self.complexes, self.reactions, 
                    self._reactions_consuming, self.is_fast)
        return self._SCCs

    def SCC_containing(self, cplx):
        """ Map each complex to the SCC which contains the complex. 
        Each complex should be in 1 SCC.
        """
        if self._scc_containing is None:
            self._scc_containing = {c: scc for scc in self.SCCs for c in scc}
        return self._scc_containing[cplx]

    def calculate_reaction_decay_probabilities(self, rxn, fate, combinations=None):
        """
        Calculate the decay probability of a reaction to a particular fate.

        Each combination (`fates`) that sums to `fate` constitutes a possible
        way this reaction can decay to `fate`, and therefore contributes to 
        P(r -> fate). This contribution is the joint probability that each
        product `d` of r decays to the corresponding fate `f` in `fates`.
        """
        if combinations is None:
            combinations = cartesian_product(list(map(self.get_fates, rxn.products)))

        # NOTE: introdcued as a consequence of a bug discovered in the cooperative
        # binding case. Before you calculate (add to) reaction decay
        # probabilities, make sure that the dictionary is empty.
        self._reaction_decay_probabilities = collections.defaultdict(float)

        for fates in combinations :
            if tuple_sum_sort(fates) != fate:
                continue

            # NOTE: this code is equivalent to the code below, but seems to be
            # a little less accurate!?

            #p = 1
            #for (d, f) in zip(rxn.products, fates):
            #    p *= self._decay_probabilities[(d, f)]
            #self._reaction_decay_probabilities[(rxn, fate)] += p

            self._reaction_decay_probabilities[(rxn, fate)] += times(
                    self._decay_probabilities[(d, f)] for (d, f) in zip(rxn.products, fates))

    def get_fates(self, cplx):
        if cplx not in self._complex_fates :
            self.compute_fates(self.SCC_containing(cplx))
        return self._complex_fates[cplx]

    def compute_fates(self, scc):
        """
        Processes a single SCC neighborhood, generating resting set multisets
        for each complex, and storing the mappings in the outer-scope
        `complex_fates` dict

        """
        complex_fates = self._complex_fates

        if scc[0] in complex_fates:
            # Dirty check to see if the scc has been processed before.
            return 

        # Convert to a set for fast lookup
        scc_set = frozenset(scc)

        outgoing_reactions = []
        for c in scc:
            for r in self.reactions_consuming(c):
                if self.is_fast(r) and is_outgoing(r, scc_set) :
                    outgoing_reactions.append(r)

        # If this SCC is a resting set:
        if len(outgoing_reactions) == 0 :
            # build new resting set
            try :
                resting_set = PepperMacrostate(scc)
            except DSDDuplicationError as e:
                resting_set = e.existing

            self._set_to_fate[scc_set] = resting_set

            # calculate stationary distribution
            self._stationary_distributions[resting_set] = self.get_stationary_distribution(scc)

            # assign fate to each complex in the SCC
            fate = (resting_set,) # needs to be iterable..
            fates = frozenset([fate])

            for c in scc:
                if c in complex_fates:
                    raise CondensationError('complex should not be assigned yet')
                complex_fates[c] = SetOfFates(fates)

                # all complexes in this SCC decay to this SCC with probability 1,
                # e.g. P(c -> SCC) = 1 for all c in this SCC
                self._decay_probabilities[(c, fate)] = 1.0

        else :
            # Compute all possible combinations of the fates of each product of r
            reaction_fate_combinations = [cartesian_product(
                list(map(self.get_fates, rxn.products))) for rxn in outgoing_reactions]

            # Compute the fates of each of the outgoing reactions by summing
            # each element above. Returns a list, with each element
            # corresponding to the frozenset of fates for one reaction.
            reaction_fates = [sorted(map(tuple_sum_sort, combination))
                              for combination in reaction_fate_combinations]
            
            # note that these two are equivalent; the intermediate
            # reaction_fate_combinations is only calculated for rates
            assert reaction_fates == [sorted(cartesian_sum(
                list(map(self.get_fates, rxn.products)))) for rxn in outgoing_reactions]

            # calculate the exit probabilities
            self._exit_probabilities[scc_set] = self.get_exit_probabilities(scc)

            # calculate the probability of each fate:
            # for each outgoing reaction, `r`
            for (i, r) in enumerate(outgoing_reactions):
                # for each possible `fate` of that reaction `r`
                for fate in reaction_fates[i]:
                    # NOTE: this is weird, first you calculate the reaction decay
                    # probabilities from decay probabilities, then calculate decay
                    # probabilities from reaction decay probabilities?
                    fate = tuple(sorted(fate))
                    # calculate the probability that the products of `r`
                    # decay to `fate`, e.g. P(r -> fate)

                    # iterate over all combinations of the fates of `r`
                    # that sum to equal `fate`
                    combinations = reaction_fate_combinations[i]
                    self.calculate_reaction_decay_probabilities(r, fate, combinations)

                    # the decay probability will be different for different
                    # complexes in the SCC
                    for c in scc:
                        # P(x decays to F) = P(SCC exits via r | SCC was entered via c) 
                        #   * P(r decays to F)
                        self._decay_probabilities[(c, fate)] += self._exit_probabilities[scc_set][
                                (c, r)] * self._reaction_decay_probabilities[(r, fate)]

            # The set of fates for the complexes in this SCC is the union of
            # the fates for all outgoing reactions.

            # Note that frozenset().union(*X) === X[0] U X[1] U X[2] ...
            # where X[i] is a frozenset and a U b represents the union of
            # frozensets a and b
            set_of_fates = frozenset().union(*reaction_fates)

            for c in scc:
                if c in complex_fates:
                    raise CondensationError('complex should not be assigned yet')
                complex_fates[c] = SetOfFates(set_of_fates)

    # Main function!
    def condense(self):
        """
        Reaction condensation.
        """
        if self._condensed_reactions is not None:
            raise CondensationError('condensation called twice')
        else :
            self._condensed_reactions = set()

        enum = self.enumerator

        for scc in self.SCCs:
            self.compute_fates(scc)

        for reaction in enum.reactions :
            # NOTE: Fast reactions have been handled by the fates
            if self.is_fast(reaction):
                continue

            # Filter reactions with no products/reactants
            if len(reaction.reactants) == 0 and len(reaction.products) == 0:
                raise Exception('before we continue, see why that happens...')
                #continue

            # Get the corresponding fates (resting sets)
            reactant_fates = list(map(self.get_fates, reaction.reactants))
            product_fates = list(map(self.get_fates, reaction.products))

            # Get all combinations of reactant and product fates,
            new_reactant_combinations = cartesian_sum(reactant_fates)
            new_product_combinations = cartesian_sum(product_fates)

            # Make sure each reactant is a resting set: 
            #   F(c) is a singleton for each reactant c
            for e, Fc in enumerate(reactant_fates):
                if not Fc.is_singleton() :
                    log.error("Cannot condense reaction {}: ".format(reaction) + \
                            "reactant {} has multiple fates: F({}) = {}".format(
                                str(reaction.reactants[e]), str(reaction.reactants[e]), str(Fc)))
                    raise CondensationError()

            # Take all combinations of reactants and products
            new_reactant_product_combinations = it.product(
                    new_reactant_combinations, new_product_combinations)

            # Generate new reactions with each non-trivial combination
            for (reactants, products) in new_reactant_product_combinations:
                reactants = sorted(reactants)
                products = sorted(products)

                # skip the trivial case
                if reactants == products :
                    continue

                try :
                    reaction = PepperReaction(reactants, products, rtype='condensed')
                except DSDDuplicationError as e:
                    log.debug('Duplicating PepperReaction: {}'.format(e.existing))
                    reaction = e.existing

                reaction.const = self.get_condensed_rate(reaction)
                self._condensed_reactions.add(reaction)
        return

    def get_condensed_rate(self, rxn):
        """

        Condensed rate from P + Q -> Fate
        
        Sum_{i,j,l} (p_i * q_j * K_{ij}^l * d^l(Fate)

        """
        reactants = tuple(rxn.reactants)
        products = tuple(rxn.products)

        # calculate_reaction_decay_probabilities
        # reaction_decay_probabilities

        # calculate reaction rate by summing over all representative detailed reactions
        detailed_reactions_consuming = set(
                [r for reactant in reactants 
                    for c in reactant.complexes 
                    for r in self.reactions_consuming(c)])

        reaction_rate = 0.0
        for r in detailed_reactions_consuming:
            # Calculate the probability that this detailed reaction decays to the products
            # P(r -> F), where F = products.
            self.calculate_reaction_decay_probabilities(r, products)

            # probability that the products of this detailed reaction decay into fates
            # that yield the condensed reaction
            product_probability = self._reaction_decay_probabilities[(r, products)]
            assert product_probability >= 0
            assert product_probability <= 1.000001

            # probability that the resting sets comprising the reactants of
            # the condensed reaction will be in the right configuration to
            # perform the detailed reaction

            #reactant_probabilities = 1
            #for a in r.reactants:
            #    fset = frozenset(self.SCC_containing(a))
            #    rset = self._set_to_fate[fset]
            #    reactant_probabilities *= self._stationary_distributions[rset][a]

            reactant_probabilities = times(self._stationary_distributions[self._set_to_fate[
                frozenset(self.SCC_containing(a))]][a] for a in r.reactants)
            assert reactant_probabilities >= 0
            assert reactant_probabilities <= 1.000001

            # rate of the detailed reaction
            k = r.const
            assert k > 0

            # overall contribution of detailed reaction r to rate of the condensed reaction 
            reaction_rate += reactant_probabilities * k * product_probability

            if isinstance(reaction_rate, complex):
                if reaction_rate.imag > 0:
                    log.warn("Detailed reaction {} contributes a complex rate of {} + {} " + \
                         " to condensed reaction {}.".format(r, reaction_rate.real, reaction_rate.imag, rxn))

        return reaction_rate

    def get_stationary_distribution(self, scc):
        """
        Take a strongly connected component and calculate the stationary distribution.

        Returns:
            [:obj:`dict()`]: Stationary distributions: dict['cplx'] = sdist (flt)
        """
        scc_set = frozenset(scc)
        scc_list = sorted(scc)
        L = len(scc)
    
        # assign a numerical index to each complex
        complex_indices = {c: i for (i, c) in enumerate(scc_list)}
    
        # find all reactions between complexes in this SCC (non-outgoing 1,1 reactions)
        reactions = [r for c in scc for r in self.reactions_consuming(c)
                     if (r.arity == (1, 1) and not is_outgoing(r, scc_set))]
    
        # T is the transition rate matrix, defined as follows:
        # T_{ij} = { rate(j -> i)       if  i != j
        #          { - sum_{j'} T_{j'i} if  i == j
        T = np.zeros((L, L))
    
        # add transition rates for each reaction to T
        for r in reactions:
            assert len(r.reactants) == 1
            assert len(r.products) == 1
            # r : a -> b
            # T_{b,a} = rate(r : a -> b)
            a = r.reactants[0]
            b = r.products[0]
            T[complex_indices[b]][complex_indices[a]] = r.const
    
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
        epsilon = 1e-5
        i = np.argmin(np.abs(w))
        if abs(w[i]) > epsilon:
            log.warn(
                ("Bad stationary distribution for resting set transition matrix. " +
                 "Eigenvalue found %f has magnitude greater than epsilon = %f. " +
                 "Markov chain may be periodic, or epsilon may be too high. Eigenvalues: %s") %
                (w(i), epsilon, str(w)))
        s = v[:, i]
    
        # check that the stationary distribution is good
        if not ((s >= 0).all() or (s <= 0).all()) : 
            log.error('Stationary distribution of resting set complex should not be an eigenvector of mixed sign. Condensed reaction rates may be incorrect.')

        s = s / np.sum(s)
        if not (abs(np.sum(s) - 1) < epsilon) :
            log.error('Stationary distribution of resting set complex should sum to 1 after normalization. Condensed reaction rates may be incorrect.')
    
        # return dict mapping complexes to stationary probabilities
        return {c: s[i] for (c, i) in complex_indices.items()}
    
    def get_exit_probabilities(self, scc):
        """
        """
        log.debug("Exit probabilities: {}".format([x.name for x in scc]))
        # build set and list of elements in SCC; assign a numerical index to each complex
        scc_set = frozenset(scc)
        scc_list = sorted(scc)
        complex_indices = {c: i for (i, c) in enumerate(scc_list)}
    
        # find all fast reactions involving complexes in this SCC
        reactions = [r for c in scc for r in self.reactions_consuming(c) if self.is_fast(r)]
    
        # sort reactions into internal and outgoing; assign a numerical index to each reaction
        r_outgoing = sorted([r for r in reactions if is_outgoing(r, scc_set)])
        r_internal = sorted([r for r in reactions if not is_outgoing(r, scc_set)])
        exit_indices = {r: i for (i, r) in enumerate(r_outgoing)}
    
        # L = # of complexes in SCC
        L = len(scc)
    
        # e = # of exit pathways
        e = len(r_outgoing)
    
        # add transition rates for each internal reaction
        T = np.zeros((L, L))
        for r in r_internal:
            assert len(r.reactants) == 1
            assert len(r.products) == 1
            a = r.reactants[0]
            b = r.products[0]
            T[complex_indices[a]][complex_indices[b]] = r.const

        # add transition rates for each outgoing reaction
        Te = np.zeros((L, e))
        for r in r_outgoing:
            a = r.reactants[0]
            Te[complex_indices[a]][exit_indices[r]] = r.const
    
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
        log.debug("P:\n{}".format(P))
    
        # extract the interior transition probabilities (Q_{LxL})
        Q = P[:, 0:L]
        log.debug("Q:\n{}".format(Q))
    
    
        # extract the exit probabilities (R_{Lxe})
        R = P[:, L:L + e]
        log.debug("R:\n{}".format(R))

    
        # calculate the fundamental matrix (N = (I_L - Q)^-1)
        N = np.linalg.inv(np.eye(L) - Q)
        log.debug("N:\n{}".format(N))

        # make sure all elements of fundamental matrix are >= 0
        if not (N >= 0).all() :  # --- commented out by EW (temporarily)
            log.error('Negative elements in fundamental matrix. Condensed reaction rates may be incorrect.')
    
        # calculate the absorption matrix (B = NR)
        B = np.dot(N, R)

        # --- added by EW as a weaker surrugate for the above, when necessary
        # assert (B >= 0).all()
    
        # return dict mapping tuples of (incoming complex, outgoing reaction)
        # to exit probabilities
        return {(c, r): B[i, j] for (c, i) in complex_indices.items()
                for (r, j) in exit_indices.items()}

class ReactionGraph(PepperCondensation):
    def __init__(self, *kargs, **kwargs):
        print('''# WARNING: peppercorn-v0.6: using depricated object name:''')
        print('''# Please rename peppercornenumerator.condense.ReactionGraph to peppercornenumerator.condense.PepperCondensation.  ''')
        super(ReactionGraph, self).__init__(*kargs, **kwargs)


class SetOfFates(object):
    """
    Corresponds to the set of multisets that can be reached from a given
    complex. Thin wrapper around the frozenset, mostly for prettier printing.
    """

    def __init__(self, states):
        # self._states = frozenset(states)
        self._states = frozenset([tuple(sorted(s)) for s in states])
        self._hash = None

    @property
    def states(self):
        return self._states

    def is_singleton(self):
        x = list(self._states)[0]
        return (len(self._states) == 1) and (len(x) == 1)

    def get_singleton(self):
        x = list(self._states)[0]
        return x[0]

    def __iter__(self):
        return iter(self.states)

    def __hash__(self):
        if(self._hash is None):
            self._hash = hash(self._states)
        return self._hash

    def __eq__(self, other):
        return sorted(list(self.states)) == sorted(list(other.states))

    def __str__(self):
        return "{ " + ", ".join(["{" + ", ".join(map(str, inner)) +
                                 "}" for inner in self.states]) + " }"

    def __len__(self):
        return len(self._states)

    def __repr__(self):
        return "SetOfFates([ " + ", ".join(
            ["(" + ", ".join(map(repr, inner)) + ")" for inner in self.states]) + " ])"


def is_outgoing(reaction, SCC_set):
    """
    Determines if the passed PepperReaction leads out of the passed SCC
    """
    if SCC_set is not set:
        SCC_set = set(SCC_set)

    # is subset of
    return not (set(reaction.products) <= SCC_set)


def times(iterable):
    """
    Multiplies each element of the given `iterable`, returns the product
    """
    return reduce(operator.mul, iterable)


def tuple_sum(iterable):
    """
    Given an iterable of tuples, concatenates each of the tuples
    and returns one big tuple. Precisely, for some list of tuples
    P = [p1, p2, p3 ...]::

        tuple_sum(P) = p1 + p2 + p3 + ...

    where + represents the concatenation of two tuples, e.g.
    (a1, a2) + (a3, a4) = (a1, a2, a3, a4).

    Example::

        tuple_sum([(1,2,3),(4,),(5,6)]) == (1, 2, 3, 4, 5, 6)

    """
    return reduce(operator.concat, tuple(iterable), tuple())


def tuple_sum_sort(iterable):
    return tuple_sort(tuple_sum(iterable))


def tuple_sort(tup):
    return tuple(sorted(tup))


def cartesian_product(iterable):
    """
    Gives the cartesian product of the passed `iterable`
    """
    # return map(tuple_sort, it.product(*iterable))
    return [tuple(x) for x in it.product(*iterable)]


def cartesian_sum(iterable):
    """
    Takes the cartesian sum over a set of multisets (counters). We define the
    cartesian sum X [+] Y between two sets X and Y as::

        X [+] Y = {x + y : x is in X and y is in Y}

    where + represents the sum of the multisets.

    Equivalently, in terms of the cartesian product X x Y, we can write X [+] Y as::

        X [+] Y = { (tuple_sum_{a in p} a) : p in X x Y }
    """
    # return map(tuple_sum,it.product(*iterable))
    return [tuple_sort(tuple_sum(x)) for x in cartesian_product(iterable)]


def get_reactions_consuming(complexes, reactions):
    """
    Returns a dict mapping each passed complex C to a list of reactions for which
    C is a reactant.
    """
    return dict((c, [r for r in reactions if (c in r.reactants)])
                for c in complexes)


def get_reactions_producing(complexes, reactions):
    """
    Returns a dict mapping each passed complex C to a list of reactions for which
    C is a product.
    """
    return dict((c, [r for r in reactions if (c in r.products)])
                for c in complexes)


def tarjans(complexes, reactions, reactions_consuming, is_fast):
    """
    Implementation of tarjans algorithm which segments SCCs of fast reactions
    """

    def strongconnect(cplx, index):
        """
        Traverse a complex and all outgoing reactions
        """
        cplx._index = index
        cplx._lowlink = index
        index += 1
        S.append(cplx)

        for reaction in reactions_consuming[cplx]:
            if reaction.arity == (1, 1) and is_fast(reaction):
                for product in reaction.products:

                    # Product hasn't been traversed; recurse
                    if product._index is None :
                        index = strongconnect(product, index)
                        cplx._lowlink = min(product._lowlink, cplx._lowlink)

                    # Product is in the current neighborhood
                    elif product in S:
                        cplx._lowlink = min(product._index, cplx._lowlink)

        if cplx._lowlink == cplx._index:
            scc = []
            while True:
                next = S.pop()
                scc.append(next)
                if next == cplx:
                    break
            SCCs.append(scc)
        return index

    index = 0
    S = []
    SCCs = []

    # Just to make sure we don't try to access an undeclared property, initialize
    # all indicies to None
    for complex in complexes:
        complex._index = None

    # Outer loop ensures that no complexes get left behind by the DFS
    for complex in complexes:
        if(complex._index is None):
            index = strongconnect(complex, index)

    return SCCs

# NOTE: temporary -- mainly a copy of PepperCondensation.get_stationary_distribution(scc)
def stationary_distribution(T, nodes = None):
    """
    Take a strongly connected component and calculate the stationary distribution.

    Args:
        T (numpy.matrix): a rate matrix 
        nodes (list, optional): A list of objects for which stationary distribution is to
            be determined.

    Returns:
        [:obj:`dict()`]: Stationary distributions: dict['cplx'] = sdist (flt)
    """
    L = len(T)
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
    epsilon = 1e-5
    i = np.argmin(np.abs(w))
    if abs(w[i]) > epsilon:
        log.warn(
            ("Bad stationary distribution for resting set transition matrix. " +
             "Eigenvalue found %f has magnitude greater than epsilon = %f. " +
             "Markov chain may be periodic, or epsilon may be too high. Eigenvalues: %s") %
            (w(i), epsilon, str(w)))
    s = v[:, i]

    # check that the stationary distribution is good
    if not ((s >= 0).all() or (s <= 0).all()) : 
        #for x,y in zip(s, nodes):
        #    print y, '"', y.kernel_string, x
        log.error('Stationary distribution of resting set complex should not be an eigenvector of mixed sign. Condensed reaction rates may be incorrect.')

    s = s / np.sum(s)
    if not (abs(np.sum(s) - 1) < epsilon) :
        log.error('Stationary distribution of resting set complex should sum to 1 after normalization. Condensed reaction rates may be incorrect.')

    return s
    ## return dict mapping complexes to stationary probabilities
    #return {c: s[i] for (c, i) in complex_indices.iteritems()}
    
