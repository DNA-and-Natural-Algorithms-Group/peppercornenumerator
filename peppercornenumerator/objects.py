#
#  peppercornenumerator/objects.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

from itertools import chain

from dsdobjects import SingletonError, clear_singletons, show_singletons
from dsdobjects.base_classes import ObjectInitError, DomainS, ComplexS, MacrostateS, ReactionS

def clear_memory():
    # This is for unittests only!  The memory of a Singleton clears
    # automatically if there is no hard reference to the object.
    clear_singletons(PepperReaction)
    clear_singletons(PepperMacrostate)
    clear_singletons(PepperComplex)
    clear_singletons(PepperDomain)

def show_memory():
    # This is for unittests only!  The memory of a Singleton clears
    # automatically if there is no hard reference to the object.
    for x in chain(show_singletons(PepperReaction),
                   show_singletons(PepperMacrostate),
                   show_singletons(PepperComplex),
                   show_singletons(PepperDomain)):
        print(x)

class PepperDomain(DomainS):
    def __init__(self, *args, **kwargs):
        super(PepperDomain, self).__init__(*args, **kwargs)
        self.sequence = None

    def can_pair(self, other):
        """
        Returns True if this domain is complementary to the argument.
        """
        return self is ~other

class PepperComplex(ComplexS):
    PREFIX = 'e'
    @property
    def available_domains(self):
        ad = []
        for x, y in self.exterior_domains:
            ad.append((self.get_domain((x,y)), x, y))
        return ad

class PepperMacrostate(MacrostateS):
    def __init__(self, *kargs, **kwargs):
        super(PepperMacrostate, self).__init__(*kargs, **kwargs)
        self._internal_reactions = set()
        self._exit_reactions = set()

        # all reactions consuming a particular species
        #self._reactions_consuming = dict()
        self._stationary_distribution = []
        self._exit_probabilities = []

    @property
    def is_resting(self):
        """True iff it is a resting macrostate. """
        return len(self._exit_reactions) == 0

    @property
    def is_transient(self):
        """True iff it is a transient macrostate. """
        return not self.is_resting

    def add_reaction(self, rxn):
        """Adds a reaction to the macrostate. """
        assert len(rxn.reactants) == 1
        assert rxn.reactants[0] in self._complexes
        is_internal = rxn.products[0] in self._complexes
        if is_internal:
            assert len(rxn.products) == 1
            self._internal_reactions.add(rxn)
        else :
            self._exit_reactions.add(rxn)

        self._stationary_distribution = []
        self._exit_probabilities = []

    def get_exit_probabilities(self, warnings=True):
        """
        """
        # build set and list of elements in SCC; assign a numerical index to each complex
        scc_set = frozenset(self._complexes)
        complex_indices = {c: i for (i, c) in enumerate(self._complexes)}
    
        # sort reactions into internal and outgoing; assign a numerical index to each reaction
        r_outgoing = self._exit_reactions
        r_internal = self._internal_reactions
        exit_indices = {r: i for (i, r) in enumerate(r_outgoing)}
    
        # L = # of complexes in SCC
        L = len(self._complexes)
    
        # e = # of exit pathways
        e = len(r_outgoing)
    
        # add transition rates for each internal reaction
        T = np.zeros((L, L))
        for r in r_internal:
            assert len(r.reactants) == 1
            assert len(r.products) == 1
            a = r.reactants[0]
            b = r.products[0]
            T[complex_indices[a]][complex_indices[b]] = r.rate

        # add transition rates for each outgoing reaction
        Te = np.zeros((L, e))
        for r in r_outgoing:
            a = r.reactants[0]
            Te[complex_indices[a]][exit_indices[r]] = r.rate
    
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
    
        # extract the interior transition probabilities (Q_{LxL})
        Q = P[:, 0:L]
    
        # extract the exit probabilities (R_{Lxe})
        R = P[:, L:L + e]
    
        # calculate the fundamental matrix (N = (I_L - Q)^-1)
        N = np.linalg.inv(np.eye(L) - Q)

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


    def get_stationary_distribution(self, warnings=True):
        """
        Take a strongly connected component and calculate the stationary distribution.

        Args:
            T (numpy.matrix): a rate matrix 
            nodes (list, optional): A list of objects for which stationary distribution is to
                be determined.

        Returns:
            [:obj:`dict()`]: Stationary distributions: dict['cplx'] = sdist (flt)
        """
        if not self._stationary_distribution :

            reactions = self._internal_reactions

            # Initialize a Transition Matrix where numerical index corresponds to each complex
            L = len(self._complexes)
            indices = {c: i for (i, c) in enumerate(self._complexes)}
            T = np.zeros((L, L))

            for rxn in reactions:
                # r : a -> b
                # T_{b,a} = rate(r : a -> b)
                #if rxn.reactants[0] in indices and rxn.products[0] in indices:
                a = rxn.reactants[0]
                b = rxn.products[0]
                T[indices[b]][indices[a]] = rxn.rate

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
            if warnings and not ((s >= 0).all() or (s <= 0).all()) : 
                #for cl,sd in zip(self._complexes, s):
                #    print(cl, '=', cl.kernel_string, sd)
                #for rxn in reactions:
                #    print(rxn, rxn.rate)
                log.error('Stationary distribution of resting set complex' +
                        'should not be an eigenvector of mixed sign.')

            self._stationary_distribution = s / np.sum(s)

            if not (abs(np.sum(self._stationary_distribution) - 1) < epsilon) :
                log.error('Stationary distribution of resting set complex' +
                        'should sum to 1 after normalization. Condensed reaction' +
                        'rates may be incorrect.')

        return list(zip(self._complexes, self._stationary_distribution))

class PepperReaction(ReactionS):
    RTYPES = set(['condensed', 'open', 'bind11', 'bind21', 'branch-3way', 'branch-4way'])

    @property
    def const(self):
        log.warning('deprecated property PepperReaction.const')
        return self.rate_constant[0]

    @const.setter
    def const(self, value):
        log.warning('deprecated property PepperReaction.const')
        print(value)
        self.rate_constant = value
        print(self.rate_constant[0])


class Loop:
    """ A (part of a) single open or closed loop region.

    Args: 
        loop [(domain, paired_loc, domain_loc), ...]: A list of domains and
            their structure: 
                [0] A Domain Object, 
                [1] its binding partner position (or None), 
                [2] the domain position.
            None marks the 3' or 5' end of a strand.
    """
    def __init__(self, loop):
        self._parts = loop

        bases = 0
        stems = 0
        is_open = False
        stem_list = set()
        for step in loop:
            if step is None:
                if is_open:
                    raise ObjectInitError('Double strand break in Loop Object.')
                is_open = True
            else:
                (dom, ploc, dloc) = step
                if ploc is None:
                    bases += len(dom)
                elif ploc in stem_list:
                    raise ObjectInitError('Double stem count in Loop Object.')
                else:
                    stems += 1
                    stem_list.add(dloc)
        # cache properties
        self._is_open = is_open
        self._bases = bases
        self._stems = stems

    @property
    def domains(self):
        return (part[0] if part is not None else None for part in self._parts)

    @property
    def pair_locs(self):
        return (part[1] if part is not None else None for part in self._parts)

    @property
    def domain_locs(self):
        return (part[2] if part is not None else None for part in self._parts)

    @property
    def parts(self):
        """ list: (dom, struct, loc) tuples associated with this loop.  """
        return self._parts[:]

    @property
    def stems(self):
        """ int: the number of stems. """
        return self._stems

    @property
    def bases(self):
        """ int: the number of bases. """
        return self._bases

    @property
    def is_open(self):
        """ bool: True if there is a strand break in the loop. """
        return self._is_open

    @property
    def dlength(self):
        """ int: the sum over all involved domain lengths. """
        assert None not in self.domains
        return sum(len(d) for d in self.domains)

    @property
    def llength(self):
        """ flt: the *linker length* expressed in number of nucleotides. """
        if self.is_open:
            return float('inf')
        else:
            # approximate number of single-stranded nucleotides to span a stem
            stemL = 2.0 / 0.43
            return 1 + self.bases + self.stems + (stemL * self.stems)

    def __repr__(self):
        return f"Loop({self.parts})"

