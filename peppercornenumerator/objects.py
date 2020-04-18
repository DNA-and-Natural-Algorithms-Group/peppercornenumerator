#
#  peppercornenumerator/objects.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

import numpy as np
from collections import namedtuple

from dsdobjects import clear_memory
from dsdobjects import DSDObjectsError, DSDDuplicationError
from dsdobjects.core import DL_Domain, DSD_Complex, DSD_Reaction, DSD_Macrostate
from dsdobjects.utils import split_complex, convert_units
from dsdobjects.utils import make_pair_table, pair_table_to_dot_bracket 

class PepperObjectsError(Exception):
    """
    pepperobjects error class.
    """

    def __init__(self, message, *kargs):
        if kargs:
            self.message = "{} [{}]".format(message, ', '.join(map(str,kargs)))
        else :
            self.message = message
        super(PepperObjectsError, self).__init__(self.message)

class PepperDomain(DL_Domain):
    """
    Represents a single domain. We allow several options for specifying domain
    properties. Domains might have an explicit integer (bp) length, or may be
    designated as short or long. If the latter method is used, the code will use
    the relevant constant as the integer domain length.
    """

    def __new__(cls, name, dtype=None, length=None):
        # The new method returns the present instance of an object, if it exists
        self = DL_Domain.__new__(cls)
        try:
            super(PepperDomain, self).__init__(name, dtype, length)
        except DSDDuplicationError as e :
            other = e.existing
            if dtype and (other.dtype != dtype) :
                raise PepperObjectsError('Conflicting dtype assignments for {}: "{}" vs. "{}"'.format(
                    name, dtype, other.dtype))
            elif length and (other.length != length) :
                raise PepperObjectsError('Conflicting length assignments for {}: "{}" vs. "{}"'.format(
                    name, length, other.length))
            return e.existing
        self._nucleotides = None
        return self

    def __init__(self, name, dtype=None, length=None):
        # Remove default initialziation to get __new__ to work
        pass

    @property
    def nucleotides(self):
        return self._nucleotides

    @nucleotides.setter
    def nucleotides(self, value):
        self._nucleotides = value

    @property
    def complement(self):
        # If we initialize the complement, we need to know the class.
        if self._complement is None:
            cname = self._name[:-1] if self.is_complement else self._name + '*'
            if cname in DL_Domain.MEMORY:
                self._complement = DL_Domain.MEMORY[cname]
            else :
                self._complement = PepperDomain(cname, self.dtype, self.length)
        return self._complement

    def can_pair(self, other):
        """
        Returns True if this domain is complementary to the argument.
        """
        return self == ~other

    @property
    def identity(self):
        """
        Returns the identity of this domain, which is its name without a
        complement specifier (i.e. A and A* both have identity A).
        """
        return self._name[:-1] if self._name[-1] == '*' else self._name

    @property
    def is_complement(self):
        """
        Returns true if this domain is a complement (e.g. A* rather than A),
        false otherwise.
        """
        return self._name[-1:] == '*'

class PepperComplex(DSD_Complex):
    """
    Peppercorn complex object. 

    Overwrites some functions with new names, adds some convenient stuff..
    """

    PREFIX = 'e'
    CONCENTRATION = namedtuple('concentration', 'mode value unit')

    @staticmethod
    def clear_memory(memory=True, names=True, ids=True):
        if memory:
            DSD_Complex.MEMORY = dict()
        if names:
            DSD_Complex.NAMES = dict()
        if ids:
            DSD_Complex.ID = dict()

    def __init__(self, sequence, structure, name='', prefix='', memorycheck=True):
        try :
            if not prefix :
                prefix = PepperComplex.PREFIX
            super(PepperComplex, self).__init__(sequence, structure, name, prefix, memorycheck)
        except DSDObjectsError :
            backup = 'enum' if prefix != 'enum' else 'pepper'
            super(PepperComplex, self).__init__(sequence, structure, name, backup, memorycheck)
            log.warning('Complex name existed, prefix has been changed to: {}'.format(backup))
        
        # Peppercorn IO:
        self._concentration = None
        self._elevation = None
        assert self.is_domainlevel_complement

    @property
    def concentration(self):
        if self._concentration is not None:
            return self._concentration
        return None

    @concentration.setter
    def concentration(self, trip):
        if trip is None:
            self._concentration = None
        else:
            (mode, value, unit) = trip
            assert isinstance(value, (int, float))
            self._concentration = PepperComplex.CONCENTRATION(mode, value, unit)

    def concentrationformat(self, out):
        mod = self._concentration.mode
        val = self._concentration.value
        uni = self._concentration.unit
        val = convert_units(val, uni, out) 
        return PepperComplex.CONCENTRATION(mod, val, out)

    @property
    def pair_table(self):
        return super(PepperComplex, self).pair_table
    
    @pair_table.setter
    def pair_table(self, pt):
        self._pair_table = pt

    def full_string(self):
        return "Complex(%s): %s %s" % (
            self.name, str(self.strands), str(self.structure))

    def get_structure(self, loc):
        return self.get_paired_loc(loc)

    def triple(self, *loc):
        # overwrite standard func
        return (self.get_domain(loc), self.get_paired_loc(loc), loc)

    @property
    def available_domains(self):
        ad = []
        for (x,y) in self.exterior_domains:
            ad.append((self.get_domain((x,y)), x, y))
        return ad

    @property
    def pk_domains(self):
        pd = []
        for (x,y) in self.enclosed_domains:
            pd.append((self.get_domain((x,y)), x, y))
        return pd


    def rotate_location(self, loc, n=None):
        return self.rotate_pairtable_loc(loc, n)

    def split(self):
        """ Split PepperComplex into disconneted components.
        """
        if self.is_connected:
            return [self]
        else :
            ps = self.lol_sequence
            pt = self.pair_table
            parts = split_complex(ps, pt)
            cplxs = []
            # assign new_complexes
            for (se,ss) in parts:
                try:
                    cplxs.append(PepperComplex(se, ss))
                except DSDDuplicationError as e:
                    cplxs.append(e.existing)
            return sorted(cplxs)

class PepperReaction(DSD_Reaction):
    RTYPES = set(['condensed', 'open', 'bind11', 'bind21', 'branch-3way', 'branch-4way'])

    def __init__(self, *kargs, **kwargs):
    #def __init__(self, reactants, products, rtype=None, rate=None, memorycheck=True):
        super(PepperReaction, self).__init__(*kargs, **kwargs)
        if self._rtype not in PepperReaction.RTYPES:
            try:
                del DSD_Reaction.MEMORY[self.canonical_form]
            except KeyError:
                pass
            raise DSDObjectsError('Reaction type not supported! ' + 
            'Set supported reaction types using PepperReaction.RTYPES')

        self._reverse_reaction = None

        # Store the 4 relevant Loop() objects (reactant based):
        #   * initial-locus
        #   * target-locus
        #   * x-linker
        #   * y-linker
        self.meta = None
        self.rotations = None

    @property
    def reverse_reaction(self):
        return self._reverse_reaction

    @reverse_reaction.setter
    def reverse_reaction(self, rxn):
        def rev_rtype(rtype, arity):
            """Returns the reaction type of a corresponding reverse reaction. """
            if rtype == 'open' and arity == (1,1):
                return 'bind11'
            elif rtype == 'open' and arity == (1,2):
                return 'bind21'
            elif rtype == 'bind11' or rtype == 'bind21':
                return 'open'
            else:
                return rtype

        if rxn is not False:
            assert rxn.rtype == rev_rtype(self.rtype, self.arity)
        self._reverse_reaction = rxn

    def full_string(self, molarity='M', time='s'):
        """Prints the reaction in PIL format.
        Reaction objects *always* specify rate in /M and /s.  """

        if self.rate :
            newunits = [molarity] * (self.arity[0] - 1) + [time]
            newrate = self.rateformat(newunits)
            rate = newrate.constant
            assert newunits == newrate.units
            units = ''.join(map('/{}'.format, newrate.units))
        else :
            rate = float('nan')
            units = ''

        if self.rtype :
            return '[{:14s} = {:12g} {:4s} ] {} -> {}'.format(self.rtype, rate, units,
                    " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))
        else :
            return '[{:12g} {:4s} ] {} -> {}'.format(rate, units,
                    " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))

class PepperMacrostate(DSD_Macrostate):
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

