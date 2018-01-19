
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
import numpy as np

from dsdobjects import clear_memory
from dsdobjects import DL_Domain, DSD_Complex, DSD_Reaction, DSD_RestingSet
from dsdobjects import DSDObjectsError, DSDDuplicationError
from dsdobjects.utils import split_complex 

# not needed here, but passing it on...
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
        except DSDDuplicationError, e :
            other = e.existing
            if dtype and (other.dtype != dtype) :
                raise PepperObjectsError('Conflicting dtype assignments for {}: "{}" vs. "{}"'.format(
                    name, dtype, other.dtype))
            elif length and (other.length != length) :
                raise PepperObjectsError('Conflicting length assignments for {}: "{}" vs. "{}"'.format(
                    name, length, other.length))
            return e.existing
        return self

    def __init__(self, name, dtype=None, length=None):
        # Remove default initialziation to get __new__ to work
        pass

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
            logging.warning('Complex name existed, prefix has been changed to: {}'.format(backup))

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

    def get_strand(self, loc):
        """
        Returns the strand at the given index in this complex
        """
        if(loc is not None):
            return self._strands[loc]
        return None

    @property
    def available_domains(self):
        ad = []
        for (x,y) in self.exterior_domains:
            ad.append((self.get_domain((x,y)), x, y))
        return ad

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
                except DSDDuplicationError, e:
                    cplxs.append(e.existing)
            return sorted(cplxs)

    def __cmp__(self, other):
        """
        Two complexes are compared on the basis of their complexes
        """
        return cmp(self.canonical_form, other.canonical_form)

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

    def __cmp__(self, other):
        """
        ReactionPathway objects are sorted by their canonical form.
        """
        return cmp(self.canonical_form, other.canonical_form)

    @property
    def full_string(self):
        """prints the formal chemical reaction."""
        if self.rtype :
            return '[{:14s} = {:12g} {:4s} ] {} -> {}'.format(self.rtype, self.rate, self.rateunits,
                    " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))
        else :
            return '[{:12g} {:4s} ] {} -> {}'.format(self.rate, self.rateunits,
                    " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))

class PepperMacrostate(DSD_RestingSet):
    def __init__(self, *kargs, **kwargs):
        super(PepperMacrostate, self).__init__(*kargs, **kwargs)
        self._internal_reactions = set()
        self._exit_reactions = set()

    def __str__(self):
        return self.name

    def __cmp__(self, other):
        """
        Two resting sets are compared on the basis of their complexes
        """
        return cmp(self.canonical_form, other.canonical_form)

    def __len__(self):
        """
        The number of species in a resting set
        """
        return len(self._complexes)


    def add_reaction(self, rxn):
        assert len(rxn.reactants) == 1
        assert rxn.reactants[0] in self._complexes
        is_internal = rxn.products[0] in self._complexes
        if is_internal:
            assert len(rxn.products) == 1
            self._internal_reactions.add(rxn)
        else :
            self._exit_reactions.add(rxn)

    #def reactions_within(self):
    #    raise NotImplementedError

    #def reactions_outgoing(self):
    #    # Should be empty!
    #    raise NotImplementedError
    
    @property
    def is_resting(self):
        return len(self._exit_reactions) == 0

    @property
    def is_transient(self):
        return not self.is_resting

    def is_exit(self, cplx):
        return any(map(lambda r: cplx in r.reactants, self._exit_reactions))

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
        for i in xrange(L):
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
            logging.warn(
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
            logging.error('Stationary distribution of resting set complex' +
                    'should not be an eigenvector of mixed sign.')

        s = s / np.sum(s)
        if not (abs(np.sum(s) - 1) < epsilon) :
            logging.error('Stationary distribution of resting set complex' +
                    'should sum to 1 after normalization. Condensed reaction' +
                    'rates may be incorrect.')

        return zip(self._complexes, s)

