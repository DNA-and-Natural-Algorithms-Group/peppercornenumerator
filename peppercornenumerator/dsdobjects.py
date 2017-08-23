# -*- coding: utf-8 -*-
#
# dsdobjects.py
#
# copy and modify together with test_dsdobjects.py.
#
# Written by Stefan Badelt (badelt@caltech.edu)
#
#   This file contains adapted code from the "DNAObjects" Python package coded
#   by Joseph Berleant and Joseph Schaeffer, the "Peppercornenumerator" Python 
#   packate coded by Casey Grun, Kathrik Sarma, Erik Winfree and Brian Wolfe and 
#   the "Nuskell" DNA strand displacement compiler coded by Seung Woo Shin.
#
# Distributed under the MIT License, use at your own risk.
#

"""dsdobjects.py

Base-classes and utility functions for domain-level strand displacement programming.

"""

import networkx as nx
from collections import Counter
from nuskell.parser import parse_pil_file

class DSDObjectsError(Exception):
    """
    dna_objects.py error class.
    """

    def __init__(self, message, *kargs):
        if kargs:
            self.message = "{} [{}]".format(message, ', '.join(map(str,kargs)))
        else :
            self.message = message
        self._solution = None
        self._rotations = None
        super(DSDObjectsError, self).__init__(self.message)

    @property
    def solution(self):
        return self._solution

    @solution.setter
    def solution(self, value):
        self._solution = value

    @property
    def rotations(self):
        return self._rotations

    @rotations.setter
    def rotations(self, value):
        self._rotations = value


def reset_names(): # TODO: find a better name
    """
    Reset all names, counters of dnaobjects.
    """
    DSD_Complex.id_counter = 0
    DSD_Complex.dictionary = dict()
    DSD_Complex.canonical_forms = dict()
    DSD_Domain.id_counter = 0
    DSD_Domain.dictionary = dict()
    DSD_Reaction.dictionary = dict()
    DSD_Reaction.rtypes = set([None])

def make_pair_table(ss, strand_break='+', ignore=set('.')):
    """Return a secondary struture in form of pair table:

    Args:
      ss (str or list): secondary structure in dot-bracket format.
      cut (str, optional): character that defines a cut-point. Defaults to '+'.
      ignore (set, optional): a list of characters that are ignored. Defaults to ['.']

    Example:
                             0,0  0,1  0,2   0,3   0,4     1,0   1,1
             "...((+))" -> [[None,None,None,(1,0),(1,1)],[(0,4),(0,3)]]

    Raises:
       DSDObjectsError: Too many closing parenthesis ')' in secondary structure.
       DSDObjectsError: Too many opening parenthesis '(' in secondary structure.
       DSDObjectsError: Unexpected character in sequence: "{}".

    Returns:
      [list]: A pair-table as list of lists.
    """

    assert len(strand_break) == 1
    assert '.' in ignore

    # Return value
    pair_table = []

    # Helpers
    stack = []
    strand_index = 0
    domain_index = 0

    strand = []
    pair_table.append(strand)
    for char in ss:
        if char == strand_break :
            strand_index += 1
            domain_index = 0
            strand = []
            pair_table.append(strand)
            continue

        if char == '(' :
            strand.append(None)
            stack.append((strand_index, domain_index))
            domain_index += 1
 
        elif char == ")" :
            try:
                loc = stack.pop()
            except IndexError as e:
                raise DSDObjectsError("Too few closing parenthesis ')' in secondary structure.")
            strand.append(loc)
            pair_table[loc[0]][loc[1]] = (strand_index, domain_index)
            domain_index += 1

        elif char in set(ignore): # unpaired
            strand.append(None)
            domain_index += 1

        else :
            raise DSDObjectsError("Unexpected character in sequence: '{}'.".format(char))

    if len(stack) > 0 :
        raise DSDObjectsError("Too few opening parenthesis '(' in secondary structure.")

    return pair_table 

def pair_table_to_dot_bracket(pt, strand_break='+'):
    """
    """
    assert len(strand_break) == 1
    out = []
    for si, strand in enumerate(pt):
        for di, pair in enumerate(strand):
            if pair is None:
                out.append('.')
            else :
                locus = (si, di)
                if locus < pair :
                    out.append('(')
                else :
                    out.append(')')
        out.append(strand_break)
    
    return out[:-1] #remove last '+'

def make_lol_sequence(seq):
    indices = [-1] + [i for i, x in enumerate(seq) if x == "+"]
    indices.append(len(seq))
    return [seq[indices[i - 1] + 1 : indices[i]] for i in range(1, len(indices))]


def make_loop_index(ptable):
    """
    number loops and assign each position its loop-number
    handy for checking which pairs can be added
    """
    loop_index = []
    exterior = set()

    stack = []
    (cl, nl) = (0, 0)

    for si, strand in enumerate(ptable):
        loop = []
        loop_index.append(loop)
        for di, pair in enumerate(strand):
            loc = (si, di)
            if pair is None:
                pass
            elif loc < pair : # '('
                nl += 1
                cl = nl
                stack.append(loc)
            loop.append(cl)
            if pair and pair < loc : # ')'
                _ = stack.pop()
                try :
                    ploc = stack[-1]
                    cl = loop_index[ploc[0]][ploc[1]]
                except IndexError:
                    cl = 0
        # strand break
        if cl in exterior :
            raise DSDObjectsError('Complexes not connected.')
        else :
            exterior.add(cl)
    # ptable end
    return loop_index, exterior

def split_complex(lol_seq, ptable):
    """
    number loops and assign each position its loop-number
    handy for checking which pairs can be added

    NOTE: modifies its arguments!
    """
    loop_index = []
    exterior = set()

    complexes = dict()
    splice = []

    stack = []
    (cl, nl) = (0, 0)

    # Identify exterior loops using the loop index function,
    # rewrite ptable to contain characters instead of paired locus.
    for si, strand in enumerate(ptable):
        loop = []
        loop_index.append(loop)
        for di, pair in enumerate(strand):
            loc = (si, di)
            if pair is None:
                ptable[si][di] = '.'
            elif loc < pair : # '('
                nl += 1
                cl = nl
                stack.append(loc)
                ptable[si][di] = '('
            loop.append(cl)
            if pair and pair < loc : # ')'
                _ = stack.pop()
                try :
                    ploc = stack[-1]
                    cl = loop_index[ploc[0]][ploc[1]]
                except IndexError:
                    cl = 0
                ptable[si][di] = ')'

        # store the strand_index where a given loop-cut starts
        if cl in complexes :
            start = complexes[cl]
            end = si
            splice.append((start,end))
            complexes[cl] = si
        else :
            complexes[cl] = si
    # for ptable end

    parts = []
    for (s,e) in splice:
        new_seq = lol_seq[s+1:e+1]
        new_pt = ptable[s+1:e+1]

        lol_seq[s+1:e+1] = [[] for x in range(s+1,e+1)]
        ptable[s+1:e+1] = [[] for x in range(s+1,e+1)]

        # change format from lol to regular list format
        new_seq=reduce(lambda a,b:a+['+']+b, new_seq)
        new_pt=reduce(lambda a,b:a+['+']+b, new_pt)
        parts.append((new_seq, new_pt))

    # append the remainder (or original) sequence
    new_seq = lol_seq[:]
    new_pt = ptable[:]
    new_seq=reduce(lambda a,b:a+['+']+b if b else a+b,new_seq)
    new_pt=reduce(lambda a,b:a+['+']+b if b else a+b,new_pt)
    parts.append((new_seq, new_pt))
    return parts

def _find(l, key):
    for i in range(len(l)):
        if l[i] == key:
            return i
    return None

def _flatten(l):
    if l == []:
        return l
    if isinstance(l[0], list):
        return _flatten(l[0]) + _flatten(l[1:])
    return l[:1] + _flatten(l[1:])

class SequenceConstraint(object):
    """ A nucleic acid sequence constraint in IUPAC format.
    """
    def __init__(self, sequence, molecule='DNA'):
        assert molecule == 'DNA' or molecule == 'RNA'
        self.ToU = 'T' if molecule == 'DNA' else 'U'

        self._molecule = molecule
        self._sequence = list(sequence)

    @property
    def constraint(self):
        return ''.join(self._sequence)

    @property
    def complement(self):
        return ''.join(map(self._iupac_complement, self._sequence))

    @property
    def wc_complement(self):
        return ''.join(map(self._wc_complement, self._sequence))

    @property
    def reverse_complement(self):
        return ''.join(map(self._iupac_complement, self._sequence[::-1]))

    @property
    def reverse_wc_complement(self):
        return ''.join(map(self._wc_complement, self._sequence[::-1]))

    def add_constraint(self, con):
        """Apply nucleotide sequence constraint and check if it is compatible.

        Args:
          con (str): A new sequence constraint.

        Raises:
          DSDObjectsError: Sequence constraints differ in length.
          DSDObjectsError: Sequence constraints cannot be satisfied.
        """
        con = list(con) # just to make sure...

        if len(self._sequence) != len(con):
            raise DSDObjectsError("Sequence constraints differ in length.")

        new = self._merge_constraints(self._sequence, con)

        if '' in new:
            raise DSDObjectsError("Sequence constraints cannot be satisfied.")
        else:
            self._sequence = new

    def __add__(self, other):
        assert isinstance(other, SequenceConstraint)
        assert self._molecule == other._molecule
        new = SequenceConstraint(self.constraint, self._molecule)
        new.add_constraint(other.constraint)
        return new

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __eq__(self, other):
        return self.constraint == other.constraint

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return self.constraint

    def __len__(self):
        return len(self._sequence)

    def __invert__(self):
        self._sequence = list(self.complement)
        return self

    def _merge_constraints(self, con, con2):
        """Return a new list of unified constraints. """
        return map(self._iupac_union, zip(con, con2))

    def _iupac_to_bases(self, nuc):
        # This is more of educational purpose, use _iupac_bin and _bin_iupac
        T = self.ToU
        iupac_dict = {
            'A': 'A',
            'C': 'C',
            'G': 'G',
             T :  T,
            'R': 'AG',      # purine
            'Y': 'C' + T,   # pyrimidine
            'S': 'CG',      # strong
            'M': 'AC',
            'W': 'A' + T,   # weak
            'K': 'G' + T,
            'V': 'ACG',     # not T
            'H': 'AC' + T,  # not G
            'D': 'AG' + T,  # not C
            'B': 'CG' + T,  # not A
            'N': 'ACG' + T}
        return iupac_dict[nuc]

    def _iupac_complement(self, nuc):
        T = self.ToU
        neighbor_dict = {   # ACGT => ACGT
            'A': T,         # 1000 => 0001
            'C': 'G',       # 0100 => 0010
            'G': 'Y',       # 0010 => 0101
             T : 'R',       # 0001 => 1010
            'T': 'R',       # 0001 => 1010
            'R': 'Y',       # 1010 => 0101
            'Y': 'R',       # 0101 => 1010
            'S': 'B',       # 0110 => 0111
            'M': 'K',       # 1100 => 0011
            'W': 'D',       # 1001 => 1011
            'K': 'N',       # 0011 => 1111
            'V': 'B',       # 1110 => 0111
            'H': 'D',       # 1101 => 1011
            'D': 'N',       # 1011 => 1111
            'B': 'N',       # 0111 => 1111
            'N': 'N'}       # 1111 => 1111
        return neighbor_dict[nuc]

    def _wc_complement(self, nuc):
        T = self.ToU
        neighbor_dict = {   # ACGT => ACGT
            'A': T,         # 1000 => 0001
            'C': 'G',       # 0100 => 0010
            'G': 'C',       # 0010 => 0100
             T : 'A',       # 0001 => 1000
            'N': 'N'}       # 1111 => 1111
        return neighbor_dict[nuc]

    def _iupac_union(self, nucs):
        # Return the maximal common constraint
        u = 'N'
        for n in nucs:
            u = self._bin_iupac(self._iupac_bin(u) & self._iupac_bin(n))
        return u

    def _iupac_bin(self, nuc):
        T = self.ToU
        iupac_bin_dict = {  # ACGT
            'A': 8,         # 1000,
            'C': 4,         # 0100,
            'G': 2,         # 0010,
             T : 1,         # 0001,
            'R': 10,        # 1010,  # purine
            'Y': 5,         # 0101,  # pyrimidine
            'S': 6,         # 0110,
            'M': 12,        # 1100,
            'W': 9,         # 1001,
            'K': 3,         # 0011,
            'V': 14,        # 1110,  # not T
            'H': 13,        # 1101,  # not G
            'D': 11,        # 1011,  # not C
            'B': 7,         # 0111,  # not A
            'N': 15}        # 1111,
        return iupac_bin_dict[nuc]

    def _bin_iupac(self, nuc):
        T = self.ToU
        bin_iupac_dict = [  # ACGT
            '',             # 0000  0
            T,              # 0001  1
            'G',            # 0010  2
            'K',            # 0011  3
            'C',            # 0100  4
            'Y',            # 0101  5
            'S',            # 0110  6
            'B',            # 0111  7
            'A',            # 1000  8
            'W',            # 1001  9
            'R',            # 1010 10
            'D',            # 1011 11
            'M',            # 1100 12
            'H',            # 1101 13
            'V',            # 1110 14
            'N']            # 1111 15
        return bin_iupac_dict[nuc]

class DSD_Domain(object):
    """Nucleic acid domain sequence.

    A domain is a sequence of consecutive nucleotides (or domains). Domains can
    form different secondary structures, however, they always fold as one entity,
    just like single nucleotides.

    Implementation:
      - Domains with the same name have the same nucleotide sequence.
      - Domains with the same nucleotide sequence have the same name.
      - Every DSD_Domain() has exactly one ComplementaryDomain().

    Not implemented:
      - A domain cannot be complementary to multiple domains, e.g. due to
        suboptimal hybridizations or wobble-base-pairs.

    Globals:
      id_counter (int): An automatically assigned unique descriptor that allows
        for accessing/modifying every occurrence in the System.

    Args:
      sequence (list): A list of nucleotides or Domain() objects.
      name (str, optional): Name of this domain. If not specified, an automatic
        name is generated, defaults to ''
      prefix (str, optional): A prefix for automated naming of Domains. Has to be
        set if no name is provided. Defaults to 'd' for 'domain'. Usually, the
        prefix 't' is used for 'toehold domains', and 'h' for 'histoy domains'.

    Raises:
      DSDObjectsError: Domain name must not end with '*'!
      DSDObjectsError: Domain name must not end with a digit!
      DSDObjectsError: Domain prefix must not be empty!
      DSDObjectsError: Domain prefix must not end with a digit!
      DSDObjectsError: Must pass a sequence of constraints or domain objects.

    Note:
      The Domain function does not allow names or prefixes ending with digits or
      '*'. Proximal digits are reseved to name domains automatically using their
      unique domain-IDs. A mixture of both naming modes is forbidden.  Beware,
      that it is still possible to initialize two different domains with the same
      name, but that the information that these domain are different is lost when
      writing them to a file.
    """

    id_counter = 0
    dictionary = dict()

    def __init__(self, sequence, name='', prefix='d', is_complement=False):
        # Assign name #
        if name:
            assert is_complement == (name[-1] == '*')
            self._name = name
        else:
            if prefix == '':
                raise DSDObjectsError('DSD_Domain prefix must not be empty!')
            if prefix[-1].isdigit():
                raise DSDObjectsError('DSD_Domain prefix must not end with a digit!')
            self._name = prefix + str(DSD_Domain.id_counter)
            if is_complement: self._name += '*'
            DSD_Domain.id_counter += 1

        if self._name in DSD_Domain.dictionary:
            raise DSDObjectsError('Duplicate DSD_Domain name!', self._name)

        DSD_Domain.dictionary[self._name] = self
        self._identity = self._name if not is_complement else self._name[:-1]

        # Assign domain sequence
        if not sequence or not isinstance(sequence, list):
            raise DSDObjectsError("Must pass a sequence of constraints or domain objects.")
        elif all(isinstance(s, str) for s in sequence):
            self._sequence = sequence
        elif all(isinstance(s, DSD_Domain) for s in sequence) :
            self._sequence = sequence
        elif all(isinstance(s, (str, DSD_Domain)) for s in sequence):
            raise DSDObjectsError(
                'Mixing of composite and non-composite domains is not supported.')
        else :
            raise DSDObjectsError('Invalid sequence argument.')

        # Resolve sequence constraints for complementary domains:
        if is_complement and self._identity in DSD_Domain.dictionary :
            other = DSD_Domain.dictionary[self._identity]
            #self.update_constraints(other.sequence[::-1])
        elif not is_complement and self._identity + '*' in DSD_Domain.dictionary :
            other = DSD_Domain.dictionary[self._identity + '*']
            #self.update_constraints(other.sequence[::-1])

    @property
    def name(self):
        """str: name of the domain."""
        return self._name

    @property
    def sequence(self):
        """list: DSD_Domain sequence."""
        return self._sequence

    @property
    def sequence_string(self):
        """list: DSD_Domain sequence."""
        return ''.join(map(str,self._sequence))

    @property
    def length(self):
        """int: length of the sequence."""
        return len(self.sequence)

    @property
    def base_length(self):
        """int: length of the nucleotide sequence."""
        return len(self.base_sequence)

    @property
    def base_sequence(self):
        """Breaks the domain into non-composite domains.

        Raises:
          NotImplementedError: Mixed composite and non-composite domain sequence.

        Returns:
          [list]: non-composite domain sequence.
        """
        if all(isinstance(s, str) for s in self._sequence):
            return self._sequence
        else:
            if all(isinstance(s, DSD_Domain) for s in self._sequence):
                return list(''.join(
                    map(lambda x: ''.join(x.base_sequence), self._sequence)))
            else:
                raise NotImplementedError('''Mixed composite and non-composite domain
            sequence.''')

    def update_constraints(self, con):
        """Apply nucleotide sequence constraint and check if it is compatible.

        Args:
          con (list): A new sequence constraint.

        Raises:
          DSDObjectsError: Must pass constraints as list.
          DSDObjectsError: Constraints have different length.
          DSDObjectsError: Constraints cannot be satisfied.
          NotImplementedError: Cannot update constraints on composite domains.

        """

        if not con or not isinstance(con, list):
            raise DSDObjectsError("Must pass a constraints as a list.")

        # Implement this when needed!
        if not all(isinstance(s, str) for s in self.sequence + con):
            raise NotImplementedError('Cannot update constraints on composite domains.')

        con1 = SequenceConstraint(''.join(self.sequence))
        con2 = SequenceConstraint(''.join(con))

        self._sequence = list(str(con1 + con2))

    def get_ComplementDomain(self, compseq):
        """This function returns a complementary domain.

        Args:
          compseq (list): list of IUPAC nucleic acid sequence constraints.

        Note:
          To simply return the complement, use the '~' operator. This function
          returns a DSD_Domain(), if you inherit from DSD_Domain, you might
          want to overwrite this function to return your object.

        Returns:
          [dsdobjects.DSD_Domain()]: A domain complementary to self
        """
        if not all(isinstance(s, str) for s in self.sequence + compseq):
            raise NotImplementedError('Cannot initialize composite DSD_ComplementDomain.')

        # Make sure sequence constraint is possible
        con1 = SequenceConstraint(''.join(self.sequence[::-1]))
        con2 = SequenceConstraint(''.join(compseq))
        # Invert con1 and merge with con2
        compseq = list(str(~con1 + con2))

        if self.is_complement:
            if self._identity in DSD_Domain.dictionary :
                other = DSD_Domain.dictionary[self._identity]
                other.update_constraints(compseq)
                #self.update_constraints(other.sequence[::-1])
            else :
                other = DSD_Domain(compseq, self._identity, is_complement = False)
        else :
            if self._identity + '*' in DSD_Domain.dictionary :
                other = DSD_Domain.dictionary[self._identity + '*']
                other.update_constraints(compseq)
                #self.update_constraints(other.sequence[::-1])
            else :
                other = DSD_Domain(compseq, self._identity + '*', is_complement = True)

        return other

    @property
    def is_complement(self):
        """ Checks if this domain is a complementary domain.

        Returns:
          [bool]: True if the domain is a complement domain, False othwerwise. """
        return self._name[-1:] == '*'

    def __len__(self):
        return len(self.sequence)

    def __eq__(self, other):
        """
        Domains are equal if they have the same name.
        """
        if isinstance(self, DSD_Domain) and isinstance(other, DSD_Domain):
            return self.name == other.name
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        """Return the name of this domain. """
        return self._name

    def __invert__(self):
        """ Return a complement of this Domain.

        Note:
          This function will not automatically make a DSD_ComplementDomain, as it is
          unclear how constraints should be handled in this case. Use
          get_ComplementDomain() first to initialize the complementary domain.

        Returns:
          [nuskell.objects.DSD_ComplementDomain()]
        """
        complement = self._identity if self.is_complement else self._identity + '*'
        try:
            return DSD_Domain.dictionary[complement]
        except IndexError:
            raise DSDObjectsError('Initialize complementary domain before you of "~" operator.')

class DSD_Complex(object):
    """A sequence and structure pair.

    Sequence and structure can be specified on the domain or on the nucleotide
    level, but they have to be of the same length. Although not yet implemented,
    one may define special characters for secondary structures that are more
    diverse than a regual dot-parens string, e.g. 'h' = '(', '.', ')'

    Global Variables:
      id_counter: the next automatically assigned complex-ID in a system.

    Args:
      sequence (list): A list of sequence constraints on this domain.
      structure (list): A list of dot-parens characters corresponding to the
        sequence.
      name (str, optional): Name of this domain. If not specified, an automatic
        name is generated.
      prefix (str): A prefix for automatic naming of Domains. Defaults to 'cplx'.

    """

    # Automatic name assignments
    id_counter = 0

    # Complex and complex names sanity checks:
    dictionary = dict()
    canonical_forms = dict()

    def __init__(self, sequence, structure, name='', prefix='cplx'):
        # Assign name
        if name:
            self._name = name
        else:
            if prefix == '':
                raise DSDObjectsError('DSD_Complex prefix must not be empty!')
            if prefix[-1].isdigit():
                raise DSDObjectsError('DSD_Complex prefix must not end with a digit!')
            self._name = prefix + str(DSD_Complex.id_counter)
            DSD_Complex.id_counter += 1

        if self._name in DSD_Complex.dictionary:
            raise DSDObjectsError('Duplicate DSD_Complex name!', self._name)

        DSD_Complex.dictionary[self._name] = self

        if not sequence:
            raise DSDObjectsError(
                'DSD_Complex() requires Sequence and Structure Argument')

        if len(sequence) != len(structure):
            raise DSDObjectsError(
                "DSD_Complex() sequence and structure must have same length")

        # should remain in same order as initialized
        self._sequence = sequence
        self._structure = structure
        self._canonical_form = None
        self._rotations = None

        # Initialized on demand:
        self._pair_table = None
        self._strands = None
        self._domains = None

        if self.canonical_form in DSD_Complex.canonical_forms:
            realname = DSD_Complex.canonical_forms[self.canonical_form]
            other = DSD_Complex.dictionary[realname]

            error = DSDObjectsError('Duplicate Complex specification:', self.name, realname)
            error.solution = realname
            error.rotations = self._rotations - other._rotations
            raise error
        else :
            DSD_Complex.canonical_forms[self.canonical_form] = self.name


        # rotate strands into a canonical form.
        # two complexes are equal if they have equal canonial form
        self._exterior_domains = None

    @property
    def name(self):
        """str: name of the complex object. """
        return self._name

    @name.setter
    def name(self, name):
        if name in DSD_Complex.dictionary:
            raise DSDObjectsError('name reserved for different complex')
        if self._name in DSD_Complex.dictionary:
            del DSD_Complex.dictionary[self._name]
        self._name = name
        DSD_Complex.dictionary[self._name] = self

    @property
    def sequence(self):
        """list: sequence the complex object. """
        #TODO: return copy? [:]
        return self._sequence

    @property
    def structure(self):
        """list: the complex structure. """
        #TODO: return copy? [:]
        return self._structure

    @property
    def strands(self):
        if not self._strands:
            self._strands = sorted('_'.join(tuple(map(str,self.sequence))).split('_+_'))
        return self._strands

    @property
    def domains(self):
        if not self._domains:
            self._domains = set(self.sequence)
            if '+' in self._domains: self._domains.remove('+')
            self._domains = sorted(self._domains, key = lambda x: x.name)
        return self._domains

    @property
    def kernel_string(self):
        """str: print sequence and structure in `kernel` notation. """
        seq = self.sequence
        sst = self.structure
        knl = ''
        for i in range(len(seq)):
            if sst[i] == '+':
                knl += str(sst[i]) + ' '
            elif sst[i] == ')':
                knl += str(sst[i]) + ' '
            elif sst[i] == '(':
                knl += str(seq[i]) + str(sst[i]) + ' '
            else:
                knl += str(seq[i]) + ' '
        return knl

    @property
    def canonical_form(self):
        """ Sort sequence & structure in a unique way for each complex.

        Complexes are equal if they have the same canonical form:
            1) sort by sequence string
            2) sort by structure string

        Sets a variable self._rotations.
        """
        def string(s):
            return '_'.join(map(str,s))

        if not self._canonical_form:
            all_variants = dict()
            for e, new in enumerate(self.rotate, 1):
                canon = tuple((string(self.sequence), string(self.structure)))
                if canon not in all_variants:
                    all_variants[canon] = e
            self._canonical_form = sorted(list(all_variants.keys()), key = lambda x:(x[0],x[1]))[0]
            # Invert rotations to be compatible with 'wrap'
            self._rotations = abs(all_variants[self._canonical_form] - len(self.strands))
        return self._canonical_form

    def rotate_pairtable_loc(self, loc, n=None):
        def wrap(x, m):
            """
            Mathematical modulo; wraps x so that 0 <= wrap(x,m) < m. x can be negative.
            """
            return (x % m + m) % m

        if n is None:
            n = len(self.strands)
        return (wrap(loc[0] + n, len(self.strands)), loc[1])

    # Move set
    @property
    def pair_table(self):
        # Make a new pair_table every time, it might get modified.
        return make_pair_table(self.structure)

    def get_domain(self, loc):
        return self.lol_sequence[loc[0]][loc[1]]
    
    def get_paired_loc(self, loc):
        if not self._pair_table:
            self._pair_table = make_pair_table(self.structure)
        return self.pair_table[loc[0]][loc[1]]

    def get_triple(self, *loc):
        return (loc, self.get_domain(loc), self.get_paired_loc(loc))

    @property
    def exterior_domains(self):
        if not self._exterior_domains:
            if not self._pair_table:
                self._pair_table = make_pair_table(self.structure)
            loop_index, exteriors = make_loop_index(self._pair_table)

            self._exterior_domains = []
            for si, strand in enumerate(loop_index):
                for di, domain in enumerate(strand):
                    if loop_index[si][di] in exteriors:
                        if self._pair_table[si][di] is None:
                            self._exterior_domains.append((si, di))
        return self._exterior_domains

    @property
    def lol_sequence(self):
        """ Returns sequence as a list of lists, without the '+' separator.
        
         Example:
          ['d1', 'd2', '+', 'd3'] ==> [['d1', 'd2'], ['d3']]
        """
        return make_lol_sequence(self._sequence)
        #indices = [-1] + [i for i, x in enumerate(self._sequence) if x == "+"]
        #indices.append(len(self._sequence))
        #return [self._sequence[indices[i - 1] + 1 : indices[i]] 
        #        for i in range(1, len(indices))]

    @property
    def nucleotide_sequence(self):
        """list: the complex sequence in form of a flat list of nucleotides. """
        lol = self.lol_sequence

        def my_base_sequence(seq):
            if all(isinstance(d, DSD_Domain) for d in seq):
                return map(lambda x: x.base_sequence, seq)
            else:
                return seq
        return _flatten(map(lambda x: my_base_sequence(x) + ['+'], lol))[:-1]

    @property
    def lol_structure(self):
        indices = [-1] + [i for i, x in enumerate(self._structure) if x == "+"]
        indices.append(len(self._structure))
        return [self._structure[indices[i - 1] + 1: indices[i]]
                for i in range(1, len(indices))]

    @property
    def nucleotide_structure(self):
        """list: the complex structure on nucleotide-level. """
        lol = zip(self.lol_sequence, self.lol_structure)

        def my_base_sequence(seq, sst):
            if all(isinstance(d, DSD_Domain) for d in seq):
                tups = zip(seq, sst)
                return map(lambda x_y: x_y[1] * x_y[0].base_length, tups)
            else:
                return sst
        return _flatten(map(lambda x_y1: my_base_sequence(
            x_y1[0], x_y1[1]) + ['+'], lol))[:-1]

    @property
    def rotate(self):
        """Generator function yields every rotation of the complex. """
        for i in range(len(self.lol_sequence)):
            yield self.rotate_once

    @property
    def rotate_once(self):
        """Rotate the strands within the complex and return the updated object. """
        if "+" in self._sequence:
            p = _find(self._sequence, "+")
            self._sequence = self._sequence[p + 1:] + ["+"] + self._sequence[:p]

            stack = []
            for i in range(p):
                if self._structure[i] == "(":
                    stack.append(i)
                elif self._structure[i] == ")":
                    try :
                        stack.pop()
                    except IndexError:
                        raise DSDObjectsError(
                                "Unbalanced parenthesis in secondary structure.")
            for i in stack:
                self._structure[i] = ")"

            stack = []
            for i in reversed(range(p + 1, len(self._structure))):
                if self._structure[i] == ")":
                    stack.append(i)
                elif self._structure[i] == "(":
                    try :
                        stack.pop()
                    except IndexError:
                        raise DSDObjectsError(
                                "Unbalanced parenthesis in secondary structure.")
            for i in stack:
                self._structure[i] = "("
            self._structure = self._structure[p + 1:] + ["+"] + self._structure[:p]
        self._pair_table = None
        return self

    # Sanity Checks
    @property
    def is_domainlevel_complement(self):
        """
        Determines whether the structure includes pairs only between complementary domains.
        Returns True if all paired domains are complementary, raises an Exception otherwise

        TODO: breaks if used on strings, due to the inverse operator
        """
        if not self._pair_table:
            self._pair_table = make_pair_table(self.structure)
        for si, strand in enumerate(self._pair_table):
            for di, domain in enumerate(strand):
                loc = (si,di)
                cloc = self.pair_table[si][di] 
                if cloc is None or self.get_domain(loc) == ~(self.get_domain(cloc)):
                    continue
                else :
                    raise DSDObjectsError('Domains cannot pair', 
                            self.get_domain(loc), self.get_domain(cloc))
        return True

    @property
    def is_connected(self):
        if not self._pair_table:
            self._pair_table = make_pair_table(self.structure)
        try :
            # make loop index raises disconnected error.
            make_loop_index(self._pair_table)
        except DSDObjectsError, e:
            #raise DSDObjectsError(e)
            return False
        return True

    def split_complex(self, ps, pt):
        return split_complex(ps, pt)

    def split(self):
        """ Split DSD_Complex into disconneted components.
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
                    cplxs.append(DSD_Complex(se, ss))
                except DSDObjectsError, e:
                    cplxs.append(DSD_Complex.dictionary[e.solution])
            return sorted(cplxs)

    # Object properties.
    def __str__(self):
        """str: the name of the complex. """
        return self._name

    def __eq__(self, other):
        """ Test if two complexes are equal.

        They are equal if they have the same sequence and secondary structure.
        """
        if not isinstance(self, DSD_Complex) or not isinstance(other, DSD_Complex):
            return False
        return self.canonical_form == other.canonical_form

    def __ne__(self, other):
        return not self.__eq__(other)


class DSD_Reaction(object):
    """ A reaction pathway.

    Args:
      reactants (list): A list of reactants. Reactants can be strings or
        :obj:`DSD_Complex()` objects.
      products (list): A list of products. Products can be strings or 
        :obj:`DSD_Complex()` objects.
      rtype (str, optional): Reaction type, e.g. bind21, condensed, ..
      rate (flt, optional): Reaction rate.
      name (str, optional): Name of the reaction.
      prefix (str, optional): Prefix for atomatic naming scheme.

    TODO: think about reversible reactions (invert operator, etc)
    """

    dictionary = dict()

    # Used for sanity checks, defaults to DNA branch migration reaction types.
    rtypes = set([None])

    def __init__(self, reactants, products, rtype=None, rate=None):
        self._reactants = reactants
        self._products = products
        self._rtype = rtype
        self._rate = rate
        
        # Used for __eq__ only.
        self._canonical_form = tuple(self.canonical_form)
        if self._canonical_form in DSD_Reaction.dictionary:
            raise DSDObjectsError('Duplicate reaction in system.', self._canonical_form)
        else :
            DSD_Reaction.dictionary[self._canonical_form] = self


        if rtype not in DSD_Reaction.rtypes:
            raise DSDObjectsError('Reaction type not supported! ' + 
            'Set supported reaction types using DSD_Reaction.rtypes')

        # Not Implemented, ...
        #self._normalized_form = tuple((self.normalized, rtype))

    @property
    def name(self):
        """str: name of the reaction. """
        return self._name

    @property
    def rate(self):
        """flt: reaction rate. """
        return self._rate

    @property
    def rateunits(self):
        """str: reaction rate units. """
        return "/M" * (self.arity[0] - 1) + "/s"

    @property
    def rtype(self):
        """str: *peppercorn* reaction type (bind21, condensed, ...) """
        return self._rtype

    @property
    def reactants(self):
        """list: list of reactants. """
        return self._reactants

    @property
    def products(self):
        """list: list of products. """
        return self._products

    @property
    def arity(self):
        """(int, int): number of reactants, number of products."""
        return (len(self._reactants), len(self._products))

    @property
    def kernel_string(self):
        return "{} -> {}".format("  +  ".join(r.kernel_string for r in self.reactants),
                "  +  ".join(p.kernel_string for p in self.products))

    @property
    def sorted(self):
        return "{} -> {}".format(
            " + ".join(sorted(map(str, self.reactants))), 
            " + ".join(sorted(map(str, self.products))))

    @property
    def normalized(self):
        """
        Ensures that complexes appear on only one side of the reaction by
        removing them evenly from both sides until only one side has any.
        """
        reactants = self.reactants[:]
        products = self.products[:]

        for reactant in reactants:
            while (reactant in reactants and reactant in products):
                reactants.remove(reactant)
                products.remove(reactant)
        return "{} -> {}".format(
            " + ".join(sorted(map(str, reactants))), 
            " + ".join(sorted(map(str, products))))

    @property
    def canonical_form(self):
        return self.sorted, self.rtype

    def __str__(self):
        """prints the formal chemical reaction."""
        return "{} -> {}".format(
            " + ".join(map(str, self.reactants)), " + ".join(map(str, self.products)))

    def __eq__(self, other):
        """bool: Checks if DSD_Reaction objects have the same rtype, reactants, and products. """
        if not isinstance(self, DSD_Reaction) or not isinstance(other, DSD_Reaction):
            return False
        return self._canonical_form == other._canonical_form

    def __ne__(self, other):
        return not (self == other)

#### ^^^^ end of edits ^^^^ ####
class TestTube(object):
    """A reaction network of nucleic acid complexes.

    **Description:**
      :obj:`TestTube()` objects are Nuskell's interface to enumerate and simulate
      nucleic acid systems.  Domain-level reaction networks are enumerated using
      the Python package `Peppercorn`_ and simulated using the Python package
      `crnsimulator`_. :obj:`TestTube()` objects to low-level data structures,
      e.g. to verify the equivalence between two :obj:`TestTube()` objects.

      Single or multiple :obj:`DSD_Complex()` and/or :obj:`Reaction()` objects can be
      accessed, added and removed from the system.  :obj:`TestTube()` provides
      (optional) assert statements checking if :obj:`Complex()` and
      :obj:`Domain()` instances have been duplicated, but they might be
      time-consuming for large networks.

      Built-in functions can process domain-level networks to remove strands with
      wildcard-domains (also called history-domains).  Together with the
      remainder of the molecule, a species with a wildcard forms a
      regular-expression, matching every other species in the system that differs
      only by a single domain instead of ‘?’.  If there exists a species matching
      the regular expression, then the species with the wildcard domain and every
      enumerated reaction emerging from that species is be removed from the
      system, otherwise, the wildcard domain is replaced by a regular long
      domain.

    **Structure:**
      :obj:`TestTube()` is based on networkx.MultiDiGraph(), with two types of nodes: (a)
      :obj:`DSD_Complex()` and (b) :obj:`DSD_Reaction()`.  The graph is bipartite, edges
      are directed (from reactants to prodcuts), they connect reactants to a
      reaction node and a reaction node to products.

      :obj:`TestTube()` provides an additional *concentration* attribute and *constant*
      atribute to :obj:`DSD_Complex()` nodes, as well as a rate attribute to
      :obj:`DSD_Reaction()` nodes. These attributes are accessed when writing ODE
      systems and (optionally) updated after simulations. This means a
      :obj:`TestTube()` *can* be initialized without using :obj:`DSD_Complex()` and
      :obj:`DSD_Reaction()` objects, but by using a consistent naming scheme.

    **Developers:**
      It is recommended to store all other node attributes (e.g. free energies,
      etc.) directly in the :obj:`DSD_Complex()` and :obj:`DSD_Reaction()` objects.
      :obj:`TestTube()` does not provide an I/O interface for file formats. There is a
      separate :obj:`TestTubeIO()` object explicitly to parse and write compatible file
      formats (\*.pil, \*.dom, \*.dna, etc.).

    Args:
      complexes (:obj:`dict`, optional) =  A dictionary of complex names that stores
        a tuple of [0]: the respective :obj:`DSD_Complex()` object or :obj:`None`, [1]: the
        concentration, and [2]: boolean value to specify if the concentration
        stays constant during a simulation: True = constant; False = initial
        For example: complexes['A'] = (None, 100, True) is a new
        species called 'A', which has no :obj:`DSD_Complex()` object initialized, and
        remains constant at 100 [nM] concentration.
        Note: The constant attribute defaults to False if not specified, in which
        case the concentrations specify the initial state and might be updated
        after a simulation.

      reactions <dict(), optional> =  A dictionary of reaction names that stores
        either a :obj:`DSD_Reaction()` object or a list [[reactants], [products], rate]
        Reactants and products have to be known complexes.

    Raises:
      DSDObjectsError: 'Wrong initialization of arguments for TestTube()'
      DSDObjectsError: 'Invalid Reaction format'


    .. _Peppercorn:
        http://www.github.com/DNA-and-Natural-Algorithms-Group/peppercorn
    .. _crnsimulator:
        http://www.github.com/bad-ants-fleet/crnsimulator

    """

    # A global TestTube() variable to make (time-consuming) sanity-checks when
    # separate instances are subject to boolean or arithmetic opertions.
    sanitychecks = True
    warnings = True

    def __init__(self, complexes=None, reactions=None):
        self._RG = nx.MultiDiGraph()

        if complexes:
            for name, data in complexes.items():
                if len(data) == 2 or len(data) == 3:
                    assert isinstance(data[0], DSD_Complex) or data[0] is None
                    assert isinstance(data[1], float) or data[1] is None
                    if len(data) == 3:
                        assert isinstance(data[2], bool) or data[2] is None
                else:
                    raise DSDObjectsError(
                        'Wrong initialization of arguments for TestTube()')

                cplx = data[0]
                conc = data[1]
                const = data[2] if len(data) == 3 else None
                if cplx:
                    self._RG.add_node(cplx, concentration=conc, constant=const)
                else:
                    self._RG.add_node(name, concentration=conc, constant=const)

        if reactions:
            # Also reactions have to be uniquely adressable
            for name, react in reactions.items():
                if isinstance(react, DSD_Reaction):
                    self._RG.add_node(react, rate=react.rate)
                    for r in react.reactants:
                        assert self._RG.has_node(r)
                        self._RG.add_edge(r, react)
                    for p in react.products:
                        assert self._RG.has_node(p)
                        self._RG.add_edge(react, p)
                elif isinstance(react, list):
                    assert len(react) == 3
                    self._RG.add_node(name, rate=react[2])
                    for r in react[0]:
                        assert self._RG.has_node(r)
                        self._RG.add_edge(r, name)
                    for p in react[1]:
                        assert self._RG.has_node(p)
                        self._RG.add_edge(name, p)
                else:
                    raise DSDObjectsError('Invalid Reaction format')

        self._domains = None
        self._strands = None

    @property
    def ReactionGraph(self):
        """:obj:`networkx.MultiDiGraph()`: bipartite reaction graph. """
        return self._RG

    @ReactionGraph.setter
    def ReactionGraph(self, RG):
        self._RG = RG
        self._domains = None
        self._strands = None

    @property
    def complexes(self):
        """list: a list of :obj:`DSD_Complex()` objects. """
        # TODO: This only works with Complex() objects
        return [n for n in self._RG.nodes() if isinstance(n, DSD_Complex)]

    def has_complex(self, cplx):
        return self._RG.has_node(cplx)

    @property
    def reactions(self):
        """list: a list of :obj:`DSD_Reaction()` objects. """
        # TODO: This only works with Reaction() objects
        return [n for n in self._RG.nodes() if isinstance(n, DSD_Reaction)]

    def set_complex_concentration(self, cplx, concentration, constant):
        """
        Args:
          cplx (:obj:`DSD_Complex()`): DSD_Complex.
          concentration (flt): concentration.
          constant (bool): True if the concentration is kept constant, False for
            initial concentration.
        """
        self._RG.node[cplx]['concentration'] = concentration
        self._RG.node[cplx]['constant'] = constant

    def get_complex_concentration(self, cplx):
        """flt, bool: First value is concentration, second is True if it is constant,
        False if it is variable."""
        return self._RG.node[cplx]['concentration'], self._RG.node[cplx]['constant']

    def selected_complexes(self, names):
        """list: a list of :obj:`Complex()` objects that correspond to specified names. """
        # TODO: This only works with Complex() objects
        return [n for n in self._RG.nodes() if isinstance(n, DSD_Complex)
                and n.name in names]

    def present_complexes(self, exclude=[], th=0):
        """Returns a list of :obj:`DSD_Complex()` objects with occupancy greater than a threshold.

        Args:
          exclude (list, optional): Exclude particular complex names.
          th (flt, optional): Specify the threshold to consider complexes as
            *present* in solution. Defaults to 0.
        """
        return [n for n, att in self._RG.node.items() if
                isinstance(n, DSD_Complex) and att['concentration'] > th and n.name not in exclude]

    def interpret_species(self, species, prune=True):
        """Get an interpretation dictionary.

        If a :obj:`DSD_Complex()` sequence contains a wildcard, then this function will find
        all matching complexes, and return those as interpretation.  Regex-nodes
        may have at most *one wildcard* per complex, a wildcard corresponds to
        exactly *one unpaired domain*.

        Args:
          species (list[str], optional): A list of complex names that are potential
            regular-expression complexes.
          prune (bool, optinal): True: Remove all wildcards from the network. If a
            matching complex has been found, then the regex-complex and all emerging
            complexes are removed, if no matching complex has been found, then the
            wildcard domain is replaced by a regular long domain. False: Return the
            interpretation without udating the reaction network. Defaults to True.

        Example:
          - It is possible to specify sthg like:
            A = "? a b c" | ". . . ."
            B = "a ? b + c* a*" | "( . . + . )"

          - It is not possible to specify sthg like:
            A = "? a b ?" | "( . . )"
            A = "* a b c" | "* . . ."
            A = "? a ? *" | "( . ) *"
            A = "? a ? x + z* x* f* " | "? ? ? ( + . ) ."
            A = "* a * t" | "* . * ."

        Returns:
          [:obj:`dict()`]: Interpretation of signal species: dict['A_i'] = Counter('A':1)
        """

        def patternMatch(x, y, ignore='?'):
            """Matches two complexes if they are the same, ignoring history domains.

            Note: The strand order of the second complex changes to the strand order of
            the first complex, if there is a rotation under which both complexes are
            patternMatched.

            Args:
              x (DSD_Complex()) : A nuskell DSD_Complex() object.
              y (DSD_Complex()) : A nuskell DSD_Complex() object.

            Returns: True/False
            """
            if len(x.sequence) != len(y.sequence):
                return False

            def pM_check(pMx, pMy):
                """Recursively parse the current sequences and structures.

                Args:
                  pMx [seqX,strX]: A list of two lists (sequence, structrure)
                  pMy [seqY,strY]: A list of two lists (sequence, structrure)

                Returns: True/False
                """
                if len(pMx[0]) == 0:
                    return True

                if (pMx[0][0] != ignore and pMy[0][0] != ignore) and \
                        (pMx[0][0] != pMy[0][0] or pMx[1][0] != pMy[1][0]):
                    return False
                # elif toponly and pMy[0][0] == ignore :
                #  return False
                return pM_check([pMx[0][1:], pMx[1][1:]],
                                [pMy[0][1:], pMy[1][1:]])

            pMx = [map(str, x.sequence), map(str, x.structure)]
            pMy = [map(str, y.sequence), map(str, y.structure)]
            if pM_check(pMx, pMy):
                return True
            elif '+' in map(str, x.sequence) and '+' in map(str, y.sequence):
                for yr in y.rotate:
                    pMy = [map(str, yr.sequence), map(str, yr.structure)]
                    if pM_check(pMx, pMy):
                        return True
            return False

        def get_matching_complexes(regex):
            """Find all matching complexes. """
            regseq = regex.sequence
            regstr = regex.structure
            hist = filter(lambda x: x[0] == 'h', map(str, regseq))
            if len(hist) > 1:
                raise DSDObjectsError("multiple history domains!")
            else:
                hist = hist[0]

            matching = []
            for cplx in self.complexes:
                if regex.name == cplx.name:  # found the regex complex again
                    continue
                else:
                    if patternMatch(regex, cplx, ignore=hist):
                        matching.append(cplx)
            return matching

        need_to_prune = False
        interpretation = dict()
        for fs in species:
            cplxs = [n for n in self._RG.nodes() if n.name == fs]
            if len(cplxs) == 0:
                print '====='
                print 'WARNING: No complex found with name of formal species:', fs
                print '====='
                continue
            else:
                assert len(cplxs) == 1, Warning('Duplicate complex names?')

            cplx = cplxs[0]
            # if '?' in map(str, cplx.sequence) :
            if 'h' in map(lambda d: d.name[0], [
                          d for d in cplx.sequence if d != '+']):
                matches = get_matching_complexes(cplx)
                if matches:
                    need_to_prune = True
                    for e, m in enumerate(matches, 1):
                        m.name = fs + '_' + str(e) + '_'
                        interpretation[m.name] = Counter([fs])
                    self.rm_complex(cplx, force=True)
                else:
                    # NOTE: We cannot simply remove the domain, because we would need to
                    # enumerate the network again and remove the domain everywhere! So
                    # unless we enumerate up-front with a history-pruned species, this
                    # gets us into trouble.
                    #
                    # # Remove Domain
                    # hidx = map(lambda d:d.name[0], cplx.sequence).index('h')
                    # del cplx.sequence[hidx]
                    # del cplx.structure[hidx]
                    # #print cplx.name, map(str, cplx.sequence), map(str, cplx.structure)
                    # cplx.name = fs+'_0_'
                    interpretation[cplx.name] = Counter([fs])
            else:
                interpretation[cplx.name] = Counter([fs])

        if prune and need_to_prune:
            # Get rid of all reactions with history wildcards. Start with a set
            # of produce molecules and see what species emerge from reactions
            # consuming these molecules.
            # Alternative: enumerate again using history-replaced species.
            rxns = self.reactions
            [prev, total] = [set(), set(interpretation.keys() + map(str,
                                                                    self.present_complexes(exclude=map(str, species))))]
            while prev != total:
                prev = set(list(total))
                for rxn in rxns:
                    self.rm_reaction(rxn)
                    r = map(str, rxn.reactants)
                    p = map(str, rxn.products)
                    if set(r).intersection(total) == set(r):
                        total = total.union(set(p))
            map(self.add_reaction, filter(
                lambda x: set(map(str, x.reactants)).intersection(
                    total) == set(map(str, x.reactants)),
                rxns))

            # Now remove all the left-over complexes from the graph.
            all_nodes = set(self.complexes)
            assert set(map(str, all_nodes)).issuperset(total)
            total = set([n for n in self.complexes if str(n) in total])
            remove = all_nodes.difference(total)
            self._RG.remove_nodes_from(remove)

        return interpretation

    def add_complex(self, cplx, xxx_todo_changeme, sanitycheck=True):
        """Add a complex to the TestTube.

        Args:
          cplx (:obj:`DSD_Complex()`): The complex object.
          (conc, const) (flt, bool): Concentration and True/False for constant or
            initial concentrations.
          sanitycheck (bool): True: Check if complex exists under a different name.
            This can be time consuming. Defaults to True.

        Note:
          A new complex resets TestTube.domains and TestTube.strands

        """
        (conc, const) = xxx_todo_changeme
        if not isinstance(cplx, DSD_Complex):
            # TODO: do the formal stuff later
            raise NotImplementedError

        if self._RG.has_node(cplx):
            if conc is not None:
                if self._RG.node[cplx]['concentration'] is None:
                    self._RG.node[cplx]['concentration'] = conc
                else:
                    assert self._RG.node[cplx]['concentration'] == conc, \
                        Warning("Conflicting complex concentrations")
            if const is not None:
                if self._RG.node[cplx]['constant'] is None:
                    self._RG.node[cplx]['constant'] = const
                else:
                    assert self._RG.node[cplx]['constant'] == const, \
                        Warning("Conflicting complex concentrations")
        else:
            # NOTE: This might become inefficient at some point, but it has been
            # introduced to overcome issues with some translation schemes that
            # produce the same fuel strand multiple times.
            # TODO: not only it is inefficient, it doesn't check for strand
            # ordering!
            if sanitycheck and (cplx.sequence, cplx.structure) in map(
                    lambda x: (x.sequence, x.structure), self._RG.nodes()):
                if TestTube.warnings:
                    print 'WARNING: One complex, one name! Skipping complex:', cplx.name, \
                        map(str, cplx.sequence), cplx.structure
            else:
                self._RG.add_node(cplx, concentration=conc, constant=const)
                self._domains = None
                self._strands = None

    def rm_complex(self, cplx, force=False):
        """Remove a Complex from the TestTube.

        Args:
          cplx (:obj:`DSD_Complex()`): The complex object.
          force (bool): True: remove complex and all reactions it is engaged in.
            False: Raise an Error if complex is engaged in a reaction.

        Raises:
          RuntimeError: Cannot remove a complex engaged in reactions.

        Note:
          Removing a complex resets TestTube.domains and TestTube.strands
        """
        if self._RG.has_node(cplx):
            if force:
                for (r, c) in self._RG.in_edges(cplx):
                    assert isinstance(r, DSD_Reaction)
                    self.rm_reaction(r)
                for (c, r) in self._RG.out_edges(cplx):
                    assert isinstance(r, DSD_Reaction)
                    self.rm_reaction(r)
            elif (self._RG.in_edges(cplx) or self._RG.out_edges(cplx)):
                raise DSDObjectsError(
                    "Cannot remove a complex engaged in reactions.")
            self._RG.remove_node(cplx)
            self._domains = None
            self._strands = None

    def add_reaction(self, react, sanitycheck=True):
        """Add a reaction to the TestTube.

        Args:
          react (:obj:`Reaction()`): The *irreversible* reaction to be added.
          sanitycheck (bool): True: Check if reaction exists under a different name.
            This can be time consuming. Defaults to True.
        """

        if not isinstance(react, DSD_Reaction):
            # TODO: do the formal stuff once we need it...
            raise NotImplementedError

        if self._RG.has_node(react):
            assert self._RG.node[react]['rate'] == react.rate
        else:
            # NOTE: This might become inefficient at some point, but there might be
            # cases where reactions are duplicated, so we check if the very same
            # reaction exists as a different node:
            if sanitycheck and filter(lambda x:
                                      (set(x.reactants) == set(react.reactants) and
                                       set(x.products) == set(react.products) and x.name == react.name),
                                      self.reactions):
                print 'WARNING: One reaction, one name! Skipping reaction:', react.kernel
            else:
                self._RG.add_node(react, rate=react.rate)
                for r in react.reactants:
                    assert self._RG.has_node(r)
                    self._RG.add_edge(r, react)
                for p in react.products:
                    assert self._RG.has_node(p)
                    self._RG.add_edge(react, p)

    def rm_reaction(self, react):
        """Remove a reaction from the TestTube.

        Args:
          react (:obj:`DSD_Reaction()`): The reaction object to be removed.
        """
        if self._RG.has_node(react):
            self._RG.remove_node(react)

    @property
    def strands(self):
        """Return a dictionary of strands present in the TestTube.

        A strand is a nucleic-acid molecule connected by a single covalent
        backbone. Strands are named automatically, and their names may change
        whenever a new Complex is added to the TestTube.

        Returns:
          [:obj:`dict()`]: strands[strand_1] = [Domain(X), Domain(Y), Domain(Z)]
        """
        if not self._strands:
            count = 0
            self._strands = dict()
            self._strand_names = dict()
            for cplx in self.complexes:
                for s in cplx.lol_sequence:
                    strand = tuple(map(str, s))
                    if strand not in self._strand_names:
                        name = 'strand_{}'.format(count)
                        count += 1
                        self._strand_names[strand] = name
                        self._strands[name] = s
        return self._strands

    @property
    def domains(self):
        """Return a dictionary of Domain Objects present in the TestTube.

        Returns:
          [:obj:`dict()`]: domains[Domain.name] = Domain
        """
        if not self._domains:
            self._domains = dict()
            for cplx in self.complexes:
                for d in cplx.sequence:
                    if d == '+':
                        continue
                    if d.name in self._domains:
                        assert self._domains[d.name] is d
                    else:
                        self._domains[d.name] = d
        return self._domains

    def enumerate_reactions(
            self, args=None, condensed=True, rename=None, prefix='e'):
        """Enumerate reactions using the *peppercorn* enumerator.
        Args:
          args(:obj:`argparse.ArgumentParser()`, optional): Arguments for *peppercorn*.
          condensed (bool, optional): Udate the reaction graph using *condensed* format.
        """
        from nuskell.enumeration import TestTubePeppercornIO
        TestTubePeppercornIO.condensed = condensed
        interface = TestTubePeppercornIO(testtube=self, enumerator=None,
                                         pargs=args, rename=rename, prefix=prefix)
        interface.enumerate()
        self.ReactionGraph = nx.compose(
            self.ReactionGraph,
            interface.testtube.ReactionGraph)

        # TODO: do we really need to reset domains and strands??
        self._domains = None
        self._strands = None

    def simulate_crn(self, odename, sorted_vars=None, unit='M'):
        raise NotImplementedError('removed crnsimulator dependency')

    def __add__(self, other):
        assert isinstance(other, TestTube)
        combined = TestTube()

        # global TestTube() variable
        if not TestTube.sanitychecks:
            if TestTube.warnings:
                print Warning('TestTube() - sanity checks turned off!')
            combined.ReactionGraph = nx.compose(
                self.ReactionGraph, other.ReactionGraph)

        elif len(other.complexes) > len(self.complexes):
            combined.ReactionGraph.add_nodes_from(
                other.ReactionGraph.nodes(data=True))
            combined.ReactionGraph.add_edges_from(
                other.ReactionGraph.edges(data=True))
            map(lambda c: combined.add_complex(c, self.get_complex_concentration(c),
                                               sanitycheck=True), self.complexes)
            map(lambda r: combined.add_reaction(
                r, sanitycheck=True), self.reactions)
        else:
            combined.ReactionGraph.add_nodes_from(
                self.ReactionGraph.nodes(data=True))
            combined.ReactionGraph.add_edges_from(
                self.ReactionGraph.edges(data=True))
            map(lambda c: combined.add_complex(c, other.get_complex_concentration(c),
                                               sanitycheck=True), other.complexes)
            map(lambda r: combined.add_reaction(
                r, sanitycheck=True), other.reactions)
        return combined

    def __radd__(self, other):
        # Reverse add is used for: sum([Testtube1, Testtube2, ...])
        if other == 0:
            return self
        else:
            return self.__add__(other)

class TestTubeIO(object):
    # Requires parser modules
    """A wrapper class to handle I/O of TestTube objects.

    Args:
      ttube (obj:`TestTube()`): A :obj:`TestTube()` object that should be initialized
        or written to a text format.

    """

    def __init__(self, ttube):
        assert isinstance(ttube, TestTube)
        self._testtube = ttube

    @property
    def testtube(self):
        """:obj:`TestTube()` property."""
        return self._testtube

    def write_pil_kernel(self, pil, unit='M', crn=None, ts=None):
        """Write the contents of :obj:`TestTube()` into a PIL file -- KERNEL notation).

        Args:
          pil (filehandle): A filehandle that the output is written to.
          unit (str, optional): Specify a unit of concentrations (M, mM, uM, nM, pM).
          crn (list[list], optional): a nuskell-style CRN expression
          ts (str, optional): name of the translation scheme

        Example:
          length d1 = 6
          length d2 = 4
          length h3 = 1
          cplx1 = h3 d1( d2( + )) @ initial 10 nM
        """
        pil.write("# File autogenerated by nuskell. ")

        if ts:
            pil.write("\n# - Translation Scheme: {}".format(ts))
        if crn:
            pil.write("\n# - Input CRN: \n")
            for rxn in crn:
                assert len(rxn) == 3
                if len(rxn[2]) == 2:
                    pil.write("#    {} <=> {}\n".format(
                        ' + '.join(rxn[0]), ' + '.join(rxn[1])))
                else:
                    pil.write("#    {} -> {}\n".format(
                        ' + '.join(rxn[0]), ' + '.join(rxn[1])))
        pil.write("#\n\n".format(crn))

        domains = self._testtube.domains

        def adjust_conc(conc, unit):
            units = ['M', 'mM', 'uM', 'nM', 'pM']
            # 0,  3,   6,   9,   12
            assert unit in units
            mult = units.index(unit) * 3
            return conc * (10**mult), unit

        # Print Domains
        pil.write("# Domain Specifications\n")
        for k, v in sorted(domains.items(), key=lambda x: x[1].name):
            if v.name[-1] == '*':
                continue
            pil.write("length {:s} = {:d}\n".format(v.name, len(v)))
            #pil.write("sequence {:s} = {:s}\n".format(v.name, ''.join(v.sequence)))

        pil.write("\n# Complex Specifications\n")

        # Print Complexes
        for cplx in sorted(self._testtube.complexes, key=lambda x: str(x)):
            pil.write("{:s} = ".format(cplx.name))
            seq = cplx.sequence
            sst = cplx.structure
            for i in range(len(seq)):
                if sst[i] == '+':
                    pil.write("{:s} ".format(str(sst[i])))
                elif sst[i] == ')':
                    pil.write("{:s} ".format(str(sst[i])))
                elif sst[i] == '(':
                    pil.write("{:s} ".format(str(seq[i]) + str(sst[i])))
                else:
                    pil.write("{:s} ".format(str(seq[i])))

            conc, const = self._testtube.get_complex_concentration(cplx)
            if const is True:
                pil.write(" @ constant {} {}".format(*adjust_conc(conc, unit)))
            elif const is False:
                pil.write(" @ initial {} {}".format(*adjust_conc(conc, unit)))
            pil.write(" \n")

    def load_pil_kernel(self, pilfile):
        """Parses a file written in PIL - KERNEL notation! """
        ppil = parse_pil_file(pilfile)

        def resolve_loops(loop):
            """ Return a sequence, structure pair from kernel format with parenthesis. """
            sequen = []
            struct = []
            for dom in loop:
                if isinstance(dom, str):
                    sequen.append(dom)
                    if dom == '+':
                        struct.append('+')
                    else:
                        struct.append('.')
                elif isinstance(dom, list):
                    struct[-1] = '('
                    old = sequen[-1]
                    se, ss = resolve_loops(dom)
                    sequen.extend(se)
                    struct.extend(ss)
                    sequen.append(old + '*' if old[-1] != '*' else old[:-1])
                    struct.append(')')
            return sequen, struct

        domains = []
        for line in ppil:
            if line[0] == 'domain':
                domains.append(
                    DSD_Domain(name=line[1], sequence=list('N' * int(line[2]))))
            elif line[0] == 'complex':
                name = line[1]
                sequence, structure = resolve_loops(line[2])
                constant, concentration = None, float('inf')
                if len(line) > 3:
                    i, c, u = line[3]
                    constant = (i == 'constant')
                    if u == 'M':
                        concentration = float(c)
                    elif u == 'mM':
                        concentration = float(c) * 1e-3
                    elif u == 'uM':
                        concentration = float(c) * 1e-6
                    elif u == 'nM':
                        concentration = float(c) * 1e-9
                    elif u == 'pM':
                        concentration = float(c) * 1e-12
                    else:
                        raise DSDObjectsError(
                            'unknown unit for concentrations specified.')

                for e in range(len(sequence)):
                    d = sequence[e]
                    if d == '+':
                        continue
                    if d[-1] == '*':
                        dname = d[:-1]
                        dom = filter(lambda x: x.name == dname, domains)
                        if len(dom) < 1:
                            raise DSDObjectsError(
                                'Missing domain specification', d)
                        elif len(dom) > 1:
                            raise DSDObjectsError(
                                'Conflicting matches for domain specification', d)
                        sequence[e] = dom[0].get_ComplementDomain(
                            list('R' * dom[0].length))

                    else:
                        dname = d
                        dom = filter(lambda x: x.name == dname, domains)
                        if len(dom) < 1:
                            raise DSDObjectsError(
                                'Missing domain specification', d)
                        elif len(dom) > 1:
                            raise DSDObjectsError(
                                'Conflicting matches for domain specification', d)
                        sequence[e] = dom[0]

                self._testtube.add_complex(
                    DSD_Complex(
                        sequence=sequence,
                        structure=structure,
                        name=name),
                    (concentration,
                     constant))
            else:
                raise NotImplementedError(
                    'Weird expression returned from pil_parser!')

        return self._testtube

    def write_dnafile(self, fh, signals=[], crn=None, ts=None):
        """ Write a TestTube Object into VisualDSD \*.dna format.

        Note:
          This function assumes that toehold domains are named starting with a 't',
          history domains start with a 'h' and anti-sense domains end with '*'.

        Args:
          fh (filehandle): The function prints to this filehandle.
          signals (list[str], optional): A list of signal species.
          crn (list[list], optional): a nuskell-style CRN expression
          ts (str, optional): name of the translation scheme
        """
        fh.write("(* File autogenerated by nuskell. ")

        if ts:
            fh.write("\n - Translation Scheme: {}".format(ts))
        if crn:
            fh.write("\n - Input CRN: \n")
            for rxn in crn:
                assert len(rxn) == 3
                if len(rxn[2]) == 2:
                    fh.write("    {} <=> {}\n".format(
                        ' + '.join(rxn[0]), ' + '.join(rxn[1])))
                else:
                    fh.write("    {} -> {}\n".format(
                        ' + '.join(rxn[0]), ' + '.join(rxn[1])))

        fh.write("*)\n\n".format(crn))

        fh.write("def Fuel = 20\n")
        fh.write("def Signal = 5\n\n")

        first = True
        for cplx in sorted(self._testtube.complexes, key=lambda x: x.name):

            if first:
                fh.write('( ')
                first = False
            else:
                fh.write('| ')

            if cplx.name in signals:
                fh.write("Signal * ")
            else:
                fh.write("constant Fuel * ")

            name = cplx.name
            sequ = cplx.sequence
            stru = cplx.structure

            ptab = pair_table(stru)

            dnaexpr = [[]]
            pos = 0
            for e, d in enumerate(ptab):
                if d == '+':
                    flag = 'top' if flag == 'bound' else flag
                    expr = 'cut'

                elif d == -1:
                    toe = '^' if sequ[e].name[0] == 't' else ''
                    if sequ[e].name[-1] == '*':
                        flag = 'bottom'
                        expr = sequ[e].name[:-1] + toe + '*'
                    elif sequ[e].name[0] == 'h':
                        flag = 'top'
                        expr = '_'
                    else:
                        flag = 'top'
                        expr = sequ[e].name + toe

                elif d > e:  # '('
                    flag = 'bound'
                    toe = '^' if sequ[e].name[0] == 't' else ''
                    if sequ[e].name[-1] == '*':
                        expr = sequ[e].name[:-1] + toe + '*'
                    elif sequ[e].name[0] == 'h':
                        raise DSDObjectsError(
                            'Unexpected bound history domain.')
                    else:
                        expr = sequ[e].name + toe

                    dnaexpr.append([])
                    pos += 1

                elif d < e:  # ')'
                    flag = 'bottom'
                    expr = None
                    pos -= 1
                    if pos < 0:
                        raise DSDObjectsError('too many closing base-pairs')
                    continue
                else:
                    raise DSDObjectsError('strange case:', e, d)

                if dnaexpr[pos] == []:
                    dnaexpr[pos] = [[flag, expr]]
                else:
                    dnaexpr[pos].append([flag, expr])

            # decode dnaexpr
            dnaflat = []
            for d in dnaexpr:
                for dd in d:
                    dnaflat.append(dd)

            # PRINT TO FILE
            close = None
            for e, d in enumerate(dnaflat):
                if d[1] == 'cut':
                    fh.write(close)
                    close = None
                    if e == len(dnaflat) - 1:
                        continue
                    if d[0] == 'bottom':
                        fh.write('::')
                    else:
                        fh.write(':')
                    continue

                if d[0] == 'bottom':
                    if close is None:
                        fh.write('{')
                        close = '}'
                    elif close == ']' or close == '>':
                        fh.write('{}{'.format(close))
                        close = '}'

                if d[0] == 'bound':
                    if close is None:
                        fh.write('[')
                        close = ']'
                    elif close == '}' or close == '>':
                        fh.write('{}['.format(close))
                        close = ']'

                if d[0] == 'top':
                    if close is None:
                        fh.write('<')
                        close = '>'
                    elif close == '}' or close == ']':
                        fh.write('{}<'.format(close))
                        close = '>'
                fh.write(" {} ".format(d[1]))
            if close:
                fh.write("{} (* {} *)\n".format(close, name))
            else:
                fh.write(" (* {} *)\n".format(name))

        fh.write(")\n")

def pair_table(ss, chars=['.']):
    print DeprecationWarning('old format')
    """Return a secondary struture in form of pair table:

    Args:
      ss (str): secondary structure in dot-bracket format
      chars (list, optional): a list of characters that are ignored. Defaults to
        ['.']

    Example:
       ((..)). => [5,4,-1,-1,1,0,-1]

    Raises:
       DSDObjectsError: Too many closing brackets in secondary structure.
       DSDObjectsError: Too many opening brackets in secondary structure.
       DSDObjectsError: Unexpected character in sequence: "{}"

    Returns:
      [list]: A pair-table
    """
    stack = []

    pt = [-1] * len(ss)

    for i, char in enumerate(ss):
        if (char == '('):
            stack.append(i)
        elif (char == ')'):
            try:
                j = stack.pop()
            except IndexError as e:
                raise DSDObjectsError(
                    "Too many closing brackets in secondary structure")
            pt[i] = j
            pt[j] = i
        elif (char == '+'):
            pt[i] = '+'
        elif (char not in set(chars)):
            raise DSDObjectsError(
                "Unexpected character in sequence: '" + char + "'")

    if stack != []:
        raise DSDObjectsError(
            "Too many opening brackets in secondary structure")
    return pt

