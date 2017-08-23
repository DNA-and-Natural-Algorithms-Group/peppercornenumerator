#
#  utils.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/18/10.
#  Modifications by Casey Grun and Erik Winfree 8/15/2014.#

import re
import sys
import copy
import logging
from math import log10
from dsdobjects import DSD_Domain, DSD_Complex, SequenceConstraint, DSDObjectsError

SHORT_DOMAIN_LENGTH = 6
LONG_DOMAIN_LENGTH = 12

class PeppercornUsageError(Exception):
    """Error class to catch usage errors."""

    def __init__(self, msg, val=None):
        self.message = msg
        if val :
            self.message += " ({})".format(val)
        super(PeppercornUsageError, self).__init__(self.message) 

class colors:
    RED = '\033[91m'
    YELLOW = '\033[93m'
    GREEN = '\033[92m'
    BLUE = '\033[94m'
    PINK = '\033[95m'
    CYAN = '\033[96m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    colors = [RED, YELLOW, GREEN, CYAN, BLUE, PINK]

    @staticmethod
    def color(string):
        pass

    @staticmethod
    def legend(keys=None):
        if keys is None:
            l = enumerate(colors.colors)
        else:
            l = zip(keys, colors.colors)
        return "\n".join([(c + str(i) + colors.ENDC) for i, c in l])


def wrap(x, m):
    """
    Mathematical modulo; wraps x so that 0 <= wrap(x,m) < m. x can be negative.
    """
    return (x % m + m) % m


def natural_sort(l):
    """
    Sorts a collection in the order humans would expect. Implementation from
    http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    """
    def convert(text): return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): return [convert(c)
                                   for c in re.split('([0-9]+)', str(key))]
    return sorted(l, key=alphanum_key)


def find(f, seq, default=None):
    """
    Return first item in sequence where f(item) == True.
    """
    for item in seq:
        if f(item):
            return item
    return default


def warning(message):
    # from termcolor import colored, cprint
    # cprint("Warning: " + message, 'yellow')
    logging.warning(message)


def error(message):
    logging.error(message)
    sys.exit(1)


def wait_for_input(message="[Press Enter to continue...]"):
    raw_input(message)
    print ""


def parse_parameters(parameters):
    match = re.match(r"\s*\[([^\]]+)\]\s*", parameters)
    if match is not None:
        parameters, = match.groups()

    params = {'concentration': None}

    # parse parameters into <list of targets> @ <list of conditions>
    parameters = parameters.split("@")
    if len(parameters) > 1:
        targets, conditions = parameters

        # split conditions into comma-separated list, discard whitespace
        conditions = conditions.split(",")
        for condition in conditions:
            condition = condition.strip()

            # try to parse a concentration
            concentration = parse_concentration(condition)
            if concentration is not None:
                params['concentration'] = concentration
    return params


exp_to_si_prefix = {9: 'G', 6: 'M', 3: 'k', 0: '', -3: 'm', -
                    6: 'u', -9: 'n', -12: 'p', -15: 'f', -18: 'a', -21: 'z'}

si_prefix_to_exp = dict((pre, 10**exp)
                        for (exp, pre) in exp_to_si_prefix.iteritems())


def parse_concentration(condition):
    parts = re.match(
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(p|n|u|m|d|)M",
        condition)

    if parts is not None:
        conc, unit = parts.groups()
        base = si_prefix_to_exp[unit]
        concentration = float(conc) * base
        return concentration
    return None


def format_si(n):
    try:
        x = int(log10(n) // 3) * 3
        u = exp_to_si_prefix[x]
    except (ValueError, KeyError):
        x = 0
    return (n / 10**x), exp_to_si_prefix[x]


def resolve_length(length):
    if (isinstance(length, type(0))):
        return length
    elif (length == "short"):
        return SHORT_DOMAIN_LENGTH
    elif (length == "long"):
        return LONG_DOMAIN_LENGTH


def parse_basewise_dot_paren(structure_line, strands):
    parts = [x.strip() for x in structure_line.split("+")]
    assert len(parts) == len(strands), "Structure '%s' has %d parts, but corresponds to %d strands" % (
        structure_line, len(parts), len(strands))

    segment_struct = ["" for s in strands]
    for (i, s) in enumerate(strands):
        strand_part = parts[i]
        for d in s.domains:
            assert len(strand_part) >= len(
                d), "Not enough characters for domain %s" % str(d)
            domain_part, strand_part = strand_part[:len(
                d)], strand_part[len(d):]

            assert all(
                c == domain_part[0] for c in domain_part), "Not all parts of structure for %s are the same" % str(d)
            segment_struct[i] += domain_part[0]

    return parse_dot_paren("+".join(segment_struct))


def parse_dot_paren(structure_line):
    """
    Parses a dot-parenthesis structure into the list of lists representeation
    used elsewhere in the enumerator.

    Example::

                             0,0  0,1  0,2   0,3   0,4     1,0   1,1
             "...((+))" -> [[None,None,None,(1,0),(1,1)],[(0,4),(0,3)]]

    """
    complex_structure = []
    dot_paren_stack = []
    strand_index = 0
    domain_index = 0
    curr_strand = []
    complex_structure.append(curr_strand)
    for part in structure_line:
        # stand break
        if (part == "+"):
            strand_index += 1
            domain_index = 0
            curr_strand = []
            complex_structure.append(curr_strand)
            continue

        # unpaired
        if (part == "."):
            curr_strand.append(None)
            domain_index += 1

        # paired to later domain
        elif (part == "("):
            curr_strand.append(None)
            dot_paren_stack.append((strand_index, domain_index))
            domain_index += 1

        # paired to earlier domain
        elif (part == ")"):
            loc = dot_paren_stack.pop()
            curr_strand.append(loc)
            complex_structure[loc[0]][loc[1]] = (strand_index, domain_index)
            domain_index += 1
    return complex_structure


def index_parts(enum):
    """
    Testing tool. Accepts an enumerator, produces a tuple
    (domains, strands, complexes) where each element is a dict mapping names
    of those objects to the objects in the enumerator. For instance, domains
    maps the name of each domain in the enumerator to the Domain object
    """
    domains = {}
    strands = {}
    complexes = {}

    for domain in enum.domains:
        domains[domain.name] = domain

    for strand in enum.strands:
        strands[strand.name] = strand

    for complex in enum.initial_complexes:
        complexes[complex.name] = complex

    return (domains, strands, complexes)


class Loop(object):
    """
    Represents (possibly a part) of a single open or closed loop

    Args: 
        loop [(Domain, bound_to, position),..]: Is a list of descriptors for
            involved domains.  In particular: [0] A Domain Object, [1] its
            binding partner position (or None), [2] the domain position.
            None marks the 3' or 5' end of a strand.
    
    TODO: Loop == Loop is only true if the two loops were initialized with the
    same domain ordering.

    """

    def __init__(self, loop):
        self._parts = loop

        is_open = False
        bases = 0
        stems = 0

        # Sanity check to alert that a stem has been specified with both
        # complementary domains.
        stem_list = set()

        # calculate stems, bases, and is_open   #loop re-written by EW
        for step in loop:
            if step is None:
                is_open = True
            else:
                (dom, struct, loc) = step
                if struct is None:
                    bases += len(dom)
                elif struct in stem_list:
                    raise PeppercornUsageError('Double stem count in Loop() Object.')
                else :
                    stems += 1
                    stem_list.add(loc)

        # update cached properties
        self._is_open = is_open
        self._bases = bases
        self._stems = stems

    @property
    def locs(self):
        return (part[2] if part is not None else None for part in self._parts)

    @property
    def revlocs(self):
        return (part[2] if part is not None else None for part in self._parts[::-1])

    @property
    def domains(self):
        return (part[0] if part is not None else None for part in self._parts)

    @property
    def structures(self):
        return (part[1] if part is not None else None for part in self._parts)

    @property
    def structs(self):
        raise DeprecationWarning('SB: Use of Loop.structs has been replaced by Loop.structures.')
        return self.structures

    @property
    def parts(self):
        """
        Gives the list of (dom, struct, loc) tuples associated with this loop
        """
        return self._parts[:]

    @property
    def stems(self):
        """
        Gives the number of stems in the loop
        """
        return self._stems

    @property
    def bases(self):
        """
        Gives the number of bases in the loop
        """
        return self._bases

    @property
    def is_open(self):
        """
        True if the loop is an open loop, else false
        """
        return self._is_open

    def __contains__(self, item):
        if isinstance(item, Domain):
            return item in self.domains
        elif isinstance(item, tuple):
            return item in self.locs

    def __str__(self):
        return ' '.join([(str(s) if s is not None else '+') for s in self.domains])

    def __repr__(self):
        return "Loop(%s)" % list(self.domains)

    def __len__(self):
        return sum(len(dom) for dom in self.domains)

    def __eq__(self, other):
        #TODO print Warning('SB: not obvious!!')
        return self._parts == other._parts

#class Domain(object):
#    """
#    Represents a single domain. We allow several options for specifying domain
#    properties. Domains might have an explicit integer (bp) length, or may be
#    designated as short or long. If the latter method is used, the code will use
#    the relevant constant as the integer domain length.
#    """
#
#    def __init__(self, name, length, is_complement=False, sequence=None):
#        """
#        Default constructor. Takes a domain name, a length (positive integer or
#        "short" or "long"), and optionally a base sequence.
#        """
#        assert isinstance(length, int) or length in ['long', 'short']
#        self._name = name
#        self._length = length
#        self._sequence = None
#        self._complement_sequence = None
#        self._is_complement = is_complement
#        if (sequence):
#            self._sequence = sequence.upper()
#            self._complement_sequence = self._sequence
#            self._complement_sequence = list(self._complement_sequence[::-1])
#            for i, l in enumerate(self._complement_sequence):
#                if (l == 'A'):
#                    self._complement_sequence[i] = 'T'
#                elif (l == 'T'):
#                    self._complement_sequence[i] = 'A'
#                elif (l == 'C'):
#                    self._complement_sequence[i] = 'G'
#                elif (l == 'G'):
#                    self._complement_sequence[i] = 'C'
#            self._complement_sequence = ''.join(self._complement_sequence)
#
#    def __repr__(self):
#        return "Domain(%s)" % (self.name)
#
#    def __str__(self):
#        return self.name
#
#    def __cmp__(self, other):
#        out = cmp(self.name, other.name)
#        if (out != 0):
#            return out
#
#        out = self.length - other.length
#        if (out != 0):
#            return out
#
#        if (self.is_complement == other.is_complement):
#            return 0
#        elif self.is_complement:
#            return 1
#        else:
#            return 0
#
#    def __len__(self):
#        return self.length
#
#    def __hash__(self):
#        return hash((self.name, self.length, self.is_complement))
#
#    def can_pair(self, other):
#        """
#        Returns True if this domain is complementary to the argument.
#        """
#        return ((self.identity == other.identity) and
#                (self.is_complement or other.is_complement) and
#                (not (self.is_complement and other.is_complement)))
#
#    @property
#    def identity(self):
#        """
#        Returns the identity of this domain, which is its name without a
#        complement specifier (i.e. A and A* both have identity A).
#        """
#        return self._name
#
#    @property
#    def name(self):
#        """
#        The name of this domain.
#        """
#        return (self._name + ("" if not self.is_complement else "*"))
#
#    @property
#    def length(self):
#        """
#        The length of this domain. Either uses the integer length previously
#        specified or the constant associated with "short" and "long" domains as
#        appropriate.
#        """
#        return resolve_length(self._length)
#
#    @property
#    def sequence(self):
#        """
#        This domain's sequence. May be None if no sequence was specified.
#        """
#        if (not self.is_complement):
#            return self._sequence
#        else:
#            return self._complement_sequence
#
#    @property
#    def is_complement(self):
#        """
#        Returns true if this domain is a complement (e.g. A* rather than A),
#        false otherwise.
#        """
#        return self._is_complement

class PepperDomain(DSD_Domain):
    """
    Represents a single domain. We allow several options for specifying domain
    properties. Domains might have an explicit integer (bp) length, or may be
    designated as short or long. If the latter method is used, the code will use
    the relevant constant as the integer domain length.
    """
    def __init__(self, *kargs, **kwargs):
        super(PepperDomain, self).__init__(*kargs, **kwargs)

    def __repr__(self):
        return "%s" % (self.name)

    def __cmp__(self, other):
        out = cmp(self.name, other.name)
        if (out != 0):
            return out

        out = self.length - other.length
        if (out != 0):
            return out

        if (self.is_complement == other.is_complement):
            return 0
        elif self.is_complement:
            return 1
        else:
            return 0

    def __hash__(self):
        return hash((self.name, self.length, self.is_complement))

    # Overwrite original
    def get_ComplementDomain(self, compseq):
        """This function returns a complementary domain.

        Args:
          compseq (list): list of IUPAC nucleic acid sequence constraints.

        Note:
          To simply return the complement, use the '~' operator.

        Returns:
          [PepperDomain()]
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
                other = PepperDomain(compseq, self._identity, is_complement = False)
        else :
            if self._identity + '*' in DSD_Domain.dictionary :
                other = DSD_Domain.dictionary[self._identity + '*']
                other.update_constraints(compseq)
                #self.update_constraints(other.sequence[::-1])
            else :
                other = PepperDomain(compseq, self._identity + '*', is_complement = True)

        return other

    def can_pair(self, other):
        """
        Returns True if this domain is complementary to the argument.
        """
        try :
            return self == ~other
        except KeyError, e:
            return ((not self.is_complement == other.is_complement) and 
                    (self.identity == other.identity))
            #raise DSDObjectsError('Complementary domain has not been initialized.')

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


class Strand(object):
    """
    Represents a strand---an ordered sequence of domains.
    """

    def __init__(self, name, domains):
        """
        Constructor for Strand objects. Takes a strand name and the ordered list
        of Domain objects from 5' to 3'.

        :param string name:
        :param list domains: List of :py:class:Domain objects
        """
        self._name = name
        self._domains = domains
        self._hash = None  # assigned lazily by __hash__

    def __eq__(self, other):
        """
        Strands are compared on the basis of their name and their :py:meth:`domains`
        """
        return (self.name == other.name) and (self.domains == other.domains)

    def __hash__(self):
        if(self._hash is None):
            self._hash = hash((self.name, tuple(self.domains)))

        return self._hash

    def __cmp__(self, other):
        """
        Strands are compared on the basis of their names
        """
        return cmp(self.name, other.name)

    def __len__(self):
        """
        The length of a strand is equal to the number of :py:meth:`domains` in the strand
        """
        return len(self.domains)

    def __repr__(self):
        return "Strand(%s)" % (self.name)

    def __str__(self):
        return self.name

    @property
    def name(self):
        """
        Gives the name of the strand
        """
        return self._name

    @property
    def length(self):
        """
        Returns the number of domains in this strand.
        """
        return len(self._domains)

    @property
    def domains(self):
        """
        Gives the list of domains of a strand
        """
        return self._domains

class PepperComplex(DSD_Complex):
    def __init__(self, *kargs, **kwargs):
        super(PepperComplex, self).__init__(*kargs, **kwargs)

    def __hash__(self):
        return hash(self.canonical_form)

    #def __cmp__(self):
    #    raise NotImplementedError

    def __repr__(self):
        return "%s" % (self.name)

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
        """ Split DSD_Complex into disconneted components.
        """
        if self.is_connected:
            return [self]
        else :
            ps = self.lol_sequence
            pt = self.pair_table
            parts = self.split_complex(ps, pt)
            cplxs = []
            # assign new_complexes
            for (se,ss) in parts:
                try:
                    cplxs.append(PepperComplex(se, ss))
                except DSDObjectsError, e:
                    cplxs.append(PepperComplex.dictionary[e.solution])
            return sorted(cplxs)


class Complex(object):
    pass

class RestingState(object):
    """
    Represents a resting state, which is a collection of complexes.
    """

    def __init__(self, name, complexes):
        """
        Constructor for RestingState objects. Takes a (unique) name and list of
        complexes.
        """
        complexes.sort()
        self._complexes = complexes
        self._name = name
        self._canonical = find(lambda s: not str(
            s).isdigit(), sorted(complexes), complexes[0])
        self._hash = None

    @property
    def name(self):
        """
        Gives the name of the resting state
        """
        return self._name

    @property
    def complexes(self):
        """
        Gives a list of complexes in the resting state
        """
        return self._complexes[:]

    @property
    def canonical_name(self):
        """
        Gives the canonical name of the resting state, chosen by the lexicographically lowest
        name of a complex in the resting state.
        """
        return str(self._canonical)

    @property
    def canonical(self):
        """
        See ``canonical_name``.
        """
        return self._canonical

    def kernel_string(self):
        return self.canonical.kernel_string()

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(tuple(self.complexes))
        return self._hash

    def __eq__(self, other):
        """
        Two resting states are equal if their complexes are equal
        """
        return (self.complexes == other.complexes)

    def __cmp__(self, other):
        """
        Two resting states are compared on the basis of their complexes
        """
        return cmp(self.complexes, other.complexes)

    def __str__(self):
        return self.canonical_name

    def __repr__(self):
        return "RestingState(\"%s\", %s)" % (self.name, str(self.complexes))
