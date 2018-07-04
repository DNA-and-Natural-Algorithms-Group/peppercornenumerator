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

SHORT_DOMAIN_LENGTH = 6
LONG_DOMAIN_LENGTH = 12

class PeppercornUsageError(Exception):
    """Error class to catch usage errors."""

    def __init__(self, msg, val=None):
        self.message = msg
        if val :
            self.message += " ({})".format(val)
        super(PeppercornUsageError, self).__init__(self.message) 


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
    def convert(text): 
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): 
        return [convert(c) for c in re.split('([0-9]+)', str(key))]

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
    def domains(self):
        return (part[0] if part is not None else None for part in self._parts)

    @property
    def structures(self):
        return (part[1] if part is not None else None for part in self._parts)

    @property
    def structs(self):
        raise DeprecationWarning('Use of Loop.structs has been replaced by Loop.structures.')
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

