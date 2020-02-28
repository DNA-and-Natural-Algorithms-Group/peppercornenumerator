#
#  peppercornenumerator/utils.py
#  EnumeratorProject
#
from __future__ import absolute_import, print_function, division

import logging
log = logging.getLogger(__name__)

import re
import sys

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
        return "Loop({:s})".format(list(self.domains))

    def __len__(self):
        return sum(len(dom) for dom in self.domains)

    def __eq__(self, other):
        return self._parts == other._parts

