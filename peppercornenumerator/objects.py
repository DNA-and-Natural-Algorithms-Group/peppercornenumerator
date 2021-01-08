#
#  peppercornenumerator/objects.py
#  EnumeratorProject
#
import gc
import logging
log = logging.getLogger(__name__)

from itertools import chain

from dsdobjects import SingletonError, clear_singletons, show_singletons
from dsdobjects.base_classes import ObjectInitError, DomainS, StrandS, ComplexS, MacrostateS, ReactionS

def clear_memory():
    # This is for unittests only!  The memory of a Singleton clears
    # automatically if there is no hard reference to the object.
    gc.collect()
    if list(show_memory()):
        log.warning('Could not collect singleton objects, trying other means.')
    clear_singletons(PepperDomain)
    clear_singletons(PepperStrand)
    clear_singletons(PepperComplex)
    clear_singletons(PepperMacrostate)
    clear_singletons(PepperReaction)

def show_memory():
    # This is for unittests only!  The memory of a Singleton clears
    # automatically if there is no hard reference to the object.
    for x in chain(show_singletons(PepperDomain),
                   show_singletons(PepperStrand),
                   show_singletons(PepperComplex),
                   show_singletons(PepperMacrostate),
                   show_singletons(PepperReaction)):
        yield x

class PepperDomain(DomainS):
    pass

class PepperStrand(StrandS):
    pass

class PepperComplex(ComplexS):
    PREFIX = 'e'
    @property
    def available_domains(self):
        ad = []
        for x, y in self.exterior_domains:
            ad.append((self.get_domain((x, y)), x, y))
        return ad

class PepperMacrostate(MacrostateS):
    pass

class PepperReaction(ReactionS):
    RTYPES = set(['condensed', 'open', 'bind11', 'bind21', 'branch-3way', 'branch-4way'])

    @property
    def const(self):
        log.warning('deprecated property PepperReaction.const')
        return self.rate_constant[0]

    @const.setter
    def const(self, value):
        log.warning('deprecated property PepperReaction.const')
        self.rate_constant = value

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
    def domain_locs(self):
        return (part[2] if part is not None else None for part in self._parts)

    @property
    def pair_locs(self):
        return (part[1] if part is not None else None for part in self._parts)

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

