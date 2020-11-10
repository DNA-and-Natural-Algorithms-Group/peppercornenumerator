#
#  peppercornenumerator/enumerator.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

from .utils import PeppercornUsageError, tarjans
from .input import read_pil, read_seesaw
from .output import write_pil, write_sbml
from .objects import (SingletonError, PepperMacrostate, clear_memory)
from .condense import PepperCondensation
from .reactions import (bind11, bind21, open1N, 
                        branch_3way, branch_4way,
                        find_on_loop, filter_bind11)
from .ratemodel import opening_rate

"""
WASTEFUL (bool): A global parameter to erase singleton object memory on the
beginning of every call of "enumerate_pil" or "enumerate_ssw". This should not
be necessary if you make sure that all "old" peppercornenumerator.objects are
deleted.  Strangely, some dependencies (in particluar matplotlib.pyplot)
prevent the garbage collection of "old" objects, and it is necessary to use
this parameter once such a depencency has been imported. I suspect wasteful
behavior to be the cause of occasional segfaults when enumerating *many*
different systems in order to do some meta analysis. More details in 
tests/test_wasteful.py
"""
WASTEFUL = False

UNI_REACTIONS = [bind11, open1N, branch_3way, branch_4way]
BI_REACTIONS = [bind21]

class PolymerizationError(Exception):
    pass

class Enumerator:
    """
    The main enumerator object, collecting all the parameters to enumerate the
    Chemical Reaction Network (CRN) from an initial set of domain-level
    complexes.  Enumerators have immutable starting conditions. 

    Args:
        initial_complexes(list[PepperComplex()]): The list of complexes
            presenting a starting condition for the enumeration process.
        initial_reactions(list[PepperReaction()], optional): An optional list
            of Reactions to be added to the enumeration results. These 
            reactions will be considered when finding SCCs and when condensing
            the CRN after enumeration.
    """
    def __init__(self, initial_complexes, initial_reactions = None, 
                 named_complexes = None):
        # Sanity checks for input.
        for cplx in initial_complexes:
            if not cplx.is_connected:
                raise PeppercornUsageError(f'Initial complex is not connected: {cplx}')
        # A set of initially present complexes.
        self._initial_complexes = set(initial_complexes)
        # A set of all detailed reactions before and after enumeration.
        self._initial_reactions = set(initial_reactions) if initial_reactions else set()
        # A list of "known" complexes: 
        #   - expected to show up during enumeration.
        #   - used to prioritize macrostate representatives.
        #   - preserving names for output (could be just strand names!).
        self._named_complexes = set(named_complexes) if named_complexes else set()

        self.representatives = sorted(self._initial_complexes | self._named_complexes)
        
        # Enumeration results
        self._reactions = set(initial_reactions) if initial_reactions else set()
        self._complexes = None
        self._resting_complexes = None
        self._transient_complexes = None
        self._resting_macrostates = None

        # Condensation
        self.condensation = None
        self._detailed_reactions = None
        self._condensed_reactions = None

        # Polymerization settings to prevent infinite looping
        self._max_complex_count = max(200, len(self._initial_complexes))
        self._max_reaction_count = max(1000, len(self._reactions))

        # Settings for reaction enumeration. 
        self.max_helix = True
        self.reject_remote = False
        self._release_11 = 7
        self._release_12 = 7
        self._max_complex_size = 6

        # Rate-dependent enumeration settings:
        #
        # Set separation of timescales for *unimolecular* reactions.
        #
        #  ignore-reaction | slow-reaction | fast-reaction
        # -----------------|---------------|---------------> rate
        #                 k_slow          k_fast
        #
        # Default: OFF = all unimolecular reactions are fast or negligible
        self._k_slow = 0 
        self._k_fast = 0

        # Rate parameter adjustments
        self.dG_bp = -1.7 # adjust bimolecular binding strength

        # Debugging features
        self.DFS = True
        self.interactive = False
        self.interruptible = False

    @property
    def initial_complexes(self):
        """ Complexes present in the system's initial configuration. """
        return set(self._initial_complexes)

    @initial_complexes.setter
    def initial_complexes(self, value):
        raise PeppercornUsageError('Initial complexes must be provided at initialization!')

    @property
    def initial_reactions(self):
        """ Reactions to be included in the system's enumeration. """
        return set(self._reactions)

    @initial_reactions.setter
    def initial_reactions(self, value):
        # NOTE: there is no need to enforce this.
        raise PeppercornUsageError('Initial reactions must be provided at initialization!')

    @property
    def named_complexes(self):
        """ Known complexes which do not necessarily take place in enumeration. """
        return set(self._named_complexes)

    @named_complexes.setter
    def named_complexes(self, value):
        # NOTE: there is no need to enforce this.
        raise PeppercornUsageError('Named complexes must be provided at initialization!')












    @property
    def max_complex_size(self):
        return self._max_complex_size

    @max_complex_size.setter
    def max_complex_size(self, value):
        assert isinstance(value, int)
        if not all(cx.size <= value for cx in self._initial_complexes):
            msize = max([c.size for c in self._initial_complexes])
            raise PeppercornUsageError(
                    'Maximum complex size must include all input complexes.',
                    'Set to at least: {}'.format(msize))
        self._max_complex_size = value

    @property
    def max_reaction_count(self):
        return self._max_reaction_count

    @max_reaction_count.setter
    def max_reaction_count(self, value):
        assert isinstance(value, int) and value >= len(self._reactions)
        self._max_reaction_count = value

    @property
    def max_complex_count(self):
        return self._max_complex_count

    @max_complex_count.setter
    def max_complex_count(self, value):
        assert isinstance(value, int) and value >= len(self._initial_complexes)
        self._max_complex_count = value

    @property
    def k_slow(self):
        return self._k_slow

    @k_slow.setter
    def k_slow(self, value):
        if value > 0:
            (rc, k_rc) = (0, None)
            while True:
                rc += 1
                k_rc = opening_rate(rc, dG_bp = self.dG_bp)
                if k_rc < value: break
            self.release_cutoff = max(rc, self.release_cutoff)

        if 0 < self.k_fast < value :
            raise PeppercornUsageError('k-slow must not be bigger than k-fast.')

        self._k_slow = value

    @property
    def k_fast(self):
        return self._k_fast

    @k_fast.setter
    def k_fast(self, value):
        if 0 < value < self.k_slow:
            raise PeppercornUsageError('k-fast must not be smaller than k-slow.')
        self._k_fast = value

    @property
    def release_cutoff(self):
        if self._release_11 != self._release_12 :
            raise PeppercornUsageError('Ambiguous release cutoff request.')
        return self._release_11

    @release_cutoff.setter
    def release_cutoff(self, value):
        assert isinstance(value, int) and value >= 0
        self._release_11 = value
        self._release_12 = value

    @property
    def release_cutoff_1_1(self):
        return self._release_11

    @release_cutoff_1_1.setter
    def release_cutoff_1_1(self, value):
        assert isinstance(value, int) and value >= 0
        self._release_11 = value

    @property
    def release_cutoff_1_2(self):
        return self._release_12

    @release_cutoff_1_2.setter
    def release_cutoff_1_2(self, value):
        assert isinstance(value, int) and value >= 0
        self._release_12 = value

    @property
    def domains(self):
        domains = set()
        for cplx in self._initial_complexes:
            [domains.add(d) for d in cplx.domains]
        return list(domains)

    @property
    def reactions(self):
        """ List of detailed reactions before or after enumeration. """
        return list(self._reactions)

    @property
    def detailed_reactions(self):
        return list(self._reactions)

    @property
    def condensed_reactions(self):
        if self.condensation is None:
            self.condense()
        return self.condensation.condensed_reactions

    @property
    def resting_macrostates(self):
        """
        List of resting states enumerated. :py:meth:`.enumerate` must be
        called before access.
        """
        if self._resting_macrostates is None:
            raise PeppercornUsageError("enumerate not yet called!")
        return self._resting_macrostates[:]

    @property
    def complexes(self):
        """
        List of complexes enumerated. :py:meth:`.enumerate` must be
        called before access.
        """
        if self._complexes is None:
            raise PeppercornUsageError("enumerate not yet called!")
        return self._complexes[:]

    @property
    def resting_complexes(self):
        """
        List of complexes enumerated that are within resting states.
        :py:meth:`.enumerate` must be called before access.
        """
        if self._resting_complexes is None:
            raise PeppercornUsageError("enumerate not yet called!")
        return self._resting_complexes[:]

    @property
    def transient_complexes(self):
        """
        List of complexes enumerated that are not within resting sets (e.g.
        complexes which are transient). :py:meth:`.enumerate` must be
        called before access.
        """
        if self._transient_complexes is None:
            raise PeppercornUsageError("enumerate not yet called!")
        return self._transient_complexes[:]

    def __eq__(self, other):
        raise NotImplementedError('No notion of equality implemented for Enumerator objects.')

    def __ne__(self, other):
        return not (self == other)

    def dry_run(self):
        """
        Make it look like you've enumerated, but actually do nothing...
        """
        if self._complexes:
            raise PeppercornUsageError('Cannot call dry-run after enumeration!')

        rxs = [r for r in list(self._reactions) if len(r.reactants) == 1]
        info = segment_neighborhood(self.initial_complexes, rxs, represent = self.representatives)

        self._complexes = list(self.initial_complexes)
        self._resting_complexes = info['resting_complexes']
        self._transient_complexes = info['transient_complexes']
        self._resting_macrostates = info['resting_macrostates']

    def enumerate(self):
        """ Enumerate the reaction network from initial complexes.

        The results of this method are accessible via the following properties
        of this class:
            - complexes
            - resting_complexes
            - transient_complexes
            - resting_macrostates
            - reactions
            - detailed_reactions
        """

        # List E contains enumerated resting complexes. Every time a new
        # complex is added (from S), all cross reactions with other resting
        # complexes are enumerated. Complexes remain in this list throughout
        # function execution, the products of slow reactions go into list B.
        self._E = []

        # List S contains newly determined resting complexes after all fast
        # reactions have been enumerated. They will be moved to E and thereby
        # tested for cross reactions with set E. 
        self._S = []

        # List T contains newly determined transient states after fast-reaction
        # enumeration. These complexes will remain in this list throughout
        # function execution.
        self._T = []

        # List N contains the neighborhood of some initial complexes with all
        # fast reactions enumerated. Complexes in N have not yet been
        # characterized as transient or resting sets.
        self._N = []

        # List F contains components of the current 'neighborhood' which have
        # not yet been reactants for potential fast reactions.  They will be
        # moved to N once they were enumerated.
        self._F = []

        # List B contains initial complexes, or products of bimolecular
        # reactions that have had no reactions enumerated yet. They will be
        # moved to F for their fast neighborhood to be enumerated.
        self._B = list(self._initial_complexes)

        def do_enumerate():
            log.debug("Enumerating fast reactions from initial complexes.")
            while len(self._B) > 0:
                source = self._B.pop()
                self.process_fast_neighborhood(source)

            log.debug("Enumerating slow reactions between resting complexes.")
            while len(self._S) > 0:
                element = self._S.pop() if self.DFS else self._S.pop(0)
                log.debug(f"Finding slow reactions from '{element}'. " + \
                          f"({len(self._S)} remaining in S).")
                slow_reactions = set(self.get_slow_reactions(element))
                self._B = list(self.get_new_products(slow_reactions))
                log.debug(f"Generated {len(slow_reactions)} new slow reactions " + \
                          f"with {len(self._B)} new products.")

                self._E.append(element)
                self._reactions |= slow_reactions
                if self.interactive:
                    # Display new reactions in interactive mode
                    self.reactions_interactive(element, slow_reactions, 'slow')

                while len(self._B) > 0:
                    # Check whether too many complexes have been generated
                    if len(self._E) + len(self._T) + len(self._S) > self._max_complex_count:
                        raise PolymerizationError("Too many complexes enumerated! ({:d})".format(
                                len(self._E) + len(self._T) + len(self._S)))
                    # Check whether too many reactions have been generated
                    if len(self._reactions) > self._max_reaction_count:
                        raise PolymerizationError("Too many reactions enumerated! ({:d})".format(
                                len(self._reactions)))
                    # Generate a neighborhood from `source`
                    source = self._B.pop()
                    self.process_fast_neighborhood(source)

        self._complexes = []
        self._resting_macrostates = []
        def finish(premature = False):
            # Will be called once enumeration halts, either because it's finished or
            # because too many complexes/reactions have been enumerated

            # copy E and T into complexes
            self._complexes += (self._E)
            self._complexes += (self._T)

            # preserve resting and transient complexes separately
            self._transient_complexes = self._T
            self._resting_complexes = self._E

            # If we're bailing because of too many reactions or complexes, search
            # self._reactions to ensure there are no reactions which contain
            # products that never made it into self._complexes ...
            if premature:
                self._resting_complexes += self._S
                self._complexes += self._S
                complexes = set(self._complexes)
                rm_reactions = []
                for reaction in self.reactions:
                    #NOTE: This should only matter in the Ctrl-C case, right?
                    reaction_ok = all((prod in complexes) for prod in reaction.products) and \
                                  all((reac in complexes) for reac in reaction.reactants)
                    if reaction_ok:
                        pass
                    else :
                        rm_reactions.append(reaction)
                self._reactions -= set(rm_reactions)

        if self.interruptible:
            try:
                do_enumerate()
                finish()
            except KeyboardInterrupt:
                log.warning("Interrupted. Gracefully exiting ...")
                finish(premature = True)
            except PolymerizationError as err:
                log.exception(err)
                log.error("Polymerization error. Gracefully exiting ...")
                finish(premature = True)
        else:
            do_enumerate()
            finish()

    def condense(self):
        self.condensation = PepperCondensation(self)
        self.condensation.condense()

    def to_pil(self, filename = None, **kwargs):
        if filename:
            with open(filename, 'w') as pil:
                write_pil(self, fh=pil, **kwargs)
            return ''
        else:
            return write_pil(self, **kwargs)

    def to_sbml(self, filename = None, **kwargs):
        if filename:
            with open(filename, 'w') as sbml:
                write_sbml(self, fh = sbml, **kwargs)
            return ''
        else:
            return write_sbml(self, **kwargs)

    def reactions_interactive(self, root, reactions, rtype='fast'):
        """
        Prints the passed reactions as a kernel string, then waits for keyboard
        input before continuing.
        """
        print("{} = {} ({})".format(root.name, root.kernel_string, rtype))
        if len(reactions) == 0:
            print("(No {} reactions)".format(rtype))
        for r in reactions:
            print("{} ({}: {})".format(r, r.rtype, r.const))
            print(" {}".format(r.kernel_string))
        input("\n[Press Enter to continue...]")

    def process_fast_neighborhood(self, source):
        """ Enumerate neighborhood of fast reactions.

        Takes a single complex, generates the 'neighborhood' of complexes
        reachable from that complex through fast reactions, classifies these
        complexes as transient or resting complexes, and modifies the global
        lists and list of reactions accordingly.

        Args:
            source (:obj:`PepperComplex`): Initial complex to generate a
                neighborhood of fast reactions.
        """
        log.debug(f"Processing neighborhood of '{source}'")
        assert len(self._N) == 0
        assert len(self._F) == 0
        self._F = [source]
        rxns_N = [] # reactions which are part of the current neighborhood.
        interrupted = False
        try:
            # Find fast reactions from `element`
            while len(self._F) > 0 :
                element = self._F.pop()
                log.debug(f"Finding fast reactions from '{element}'. " + \
                          f"({len(self._F)} remaining in F).")
                reactions = list(self.get_fast_reactions(element))
                nproducts = self.get_new_products(reactions)
                log.debug(f"Generated {len(reactions)} new fast reactions " + \
                          f"with {len(nproducts)} new products.")
                rxns_N += reactions
                self._F += list(nproducts)
                self._N.append(element)
                if self.interactive:
                    # Display new reactions in interactive mode
                    self.reactions_interactive(element, reactions, 'fast')
        except KeyboardInterrupt:
            log.warning(f"Exiting neighborhood of '{source}' prematurely ...")
            if self.interruptible:
                interrupted = True
            else:
                raise KeyboardInterrupt
        log.debug(f"Segmenting neighborhood of {source}: " + \
                  f"{len(self._N)} complexes and {len(rxns_N)} reactions.")
        sn = segment_neighborhood(self._N, rxns_N, represent = self.representatives)
        self._S += sn['resting_complexes']
        self._T += sn['transient_complexes']
        self._resting_macrostates += sn['resting_macrostates']
        self._reactions |= set(rxns_N)
        log.debug(f"Of {len(self._N)} complexes: " + \
                  "{:d} are transient, {:d} are resting.".format(
                    len(sn['transient_complexes']), 
                    len(sn['resting_complexes'])))
        log.debug("in a total of {:d} resting macrostates".format(
                    len(sn['resting_macrostates'])))
        self._N = [] # Reset neighborhood
        if interrupted:
            raise KeyboardInterrupt
        return

    def get_uni_reactions(self, cplx):
        maxsize = self.max_complex_size
        for move in UNI_REACTIONS:
            if move.__name__ == 'bind11':
                reactions = move(cplx, 
                                 max_helix = self.max_helix)
            elif move.__name__ == 'open1N':
                reactions = move(cplx, 
                                 max_helix = self.max_helix, 
                                 release_11 = self.release_cutoff_1_1,
                                 release_1N = self.release_cutoff_1_2,
                                 dG_bp = self.dG_bp)
            else:
                reactions = move(cplx, 
                                 max_helix = self.max_helix, 
                                 remote = not self.reject_remote)
            for rxn in reactions:
                if maxsize and any(p.size > maxsize for p in rxn.products):
                    log.warning(f"Ignoring unimolecular reaction {rxn}: " + \
                            "product complex size {:d} > --max-complex-size {:d}.".format(
                            max([p.size for p in rxn.products]), maxsize))
                    continue
                yield rxn
        return

    def get_bi_reactions(self, cplx1, cplx2):
        maxsize = self.max_complex_size
        for move in BI_REACTIONS:
            reactions = move(cplx1, cplx2, max_helix = self.max_helix)
            for rxn in reactions:
                if maxsize and any(p.size > maxsize for p in rxn.products):
                    log.warning(f"Ignoring bimolecular reaction {rxn}: " + \
                            "product complex size {:d} > --max-complex-size {:d}.".format(
                                max([p.size for p in rxn.products]), maxsize))
                    continue
                yield rxn
        return

    def get_fast_reactions(self, cplx):
        """ set: Fast reactions using complex as a reactant. """
        for rxn in self.get_uni_reactions(cplx):
            # Do unimolecular reactions that are (mostly) fast.
            if rxn.rate_constant[0] >= max(self._k_slow, self._k_fast):
                yield rxn
        return

    def get_slow_reactions(self, cplx):
        """ set: Slow reactions using complex as a reactant. """
        if self.k_fast > self.k_slow:
            # Do unimolecular reactions that are sometimes slow.
            for rxn in self.get_uni_reactions(cplx):
                if self.k_slow <= rxn.rate_constant[0] < self.k_fast:
                    log.info(f'Found unimolecular slow reaction: {rxn}')
                    yield rxn

        # Do bimolecular reactions that are always slow
        for rxn in self.get_bi_reactions(cplx, cplx):
            yield rxn
        for cplx2 in self._E:
            for rxn in self.get_bi_reactions(cplx, cplx2):
                yield rxn
        return

    def get_new_products(self, reactions):
        """ list: all products in the list of reactions. """
        B = set(self._B)
        ESTNF = set(self._E + self._S + self._T + self._N + self._F)
        new_products = set()
        for rxn in reactions:
            for product in rxn.products:
                # This must have been checked earlier ...
                assert product.size <= self.max_complex_size
                if product in B:
                    # If the product is in list B, then we need to remove it from
                    # that list so that it can be enumerated for self-interactions
                    # as part of this neighborhood
                    self._B.remove(product)
                    B.remove(product)
                if product not in ESTNF:
                    # The product has not been enumerated ...
                    new_products.add(product)
        assert (ESTNF - new_products) == ESTNF
        assert (ESTNF - B) == ESTNF
        return new_products

def enumerate_pil(pilstring, 
        is_file = False, 
        enumfile = None,
        detailed = True, 
        condensed = False, 
        enumconc = 'nM', 
        **kwargs):
    """ A wrapper function to directly enumerate a pilstring or file.

    Args:
        pilstring (str): Either a full pil string or a file name containing that string.
        is_file (bool, optional): True if pilstring is a filename, False otherwise. 
            Defaults to False.
        enumfile (str, optional): A filename for the enumerated output. The function 
            returns the output string if no enumfile is specified. Defaults to None.
        detailed (bool, optional): Returns the detailed reaction network. Defaults to True.
        condensed (bool, optional): Returns the condensed reaction network. Defaults to False.
        enumconc (str, optional): Chose concentation units for output: 'M', 'mM', 'uM', 'nM', ...
        **kwargs: Attributes of the peppercornenumerator.Enumerator object.

    Returns:
        Enumerator-object, Outputstring
    """
    global WASTEFUL
    if WASTEFUL: clear_memory()

    cxs, rxns = read_pil(pilstring, is_file)
    cplxs = list(cxs.values())
    init_cplxs = [x for x in cxs.values() if x.concentration is None or x.concentration[1] != 0]
    enum = Enumerator(init_cplxs, rxns, named_complexes = cplxs)
    # set kwargs parameters
    for k, w in kwargs.items():
        if hasattr(enum, k):
            setattr(enum, k, w)
        else:
            raise PeppercornUsageError('No Enumerator attribute called: {}'.format(k))
    enum.enumerate()
    outstring = enum.to_pil(enumfile, detailed = detailed, condensed = condensed, molarity = enumconc)
    return enum, outstring

def enumerate_ssw(sswstring, 
        is_file = False, 
        ssw_expl = False,
        ssw_conc = 100e-9,
        ssw_rxns = 'T25-utbr-leak-reduced',
        dry_run = True,
        enumfile = None,
        detailed = True, 
        condensed = False, 
        enumconc = 'nM', 
        **kwargs):
    """ A wrapper function to directly enumerate a seesaw string or file.

    Args:
        sswstring (str): Either a full seesaw string or a file name containing that string.
        is_file (bool, optional): True if sswstring is a filename, False otherwise. 
            Defaults to False.

        ssw_expl (bool, optional): Explict domain-level specification? Defaults to False.
        ssw_conc (str, optional): 1x concentation in 'M'. Defaulst to 1e-7.
        ssw_rxns (str, optional): Specify the kind of seesaw reactions. Defaults to 
            T25-utbr-leak-reduced.
        dry_run (bool, optional): True if peppercorn enumeration is desired. Defaults to False.
        enumfile (str, optional): A filename for the output. The function 
            returns the output string if no enumfile is specified. Defaults to None.
        detailed (bool, optional): Returns the detailed reaction network. Defaults to True.
        condensed (bool, optional): Returns the condensed reaction network. Defaults to False.
        enumconc (str, optional): Chose concentation units for output: 'M', 'mM', 'uM', 'nM', ...
        **kwargs: Attributes of the peppercornenumerator.Enumerator object.

    Returns:
        Enumerator-object, Outputstring
    """
    global WASTEFUL
    if WASTEFUL: clear_memory()

    utbr = kwargs.get('utbr_species', True)
    cxs, rxns = read_seesaw(sswstring, is_file, 
            explicit = ssw_expl,
            conc = ssw_conc,
            reactions = ssw_rxns)
    enum = Enumerator(list(cxs.values()), rxns)
    # set kwargs parameters
    for k, w in kwargs.items():
        if hasattr(enum, k):
            setattr(enum, k, w)
        else:
            raise PeppercornUsageError('No Enumerator attribute called: {}'.format(k))
    if dry_run:
        enum.dry_run()
    else:
        enum.enumerate()
    outstring = enum.to_pil(enumfile, detailed = detailed, condensed = condensed, molarity = enumconc)
    return enum, outstring

def segment_neighborhood(complexes, reactions, represent = None):
    """ Segmentation of a neighborhood of fast reactions.

    Note that all reactants of the specified reactions must be contained in the
    complexes. Since only the specified complexes are of interest, any other
    complexes that can be produced by the given reactions are not contained
    in the output. 

    Args:
        complexes (list[:obj:`PepperComplex`]): 
        reactions (list[:obj:`PepperReaction`]):
        represent (bool, optional):

    Returns:
        dict: 'resting_complexes', 'transient_complexes', 'resting_macrostates'.
    """
    cplxs = set(complexes)
    total = list(cplxs.union(*[rxn.products for rxn in reactions]))

    # filters reaction products such that there are only species from within
    # complexes. This is ok, because it must not change the assignment of SCCs!
    products = {k: set() for k in cplxs}
    for rxn in reactions:
        assert len(rxn.reactants) == 1
        products[rxn.reactants[0]] |= set().union((p for p in rxn.products if p in cplxs))
    SCCs = tarjans(cplxs, products)

    rxns_consuming = {k: [r for r in reactions if (k in r.reactants)] for k in total}

    resting_macrostates = []
    resting_complexes = []
    transient_complexes = []
    for scc in SCCs:
        if represent:
            for cx in scc:
                if cx in represent:
                    rcx = cx.name
                    break
            else:
                rcx = None
        else:
            rcx = None
        try:
            ms = PepperMacrostate(scc[:], name = rcx)
        except SingletonError as err:
            ms = err.existing
        for c in scc:
            for rxn in rxns_consuming[c]:
                assert rxn.reactants[0] in ms._complexes
                if rxn.products[0] not in ms._complexes:
                    # an exit reaction!
                    transient_complexes += (scc)
                    break
                assert len(rxn.products) == 1
            else:
                continue
            break
        else:
            resting_macrostates.append(ms)
            resting_complexes += (scc)
    return {
        'resting_macrostates': resting_macrostates,
        'resting_complexes': resting_complexes,
        'transient_complexes': transient_complexes
    }

