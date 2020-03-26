#
#  peppercornenumerator/enumerator.py
#  EnumeratorProject
#
from __future__ import absolute_import, print_function, division
from builtins import input

import logging
log = logging.getLogger(__name__)

import sys
import math

from peppercornenumerator.condense import PepperCondensation
from peppercornenumerator.objects import PepperMacrostate, PepperComplex
from peppercornenumerator.objects import DSDDuplicationError, DSDObjectsError
from peppercornenumerator.utils import PeppercornUsageError
from peppercornenumerator.input import read_pil, read_seesaw
from peppercornenumerator.output import write_pil, write_sbml
import peppercornenumerator.reactions as reactlib


# There should be better control of this limit -- only set it when
# necessary (Erik Winfree based on Chris Thachuk's advice...)
sys.setrecursionlimit(20000)


FAST_REACTIONS = {1: [reactlib.bind11, reactlib.open,
                      reactlib.branch_3way, reactlib.branch_4way]}
"""
Dictionary of reaction functions considered *fast* for a given "arity".
Keys are arities (e.g. 1 = unimolecular, 2 = bimolecular, 3 = trimolecular,
etc.), and values are lists of reaction functions. Currently, only
unimolecular fast reactions (arity = 1) are supported.
"""

SLOW_REACTIONS = {
    1: [],
    2: [reactlib.bind21]
}
"""
A dictionary of reaction functions considered *slow* for a given "arity".
Keys are arities (e.g. 1 = unimolecular, 2 = bimolecular, 3 = trimolecular,
etc.), and values are lists of reaction functions. Currently, only
bimolecular reactions (arity = 2) are by default slow. However, when using
k-fast, also unimolecular reactions can fall into the *slow* regime.
"""

class PolymerizationError(Exception):
    """Error class to catch polymerization."""

    def __init__(self, msg, val=None):
        self.message = msg
        if val :
            self.message += " ({})".format(val)
        super(PolymerizationError, self).__init__(self.message) 

class Enumerator(object):
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

    def __init__(self, initial_complexes, initial_reactions=None):
        # System initialization
        self._initial_complexes = initial_complexes
        assert isinstance(initial_complexes, list)
        for cplx in initial_complexes:
            if not cplx.is_connected:
                raise PeppercornUsageError('Initial complex is not connected: {}'.format(cplx))

        # list containing all detailed reactions after enumeration
        if initial_reactions:
            self._reactions = set(initial_reactions)
        else :
            self._reactions = set()

        # A list of "known" complexes, 
        #   e.g. expected to show up during enumeration.
        # This is useful when chosing representatives of resting macrostates.
        self._named_complexes = list(PepperComplex.MEMORY.values())
        
        # Enumeration results
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

        # Not officially supported parameters
        self.local_elevation = False
        self.p_min = 0 # minimum complex steady state probability to engage in slow reaction.

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
                k_rc = reactlib.opening_rate(rc, dG_bp = self.dG_bp)
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
    def remote_migration(self):
        print(DeprecationWarning('Peppercorn>=v0.7: use Enumerator.reject_remote instead of Enumerator.remote_migration.'))
        return not self.reject_remote

    @remote_migration.setter
    def remote_migration(self, remote):
        assert isinstance(remote, bool)
        print(DeprecationWarning('Peppercorn>=v0.7: use Enumerator.reject_remote instead of Enumerator.remote_migration.'))
        self.reject_remote = not remote

    @property
    def max_helix_migration(self):
        print(DeprecationWarning('Peppercorn>=v0.7: use Enumerator.max_helix instead of Enumerator.max_helix_migration.'))
        return self.max_helix

    @max_helix_migration.setter
    def max_helix_migration(self, max_helix):
        print(DeprecationWarning('Peppercorn>=v0.7: use Enumerator.max_helix instead of Enumerator.max_helix_migration.'))
        self.max_helix = max_helix

    @property
    def initial_complexes(self):
        """
        Complexes present in the system's initial configuration
        """
        return self._initial_complexes[:]

    @property
    def domains(self):
        domains = set()
        for cplx in self._initial_complexes:
            [domains.add(d) for d in cplx.domains]
        return list(domains)

    @property
    def reactions(self):
        """
        List of reactions enumerated. :py:meth:`.enumerate` must be
        called before access.
        """
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
    def resting_sets(self):
        print("# Deprecated function: Enumerator.resting_sets. Replace with Enumerator.resting_macrostates")
        return self.resting_macrostates

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

        rep = set(self._named_complexes) if self._named_complexes else set(self._initial_complexes)
        rxs = [r for r in list(self._reactions) if len(r.reactants) == 1]
        info = segment_neighborhood(self.initial_complexes, rxs, self.p_min, represent=rep)

        self._complexes = self.initial_complexes
        self._resting_complexes = info['resting_complexes']
        self._transient_complexes = info['transient_complexes']
        self._resting_macrostates = info['resting_macrostates']

    def enumerate(self):
        """
        Generates the reaction graph consisting of all complexes reachable from
        the initial set of complexes. Produces a full list of :py:meth:`complexes`, resting
        sets, and :py:meth:`reactions, which are stored in the associated members of this
        class.
        """

        # Will be called once enumeration halts, either because it's finished or
        # because too many complexes/reactions have been enumerated
        def finish(premature=False):
            # copy E and T into complexes
            self._complexes += (self._E)
            self._complexes += (self._T)

            # preserve resting and transient complexes separately
            self._transient_complexes = self._T
            self._resting_complexes = self._E

            # If we're bailing because of too many reactions or complexes, search
            # self._reactions to ensure there are no reactions which contain
            # products that never made it into self._complexes...
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
        self._B = self.initial_complexes

        self._complexes = []
        self._resting_macrostates = []

        def do_enumerate():
            log.debug("Fast reactions from initial complexes ...")
            while len(self._B) > 0:
                # Generate a neighborhood from `source`
                source = self._B.pop()
                self.process_neighborhood(source)

            # Consider slow reactions between resting set complexes
            log.debug("Slow reactions between resting set complexes ...")
            while len(self._S) > 0:

                # Find slow reactions from 'element'
                if self.DFS:
                    element = self._S.pop()
                else:
                    element = self._S.pop(0)

                log.debug("Slow reactions from complex {:s} ({:d} remaining in S)".format(
                    str(element), len(self._S)))
                slow_reactions = self.get_slow_reactions(element)
                self._E.append(element)

                # Find the new complexes which were generated
                self._B = self.get_new_products(slow_reactions)
                self._reactions |= set(slow_reactions)
                log.debug("Generated {:d} new slow reactions".format(len(slow_reactions)))
                log.debug("Generated {:d} new products".format(len(self._B)))

                # Display new reactions in interactive mode
                if self.interactive:
                    self.reactions_interactive(element, slow_reactions, 'slow')

                # Now find all complexes reachable by fast reactions from these new complexes
                while len(self._B) > 0:

                    # Check whether too many complexes have been generated
                    if (len(self._E) + len(self._T) + len(self._S) > self._max_complex_count):
                        raise PolymerizationError("Too many complexes enumerated!", 
                                len(self._E) + len(self._T) + len(self._S))

                    # Check whether too many reactions have been generated
                    if (len(self._reactions) > self._max_reaction_count):
                        raise PolymerizationError("Too many reactions enumerated!", 
                                len(self._reactions))

                    # Generate a neighborhood from `source`
                    source = self._B.pop()
                    self.process_neighborhood(source)

        if self.interruptible:
            try:
                do_enumerate()
                finish()
            except KeyboardInterrupt:
                log.warning("Interrupted; gracefully exiting...")
                finish(premature=True)
            except PolymerizationError as err:
                log.exception(err)
                log.error("Polymerization error; gracefully exiting...")
                finish(premature=True)
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
        if len(reactions) is 0:
            print("(No {} reactions)".format(rtype))
        for r in reactions:
            print("{} ({}: {})".format(r, r.rtype, r.const))
            print(" {}".format(r.kernel_string))
        input("\n[Press Enter to continue...]")

    def process_neighborhood(self, source):
        """ Enumerate neighborhood of fast reactions.

        Takes a single complex, generates the 'neighborhood' of complexes
        reachable from that complex through fast reactions, classifies these
        complexes as transient or resting complexes, and modifies the global
        lists and list of reactions accordingly.

        Args:
            source (:obj:`PepperComplex`): Initial complex to generate a
                neighborhood of fast reactions.
        """

        assert len(self._F) == 0
        self._F = [source]

        # N_reactions holds reactions which are part of the current neighborhood
        N_reactions = []

        log.debug("Processing neighborhood: %s" % source)

        interrupted = False
        try:
            while len(self._F) > 0 :

                # Find fast reactions from `element`
                element = self._F.pop()
                log.debug("Fast reactions from {:s}... ({:d} remaining in F)".format(
                    str(element), len(self._F)))

                # Return valid fast reactions:
                reactions = self.get_fast_reactions(element)
                   
                # Add new products to F
                new_products = self.get_new_products(reactions)
                self._F += (new_products)

                # Add new reactions to N_reactions
                N_reactions += (reactions)
                self._N.append(element)

                log.debug("Generated {:d} new fast reactions.".format(len(reactions)))
                log.debug("Generated {:d} new products.".format(len(new_products)))
                   
                # Display new reactions in interactive mode
                if self.interactive:
                    self.reactions_interactive(element, reactions, 'fast')

        except KeyboardInterrupt:
            log.warning("Exiting neighborhood %s prematurely..." % source)
            if self.interruptible:
                interrupted = True
            else :
                raise KeyboardInterrupt

        log.debug("In neighborhood %s..." % source)
        log.debug("Segmenting %d complexes and %d reactions" %
                      (len(self._N), len(N_reactions)))

        # Now segment the neighborhood into transient and resting complexes
        # by finding the strongly connected components.
        rep = set(self._named_complexes) if self._named_complexes else set(self._initial_complexes)
        segmented_neighborhood = segment_neighborhood(self._N, N_reactions, 
                                                      p_min = self.p_min,
                                                      represent = rep)

        # Resting complexes are added to S
        self._S += (segmented_neighborhood['resting_complexes'])

        # Transient complexes are added to T
        self._T += (segmented_neighborhood['transient_complexes'])

        # Resting macrostates are added to global list
        self._resting_macrostates += (segmented_neighborhood['resting_macrostates'])

        # Reactions from this neighborhood are added to the list
        self._reactions |= set(N_reactions)

        # Reset neighborhood
        log.debug("Generated {:d} new fast reactions".format(len(N_reactions)))
        log.debug("Generated {:d} new complexes: ({:d} transient, {:d} resting)".format(
            len(self._N), len(segmented_neighborhood['transient_complexes']), 
                len(segmented_neighborhood['resting_complexes'])))
        self._N = []

        log.debug("Generated {:d} resting macrostates".format(
            len(segmented_neighborhood['resting_macrostates'])))
        log.debug("Done processing neighborhood: {:s}".format(str(source)))

        if interrupted:
            raise KeyboardInterrupt

    def get_slow_reactions(self, complex):
        """
        Returns a list of slow reactions possible using complex and other
        complexes in list E as reagents.

        This only supports unimolecular and bimolecular reactions. Could be
        extended to support arbitrary reactions.
        """

        maxsize = self._max_complex_size

        reactions = []

        # Do unimolecular reactions that are always slow... not supported
        # anymore...
        assert SLOW_REACTIONS[1] == []

        # Do unimolecular reactions that are sometimes slow
        if self._k_fast > self._k_slow:
            for move in FAST_REACTIONS[1]:
                if move.__name__ == 'open':
                    move_reactions = move(complex, 
                            max_helix = self.max_helix, 
                            release_11 = self._release_11,
                            release_1N = self._release_12,
                            dG_bp = self.dG_bp)
                else :
                    move_reactions = move(complex, 
                            max_helix = self.max_helix, 
                            remote = not self.reject_remote)

                if self.local_elevation:
                    reactions += [rxn for rxn in move_reactions if self._k_slow <= local_elevation_rate(rxn) < \
                                self._k_fast]
                else :
                    reactions += [r for r in move_reactions if self._k_slow <= r.const < self._k_fast]
            for rxn in reactions:
                log.info('adding unimolecular slow reaction {}'.format(rxn.full_string()))

        # Do bimolecular reactions
        for move in SLOW_REACTIONS[2]:
            reactions += (move(complex, complex))
            for complex2 in self._E:
                reactions += (move(complex, complex2))

        # And now remove all over max-complex-size
        valid_reactions = []
        for rxn in reactions: 
            if maxsize and not all(p.size <= maxsize for p in rxn.products) :
                log.warning("Product complex size (={}) larger than --max-complex-size(={}). Ignoring slow reaction {}!".format(max([p.size for p in rxn.products]), maxsize, str(rxn)))
                continue
            valid_reactions.append(rxn)

        return valid_reactions

    def get_fast_reactions(self, cplx, rtypes = None, restrict=True):
        """Returns a list of fast reactions possible using complex as a reagent.

        Args:
            cplx: A reactant complex.
            rtypes (str, optional): Only reactions of a certain type.
            restrict (bool, optional): Use golbal parameters release_cutoff, 


        """

        maxsize = self._max_complex_size

        # Do unimolecular reactions
        reactions = []
        for move in FAST_REACTIONS[1]:
            if rtypes and move.__name__ not in rtypes:
                continue
            if move.__name__ == 'open':
                if restrict:
                    move_reactions = move(cplx, 
                        max_helix = self.max_helix, 
                        release_11 = self._release_11,
                        release_1N = self._release_12,
                        dG_bp = self.dG_bp)
                else :
                    move_reactions = move(cplx, 
                        max_helix = self.max_helix, 
                        release_11 = 0, release_1N = 0,
                        dG_bp = self.dG_bp)
            else :
                move_reactions = move(cplx, 
                        max_helix = self.max_helix, 
                        remote = not self.reject_remote)
            
            for rxn in move_reactions: 
                if maxsize and any(p.size > maxsize for p in rxn.products):
                    log.warning("Product complex size (={}) larger than --max-complex-size(={}). Ignoring fast reaction {}!".format( max([p.size for p in rxn.products]), maxsize, str(rxn)))
                    continue
                reactions.append(rxn)

        if restrict: # apply the k-fast / k-slow filter
            if self.local_elevation:
                # let's first remove reactions that cannot pass...
                reactions = [rxn for rxn in reactions if rxn.const >= self._k_slow]
                assert cplx._elevation is None
                el = self.get_local_elevation(cplx, reactions)
                cplx._elevation = el
                reactions = [rxn for rxn in reactions if local_elevation_rate(rxn, el) >= \
                                max(self._k_slow, self._k_fast)]
            else :
                reactions = [rxn for rxn in reactions if rxn.const >= max(self._k_slow, self._k_fast)]

        return reactions

    def get_new_products(self, reactions):
        """
        Checks the products in the list of reactions. Returns the new complexes in a list.

        """
        new_products = set()

        ESTNF = set(self._E + self._S + self._T + self._N + self._F)
        B = set(self._B)

        # Loop over every reaction
        for reaction in reactions:

            # Check every product of the reaction to see if it is new
            for (i, product) in enumerate(reaction.products):

                # This should have been checked earlier...
                assert product.size <= self._max_complex_size

                if product in B:
                    # If the product is in list B, then we need to remove it from
                    # that list so that it can be enumerated for self-interactions
                    # as part of this neighborhood
                    self._B.remove(product)
                    B.remove(product)

                # has not been enumerated...
                if product not in ESTNF:
                   new_products.add(product)

        assert (ESTNF - new_products) == ESTNF
        assert (ESTNF - B) == ESTNF

        return list(new_products)

    def get_local_elevation(self, cplx, reactions):
        """Calculate the local elevation of a complex.

        Local elevation is the inverse of the maximum free energy gain when
        applying a sequence of compatible moves. In other words, the energy
        difference between this complex and its gradient decent local minimum.
        Because our rate model is insufficient, this is just an approximation: 

            (1) Take all unary open/bind reactions, and their reverse reaction. 
            (2) Calculate the free energy change for each reversible reaction: 
                dG = ln(k1/k2)
            (3) Solve the max-clique problem to identify compatible downhill
                reactions. That means, reactions that can be applied
                successively leading to the same local minimum. 
            (4) Calculate the cummulative free energy change for every path
                that leads to a local minimum. The maximum value is the local
                elevation.

        """

        def rev_rtype(rtype, arity):
            """Returns the reaction type of a corresponding reverse reaction. """
            if rtype == 'open' and arity == (1,1):
                return 'bind11'
            elif rtype == 'open' and arity == (1,2):
                return 'bind21'
            elif rtype == 'bind11' or rtype == 'bind21':
                return 'open'
            else:
                raise NotImplementedError

        def elevation(rxn):
            assert rxn.arity == (1,1)
            if rxn.reverse_reaction is None:
                rr = rev_rtype(rxn.rtype, rxn.arity)
                r = self.get_fast_reactions(rxn.products[0], rtypes = [rr], restrict=False)
                r = [x for x in r if sorted(x.products) == sorted(rxn.reactants)]
                assert len(r) == 1
                if len(r) == 1:
                    rxn.reverse_reaction = r[0]
                    r[0].reverse_reaction = rxn
                else :
                    rxn.reverse_reaction = False
            
            dG = math.log(rxn.const/rxn.reverse_reaction.const)
            # downhill reaction has dG > 0
            return dG if dG > 0 else 0

        def try_move(reactant, rotations, rtype, meta):
            """
            Args:
                rectant: Is the product of a previous rection
                rotations: Rotations needed to rotate this product back in reactant form
                rtype: the type of the new reaction
                meta: the find_on_loop parameters for the new reaction
            """
            
            (invader, target, x_linker, y_linker) = meta

            # First, make sure invader and target locations map to the reactant-rotation.
            if rotations is not None and rotations != 0:
                invaders = [reactant.rotate_location(x, rotations) for x in invader.locs]
                targets = [reactant.rotate_location(x, rotations) for x in target.locs]
            else :
                invaders = list(invader.locs)
                targets = list(target.locs)

            # forbid max helix.... :-(
            assert len(invaders) == 1
            assert len(targets) == 1

            def triple(loc):
                return (reactant.get_domain(loc), reactant.get_structure(loc), loc)

            if rtype == 'bind11':
                # fore i in invaders...
                results = reactlib.find_on_loop(reactant, invaders[0], reactlib.filter_bind11)
            elif rtype == 'open':
                return reactant.get_structure(invaders[0]) == targets[0]
            else : 
                raise NotImplementedError('rtype not considered in local neighborhood: {}'.format(
                    rtype))

            targets = list(map(triple, targets))
            for [s,x,t,y] in results:
                if t == targets:
                    return True
            return False
        
        # Calculate the local elevation of a complex using only 1-1 reactions.
        compat = dict()
        for rxn in reactions:
            # exclude open and branch-migration with arity 1,2
            if rxn.arity != (1,1): continue
            if rxn.rtype in ('branch-3way', 'branch-4way'): continue
            compat[rxn] = set()
            for other in reactions:
                if rxn == other: continue
                if other.arity != (1,1): continue
                if other.rtype in ('branch-3way', 'branch-4way'): continue
                # Try to append the other reaction to this reaction
                success = try_move(rxn.products[0], rxn.rotations, other.rtype, other.meta)
                if success:
                    # We can apply other to the product of rxn
                    compat[rxn].add(other)

        # https://en.wikipedia.org/wiki/Clique_(graph_theory)
        cliques = []
        for p in sorted(compat):
            vert = p
            if compat[p] == set():
                cliques.append(set([p]))
                continue
            for l in list(compat[p]):
                edge = (p,l)
                processed = False
                for e, c in enumerate(cliques):
                    if p in c and (l in c or all(l in compat[y] for y in c)):
                        c.add(l)
                        processed = True
                if not processed:
                    cliques.append(set([p,l]))

        elevations = []
        for c in cliques:
            elevations.append(sum(map(elevation, c)))

        if elevations:
            eleven = max(elevations)
        else :
            eleven = 0

        return eleven

def enumerate_pil(pilstring, 
        is_file = False, 
        enumfile = None,
        detailed = True, 
        condensed = False, 
        enumconc = 'M', 
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
    cxs, rxns, comp = read_pil(pilstring, is_file, composite = True)

    init_cplxs = list(filter(lambda x: x._concentration is None or \
                    float(x._concentration[1]) != (0.0), cxs.values()))

    enum = Enumerator(init_cplxs, rxns)

    # set kwargs parameters
    for k, w in kwargs.items():
        if hasattr(enum, k):
            setattr(enum, k, w)
        else:
            raise PeppercornUsageError('No Enumerator attribute called: {}'.format(k))

    enum.enumerate()

    outstring = enum.to_pil(enumfile, detailed = detailed, condensed = condensed, 
                composite = comp, molarity = enumconc)
    
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
        enumconc = 'M', 
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


def segment_neighborhood(complexes, reactions, p_min=None, represent=None):
    """
    Segmentation of a potentially incomplete neighborhood. That means only the
    specified complexes are interesting, all others should not be returned.

    Beware: Complexes must contain all reactants in reactions *and* there
    must not be any incoming fast reactions into complexes, other than
    those specified in reactions. Thus, we can be sure that the SCCs found here
    are consistent with SCCs found in a previous iteration.

    Args:
        complexes(list[:obj:`PepperComplex`])
    """
    index = 0
    S = []
    SCCs = []

    total = complexes[:]
    for rxn in reactions: 
        total += rxn.products

    total = list(set(total))

    # Set up for Tarjan's algorithm
    for c in total: c._index = None

    # filters reaction products such that there are only species from within complexes
    # this is ok, because it must not change the assignment of SCCs.
    rxns_within = {k: [] for k in complexes}
    rxns_consuming = {k: [r for r in reactions if (k in r.reactants)] for k in total}
    for rxn in reactions:
        assert len(rxn.reactants) == 1
        for product in rxn.products:
            if product in complexes:
                rxns_within[rxn.reactants[0]].append(product)
        rxns_consuming[rxn.reactants[0]]

    def tarjans_scc(cplx, index):
        """
        Executes an iteration of Tarjan's algorithm (a modified DFS) starting
        at the given node.
        """
        # Set this node's tarjan numbers
        cplx._index = index
        cplx._lowlink = index
        index += 1
        S.append(cplx)

        for product in rxns_within[cplx]:
            # Product hasn't been traversed; recurse
            if product._index is None :
                index = tarjans_scc(product, index)
                cplx._lowlink = min(cplx._lowlink, product._lowlink)

            # Product is in the current neighborhood
            elif product in S:
                cplx._lowlink = min(cplx._lowlink, product._index)
    
        if cplx._lowlink == cplx._index:
            scc = []
            while True:
                next = S.pop()
                scc.append(next)
                if next == cplx:
                    break
            SCCs.append(scc)
        return index

    # We now perform Tarjan's algorithm, marking nodes as appropriate
    for cplx in complexes:
        if cplx._index is None:
            tarjans_scc(cplx, index)

    resting_macrostates = []
    transient_macrostates = []
    resting_complexes = []
    transient_complexes = []

    for scc in SCCs:
        try:
            if represent :
                rcx = [cx for cx in scc if cx in represent]

            rcx = sorted(rcx)[0] if rcx else None
            ms = PepperMacrostate(scc[:], prefix='', representative = rcx)
        except DSDDuplicationError as e:
            assert set(e.existing.complexes) == set(scc)
            ms = e.existing

        for c in scc:
            for rxn in rxns_consuming[c]:
                ms.add_reaction(rxn)

        if ms.is_transient:
            transient_complexes += (scc)
        else :
            resting_macrostates.append(ms)

            if p_min:
                for (c, s) in ms.get_stationary_distribution():
                    if s < p_min:
                        transient_complexes.append(c)
                    else :
                        resting_complexes.append(c)
            else:
                resting_complexes += (scc)

    resting_macrostates.sort()
    resting_complexes.sort()
    transient_complexes.sort()

    return {
        'resting_macrostates': resting_macrostates,
        'resting_complexes': resting_complexes,
        'transient_complexes': transient_complexes
    }

def local_elevation_rate(rxn, el=None):
    """Re-calculate a transition rate based on local elevation.

    Right now, downhill reactions are only passed through, but that
    probably needs some more thought about the consequences.

    """
    if el is None:
        assert rxn.arity[0] == 1
        cplx = rxn.reactants[0]
        if cplx._elevation is None:
            return None
        el = cplx._elevation

    if rxn.arity != (1,1) or rxn.rtype in ('branch-3way', 'branch-4way'): 
        # we don't know the reverse rate, ...
        return 1/(1+math.e**(el)) * rxn.const

    elif rxn.const < rxn.reverse_reaction.const:
        # this is a true uphill reaction...
        return 1/(1+math.e**(el)) * rxn.const
    
    # it is a downhill reaction, elevation does not matter ...
    assert rxn.arity[0] == 1
    return rxn.const
 
