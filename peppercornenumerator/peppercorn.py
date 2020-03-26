#!/usr/bin/env python
#
#  scripts/peppercorn
#  EnumeratorProject
#
from __future__ import absolute_import, division, print_function

import logging

import os
import sys
import argparse

# Import global default variables from peppercornenumerator library
import peppercornenumerator 
from peppercornenumerator import Enumerator, __version__
from peppercornenumerator.enumerator import FAST_REACTIONS
from peppercornenumerator.input import read_pil, read_seesaw, ParseException
from peppercornenumerator.reactions import branch_3way, branch_4way, opening_rate

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
            l = list(zip(keys, colors.colors))
        return "\n".join([(c + str(i) + colors.ENDC) for i, c in l])

class ColorFormatter(logging.Formatter):
    def __init__(self, msg, use_color=True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color
        self.COLORS = {
            'DEBUG': colors.CYAN,
            'INFO': colors.BLUE,
            'WARNING': colors.YELLOW,
            'ERROR': colors.RED,
            'Exception': colors.PINK,
        }
        self.RESET = colors.ENDC

    def format(self, record):
        levelname = record.levelname
        if self.use_color:
            record.levelname = self.COLORS[levelname] + \
                levelname + self.RESET
        return logging.Formatter.format(self, record)

def add_peppercorn_args(parser):
    output    = parser.add_argument_group('Peppercorn output format')
    semantics = parser.add_argument_group('Peppercorn reaction semantics')
    limits    = parser.add_argument_group('Peppercorn polymerization parameters')
    devel     = parser.add_argument_group('Peppercorn performance analysis')


    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument( '-v', '--verbose', action='count', default=0,
        help="Print logging output. (-vv increases verbosity.)")
    parser.add_argument('--logfile', default='', action='store', metavar='<str>',
        help="""Redirect logging information to a file.""")
    parser.add_argument('input_filename', default=None, nargs='?', metavar='<str>',
            help="Path to the input file.")

    output.add_argument('-o', '--output', dest='output_filename', default=None, metavar='<str>',
        help="""Redirect output to a file.""")
    output.add_argument('--sbml', default=None, metavar='<str>',
        help="""Print system as SBML (XML) file.""")
    output.add_argument('-c', '--condensed', action='store_true',
        help="Condense reactions into only resting complexes.")
    output.add_argument('-d', '--detailed', action='store_true',
        help="Print detailed reactions even if --condensed is chosen.")
    output.add_argument('--dry-run', action='store_true', 
        help="Dry run: read input, write output. Do not enumerate any reactions.")
    output.add_argument("--concentration-unit", default='M', action='store',
        choices=('M', 'mM', 'uM', 'nM', 'pM'),
        help="""Specify output concentration units for species and reaction rates.""")
    output.add_argument("--time-unit", default='s', action='store',
        choices=('s', 'm', 'h'),
        help="""Specify output time units for reaction rates.""")
    output.add_argument('--seesaw-explicit', action='store_true',
        help="""When using seesaw language input format, translate it into
        explicit form.""")
    output.add_argument('--seesaw-conc', default=100e-9, type=float, metavar='<flt>',
        help="""When using seesaw language input format, specify 1x
        concentrations [M].""")
    output.add_argument('--seesaw-reactions', default='', action='store',
        #choices=('seesaw', 'T20', 'T25', 'utbr', 'leak', 'reduced'),
        help=argparse.SUPPRESS)

    limits.add_argument('--max-complex-size', type=int, default=6, metavar='<int>',
        help="""Maximum number of strands allowed in a complex (used to prevent
        polymerization).""")
    limits.add_argument('--max-complex-count', type=int, default=200, metavar='<int>',
        help="""Maximum number of complexes that may be enumerated before the
        enumerator halts.""")
    limits.add_argument('--max-reaction-count', type=int, default=1000, metavar='<int>',
        help="Maximum number of reactions that may be enumerated before the enumerator halts.")

    semantics.add_argument('--k-slow', default=0.0, type=float, metavar='<flt>',
        help="Unimolecular reactions slower than this rate will be discarded.")
    semantics.add_argument('--k-fast', default=0.0, type=float, metavar='<flt>',
        help="Unimolecular reactions slower than this rate will be marked as slow.")
    semantics.add_argument('--p-min', default=0.0, type=float, metavar='<flt>',
        #help="""Minimal occupancy of a complex in steady state to engage in slow reactions.""")
        help=argparse.SUPPRESS)
    semantics.add_argument('--dG-bp', default = -1.7, type = float, metavar = '<flt>',
        help="""Adjust the average strength [kcal/mol] of a base-pair for toehold-binding 
        (affects only the opening rate). """)
    semantics.add_argument('--local-elevation', action='store_true', 
        #help="""Local probability threshold to accept an unfavorable reaction.""")
        help=argparse.SUPPRESS)

    semantics.add_argument('-L', '--release-cutoff', default=None, type=int, metavar='<int>',
        help="""Maximum number of bases that will be released spontaneously in
        an `open` reaction.""")
    semantics.add_argument('--release-cutoff-1-1', type=int, default=7, metavar='<int>',
        help="""Maximum number of bases that will be released spontaneously in
        an `open` reaction with one product.""")
    semantics.add_argument('--release-cutoff-1-2', type=int, default=7, metavar='<int>',
        help="""Maximum number of bases that will be released spontaneously in
        an `open` reaction with two products.""")

    semantics.add_argument('--no-max-helix', action='store_true', 
        help="Do not apply 'max helix' semantics.")
    semantics.add_argument('--ignore-branch-3way', action='store_true', 
        help="Ignore 3-way branch migration reactions during enumeration.")
    semantics.add_argument('--ignore-branch-4way', action='store_true',
        help="Ignore 4-way branch migration reactions during enumeration.")
    semantics.add_argument('--reject-remote', action='store_true', 
        help="""Discard remote toehold mediated 3-way and 4-way branch
        migration reactions.""")

    devel.add_argument('--interactive', action='store_true', 
        help="""Show new reactions after each step, confirm them in an
        interactive session.""")
    devel.add_argument('--interruptible', action='store_true', 
        help="""If polymerization errors are raised or after a keyboard
        interrupt, proceed with neighborhood segmentation and CRN
        condensation.""")
    devel.add_argument('--bfs-ish', action='store_true', dest='bfs',
        help="When searching for bimolecular reactions, look to the oldest complexes first.")
    devel.add_argument('--profile', action='store_true',
        help="Provide statistical profiling of enumeration (not condensation).")

    return

def set_handle_verbosity(h, v):
    if v == 0:
        h.setLevel(logging.WARNING)
    elif v == 1:
        h.setLevel(logging.INFO)
    elif v == 2:
        h.setLevel(logging.DEBUG)
    elif v >= 3:
        h.setLevel(logging.NOTSET)

def main():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description="""Peppercorn: Domain-level nucleic acid reaction enumerator.""")
    add_peppercorn_args(parser)
    args = parser.parse_args()
 
    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    title = "Peppercorn Domain-level Reaction Enumerator"
    logger = logging.getLogger('peppercornenumerator')
    logger.setLevel(logging.DEBUG)

    if args.logfile:
        banner = "{} {}".format(title, __version__)
        fh = logging.FileHandler(args.logfile)
        formatter = logging.Formatter('%(levelname)s - %(message)s')
        set_handle_verbosity(fh, args.verbose)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    else:
        banner = "{} {}".format(colors.BOLD + title + colors.ENDC, 
                                  colors.GREEN + __version__ + colors.ENDC)
        ch = logging.StreamHandler()
        formatter = ColorFormatter('%(levelname)s %(message)s', use_color = True)
        set_handle_verbosity(ch, args.verbose)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    logger.info(banner)

    systeminput = args.input_filename
    if not systeminput :
        if sys.stdout.isatty():
            logger.info("Reading file from STDIN, Ctrl-D to stop")
        systeminput = ''
        for l in sys.stdin:
            systeminput += l
        if args.interactive:
            loger.error("Interactive mode needs to read input from file, not STDIN.")
            raise SystemExit

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Input parsing to set initial complexes for enumeration #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    composite = None
    try :
        complexes, reactions, composite = read_pil(systeminput, 
                args.input_filename is not None, composite=True)
    except ParseException as ex_pil:
        try :
            complexes, reactions = read_seesaw(systeminput, 
                    args.input_filename is not None, 
                    conc = args.seesaw_conc, 
                    explicit = args.seesaw_explicit, 
                    reactions = args.seesaw_reactions)
        except ParseException as ex_ssw:
            logger.error('Pil-format parsing error:')
            logger.error('Cannot parse line {:5d}: "{}"'.format(ex_pil.lineno, ex_pil.line))
            logger.error('                          {} '.format(' ' * (ex_pil.col-1) + '^'))
            logger.error('SeeSaw-format parsing error:')
            logger.error('Cannot parse line {:5d}: "{}"'.format(ex_ssw.lineno, ex_ssw.line))
            logger.error('                          {} '.format(' ' * (ex_ssw.col-1) + '^'))
            raise SystemExit


    if args.dry_run:
        enum = Enumerator(list(complexes.values()), reactions)
    else:
        init_cplxs = [x for x in list(complexes.values()) if \
                x.concentration is None or float(x.concentration.value) != 0]
        enum = Enumerator(init_cplxs, reactions)

    # Log initial complexes
    logger.info("")
    logger.info("Initial complexes:")
    for c in enum.initial_complexes:
        logger.info("{}: {}".format(c, c.kernel_string))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Transfer options to enumerator object #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    logger.info("")
    logger.info("Enumeration settings:")
    enum.max_complex_size = args.max_complex_size
    logger.info("Max complex size = {}".format(enum.max_complex_size))
    enum.max_complex_count = max(args.max_complex_count, len(complexes))
    logger.info("Max complex count = {}".format(enum.max_complex_count))
    enum.max_reaction_count = max(args.max_reaction_count, len(reactions))
    logger.info("Max reaction count = {}".format(enum.max_reaction_count))
    enum.max_helix = not args.no_max_helix
    logger.info('Max-helix semantics = {}'.format(enum.max_helix))
    enum.reject_remote = args.reject_remote 
    logger.info('Reject-remote semantics = {}'.format(enum.reject_remote))
    enum.dG_bp = args.dG_bp
    logger.info('Average strength of a toehold base-pair dG_bp = {}'.format(enum.dG_bp))
    if args.ignore_branch_3way:
        if branch_3way in FAST_REACTIONS[1]:
            FAST_REACTIONS[1].remove(branch_3way)
        logger.info('No 3-way branch migration.')
    if args.ignore_branch_4way:
        if branch_4way in FAST_REACTIONS[1]:
            FAST_REACTIONS[1].remove(branch_4way)
        logger.info('No 4-way branch migration.')

    # Set either k-slow or release cutoff.
    if args.k_slow:
        if args.release_cutoff is not None:
            args.release_cutoff = None
            logger.warning('Release-cutoff overwritten by k-slow!')
        if args.release_cutoff_1_1 != args.release_cutoff_1_2:
            logger.warning('Release-cutoff (1,1) overwritten by k-slow!')
            logger.warning('Release-cutoff (1,2) overwritten by k-slow!')
        rc, k_rc = 0, None
        while True:
            rc += 1
            k_rc = opening_rate(rc)
            if k_rc < args.k_slow:
                break
        enum.release_cutoff = rc
        enum.k_slow = args.k_slow
        logger.info('Rate-dependent enumeration: k-slow = {}'.format(enum.k_slow))
        logger.info('  - corresponding release-cutoff: {} < L < {}'.format(rc-1, rc))
    else:
        if args.release_cutoff is not None:
            enum.release_cutoff = args.release_cutoff
            logger.info('Rate-independent enumeration: release cutoff L = {}'.format(
                enum.release_cutoff))
            logger.info('  - corresponding k-slow: {}'.format(opening_rate(enum.release_cutoff)))
        else:
            logger.info("Rate-independent enumeration:")
            enum.release_cutoff_1_1 = args.release_cutoff_1_1
            logger.info("  - release cutoff for reaction arity (1,1) = {}".format(
                enum.release_cutoff_1_1))
            enum.release_cutoff_1_2 = args.release_cutoff_1_2
            logger.info("  - release cutoff for reaction arity (1,2) = {}".format(
                enum.release_cutoff_1_2))
    if args.k_fast:
        enum.k_fast = args.k_fast
        logger.info('Rate-dependent enumeration: k-fast = {}'.format(enum.k_fast))

    # DEBUGGING
    enum.DFS = not args.bfs
    enum.interactive = args.interactive
    enum.interruptible = args.interruptible

    # EXPERIMENTAL
    if args.p_min != 0:
        enum.p_min = args.p_min
        logger.warning('Using experimental option: --p-min = {}'.format(enum.p_min))
    if (args.k_fast or args.k_slow) and args.local_elevation:
        logger.warning('Using experimental option: --local-elevation.')
        enum.local_elevation = True
        if enum.max_helix:
            enum.no_max_helix = False
            logger.warning('Turning off max-helix mode for local-elevation semantics.')
    elif args.local_elevation:
        logger.warning('Local-elevation and rate-independent semantics are incompatible.')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Run reaction enumeration (or not) #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    logger.info("")
    if args.profile:
        try:
            import statprof
        except ImportError as err:
            logger.warning("Python-module statprof not found, disabled Peppercorn profiling.")
            args.profile = False

    if args.dry_run:
        logger.info("Dry run (not enumerating any reactions)... ")
        enum.dry_run()
        logger.info("Done.")
    else:
        logger.info("Enumerating reactions...")
        if args.interactive:
            logger.info("Interactive mode enabled: Fast and slow reactions " + \
                        "will be printed for each complex as enumerated." + \
                        "Press ^C at any time to terminate and write accumulated" + \
                        "complexes to output.")
        if args.profile:
            statprof.start()
            try:
                enum.enumerate()
            finally:
                statprof.stop()
                statprof.display()
        else:
            enum.enumerate()
        logger.info("Done.")

    # ~~~~~~~~~~~~~~~~~~~ #
    # Handle condensation #
    # ~~~~~~~~~~~~~~~~~~~ #
    condensed = args.condensed
    detailed = (not args.condensed or args.detailed)
    if condensed:
        logger.info("Output will be condensed to remove transient complexes.")
        if args.profile:
            statprof.start()
            try:
                enum.condense()
            finally:
                statprof.stop()
                statprof.display()
        else:
            enum.condense()

    if args.sbml:
        if detailed and condensed:
            logger.error("SBML output can be detailed OR condensed, not both.")
        enum.to_sbml(args.sbml, condensed = condensed)

    output = enum.to_pil(args.output_filename, 
                         detailed=detailed, condensed=condensed, composite=composite, 
                         molarity=args.concentration_unit, time = args.time_unit)

    print(output, end='')

if __name__ == '__main__':
   main()

