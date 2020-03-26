#!/usr/bin/env python3
#
# pilsimulator: simulate PIL files using ODEs
#
# Written by Stefan Badelt (badelt@caltech.edu).
#
# Use at your own risk.
#
#
import logging
import os
import sys
import argparse

from peppercornenumerator import __version__
from peppercornenumerator.input import load_pil_crn
from crnsimulator import ReactionGraph, get_integrator
from crnsimulator.odelib_template import add_integrator_args

class SimulationSetupError(Exception):
    pass

def main():
    """Translate a CRN into a system of ODEs. Optional: Simulate ODEs on-the-fly. """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action='count', default = 0,
            help = "Print logging output. (-vv increases verbosity.)")
    parser.add_argument('--logfile', default = '', action = 'store', metavar = '<str>',
            help = "Redirect logging information to a file.")
    parser.add_argument("--force", action='store_true',
            help="Overwrite existing files")
    parser.add_argument("--dry-run", action='store_true',
            help="Do not run the simulation, only write the files.")
    parser.add_argument("-o", "--output", default='odesystem', metavar='<str>',
            help="Name of ODE library files.")
    parser.add_argument("--no-jacobian", action='store_true',
            help=argparse.SUPPRESS)
    parser.add_argument("--jacobian", action='store_true',
            help="""Symbolic calculation of Jacobi-Matrix. 
            This may generate a very large simulation file.""")
    add_integrator_args(parser)
    args = parser.parse_args()

    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    logger = logging.getLogger('crnsimulator')
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(args.logfile) if args.logfile else logging.StreamHandler()
    if args.verbose == 0:
        handler.setLevel(logging.WARNING)
    elif args.verbose == 1:
        handler.setLevel(logging.INFO)
    elif args.verbose == 2:
        handler.setLevel(logging.DEBUG)
    elif args.verbose >= 3:
        handler.setLevel(logging.NOTSET)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if args.no_jacobian:
        logging.warning('Deprecated argument: --no-jacobian, OFF by default.', DeprecationWarning)
    if args.pyplot_labels:
        warning.warning('Deprecated argument: --pyplot_labels, use --labels.', DeprecationWarning)
        args.labels = args.pyplot_labels
        args.pyplot_labels = None

    # ********************* #
    # ARGUMENT PROCESSING 1 #
    # ..................... #
    filename = args.output + \
        '.py' if args.output[-3:] != '.py' else args.output
    odename = 'odesystem'

    input_crn = sys.stdin.readlines()
    input_crn = "".join(input_crn)

    crn, species = load_pil_crn(input_crn)

    if crn == []:
        raise SystemExit('No Input CRN provided.')

    V = [] # sorted species (vertices) vector
    C = [] # corresponding concentration vector
    seen = set() # keep track of what species are covered

    # Move interesting species to the front, in their given order.
    const = []
    labels = args.pyplot_labels
    for s in labels:
        if s in seen :
            raise SimulationSetupError('Multiple occurances of {} in labels.'.format(s))
        V.append(s)
        const.append(False if species[s][0][0] == 'i' else True)
        C.append(species[s][1])
        seen.add(s)

    # Append the remaining specified species 
    for s in sorted(species):
        if s in seen : continue
        V.append(s)
        const.append(False if species[s][0][0] == 'i' else True)
        C.append(species[s][1])
        seen.add(s)

    # Split CRN into irreversible reactions
    new = []
    for [r, p, k] in crn:
        assert not (None in k)
        if len(k) == 2:
            new.append([r, p, k[0]])
            new.append([p, r, k[1]])
        else:
            new.append([r, p, k[0]])
    crn = new

    # **************** #
    # WRITE ODE SYSTEM #
    # ................ #
    if filename != 'odesystem.py' and not args.force and os.path.exists(filename):
        logger.warning(f'Reading ODE system from existing file: {filename}')
    else:
        # ******************* #
        # BUILD REACTIONGRAPH #
        # ................... #
        RG = ReactionGraph(crn)
        if len(RG.species) > len(V):
            raise SimulationSetupError('Undefined species in the CRN.')
        elif len(RG.species) < len(V):
            raise SimulationSetupError('Some species do not appear in the CRN.')

        # ********************* #
        # PRINT ODE TO TEMPLATE #
        # ..................... #
        filename, odename = RG.write_ODE_lib(sorted_vars = V, concvect = C,
                                             const = const if any(const) else None,
                                             jacobian = not args.no_jacobian, 
                                             filename = filename,
                                             odename = odename)
        logger.info(f'CRN to ODE translation successful. Wrote file: {filename}')

    # ******************* #
    # SIMULATE ODE SYSTEM #
    # ................... #
    if args.dry_run:
        logger.info('Dry run: Simulate the ODE system using:')
        logger.info(f"  python {filename} --help ")
    else:
        logger.info('Simulating the ODE system, change parameters using:')
        logger.info(f"  python {filename} --help ")

        integrate = get_integrator(filename)

        # ********************* #
        # ARGUMENT PROCESSING 2 #
        # ..................... #
        integrate(args, setlogger = True)

    return


if __name__ == '__main__':
    main()
