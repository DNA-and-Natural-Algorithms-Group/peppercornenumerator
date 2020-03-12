#!/usr/bin/env python3
#
# pilsimulator: simulate PIL files using ODEs
#
# Written by Stefan Badelt (badelt@caltech.edu).
#
# Use at your own risk.
#
#
from __future__ import absolute_import, division, print_function

import os
import sys
import argparse

from crnsimulator import ReactionGraph
from crnsimulator.odelib_template import add_integrator_args
from peppercornenumerator.input import load_pil_crn
from peppercornenumerator import __version__


def get_integrator(filename, function = 'integrate'):
    """Workaround to avoid deprecation warnings for the imp module. 

    Switch to crnsimulator.get_integrator as soon as we drop Python 2.7 support.
    """
    try: # Python 3.7
        import types
        import importlib.machinery

        loader = importlib.machinery.SourceFileLoader('myloader', filename)
        mod = types.ModuleType(loader.name)
        loader.exec_module(mod)

        return getattr(mod, function)
    except ImportError as err: # Python 2.7
        import imp
        _temp = imp.load_source('mysource', filename) # used to b odename instaed of 'mysource'
        return getattr(_temp, function)

class SimulationSetupError(Exception):
    pass

def main():
    """Translate a CRN into a system of ODEs. Optional: Simulate ODEs on-the-fly. """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pilsim = parser.add_argument_group('pilsimulator I/O')

    pilsim.add_argument('--version', action='version', 
            version='%(prog)s ' + __version__)
    pilsim.add_argument("--dry-run", action='store_true',
            help="Do not run the simulation, only write the files.")
    pilsim.add_argument("-o", "--output", default='odesystem', metavar='<str>',
            help="Name of ODE library files.")
    pilsim.add_argument("--force", action='store_true',
            help="Overwrite existing files")
    pilsim.add_argument("--no-jacobian", action='store_true',
            help="Do not write the Jacobi-Matrix. Let the solver do the work.")

    add_integrator_args(parser)

    args = parser.parse_args()


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
    labels = args.pyplot_labels
    for s in labels:
        if s in seen :
            raise SimulationSetupError('Multiple occurances of {} in labels.'.format(s))
        V.append(s)

        if species[s][0][0] != 'i': raise NotImplementedError
        C.append(species[s][1])
        seen.add(s)

    # Append the remaining specified species 
    for s in sorted(species):
        if s in seen : continue
        V.append(s)
        if species[s][0][0] != 'i':
            raise NotImplementedError('Concentrations must be given as "initial" concentrations.')
        C.append(species[s][1])
        seen.add(s)

    # Split CRN into irreversible reactions
    new = []
    for [r, p, k] in crn:
        if None in k:
            print('# Set missing rates to 1.')
            k[:] = [x if x is not None else 1 for x in k]

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
        print('# Reading ODE system from existing file:', filename)
    else:
        # ******************* #
        # BUILD REACTIONGRAPH #
        # ................... #
        RG = ReactionGraph(crn)
        if len(RG.species) > len(V):
            raise SimulationSetupError('Undefined species in the CRN')
        elif len(RG.species) < len(V):
            raise SimulationSetupError('Some species do not appear in the CRN:')

        # ********************* #
        # PRINT ODE TO TEMPLATE #
        # ..................... #
        filename, odename = RG.write_ODE_lib(sorted_vars = V, concvect = C,
                                             jacobian = not args.no_jacobian, 
                                             filename = filename,
                                             odename = odename)
        print('# Wrote ODE system:', filename)

    # ******************* #
    # SIMULATE ODE SYSTEM #
    # ................... #
    if args.dry_run:
        print('# Dry-run: Simulate the ODE system using:')
        print("#  python {} --help ".format(filename))
    else:
        print('# Simulating the ODE system. Change parameters using:')
        print("#  python {} --help ".format(filename))

        integrate = get_integrator(filename)

        # ********************* #
        # ARGUMENT PROCESSING 2 #
        # ..................... #
        integrate(args)

    return


if __name__ == '__main__':
    main()
