#!/usr/bin/env python

# TODO: 
#   1) make main function
#   2) write tests
#   3) remove old code
#   4) merge into library

from __future__ import division

import os
import sys
import pandas as pd
import subprocess as sub

from crnsimulator import parse_crn_string

from peppercornenumerator import Enumerator, __version__
from peppercornenumerator.objects import clear_memory
from peppercornenumerator.condense import PepperCondensation
from peppercornenumerator.input import read_pil, read_seesaw, ParseException
from peppercornenumerator.output import write_pil

class MissingDataError(Exception):
    pass

#todo: docstring, libraryweakness, moveit
def seesaw_model(pilstring, is_file = False, enumfile='',
        detailed = True, condensed = False, conc = 'nM', **kwargs):
    """
    # return output format:
    #   1) object
    #   2) file   -> file 
    #      string -> string 
    """
    cxs, rxns = read_seesaw(pilstring, is_file, 
            explicit = False,
            conc=kwargs['seesaw-conc'], 
            reactions = kwargs['seesaw-rxns'])

    enum = Enumerator(cxs.values(), rxns)
    enum.dry_run()

    # NOTE: library interface weakness: cannot access condensation!
    if enumfile :
        with open(enumfile, 'w') as fh :
            write_pil(enum, fh=fh, detailed = detailed, condensed = condensed, 
                    molarity = conc)
            outstring = enumfile
    else :
        outstring = write_pil(enum, fh=None, detailed = detailed, condensed = condensed, 
                molarity = conc)
    
    # do condensation, again!
    if condensed:
        enum = PepperCondensation(enum)
        enum.condense()
    return enum, outstring

def peppercorn(pilstring, is_file = False, enumfile='',
        detailed = True, condensed = False, conc = 'nM', **kwargs):
    """
    # return output format:
    #   1) object
    #   2) file   -> file 
    #      string -> string 
    """
    try:
        cxs, rxns, comp = read_pil(pilstring, is_file, composite = True)
    except ParseException as ex_pil:
        cxs, rxns = read_seesaw(pilstring, is_file, explicit = False)
        comp = []

    init_cplxs = filter(lambda x: x._concentration is None or \
                    float(x._concentration[1]) != (0.0), cxs.values())

    enum = Enumerator(init_cplxs, rxns)

    # set kwargs parameters
    for k, w in kwargs.items():
        if hasattr(enum, k):
            setattr(enum, k, w)
        else:
            raise ValueError('No peppercorn attribute called: {}'.format(k))
        # enum.release_cutoff = release_cutoff
        # enum.max_helix_migration = max_helix
        # enum.remote_migration = not reject_remote
        # enum.max_complex_size = max_complex_size
        # enum.max_complex_count = max_complex_count
        # enum.max_reaction_count = max_reaction_count
        # enum.k_fast = k_fast
        # enum.k_slow = k_slow
        # enum.ddG_bind = ddG_bind

    enum.enumerate()

    # NOTE: library interface weakness: cannot access condensation!
    if enumfile :
        with open(enumfile, 'w') as fh :
            write_pil(enum, fh=fh, detailed = detailed, condensed = condensed, 
                    composite = comp, molarity = conc)
        outstring = enumfile
    else :
        outstring = write_pil(enum, fh=None, detailed = detailed, condensed = condensed, 
                composite = comp, molarity = conc)
    
    # do condensation, again!
    if condensed:
        enum = PepperCondensation(enum)
        enum.condense()
    return enum, outstring

#cleanup
class FigureData(object):
    """Analyze a system using peppercorn and compare with experimental data.

    Examples:
        for a pilfile and arguments:
            recation x has rate y

        for a pilfile and arguments:
            for initial conditions x:
                reporter x and metric y gives tuple (t,c)
    """

    def __init__(self, name):
        self.name = name

        self.tmpdir = 'tmp/'
        self.count = 0

        # Common sets of peppercorn parameters...
        self._pepperargs = {
                'detailed' : {'conc' : 'nM'},
                'Detailed' : {'conc' : 'nM',
                             'max_complex_size' : 10, 
                             'max_complex_count' : 1000,
                             'max_reaction_count' : 5000},
                'DETAILED' : {'conc' : 'nM',
                             'max_complex_size' : 15, 
                             'max_complex_count' : 10000,
                             'max_reaction_count' : 50000},
                'condensed': {'conc' : 'nM',
                              'condensed': True},
                'Condensed': {'conc' : 'nM',
                              'condensed': True,
                              'max_complex_size' : 10, 
                              'max_complex_count' : 1000,
                              'max_reaction_count' : 5000},
                'CONDENSED': {'conc' : 'nM',
                              'condensed': True,
                              'max_complex_size' : 15, 
                              'max_complex_count' : 10000,
                              'max_reaction_count' : 50000},
                }

        self._pil_to_file = {}

        self._ratemode = False
        self._ratedata = []
        self._ratecalc = []

        self._simmode = False
        self._simdata = []
        self._simcalc = []

    @property
    def pepperargs(self):
        return self._pepperargs

    @pepperargs.setter
    def pepperargs(self, value):
        self._pepperargs.update(value)

    def make_pilfile(self, pilstring, suffix):
        name = '{}{}-{:02d}'.format(self.tmpdir, self.name, self.count)
        self.count += 1

        with open(name + suffix, 'w') as pil:
            pil.write(pilstring)

        return name 

    def get_dataframes(self):
        if self._ratemode:
            return self.get_reaction_dataframes()
        elif self._simmode:
            return self.get_system_dataframes()

    def eval(self, pepperargs='default', cmpfig=False, verbose=0):
        if self._ratemode:
            return self.eval_reactions(pepperargs, verbose = verbose)
        elif self._simmode:
            return self.eval_system(pepperargs, cmpfig, verbose = verbose)
        else:
            raise MissingDataError('Cannot find data for evaluation')

    # Reactionrate mode
    def canon_rxn(self, rxn):
        return '{} -> {}'.format(' + '.join(sorted(rxn[0])), ' + '.join(sorted(rxn[1])))

    def add_reaction_rate_setup(self, pilstring, reaction):
        # Supports only single reactions
        rxn, _ = parse_crn_string(reaction)
        reaction = self.canon_rxn(rxn[0])
        rate = float(rxn[0][2][0])

        pilname = self.make_pilfile(pilstring, suffix='-input.pil')

        self._ratedata.append((pilname, pilstring, reaction, rate))
        self._ratemode = True

    def eval_reactions(self, pepperargs='default', verbose=0):
        if not self._ratemode:
            raise MissingDataError('Cannot find ratedata for evaluation')

        if pepperargs in self._pepperargs:
            pargs = self._pepperargs[pepperargs]
        else:
            raise MissingDataError('Cannot find key "{}" in pepperargs'.format(pepperargs))

        ratecalc = []
        for name, pil, rxn, rate in self._ratedata:
            clear_memory()

            pname = name + '-input.pil'
            ename = name + '-enum.pil'
            enumOBJ, _ = peppercorn(pname, is_file=True, enumfile=ename, **pargs)
            rxns = enumOBJ.reactions #NOTE: requires new release

            prate = None
            for r in rxns:
                ed = map(str, r.reactants)
                pr = map(str, r.products)
                if self.canon_rxn([ed, pr]) == rxn:
                    prate = r.rate
                    break
            if prate is None:
                raise MissingDataError('Target reaction not found: {} not in {}'.format(rxn, name))
            ratecalc.append((pepperargs, prate))
        self._ratecalc.append((ratecalc))

    def get_reaction_dataframes(self):
        name, pil, rxn, exp = zip(*self._ratedata)

        ratecalc = [None] if not self._ratecalc else self._ratecalc

        for rcalc in ratecalc:
            pargs, calc = (None, None) if rcalc is None else zip(*rcalc)
            df = pd.DataFrame(data={
                'pilfile': name,
                'reaction': rxn,
                'literature': exp,
                'calculated': calc,
                'pepperargs': pargs},
                columns=['pilfile', 'reaction', 'literature', 'pepperargs', 'calculated'])
            yield df

    # Simulation mode
    def add_system_simulation_setup(self, pilstring, simulation, reporter, metric, (time, conc)):
        if pilstring in self._pil_to_file:
            pilname = self._pil_to_file[pilstring]
        else:
            pilname = self.make_pilfile(pilstring, suffix='-input.pil')
            self._pil_to_file[pilstring] = pilname
        self._simdata.append((pilname, pilstring, simulation, reporter, metric, time, conc))
        self._simmode = True

    def eval_system(self, pepperargs='default', cmpfig=False, verbose=0):
        if pepperargs in self._pepperargs:
            pargs = self._pepperargs[pepperargs]
            condensed = pargs.pop('condensed', False)
        else:
            raise MissingDataError('Cannot find key "{}" in pepperargs'.format(pepperargs))

        elast = None
        slast = None
        simcalc = []
        trajectories = []
        for name, pil, sim, reporter, metric, time, conc in self._simdata:
            clear_memory()
            if verbose:
                print("{}: {}".format(name, pepperargs))
                if verbose > 1:
                    print("  {}\n  {}".format(sim, pargs))

            idx = sim.find('--p0')

            pname = '{}-{}.pil'.format(name, 'input')
            ename = '{}-{}-{}.pil'.format(name, pepperargs, 'enum')

            if ename != elast :
                if 'seesaw-rxns' in pargs:
                    _ , ename = seesaw_model(pname, is_file=True, enumfile=ename, **pargs)
                else:
                    _ , ename = peppercorn(pname, is_file=True, enumfile=ename, 
                        detailed = (not condensed), condensed = condensed, **pargs)
                elast = ename

            sname = '{}-{}-{}'.format(name, pepperargs, 'simu')
            nxy = simulate_pil(ename, sname, sim.split(), force=(sname != slast))
            slast = sname

            metric = metric.split(':')

            if metric[0] == 'completion-time':
                sim_time = nxy_get_time_at_conc(nxy, reporter, conc)
                sim_conc = conc
                simcalc.append((pepperargs, sim_time, sim_conc))
            elif metric[0] == 'diagonal-crossing-time':
                tmax = float(metric[1])
                cmax = float(metric[2])
                sim_time, sim_conc = nxy_get_diagonal_points(nxy, reporter, tmax, cmax)
                simcalc.append((pepperargs, sim_time, sim_conc))
            else:
                raise NotImplementedError('Metric "{}" not supported.'.format(metric))

            if cmpfig :
                tr = nxy_get_trajectory(nxy, reporter)
                if tr is None: continue

                if trajectories :
                    for te, [time, val] in enumerate(tr):
                        if len(trajectories) > te:
                            assert trajectories[te][0] == time
                            trajectories[te].append(val)
                else :
                    trajectories = tr

        if cmpfig: 
            with open(sname + '_cmp.nxy', 'w') as cnxy:
                for td in trajectories:
                    cnxy.write("{}\n".format(' '.join(td)))

        self._simcalc.append((simcalc))
        return self.get_system_dataframe(simcalc)

    def get_system_dataframe(self, simcalc):
        name, pil, sim, rep, met, time, conc = zip(*self._simdata)
        pargs, stime, sconc = zip(*simcalc)

        # Only print initial conditions in the dataframe
        ini = map(lambda x: x[x.find('--p0'):], sim)

        df = pd.DataFrame(data={
            'pilfile': name,
            'simulation': ini,
            'reporter': rep,
            'metric': met,
            'exp-conc': conc,
            'exp-time': time,
            'sim-time': stime,
            'sim-conc': sconc,
            'pepperargs': pargs},
            columns=['pilfile', 'simulation', 'reporter', 'metric', 'exp-conc',
                'exp-time', 'pepperargs', 'sim-conc', 'sim-time'])
        return df

    def get_system_dataframes(self):
        name, pil, sim, rep, met, time, conc = zip(*self._simdata)

        simcalc = [None] if not self._simcalc else self._simcalc
        
        # Only print initial conditions in the dataframe
        ini = map(lambda x: x[x.find('--p0'):], sim)

        for scalc in simcalc:
            pargs, stime, sconc = (None, None, None) if scalc is None else zip(*scalc)
            df = pd.DataFrame(data={
                'pilfile': name,
                'simulation': ini,
                'reporter': rep,
                'metric': met,
                'exp-conc': conc,
                'exp-time': time,
                'sim-time': stime,
                'sim-conc': sconc,
                'pepperargs': pargs},
                columns=['pilfile', 'simulation', 'reporter', 'metric', 'exp-conc',
                    'exp-time', 'pepperargs', 'sim-conc', 'sim-time'])
            yield df

#done
def simulate_pil(pilfile, output, command, force = False, verbose = 0):
    """ Simulate the PIL file and return the filename storing the simulation. 

    Args:
        pilfile (str): A peppercorn output file. It is currently required to
            have the ending "_enum.pil", to emphasize that it is an enumerated
            system in pil format.
        output (str): The basename of the output files.
        command (list): The command for the pilsimulator.
        force (bool, optional): Overwrite existing files from the pilsimulator.

    Returns:
        (str): The name of a file storing the simulation in *.nxy format.
    """
    simname = output + '.py'
    if force or not os.path.exists(simname):
        simu = command + ['-o', simname, '--force']
    else :
        assert command[0] == 'pilsimulator'
        if '--no-jacobian' in command: 
            # make sure it is at the beginning...
            assert command[1] == '--no-jacobian'
            simu = ['python', simname] + command[2:]
        else :
            simu = ['python', simname] + command[1:]

    if verbose:
        print("{}".format(' '.join(simu)))

    with open(pilfile, 'r') as crn, \
            open(output + '.nxy', 'w') as nxy, \
            open(output + '.err', 'w') as err:
        proc = sub.Popen(simu, stdin=crn, stdout=nxy, stderr=err)
        proc.communicate(None)
        if proc.returncode:
            raise Exception(err)
    return output + '.nxy'

#done
def nxy_get_trajectory(nxyfile, species):
    """ Extract a specific trajectory from *.nxy file.

    Args:
        nxyfile (str): Path to a *.nxy file.
        species (str): The species of interest.

    Returns:
        (list[[str,str]]): A trajectory as lol: simulation-time & concentration.
    """
    trajectory = []
    with open(nxyfile, 'r') as nxy:
        idx = None
        for l in nxy.readlines():
            if l[0:25] == '# Initial concentrations:' : 
                data = eval(l.strip().split(': ')[1])
                for e, (sp,ini) in enumerate(data,1):
                    if sp == species:
                        idx = e
                continue
            elif l[0] == '#' : 
                continue

            if idx is None:
                raise MissingDataError('Could not find species {} in {}'.format(species, nxyfile))

            d = l.strip().split()
            trajectory.append([d[0], d[idx]])
    return trajectory

#done
def nxy_get_conc_at_time(nxyfile, species, time):
    """ Returns the concentration of a species at (or after) a specific time.
    """
    with open(nxyfile, 'r') as nxy:
        idx = None
        for l in nxy.readlines():
            if l[0:25] == '# Initial concentrations:' : 
                data = eval(l.strip().split(': ')[1])
                for e, (sp,ini) in enumerate(data,1):
                    if sp == species:
                        idx = e
                continue
            elif l[0] == '#' : 
                continue

            if idx is None:
                raise MissingDataError('Could not find species {} in {}'.format(species, nxyfile))

            d = map(float, l.strip().split())
            if d[0] >= time:
                return d[idx]
    return None

#done
def nxy_get_time_at_conc(nxyfile, species, threshold):
    """Returns first timepoint at which a concentration threshold was reached.

    Args:
        nxyfile (str): Path to a *.nxy file.
        species (str): The species of interest.
        threshold (flt): The concentration threshold (positive or negative). 
            If positive, returns the first time point with higher concentration.
            If negative, returns the first time point with lower concentration.

    Returns:
        (float): time 
    """
    [sign, th] = [-1, -threshold] if threshold < 0 else [1, threshold]

    with open(nxyfile, 'r') as nxy:
        idx = None
        for l in nxy.readlines():
            if l[0:25] == '# Initial concentrations:' : 
                data = eval(l.strip().split(': ')[1])
                for e, (sp,ini) in enumerate(data,1):
                    if sp == species:
                        idx = e
                continue
            elif l[0] == '#' : 
                continue

            if idx is None:
                raise MissingDataError('Could not find species {} in {}'.format(species, nxyfile))

            d = map(float, l.strip().split())
            if sign == 1 and d[idx] > th:
                return d[0]
            elif sign == -1 and d[idx] < th:
                return d[0]
    return None

#done
def nxy_get_diagonal_points(nxyfile, species, tmax, cmax):
    """Return time and concentration when a trajectory crosses a diagonal line.
    
    Args:
        nxyfile (str): Path to a *.nxy file.
        species (str): The species of interest.
        tmax (flt): The time of the experiment.
        cmax (flt): The maximum concentration.

    Returns:
        (float, float): time, concentration
    """

    with open(nxyfile, 'r') as nxy:
        idx = None
        for l in nxy.readlines():
            if l[0:25] == '# Initial concentrations:' : 
                data = eval(l.strip().split(': ')[1])
                for e, (sp,ini) in enumerate(data,1):
                    if sp == species:
                        idx = e
                continue
            elif l[0] == '#' : 
                continue

            if idx is None:
                raise MissingDataError('Could not find species {} in {}'.format(species, nxyfile))

            d = map(float, l.strip().split())

            if d[idx] >= cmax * (1 - d[0]/tmax):
                return d[0], d[idx]
    return None

def main():
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Import data from literature #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    from zhang2007  import data as z07
    from yin2008    import data as y08
    from zhang2009  import data as z09 # missing reaction rate data & sequence variation
    from zhang2010  import data as z10 # there is more!
    from genot2011  import data as g11
    from qian2011   import data as q11 
    from qian2011sqrt import data as q11sq # missing data
    from zhang2011  import data as z11 # Only Figure 2, and it is a bit strange...
    from kotani2017 import data as k17
    
    from zhang2009_rates import data as z09r # 3-way branch migration rates
    from dabby2013_rates import data as d13r # 4-way branch migration rates

    analysis = z07() + y08() + z09() + z10() + z11() + g11() + q11() + k17()
    #analysis = z09r() + d13r()

    for fig in analysis:
        print("\n{}:".format(fig.name))
        #fig.eval('default', verbose=1)
        for df in fig.get_dataframes():
            print(df)

    # There is more ...:
    # seelig2006.py
    # groves2015.py
    # srinivas2017.py
    # sun2018.py
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Choose data for analysis #
    # ~~~~~~~~~~~~~~~~~~~~~~~~ #
    # feedforward: z07_F3, q11, k17_F3
    # Catalysts : z07_F1, z09_F5, k17_F2
    # Auto-catalysts : z07_F4, y08_F3, k17_F4
    # missing: z10, z11c, g11

#deprecated
def analyze_rate_data(data, verbose = 0):
    """Calculate the 
    """
    
    kresults = [];
    for pn, pepperargs in data['pepperargs']:
        name = data['name']
        tmpname = 'tmp/'+data['name']
        if pn: 
            name += '-' + pn
            tmpname += '-' + pn
    
        for e, pilstring in enumerate(map(data['piltemplate'], data['pilparams'])):
            clear_memory()
            ttmpname = tmpname + '-' + str(e) if len(data['pilparams']) > 1 else tmpname

            if verbose:
                print('Analyzing', ttmpname)
    
            enumOBJ, crn = old_peppercorn(pilstring, name=ttmpname, **pepperargs)
    
            assert 'rates' in data
            if len(enumOBJ.condensed_reactions) == 1:
                rxn = enumOBJ.condensed_reactions[0]
                er = data['exp_results'][e]
                kresults.append([er, rxn.rate])
            else :
                if 'reactants' in data['rates']:
                    for rxn in enumOBJ.condensed_reactions:
                        if sorted([rs.name for rs in rxn.reactants]) \
                                == sorted(data['rates']['reactants']):
                            er = data['exp_results'][e]
                            kresults.append([er, rxn.rate])
                            break
                else:
                    raise NotImplementedError('multiple condensed reactions')
    return kresults

#deprecated
def analyze_system_data(data, cmpfig=False, verbose = 0):
    tresults = [];
    cresults = [];
    for pn, pepperargs in data['pepperargs']:
        name = data['name']
        tmpname = 'tmp/'+data['name']
        if pn: 
            name += '-' + pn
            tmpname += '-' + pn
    
        for e, pilstring in enumerate(map(data['piltemplate'], data['pilparams'])):
            clear_memory()

            ttmpname = tmpname + '-' + str(e) if len(data['pilparams']) > 1 else tmpname
            if verbose:
                print('Analyzing', ttmpname)
    
            enumOBJ, crn = old_peppercorn(pilstring, name=ttmpname, **pepperargs)
    
            assert 'simulation' in data
            for i, command in enumerate(data['simulation']):
                nxy = simulate_pil(crn, tmpname, command, force=(i==0))

                rep = data['reporter']

                if data['metric'] == 'half-completion-time':
                    t_half = data['exp-results']

                    # e = number of pilstring modifications
                    # i = number of simulations per modification
                    (et, ec) = t_half[e+i]

                    time = nxy_get_time_at_conc(nxy, rep, ec)
                    tresults.append([et, time])
                    if verbose:
                        print 'half-completion time', et, time

                elif data['metric'] == 'diagonal-crossing-time':
                    tmax = data['tmax']
                    cmax = data['cmax']
                    diag = data['exp-results']

                    # e = number of pilstring modifications
                    # i = number of simulations per modification
                    (et, ec) = diag[e+i]

                    time, conc = nxy_get_diagonal_points(nxy, rep, tmax, cmax)
                    if et < tmax/2:
                        tresults.append([et, time])
                        cresults.append([ec, conc])
                    else :
                        tresults.append([et, time])
                        cresults.append([ec, conc])
                    if verbose:
                        print 'diagonal-crossing time', et, time
                        print 'diagonal-crossing concentration', ec, conc

                elif data['metric'] == 'quadrant-crossing':
                    tmax = data['tmax']
                    cmax = data['cmax']
                    t_half = data['exp-results'][0]
                    c_half = data['exp-results'][1]
                    quad = t_half + c_half

                    # e = number of pilstring modifications
                    # i = number of simulations per modification
                    (et, ec) = quad[e+i]

                    if e+i < len(t_half):
                        time = nxy_get_time_at_conc(nxy, rep, ec)
                        tresults.append([et, time])
                        if verbose: print 'half-completion time', et, time
                    else :
                        conc = nxy_get_conc_at_time(nxy, rep, et)
                        cresults.append([ec, conc])
                        if verbose: print 'mid-experiment concentration', ec, conc
 
                else :
                    m = data['metric']
                    raise NotImplementedError('metric for comparison not found:', m)

                if cmpfig :
                    tr = nxy_get_trajectory(nxy, rep)
                    if tr is None: continue

                    if trajectories :
                        for te, [time, val] in enumerate(tr):
                            if len(trajectories) > te:
                                assert trajectories[te][0] == time
                                trajectories[te].append(val)
                    else :
                        trajectories = tr

            if cmpfig: 
                with open(tmpname + '_cmp.nxy', 'w') as cnxy:
                    for td in trajectories:
                        cnxy.write("{}\n".format(' '.join(td)))

    return tresults, cresults

#deprecated
def old_peppercorn(kernelstring, name, 
        condensed = False, 
        conc = 'nM', 
        seesawrxns='',
        release_cutoff = 8,
        max_complex_size = 10,
        max_complex_count = 10000,
        max_reaction_count = 10000,
        reject_remote = False,
        max_helix = True,
        ddG_bind = None,
        k_fast = None,
        k_slow = None):
    """ A wrapper for peppercorn.  """


    with open(name + '-input.pil', 'w') as pil:
        pil.write(kernelstring)

    try:
        complexes, reactions = read_pil(kernelstring)
    except ParseException, ex_pil:
        complexes, reactions = read_seesaw(kernelstring, explicit=False,
                reactions=seesawrxns)

    init_cplxs = filter(lambda x: x._concentration is None or \
                    float(x._concentration[1]) != (0.0), complexes.values())

    enum = Enumerator(init_cplxs, reactions)
    
    if seesawrxns:
        print "WARNING: dry run!!!"
        enum.dry_run()
        condensed = False
    else :
        enum.release_cutoff = release_cutoff
        enum.max_helix_migration = max_helix
        enum.remote_migration = not reject_remote
        enum.max_complex_size = max_complex_size
        enum.max_complex_count = max_complex_count
        enum.max_reaction_count = max_reaction_count
        if k_fast is not None:
            enum.k_fast = k_fast
        if k_slow is not None:
            enum.release_cutoff = 15
            enum.k_slow = k_slow
        if ddG_bind :
            enum.ddG_bind = ddG_bind
        enum.enumerate()
        #condensed = False

    if condensed:
        enumRG = PepperCondensation(enum)
        enumRG.condense()

    detailed = not condensed

    enumfile = name + '-enum.pil'

    with open(enumfile, 'w') as crn:
        # You can find this file in the tmp directory
        write_pil(enum, crn, detailed, condensed, molarity=conc)

    if condensed:
        return enumRG, enumfile
    else :
        return enum, enumfile

if __name__ == '__main__':
    main()

