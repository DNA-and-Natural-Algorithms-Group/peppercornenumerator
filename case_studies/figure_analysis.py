#!/usr/bin/env python

from __future__ import division, absolute_import, print_function

import os
import numpy as np
import pandas as pd
from subprocess import Popen
from crnsimulator import parse_crn_string
from pyparsing import ParseException

from peppercornenumerator import Enumerator, __version__
from peppercornenumerator.enumerator import enumerate_pil, enumerate_ssw
from peppercornenumerator.objects import clear_memory

class MissingDataError(Exception):
    pass

class FigureData(object):
    """ Produce DataFrames that compare Peppercorn's model with experimental data.

    Examples:
        for a pilfile and arguments:
            recation x has rate y

        for a pilfile and arguments:
            for initial conditions x:
                reporter x and metric y gives tuple (t,c)
    """

    def __init__(self, name, tmpdir = 'tmp/', force = True):
        self.name = name
        self.fname = name.replace(' ', '_').replace('\n', '_').replace('(', '').replace(')', '')
        self.tmpdir = tmpdir
        self.count = 0 # Assign different pil settings

        # Internal mapping of pilstring to a present file containing that string.
        self._pil_to_file = {}

        # TODO: overwrite files
        # self.force = force

        # Mapping of keyword to enumeration parameters
        self._pepperargs = {
                'detailed' : {'enumconc': 'nM'},
                'Detailed' : {'enumconc': 'nM',
                             'max_complex_size': 10, 
                             'max_complex_count': 1000,
                             'max_reaction_count': 5000},
                'DETAILED' : {'enumconc': 'nM',
                             'max_complex_size': 15, 
                             'max_complex_count': 10000,
                             'max_reaction_count': 50000},
                'condensed': {'enumconc': 'nM',
                              'condensed': True},
                'Condensed': {'enumconc': 'nM',
                              'condensed': True,
                              'max_complex_size': 10, 
                              'max_complex_count': 1000,
                              'max_reaction_count': 5000},
                'CONDENSED': {'enumconc': 'nM',
                              'condensed': True,
                              'max_complex_size': 15, 
                              'max_complex_count': 10000,
                              'max_reaction_count': 50000},
                }

        # Mapping of keyword to simulation parameters
        self._simargs = {}
        
        self._ratemode = False
        self._ratedata = []
        self._ratecalc = {}

        self._simmode = False # Internal flag to distingush mode
        self._simdata = []    # userprovided setup (except enumargs)
        self._simcalc = {}    # simcalc[pepperargs] = (time, conc)
        self.cmpfig  = {}

        self._enumerated = set() # self._enumerated[(input,enumargs)] = enumfile
        self._simulated = set() # 
        self._simexecs = set() # 

    def delete_tempfiles(self, pretend = False, verbose = True):
        for ef in self._enumerated:
            if verbose or pretend:
                print('removing', ef)
            if not pretend: os.remove(ef)
        if not pretend:
            self._enumerated = set()

    @property
    def pepperargs(self):
        return self._pepperargs

    @pepperargs.setter
    def pepperargs(self, value):
        self._pepperargs.update(value)

    def make_pilfile(self, pilstring, suffix):
        name = '{}{}-{:02d}'.format(self.tmpdir, self.fname, self.count)
        self.count += 1

        with open(name + suffix, 'w') as pil:
            pil.write(pilstring)

        return name 

    def get_dataframes(self):
        if self._ratemode:
            return self.get_reaction_dataframes()
        elif self._simmode:
            return self.get_system_dataframes()

    def eval(self, pepperargs = 'default', cmpfig = False, verbose = 0, enumprofile = False):
        if self._ratemode:
            return self.eval_reactions(pepperargs, verbose = verbose)
        elif self._simmode:
            return self.eval_system(pepperargs, cmpfig, verbose = verbose, enumprofile = enumprofile)
        else:
            raise MissingDataError('Cannot find data for evaluation')

    # Reactionrate mode
    def canon_rxn(self, rxn):
        return '{} -> {}'.format(' + '.join(sorted(rxn[0])), ' + '.join(sorted(rxn[1])))

    def add_reaction_rate_setup(self, pilstring, reaction, pilid = None):
        # Supports only single reactions
        rxn, _ = parse_crn_string(reaction)
        reaction = self.canon_rxn(rxn[0])
        rate = float(rxn[0][2][0])

        pilname = self.make_pilfile(pilstring, suffix='-input.pil')

        self._ratedata.append((pilname, pilstring, reaction, rate, pilid))
        self._ratemode = True
        assert not self._simmode

    def eval_reactions(self, pepperargs = 'default', verbose = 0):
        if not self._ratemode:
            raise MissingDataError('Cannot find ratedata for evaluation')

        if pepperargs in self._pepperargs:
            pargs = self._pepperargs[pepperargs]
        else:
            raise MissingDataError('Cannot find key "{}" in pepperargs'.format(pepperargs))

        ratecalc = []
        for name, pil, rxn, rate, pilid in self._ratedata:
            clear_memory()

            pname = name + '-input.pil'
            ename = name + '-enum.pil'
            enumOBJ, _ = enumerate_pil(pname, is_file = True, enumfile = ename, **pargs)
            rxns = enumOBJ.condensed_reactions \
                    if pargs.get('condensed', False) else enumOBJ.detailed_reactions

            prate = None
            for r in rxns:
                ed = list(map(str, r.reactants))
                pr = list(map(str, r.products))
                if self.canon_rxn([ed, pr]) == rxn:
                    prate = r.const
                    break
            if prate is None:
                raise MissingDataError(
                        'Target reaction not found: {} not in {}'.format(rxn, name))
            ratecalc.append(prate)
        self._ratecalc[pepperargs] = ratecalc

    def get_reaction_dataframes(self):
        name, pil, rxn, exp, pilid = zip(*self._ratedata)

        ratecalc = {None: None} if not self._ratecalc else self._ratecalc

        for pargs, rcalc in ratecalc.items():
            df = pd.DataFrame(data={
                'Input Filename': name,
                '(n, m)': pilid,
                'Reaction': rxn,
                'Rate (experiment)': exp,
                'Rate (calculated)': rcalc,
                'Semantics': pargs},
                columns=['Input Filename', '(n, m)', 'Reaction', 'Semantics', 
                    'Rate (calculated)', 'Rate (experiment)'])
            yield df

    # Simulation mode
    def add_system_simulation_setup(self, 
                                    pilstring, simulation, 
                                    reporter, metric, minfo, timeandconc,
                                    simargs = ''):
        (time, conc) = timeandconc
        if pilstring in self._pil_to_file:
            pilname = self._pil_to_file[pilstring]
        else:
            pilname = self.make_pilfile(pilstring, suffix='-input.pil')
            self._pil_to_file[pilstring] = pilname
        if simargs == '':
            simargs = simulation[simulation.find('p0')+3:].replace(' ','_')
        elif simargs == 'pilname':
            simargs = pilname
        self._simargs[simargs]=simulation
        self._simdata.append((pilname, pilstring, simargs, reporter, metric, minfo, time, conc))
        self._simmode = True
        assert not self._ratemode

    def eval_system(self, pepperargs = 'default', cmpfig = False, verbose = 0, enumprofile = False):
        if pepperargs in self._pepperargs:
            pargs = self._pepperargs[pepperargs].copy()
            condensed = pargs.pop('condensed', False)
        else:
            raise MissingDataError('Cannot find key "{}" in pepperargs'.format(pepperargs))

        simcalc = []
        trajectories = None
        for name, pilstring, simargs, reporter, metric, minfo, time, conc in self._simdata:
            clear_memory()

            if verbose:
                print("Evaluating {}: {}".format(name, pepperargs))

            pname = '{}-{}.pil'.format(name, 'input')
            ename = '{}-{}-{}.pil'.format(name, pepperargs, 'enum')
            sexec = '{}-{}-{}.py'.format(name, pepperargs, 'simu')
            cname = '{}-{}-{}.nxy'.format(name, pepperargs, 'cmp')
            if '/' in simargs : #NOTE hack
                sname = '{}-{}-{}'.format(name, pepperargs, 'simu')
            else :
                sname = '{}-{}-{}-{}'.format(name, pepperargs, simargs, 'simu')

            # First, enumerate
            if ename not in self._enumerated or enumprofile:
                if verbose:
                    print("Enumerating ... ")
                try:
                    enumerate_pil(pname, is_file = True, enumfile = ename, 
                        detailed = (not condensed), condensed = condensed, **pargs)
                except ParseException as err:
                    enumerate_ssw(pname, is_file = True, enumfile = ename, 
                        detailed = (not condensed), condensed = condensed, **pargs)
                self._enumerated.add(ename)
                if enumprofile:
                    return None

            # Second, simulate 
            if sname not in self._simulated:
                if verbose:
                    print("Simulating ... ")
                sim = self._simargs[simargs]
                nxy = simulate_pil(ename, sexec, sname, sim.split(), 
                        force = not (sexec in self._simexecs))
                self._simexecs.add(sexec)
                self._simulated.add(sname)
            else :
                nxy = sname + '.nxy'

            # Third, analyze from file data
            if verbose:
                print("Analyze ... ")
            if metric == 'completion-time':
                sim_time = nxy_get_time_at_conc(nxy, reporter, conc)
                sim_conc = conc
                simcalc.append((sim_time, sim_conc))
            elif metric == 'diagonal-crossing-time':
                minfo = minfo.split(':')
                tmax = float(minfo[0])
                cmi = minfo[1].split(';')
                if len(cmi) == 2:
                    cmin = float(cmi[0])
                    cmax = float(cmi[1])
                else:
                    assert len(cmi) == 1
                    cmin = None
                    cmax = float(cmi[0])
                sim_time, sim_conc = nxy_get_diagonal_points(nxy, reporter, tmax, cmax, cmin = cmin)
                simcalc.append((sim_time, sim_conc))
            else:
                raise NotImplementedError('Metric "{}" not supported.'.format(metric))

            if cmpfig and pepperargs not in self.cmpfig:
                if verbose:
                    print("Compare ... ")
                tr = nxy_get_trajectory(nxy, reporter)
                if cmpfig is True:
                    tr.rename(columns={reporter: simargs}, inplace = True)
                elif cmpfig == 'hack': # used for zhang09 ...
                    tr.rename(columns={reporter: cname}, inplace = True)

                if trajectories is None:
                    trajectories = tr
                else :
                    assert np.array_equal(trajectories['time'], tr['time'])
                    if cmpfig is True:
                        trajectories[simargs] = tr[simargs]
                    elif cmpfig == 'hack':
                        trajectories[cname] = tr[cname]

        if cmpfig and pepperargs not in self.cmpfig: 
            trajectories.to_csv(cname, sep='\t', float_format='%.9e', index=False,
                header = ['{:15s}'.format(x) for x in list(trajectories)])

            self.cmpfig[pepperargs] = cname

        self._simcalc[pepperargs] = simcalc
        return self.get_system_dataframe(pepperargs)

    def get_system_dataframe(self, pargs='default'):
        name, pil, sim, rep, met, minfo, time, conc = zip(*self._simdata)
        if pargs not in self._simcalc:
            raise MissingDataError(
                    'Cannot find key "{}" in pepperargs'.format(pargs))
        stime, sconc = zip(*self._simcalc[pargs])

        # Only print initial conditions in the dataframe
        #ini = map(lambda x: x[x.find('--p0'):], sim)

        df = pd.DataFrame(data={
            'Input Filename': name,
            'Simulation': list(map(lambda x: x.replace('_',' '), sim)),
            'Reporter': rep,
            'Metric': met,
            'Metric-values': minfo,
            'Concentration (experiment)': conc,
            'Time (experiment)': time,
            'Time (simulation)': stime,
            'Concentration (simulation)': sconc,
            'Semantics': pargs},
            columns=['Input Filename', 'Simulation', 'Reporter', 'Metric', 'Metric-values', 'Semantics',
                'Concentration (simulation)', 'Time (simulation)',
                'Concentration (experiment)', 'Time (experiment)'])
                #'exp-time', 'pepperargs', 'sim-conc', 'sim-time'])

        return df

    def get_system_dataframes(self):
        name, pil, sim, rep, met, minfo, time, conc = zip(*self._simdata)

        simcalc = {None: None} if not self._simcalc else self._simcalc

        # Only print initial conditions in the dataframe
        #ini = map(lambda x: x[x.find('--p0'):], sim)

        for pargs, scalc in simcalc.items():
            stime, sconc = (None, None) if scalc is None else zip(*scalc)
            df = pd.DataFrame(data={
                'Input Filename': name,
                'Simulation': list(map(lambda x: x.replace('_',' '), sim)),
                'Reporter': rep,
                'Metric': met,
                'Metric-values': minfo,
                'Concentration (experiment)': conc,
                'Time (experiment)': time,
                'Time (simulation)': stime,
                'Concentration (simulation)': sconc,
                'Semantics': pargs},
            columns=['Input Filename', 'Simulation', 'Reporter', 'Metric', 'Metric-values', 'Semantics',
                    'Concentration (simulation)', 'Time (simulation)',
                    'Concentration (experiment)', 'Time (experiment)'])
                    #'exp-time', 'pepperargs', 'sim-conc', 'sim-time'])
            yield df

#done
def simulate_pil(pilfile, execute, output, command, force = False, verbose = 0):
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
    simname = execute
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
        proc = Popen(simu, stdin=crn, stdout=nxy, stderr=err)
        proc.communicate()
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

    df = pd.read_csv(nxyfile, sep='\s+', comment='#')
    return df[['time', species]]

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

    df = pd.read_csv(nxyfile, sep='\s+', comment='#')
    time = df['time'].values
    traj = df[species].values
    assert len(traj)>1
    if sign == 1 :
        pts = np.where(traj > th)[0]
        return time[pts[0]] if len(pts) else None
    elif sign == -1:
        pts = np.where(traj < th)[0]
        return time[pts[0]] if len(pts) else None
    print('WARNING: Returning none: {} ', nxyfile, species)
    return float('inf')

#done
def nxy_get_diagonal_points(nxyfile, species, tmax, cmax, cmin = None):
    """Return time and concentration when a trajectory crosses a diagonal line.
    
    Args:
        nxyfile (str): Path to a *.nxy file.
        species (str): The species of interest.
        tmax (flt): The time of the experiment.
        cmax (flt): The maximum concentration.

    Returns:
        (float, float): time, concentration
    """
    df = pd.read_csv(nxyfile, sep='\s+', comment='#')

    [sign, cmax] = [-1, -cmax] if cmax < 0 else [1, cmax]

    time = df['time'].values
    traj = df[species].values
    assert len(traj) > 1

    if sign == 1:
        positions = np.where(traj >= (cmax * (1 - time/tmax)))[0]
    elif sign == -1:
        assert cmin is not None
        positions = np.where(traj <= cmin + ((cmax-cmin) * (time/tmax)))[0]

    if len(positions) : 
        elem = positions[0]
        return time[elem], traj[elem]

    print('WARNING: Returning none: {} ', nxyfile, species)
    return None, None

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
    from zhang2011  import data as z11 # Only Figure 2
    from kotani2017 import data as k17
    
    from zhang2009_rates import data as z09r # 3-way branch migration rates
    from dabby2013_rates import data as d13r # 4-way branch migration rates

    # There are more Papers ... 
    # seelig2006.py
    # groves2015.py
    # srinivas2017.py
    # sun2018.py

    analysis = z07() + y08() + z09() + z10() + z11() + g11() + q11() + k17()
    analysis = z09r() + d13r()

    for fig in analysis:
        print("\n{}:".format(fig.name))
        fig.eval('default', verbose = 0)
        for df in fig.get_dataframes():
            print(df)

   
if __name__ == '__main__':
    main()

