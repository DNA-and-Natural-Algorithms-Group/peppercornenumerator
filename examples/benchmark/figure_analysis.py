#!/usr/bin/python

import os
import sys
from math import log
import subprocess as sub
import matplotlib.pyplot as plt
import numpy as np

from peppercornenumerator import Enumerator, __version__
from peppercornenumerator.objects import clear_memory
from peppercornenumerator.condense import PepperCondensation
from peppercornenumerator.input import read_pil, read_seesaw, ParseException
from peppercornenumerator.output import write_pil

# updated
from zhang2007 import setups as z07
from yin2008 import setups as y08
from zhang2009 import setups as z09 # missing reaction rate data & sequence variation
from zhang2010 import setups as z10 # there is more!
from genot2011 import setups as g11
from qian2011 import setups as q11 # missing F3 & F4 & Supplement?
from qian2011sqrt import setups as q11sq # missing F3 & F4 & Supplement?
from kotani2017 import setups as k17

from zhang2009_rates import setups as z09r # 3-way branch migration rates
from dabby2013_rates import setups as d13r # 4-way branch migration rates
from zhang2011_cooperative import setups as z11c # Only Figure 2, and it is a bit strange...

# More:
# groves2015.py
# srinivas2017.py
# seelig2006.py

analysis = [z07(), y08(), z09(), z10(), z11c(), g11(), q11(), k17()] # all completion times
analysis = [z09r(), d13r()] # rates only
analysis = [z07(), y08(), z09(), z10(), z11c(), g11(), q11(), k17(), z09r(), d13r()] # everything

analysis = [z07(), y08(), z09(), z10(), q11(), z11c(), k17()] # all_comp.svg

#analysis = [z09r(), d13r()] # rates only
analysis = [q11()] # everything

# feedforward: z07_F3, q11, k17_F3

# Catalysts : z07_F1, z09_F5, k17_F2
# Auto-catalysts : z07_F4, y08_F3, k17_F4

# missing: z10, z11c, g11

paperdata = []
mcmax = []
for paper in analysis:
    paperdata += paper
    mcmax.append(len(paper))

def peppercorn(kernelstring, name, 
        condensed = True, 
        conc = 'nM', 
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

    with open(name + '_input.pil', 'w') as pil:
        pil.write(kernelstring)

    try:
        complexes, reactions = read_pil(kernelstring)
    except ParseException, ex_pil:
        complexes, reactions = read_seesaw(kernelstring, explicit=False)

    enum = Enumerator(complexes.values(), reactions)
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

    #condensed = True
    if condensed:
        enumRG = PepperCondensation(enum)
        enumRG.condense()

    detailed = not condensed

    enumfile = name + '_enum.pil'

    with open(enumfile, 'w') as crn:
        # You can find this file in the tmp directory
        write_pil(enum, crn, detailed, condensed, molarity=conc)

    if condensed:
        return enumRG, enumfile
    else :
        return enum, enumfile

def simulate_crn(infile, name, crnsimu):
    assert infile == name + '_enum.pil'
    # Do the simulation (catch treekin errors)
    print "{}".format(' '.join(crnsimu))
    with open(infile, 'r') as crn, \
            open(name + '.nxy', 'w') as nxy, \
            open(name + '.err', 'w') as err:
        proc = sub.Popen(crnsimu, stdin=crn, stdout=nxy, stderr=err)
        proc.communicate(None)
        if proc.returncode:
            raise Exception(err)
    return name + '.nxy'

def get_simulated_trajectory(nxyfile, species):
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

            d = l.strip().split()
            trajectory.append([d[0], d[idx]])
    return trajectory

def get_conc_at_time(nxyfile, species, threshold):
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

            d = map(float, l.strip().split())
            if d[0] >= threshold:
                return d[idx]
    return None

def get_simulated_time(nxyfile, species, threshold):
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

            d = map(float, l.strip().split())
            if sign == 1 and d[idx] > th:
                return d[0]
            elif sign == -1 and d[idx] < th:
                return d[0]
    return None

def main():
    logdata = True
    cmpfig = True

    if not os.path.exists('tmp'):
        raise SystemExit('Please make a directory called "tmp" to store temorary results')

    #fig = plt.figure()
    fig, ax1 = plt.subplots(1, 1, figsize=(3, 2.5))
    #ax1  = plt.subplot(1,1,1)
    #ax2  = plt.subplot(1,2,2)
    mycolors = ['orange'] + list('bgrcmykkkkkkkkkkk')
    mymarker = list('o*^.vph+D8,|_oooooooo')
    (mc,mm) = (0,0)

    allresults = []
    for data in paperdata:
        results = []

        trajectories = []
        for e, pilstring in enumerate(map(data['piltemplate'], data['pilparams'])):
            clear_memory()
            tmpname = 'tmp/'+data['name'] + '_' + str(e)
            print tmpname
            enumOBJ, crn = peppercorn(pilstring, name=tmpname, **data['pepperargs'])

            if 'rates' in data:
                if len(enumOBJ.condensed_reactions) == 1:
                    rxn = enumOBJ.condensed_reactions[0]
                    er = data['exp_results'][e]
                    results.append([er, rxn.rate])
                else :
                    if 'reactants' in data['rates']:
                        for rxn in enumOBJ.condensed_reactions:
                            if sorted([rs.name for rs in rxn.reactants]) == sorted(data['rates']['reactants']):
                                er = data['exp_results'][e]
                                results.append([er, rxn.rate])
                                break
                    else:
                        raise NotImplementedError('multiple condensed reactions')

            if 'simulation' in data:
                for i, command in enumerate(data['simulation']):
                    nxy = simulate_crn(crn, tmpname, command)
                    rep = data['reporter']
                    if 'off_results' in data:
                        dresults = data['exp_results'] + data['off_results']
                    else:
                        dresults = data['exp_results']
                    (et, ec) = dresults[e+i]
                    time = get_simulated_time(nxy, rep, ec)
                    conc = get_conc_at_time(nxy, rep, et)
                    print e+i, rep, et, ec, time, conc
                    if e+i < len(data['exp_results']):
                        results.append([et, time])
                    else :
                        pass
                        #results.append([ec, conc])

                    if cmpfig :
                        tr = get_simulated_trajectory(nxy, rep)

                    if trajectories :
                        for te, [time, val] in enumerate(tr):
                            if len(trajectories) > te:
                                assert trajectories[te][0] == time
                                trajectories[te].append(val)
                    else :
                        trajectories = tr
            
        with open('tmp/'+data['name'] + '_cmp.nxy', 'w') as cnxy:
            for td in trajectories:
                cnxy.write("{}\n".format(' '.join(td)))

        #assert len(results) == len(data['exp_results'])
        xs = []
        ys = []
        for e, (x,y) in enumerate(results):
            if x is None or y is None:
                print 'WARNING: Skipping data points {}'.format(e)
                continue
            if logdata:
                xs.append(log(x,10))
                ys.append(log(y,10))
            else:
                xs.append(x)
                ys.append(y)
        print xs
        print ys
        ax1.scatter(xs, ys, s=10, color=mycolors[mc], marker=mymarker[mm], label=data['name'])
        mc += 1
        if mc >= mcmax[mm]:
            mc = 0
            mm += 1

    #plt.title('Remote toeholds (condensed)', y=1.08)
    if logdata:
        ax1.set_xlabel('Experimental system speed [$\log_{10}(s)$]', fontsize=8)
        ax1.set_ylabel('Peppercorn system speed [$\log_{10}(s)$]', fontsize=8)
        #ax1.set_xlabel('Experimental $\log(k)$ [$/M/s$]', fontsize=9)
        #ax1.set_ylabel('Peppercorn $\log(k)$ [$/M/s$]', fontsize=9)
        (mi,ma)=(0, 6)
        #plt.xlim(mi, ma)
        #plt.ylim(mi, ma)
        ax1.plot([mi, ma], [mi, ma], color='black')
    else:
        ax1.set_xlabel('Experimental system speed [s]', fontsize=9)
        ax1.set_ylabel('Simulated system speed [s]', fontsize=9)
        (mi,ma) = (-10, 600)
        plt.xlim(mi, ma)
        plt.ylim(mi, ma)
        plt.plot([mi, ma], [mi, ma], color='black')

    #lgd = plt.legend(bbox_to_anchor=(1.40, 1.0), loc='upper right', borderaxespad=0.)
    lgd = plt.legend(loc='upper left', borderaxespad=1., fontsize=7)
    plt.xticks(np.arange(mi, ma+1, step=1))
    plt.yticks(np.arange(mi, ma+1, step=1))
    #plt.legend(loc='lower right');
    #plt.legend(loc='lower right');
    
    #pfile = 'all_systems.svg'
    #plt.savefig(pfile, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.savefig('test.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')

    return

if __name__ == '__main__':
   main()

