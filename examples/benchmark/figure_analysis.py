#!/usr/bin/python

import os
import sys
from math import log
import subprocess as sub
import matplotlib.pyplot as plt

from peppercornenumerator import Enumerator, __version__
from peppercornenumerator.objects import clear_memory
from peppercornenumerator.condense import PepperCondensation
from peppercornenumerator.input import read_kernel, read_pil
from peppercornenumerator.output import write_kernel

# updated
from zhang2007 import setups as z07
from zhang2009 import setups as z09 # missing reaction rate data & sequence variation
from zhang2010 import setups as z10 # there is more!
from genot2011 import setups as g11
from qian2011 import setups as q11 # missing F3 & F4 & Supplement?
from kotani2017 import setups as k17

from zhang2009_rates import setups as z09r # 3-way branch migration rates
from dabby2013_rates import setups as d13r # 4-way branch migration rates

# Troublemakers:
# from yin2008 import setups as y08
# from zhang2010_cooperative import setups as z10c

# More:
# groves2015.py
# srinivas2017.py
# seelig2006.py

analysis = [z07(), z09(), z10(), g11(), q11(), k17()] # all completion times
analysis = [z09r(), d13r()] # rates only
analysis = [z07(), z09(), z10(), g11(), q11(), k17(), z09r(), d13r()] # everything

paperdata = []
mcmax = []
for paper in analysis:
    paperdata += paper
    mcmax.append(len(paper))

def peppercorn(kernelstring, name, 
        condensed=True, 
        conc='nM', 
        max_complex_size = 10,
        release_cutoff = 13,
        k_fast=0):
    """ A wrapper for peppercorn.
    """
    #print kernelstring
    complexes, reactions = read_pil(kernelstring)
    enum = Enumerator(complexes.values(), reactions)
    enum.release_cutoff = release_cutoff
    enum.max_complex_size = max_complex_size
    enum.k_fast = k_fast
    enum.enumerate()

    #condensed = False
    if condensed:
        enumRG = PepperCondensation(enum)
        enumRG.condense()

    detailed = not condensed

    #write_kernel(enum, sys.stdout, detailed, condensed)
    with open(name + '.crn', 'w') as crn:
        write_kernel(enum, crn, detailed, condensed, molarity=conc)

    if condensed:
        return enumRG, name + '.crn'
    else :
        return enum, name + '.crn'

def simulate_crn(infile, name, crnsimu):
    assert infile == name + '.crn'
    # Do the simulation (catch treekin errors)
    print "{}".format(' '.join(crnsimu))
    with open(infile, 'r') as crn, \
            open(name + '.nxy', 'w') as nxy, \
            open(name + '.err', 'w') as err:
        proc = sub.Popen(crnsimu, stdin=crn, stdout=nxy, stderr=err)
        proc.communicate(None)
        if proc.returncode:
            raise Exception
    return name + '.nxy'

def get_simulated_time(nxyfile, species, threshold):
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
            if d[idx] > threshold:
                return d[0]
    return None

def main():
    logdata = True

    if not os.path.exists('tmp'):
        raise SystemExit('Please make a directory called "tmp" to store temorary results')

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    mycolors = list('bgrcmyk')
    mymarker = list('o*^.vph+D')
    (mc,mm) = (0,0)

    allresults = []
    for data in paperdata:
        results = []
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
                    (et, ec) = data['exp_results'][e+i]
                    time = get_simulated_time(nxy, rep, ec)
                    #print e+i, rep, et, ec, time
                    results.append([et, time])

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
        ax1.scatter(xs, ys, color=mycolors[mc], marker=mymarker[mm], label=data['name'])
        mc += 1
        if mc >= mcmax[mm]:
            mc = 0
            mm += 1

    plt.title('Peppercorn vs. experiment');
    if logdata:
        ax1.set_xlabel('Experimental system speed [$\log_{10}(s)$]', fontsize=16)
        ax1.set_ylabel('Simulated system speed [$\log_{10}(s)$]', fontsize=16)
        (mi,ma)=(0, 6)
        (mi,ma)=(-3, 7)
        #plt.xlim(mi, ma)
        #plt.ylim(mi, ma)
        plt.plot([mi, ma], [mi, ma], color='black')
    else:
        ax1.set_xlabel('Experimental system speed [s]', fontsize=16)
        ax1.set_ylabel('Simulated system speed [s]', fontsize=16)
        (mi,ma) = (-10, 600)
        plt.xlim(mi, ma)
        plt.ylim(mi, ma)
        plt.plot([mi, ma], [mi, ma], color='black')

    lgd = plt.legend(bbox_to_anchor=(1.40, 1.0), loc='upper right', borderaxespad=0.)
    #plt.legend(loc='lower right');
    #plt.legend(loc='lower right');
    
    #pfile = 'zhang2009_condensed_log.png'
    #plt.savefig(pfile, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.savefig('test.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')

    return

if __name__ == '__main__':
   main()

