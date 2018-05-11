#!/usr/bin/python

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

# Troublemakers:
# from yin2008 import setups as y08
# from zhang2010_cooperative import setups as z10c

# More:
# groves2015.py
# srinivas2017.py
# seelig2006.py

paperdata = z07() + z09() + z10() + g11() + q11() + k17()
#paperdata = q11()

mcmax = map(len, [z07(), z09(), z10(), g11(), q11(), k17()])
#mcmax = [13]

def peppercorn(kernelstring, name, condensed=True, conc='nM', crn=True, k_fast=0):
    """ A wrapper for peppercorn.
    """
    #print kernelstring
    complexes, reactions = read_pil(kernelstring)
    enum = Enumerator(complexes.values(), reactions)
    enum.release_cutoff = 13
    enum.max_complex_size = 10
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

    return None, name + '.crn'

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
            tmpname = data['name'] + '_' + str(e)
            print tmpname
            #print pilstring
            out, crn = peppercorn(pilstring, name=tmpname, crn=True, **data['pepperargs'])

            if 'rates' in data:
                raise NotImplementedError

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
        plt.xlim(mi, ma)
        plt.ylim(mi, ma)
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

