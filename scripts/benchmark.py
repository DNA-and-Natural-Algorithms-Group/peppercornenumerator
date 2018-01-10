#!/usr/bin/env python

import sys
import seaborn as sns
import matplotlib.pyplot as plt
from math import log
from peppercornenumerator import Enumerator, __version__
from peppercornenumerator.objects import PepperDomain, PepperComplex, PepperReaction, clear_memory
from peppercornenumerator.condense import ReactionGraph
from peppercornenumerator.input import read_kernel
from peppercornenumerator.output import write_kernel

# INPUT experiment for 3-way strand displacement reactions
def k3way_kernel(c,d) :
    if c == 0 : # special leak case
        return """# Experiments Zhang 2009, 3-way strand displacement
        length a = 16
        length B = 19
        length b = 1
        length d = {:d} # m
        
        inv = B b
        clx = a B( b( + d* ) )
        clxO = a B( b + d* b* )
        
        i1 = a B( b + B b( + d* ) ) 
        i2 = B( b( + d* ) )
        i4 = B( b + d* b* )
        i3 = a B b
        
        # potential reactions:
        # clx -> clx0
        # clx0 -> clx
        # inv + clx0 -> i1
        # i1 -> inv + clx0
        # i1 -> i2 + i3
        """.format(m)
    else :
        return """# Experiments Zhang 2009, 3-way strand displacement
        length a = 16
        length B = 20
        length c = {:d}
        {}length d = {:d} # d
        
        inv = B c
        clx = a B( + {} c* )
        
        i1 = a B( + B c( + {} ) )
        i2 = B( c( + {} ) )
        i3 = a B
        
        # potential reactions:
        # inv + clx -> i1
        # i1 -> i2 + i3
        # i1 -> inv + clx
        """.format( c,
                '#' if d==0 else '', d,
                '' if d==0 else 'd*',
                '' if d==0 else 'd*',
                '' if d==0 else 'd*')
 
k3way_exp = [ 
       #((n,  m),  k)
        #((0,  15), 1.40), 
        ((1,  14), 8.17), 
        ((2,  13), 144), 
        ((3,  12), 1.08e3), 
        ((4,  11), 5.05e4), 
        ((5,  10), 9.64e5), 
        ((6,  9),  2.36e6), 
        ((7,  8),  3.22e6), 
        ((8,  9),  3.15e6), 
        ((9,  8),  2.77e6),
        ((10, 7),  2.83e6),
        ((15, 0),  4.78e6)]

def k3way_rates(complexes, enum):
    [k0, k1, k2] = [0.,0.,0.]
    for rxn in enum.reactions:
        if rxn.reactants == sorted([complexes['inv'], complexes['clx']]):
            # forward binding rate
            assert rxn.rtype == 'bind21'
            k0 = rxn.rate
        elif rxn.reactants == [complexes['i1']] and rxn.products == sorted([complexes['i2'], complexes['i3']]):
            # branch migration & strand displacement step
            assert rxn.rtype == 'branch-3way'
            k1 = rxn.rate
        elif rxn.reactants == [complexes['i1']] and rxn.products == sorted([complexes['inv'], complexes['clx']]):
            # toehold dissociation
            assert rxn.rtype == 'open'
            k2 = rxn.rate
    return k0*k1/(k1+k2)
 
def k3way_condensed(complexes, enumRG):
    assert len(enumRG.condensed_reactions) == 1
    return enumRG.condensed_reactions[0].rate

## from Dave Zhang's matlab script, ToeEx_rates_vs_model.m
k3wayX_exp = [
                                                              ((1, 4), 7.70),   ((1, 3),   5.48), ((1, 2),   23.5), ((1, 1),   18.9),
                                            ((2, 5), 43.6),   ((2, 4), 214.05), ((2, 3),  273.0), ((2, 2),  249.0), ((2, 1),  231.0), 
                          ((3, 6), 66.9),   ((3, 5), 215.0),  ((3, 4), 939.0),  ((3, 3),  974.0), ((3, 2),  907.0), ((3, 1),  846.0),
        ((4, 7), 131.0),  ((4, 6), 407.0),  ((4, 5), 4.25e3), ((4, 4), 2.13e4), ((4, 3), 2.41e4), ((4, 2), 2.29e4), ((4, 1), 1.97e4), 
        ((5, 7), 3.59e3), ((5, 6), 9.72e4), ((5, 5), 3.45e5), ((5, 4), 1.53e6), ((5, 3), 1.58e6), ((5, 2), 1.58e6), ((5, 1), 1.73e6),
        ((6, 7), 1.61e5), ((6, 6), 4.05e5), ((6, 5), 1.48e6), ((6, 4), 3.04e6), ((6, 3), 2.59e6), ((6, 2), 3.00e6), 
        ((7, 7), 4.7e5),  ((7, 6), 1.11e6), ((7, 5), 2.90e6), ((7, 4), 3.57e6), 
        ((8, 7), 1.94e6), ((8, 6), 2.68e6), ((8, 5), 3.14e6), ((8, 4), 3.37e6)]

def k3wayX_kernel(n,m) :
    return """# Experiments Zhang 2009, 3-way strand displacement toehold exchange
        length a = 16
        length b = {:d} # m
        length B = {:d} # 20 - m
        length c = {:d} # n
        length d = {:d} # 15 - n
        
        inv = B c
        clx = a b( B( + d* c* ) )
        
        i1 = a b( B( + B c( + d* ) ) )
        i2 = a b( B + B( c( + d* ) ) )
        i3 = B( c( + d* ) ) b*
        i4 = a b B
        
        # potential reactions:
        # inv + clx -> i1
        # i1 -> inv + clx
        # i1 -> i2
        # i2 -> i1
        # i2 -> i3 + i4
        # i3 + i4 -> i2
        """.format(m, 20-m, n, 15-n)

def k3wayX_rates(complexes, enum):
    [k0, k1, k2,k3,k4,k5] = [0.,0.,0.,0.,0.,0.]
    for rxn in enum.reactions:
        if rxn.reactants == sorted([complexes['inv'], complexes['clx']]) and rxn.products == [complexes['i1']]:
            # forward binding rate
            assert rxn.rtype == 'bind21'
            k0 = rxn.rate
        elif rxn.reactants == sorted([complexes['i3'], complexes['i4']]) and rxn.products == [complexes['i2']]:
            # reverse binding rate
            assert rxn.rtype == 'bind21'
            k1 = rxn.rate
        elif rxn.reactants == [complexes['i1']] and rxn.products == [complexes['i2']]:
            # forward branch migration 
            assert rxn.rtype == 'branch-3way'
            k2 = rxn.rate
        elif rxn.reactants == [complexes['i2']] and rxn.products == [complexes['i1']]:
            # reverse branch migration 
            assert rxn.rtype == 'branch-3way'
            k3 = rxn.rate
        elif rxn.reactants == [complexes['i2']] and rxn.products == sorted([complexes['i3'], complexes['i4']]):
            # forward dissociation
            assert rxn.rtype == 'open'
            k5 = rxn.rate
        elif rxn.reactants == [complexes['i1']] and rxn.products == sorted([complexes['inv'], complexes['clx']]):
            # reverse dissociation
            assert rxn.rtype == 'open'
            k4 = rxn.rate
        else: 
            raise Exception('missing reaction')
    return k0*(k5/(k3+k5)) / ( (k2+k4)/k2 - k3/(k3+k5) )
  
def k3wayX_condensed(complexes, enumRG):
    assert len(enumRG.condensed_reactions) == 2
    for rxn in enumRG.condensed_reactions:
        for rset in rxn.reactants:
            if complexes['inv'] in rset.complexes :
                return rxn.rate
    return None
 
# (n,m,k1_fit) from and Nadine Dabby, Caltech PhD Thesis, Table 5.2  (note m,n have consistent meaning, but order in table is swapped.)
k4way_exp = [ 
        #((0, 0), 0.034),
        ((0, 2), 0.047),
        ((2, 2), 0.10),
        ((2, 0), 0.033),
        ((4, 2), 0.93),
        ((4, 0), 0.039),
        ((0, 4), 0.97),
        ((2, 4), 56),
        ((6, 2), 490),
        ((0, 6), 58),
        ((4, 4), 770),
        ((6, 0), 5.0),
        ((2, 6), 9.4e3),
        ((4, 6), 7.0e4),
        ((6, 4), 2.8e5),
        ((6, 6), 6.9e5) ]

def k4way_kernel(n,m) :
    assert n<7 and (m<7 or m==16)
    M = 0 if (m==6 or m==16) else 6-m
    N = 6-n

    return """# Nadine Dabby experiments
        length x = 21
        {}length m = {} # m
        {}length M = {} # M
        {}length n = {} # n
        {}length N = {} # N

        # starting states

        ## x*( m* M* + N* n* )
        rep = x*( {} {} + {} {} )

        ## m x( + ) n
        clx = {} x( + ) {}

        # transient states

        ## m( x( + ) + N* x( + ) ) M* 
        int1 = {} x( + ) + {} x( + ) {} {} 

        ## x( + ) n( + ) x( + ) M*
        int2 = x( + ) {} + {} x( + ) {}

        ## x( + ) n( + N* ) x( + ) M* 
        {}int3 = x( + ) {} + {} {} x( + ) {}

        # products
        # x*( + N* ) 
        npr1 = x*( + {} ) 
        # x*( n( + ) ) 
        {}npr2 = x*( {} + {} ) 
        # x*( n( + N* ) ) 
        {}npr3 = x*( {} + {} {} ) 

        # m( x( + ) ) M* 
        mpr = {} x( + ) {} {}

        """.format(
                '#' if m==0 else '', m,
                '#' if M==0 else '', M,
                '#' if n==0 else '', n,
                '#' if N==0 else '', N,
                '' if m==0 else 'm*',
                '' if M==0 else 'M*',
                '' if N==0 else 'N*',
                '' if n==0 else 'n*',
                '' if m==0 else 'm',
                '' if n==0 else 'n',
                # int1
                '' if m==0 else 'm(',
                '' if N==0 else 'N*',
                '' if m==0 else ')',
                '' if M==0 else 'M*',
                # int2
                '' if n==0 else 'n(',
                '' if n==0 else ')',
                '' if M==0 else 'M*',
                # int3
                '#' if N==0 else '',
                '' if n==0 else 'n(',
                '' if N==0 else 'N*',
                '' if n==0 else ')',
                '' if M==0 else 'M*',
                # npr
                '' if N==0 else 'N*',
                # npr2
                '#' if N==0 else '',
                '' if n==0 else 'n(',
                '' if n==0 else ')',
                # npr3
                '#' if n==0 else '',
                '' if n==0 else 'n(',
                '' if N==0 else 'N*',
                '' if n==0 else ')',
                # mpr
                '' if m==0 else 'm(',
                '' if m==0 else ')',
                '' if M==0 else 'M*')

def k4way_rates(complexes, enum):
    if len(enum.reactions) == 3:

        # Find out wich complexes are actually specified
        rep = complexes['rep']
        clx = complexes['clx']
        if complexes['int1'] in enum.complexes :
            cint = complexes['int1'] 
        elif complexes['int2'] in enum.complexes :
            cint = complexes['int2'] 
        elif 'int3' in complexes and complexes['int3'] in enum.complexes :
            cint = complexes['int3'] 
        else :
            raise Exception
        if complexes['npr1'] in enum.complexes :
            npr = complexes['npr1'] 
        elif 'npr2' in complexes and complexes['npr2'] in enum.complexes :
            npr = complexes['npr2'] 
        elif 'npr3' in complexes and complexes['npr3'] in enum.complexes :
            npr = complexes['npr3'] 
        else :
            raise Exception
 
        mpr = complexes['mpr']

        # ok now lets get the rates we need
        for rxn in enum.reactions:
            if rxn.reactants == sorted([rep, clx]) and rxn.products == [cint]:
                k0 = rxn.rate
            elif rxn.reactants == [cint] and rxn.products == sorted([rep, clx]):
                k2 = rxn.rate
            elif rxn.reactants == [cint] and rxn.products == sorted([npr, mpr]):
                k1 = rxn.rate
            else:
                raise Exception

        return k0*k1/(k1+k2)
    else :
        assert len(enum.reactions) == 13
        
        # Find out wich complexes are actually specified
        rep = complexes['rep']
        clx = complexes['clx']

        if complexes['npr1'] in enum.complexes :
            npr = complexes['npr1'] 
        elif 'npr2' in complexes and complexes['npr2'] in enum.complexes :
            npr = complexes['npr2'] 
        elif 'npr3' in complexes and complexes['npr3'] in enum.complexes :
            npr = complexes['npr3'] 
        else :
            raise Exception
        mpr = complexes['mpr']

        # both toeholds bound, not 4-way
        LRint = filter(lambda x: x.products == sorted([npr, mpr]), enum.reactions)
        assert len(LRint) == 1
        LRint = LRint[0].reactants[0]

        # left thoehold bound and 4-way
        Lint2 = filter(lambda x: x.products == sorted([mpr]), enum.reactions)
        assert len(Lint2) == 1
        Lint2 = Lint2[0].reactants[0]
                
        # right thoehold bound and 4-way
        Rint2 = filter(lambda x: x.products == sorted([npr]), enum.reactions)
        assert len(Rint2) == 1
        Rint2 = Rint2[0].reactants[0]

        # left thoehold bound, no 4-way
        Lint = filter(lambda x: x.products == sorted([Lint2, npr]), enum.reactions)
        assert len(Lint) == 1
        Lint = Lint[0].reactants[0]

        # right thoehold bound, no 4-way
        Rint = filter(lambda x: x.products == sorted([Rint2, mpr]), enum.reactions)
        assert len(Rint) == 1
        Rint = Rint[0].reactants[0]

        # ok now lets get the rates we need
        for rxn in enum.reactions:
            if rxn.reactants == sorted([rep, clx]) and rxn.products == [Lint]:
                # binding of left toehold
                k0 = rxn.rate
            elif rxn.products == sorted([rep, clx]) and rxn.reactants == [Lint]:
                # dissociation of left toehold
                k1 = rxn.rate
            elif rxn.reactants == sorted([rep, clx]) and rxn.products == [Rint]:
                # binding of right toehold
                k2 = rxn.rate
            elif rxn.products == sorted([rep, clx]) and rxn.reactants == [Rint]:
                # dissociation of right toehold
                k3 = rxn.rate
            elif rxn.reactants == [Rint] and rxn.products == [LRint]:
                # binding of second toehold 
                k4 = rxn.rate
            elif rxn.reactants == [LRint] and rxn.products == [Rint]:
                # dissciation of second toehold 
                k5 = rxn.rate
            elif rxn.reactants == [Lint] and rxn.products == [LRint]:
                # binding of second toehold 
                k6 = rxn.rate
            elif rxn.reactants == [LRint] and rxn.products == [Lint]:
                # dissciation of second toehold 
                k7 = rxn.rate
            elif rxn.reactants == [Lint] and rxn.products == sorted([Lint2, npr]):
                # b4way of left side 
                k8 = rxn.rate
            elif rxn.reactants == [Lint2] and rxn.products == [mpr]:
                k9 = rxn.rate
            elif rxn.reactants == [Rint] and rxn.products == sorted([Rint2, mpr]):
                # b4way of right side 
                k10 = rxn.rate
            elif rxn.reactants == [Rint2] and rxn.products == [npr]:
                # reverse b4way of right side 
                k11 = rxn.rate
            elif rxn.reactants == [LRint] and rxn.products == sorted([npr, mpr]):
                # b4way of both sides
                k12 = rxn.rate
            else :
                raise Exception

        k = [ # a = up; b = down
                None, # k[0] = k11
                k6,   # k[1]
                k4,   # k[2]
                None, # k[3] = k9
                k0,   # k[4] 
                k2,   # k[5]   
                k8,   # k[6]
                k10,  # k[7]
                k12,  # k[8]   !
                k1,   # k[9]   !
                k3,   # k[10]  !
                k7,   # k[11]
                k5]   # k[12]

        #k = [
        #        None, # k[0] = k11
        #        k4,   # k[1]
        #        k6,   # k[2]
        #        None, # k[3] = k9
        #        k2,   # k[4] 
        #        k0,   # k[5]   
        #        k10,  # k[6]
        #        k8,   # k[7]
        #        k12,  # k[8]   !
        #        k3,   # k[9]   !
        #        k1,   # k[10]  !
        #        k5,   # k[11]
        #        k7]   # k[12]

        ## k1 and k2 
        a5  = k[2]/(k[10]+k[7]+k[2])
        a6  = k[1]/(k[9]+k[6]+k[1])

        ## k6 and k7 
        b5  = k[7]/(k[10]+k[7]+k[2])
        b6  = k[6]/(k[9]+k[6]+k[1])

        ## k11 and k12
        y5  = k[12]/(k[8]+k[11]+k[12])
        y6  = k[11]/(k[8]+k[11]+k[12])

        ## k[8] = k12
        b10 = k[8]/(k[8]+k[11]+k[12])
        P10 = (b10+y6*b6+y5*b5)/(1-(y6*a6+y5*a5))

        # k4 and k5 
        k_eff = (k[4]*a6+k[5]*a5)*P10 + (k[4]*b6+k[5]*b5)
        return k_eff

def k4way_condensed(complexes, enumRG):
    assert len(enumRG.condensed_reactions) == 1
    return enumRG.condensed_reactions[0].rate

def benchmark(data, template, manualrates = None, condensedrates = None, autostart=True):
    for (setup, k_exp) in data:
        clear_memory()
        kernelstring = template(*setup)
        #print kernelstring
    
        complexes, reactions = read_kernel(kernelstring)
        if autostart:
            enum = Enumerator(complexes.values(), reactions)
        else :
            #print 'selected complexres'
            enum = Enumerator([complexes['rep'], complexes['clx']], reactions)
        enum.k_slow=1e-5
        #enum.k_fast=1e-3
        enum.release_cutoff = 100
        #enum.max_helix_migration = False
        enum.enumerate()
        #write_kernel(enum, sys.stdout, condensed = False)
    
        # Extract rates for manual calculation of effective rate constant
        if manualrates:
            k_eff = manualrates(complexes, enum)
        else :
            k_eff = None
   
        enumRG = ReactionGraph(enum)
        enumRG.condense()
        #write_kernel(enum, sys.stdout, condensed = True)

        if condensedrates :
            k_con = condensedrates(complexes, enumRG)
        else :
            k_con = None
    
        #print("k_eff = {}\nk_con = {}\nk_exp = {}\n".format(k_eff, k_con, k_exp))

        assert abs(k_con - k_eff) < 0.00001
        yield (k_con, k_exp)

fig = plt.figure()
ax1 = fig.add_subplot(111)

#print [(x,y) for (x,y) in benchmark(k3way_exp, k3way_kernel, k3way_rates, k3way_condensed)]
xs = []
ys = []
for (x,y) in benchmark(k3way_exp, k3way_kernel, k3way_rates, k3way_condensed):
    xs.append(log(x,10))
    ys.append(log(y,10))

ax1.scatter(xs, ys, color='green', label='3way')

#print [(x,y) for (x,y) in benchmark(k3wayX_exp, k3wayX_kernel, k3wayX_rates, k3wayX_condensed)]
xs = []
ys = []
for (x,y) in benchmark(k3wayX_exp, k3wayX_kernel, k3wayX_rates, k3wayX_condensed):
    xs.append(log(x,10))
    ys.append(log(y,10))

ax1.scatter(xs, ys, color='red', label='3wayX')

#print [(x,y) for (x,y) in benchmark(k4way_exp, k4way_kernel, k4way_rates, k4way_condensed, False)]
xs = []
ys = []
for (x,y) in benchmark(k4way_exp, k4way_kernel, k4way_rates, k4way_condensed, False):
    xs.append(log(x,10))
    ys.append(log(y,10))

ax1.scatter(xs, ys, color='blue', label='4way')



plt.legend(loc='upper left');
plt.xlim(-3, 8)
plt.ylim(-3, 8)

pfile = 'test.pdf'
plt.savefig(pfile)

