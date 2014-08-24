
# must be run in this directory, because the link to enumerator.py is hard-coded.

# you can turn on & off the three comparison runs by editing the flags here

# Zhang & Winfree 2009, figures 3B and 4B, toehold-mediated strand displacement and toehold exchange, respectively
# TE + 12.5 mM Mg++, pH 8.0, 25 C
# (1 uM poly-T20 carrier DNA added to all samples with 1uM or less of sample DNA)
do3way = True
do3way_exchange = True

# Nadine Dabby PhD thesis, table 5.2, 
# TAE = 12.5 mM Mg++, 25 C, 
# (1 uM poly-T20 carrier DNA added to all samples with 10nM or less of sample DNA, and pipet tips pre-absorbed with carrier before use)
do4way = True

# you can turn on & off the case-by-case examination of the enumerated output PIL file here
watch = False

import subprocess

# This is a test of the command-line parsing, as well as getting number for kinetics comparisons

def CMI_enum(inputstring,maxtoesize=6,semantics='detailed'):
    """Calls enumerator via the command line interface, first writing input string as a file, then reading the output file."""

    text_file = open("tmp-run.pil", "w")
    text_file.write(inputstring)
    text_file.close()

    if semantics=='condensed':
        cmd = ['../../enumerator.py','tmp-run.pil','-c','--release-cutoff',str(maxtoesize)]
    elif semantics=='detailed':
        cmd = ['../../enumerator.py','tmp-run.pil','--release-cutoff',str(maxtoesize)]
    else:
        raise NameError( "semantics must be either 'condensed' or 'detailed'." )

    subprocess.Popen(cmd).wait()

    text_file = open ("tmp-run-enum.pil", "r")
    output=text_file.read()
    text_file.close()    

    lines = output.split('\n')

    return lines


def exchange(n,m):
    """run enumerator for Zhang & Winfree 2009, figure 4B"""

    # only if both m>0 and n>0 -- fix that  (and then it will work for figure 3B too)

    assert m<8 and n<16
    if n==0:
        return 'Leak not modeled'
    elif n==15:
        assert m==0
        sys = "length a = 16\nlength B = 20\nlength c = %d\nB c\na B( + c* )\n" % n
    elif m==0:
        sys = "length a = 16\nlength B = 20\nlength c = %d\nlength d = %d\nB c\na B( + d* c* )\n" % (n,15-n)
    elif m>0:
        sys = "length a = 16\nlength b = %d\nlength B = %d\nlength c = %d\nlength d = %d\nB c\na b(B( + d* c* ))\n" % (m,20-m,n,15-n)

    pil_enum = CMI_enum(sys,8,'detailed')
    rates = [s for s in pil_enum if len(s)>0 and s[0]=='k']
    if watch:
        for s in pil_enum:
            print s


    # trust that the enumerator always lists reactions in a consistent order!
    if len(rates)==2:  # must be irreversible toehold, detailed model  
        k_eff = float(rates[0].split()[1][1:])   # forward binding rate
    elif len(rates)==3:  # must be reversible toehold, detailed model
        k0=float(rates[0].split()[1][1:])   # forward binding rate
        k1=float(rates[1].split()[1][1:])   # branch migration & strand displacement step
        k2=float(rates[2].split()[1][1:])   # toehold dissociation
        k_eff = k0*k1/(k1+k2)
    elif len(rates)==6:  # must be reversible toehold exchange, detailed model
        k0=float(rates[0].split()[1][1:])   # forward binding rate
        k1=float(rates[1].split()[1][1:])   # reverse binding rate
        k2=float(rates[2].split()[1][1:])   # forward branch migration step
        k3=float(rates[3].split()[1][1:])   # reverse branch migration step
        k4=float(rates[4].split()[1][1:])   # invading toehold dissociation
        k5=float(rates[5].split()[1][1:])   # incumbent toehold dissociation
        k_eff = k0*(k5/(k3+k5)) / ( (k2+k4)/k2 - k3/(k3+k5) )
        
    pil_enum = CMI_enum(sys,8,'condensed')
    rates = [s for s in pil_enum if len(s)>0 and s[0]=='k']

    if len(rates)==1:  # irreversibble toehold-mediated strand displacemen
        k_con = float(rates[0].split()[1][1:])
    elif len(rates)==2:  # reversible toehold exchange   #### check by hand to make sure first one is always forward
        k_con = float(rates[0].split()[1][1:])   # forward binding rate
    # now must modify stuff below to output & compare condensed rates

    if watch:
        for s in pil_enum:
            print s
        print "Calculated k_eff = %f /M/s from detailed reactions and k_con = %f /M/s from condensed reactions." % (k_eff,k_con)
        raw_input("Press enter to continue...")  # in python 3, just input()

    return (k_eff,k_con)

# now run the simulations

# print exchange(10,0)

# gets model rates for Zhang & Winfree 2009 figure 3B.
k3way_exp = [1.40, 8.17, 144, 1.08e3, 5.05e4, 9.64e5, 2.36e6, 3.22e6, 3.15e6, 2.77e6, 2.83e6, 4.78e6]
k3way_exp = zip( range(0,11)+[15], k3way_exp )

###### set this true to run the simulations
if do3way:
    k3way = [ (n,exchange(n,0)) for n in range(1,16) ]

    print "Toehold-mediate strand displacement rate constants, c.f. Zhang & Winfree 2009, figure 3B. TE+12.5mM Mg++ at 25C.  n=toehold length."
    print " toehold :  detailed+algebra  :     condensed     :    experimental"
    i=0
    j=0
    while i<len(k3way) and j<len(k3way_exp):
        if k3way[i][0]==k3way_exp[j][0]:
            print "  n = %2d : k_eff = %10.4g : k_con = %10.4g : k_exp = %10.4g" % (k3way[i][0],k3way[i][1][0],k3way[i][1][1],k3way_exp[j][1])
            i=i+1
            j=j+1
        elif k3way[i][0]<k3way_exp[j][0]:
            print "  n = %2d : k_eff = %10.4g : k_con = %10.4g : k_exp = %10s" % (k3way[i][0],k3way[i][1][0],k3way[i][1][1],'None')
            i=i+1
        else:
            print "  n = %2d : k_eff = %10s : k_con = %10s : k_exp = %10.4g" % (k3way_exp[j][0],'None','None',k3way_exp[j][1])
            j=j+1
    raw_input("Press enter to continue...")


# now, estimate values from figure 4B  (could ask Dave Zhang for more accurate numbers)
logk3wayx_exp_visual_estimate = [(1,4,.95),(1,3,.8),(1,2,1.15),(1,1,1.1), \
  (2,5,1.8),(2,4,2.15),(2,3,2.2),(2,2,2.2),(2,1,2.15), \
  (3,6,1.9),(3,5,2.15),(3,4,3.0),(3,3,3.0),(3,2,2.95),(3,1,2.9), \
  (4,7,2.05),(4,6,2.65),(4,5,3.65),(4,4,4.1),(4,3,4.1),(4,2,4.1),(4,1,4.05), \
  (5,7,3.7),(5,6,4.9),(5,5,5.7),(5,4,6.15),(5,3,6.15),(5,2,6.15),(5,1,6.15), \
  (6,7,5.1),(6,6,5.8),(6,5,6.2),(6,4,6.4),(6,3,6.2),(6,2,6.2)]
k3wayx_exp_visual_estimate = [ (n,m,10**v) for (n,m,v) in logk3wayx_exp_visual_estimate ]

# from Dave Zhang's matlab script, ToeEx_rates_vs_model.m
k3wayx_exp = [(1,4,7.70),(1,3,5.48),(1,2,23.5),(1,1,18.9), \
  (2,5,43.6),(2,4,214.05),(2,3,273.0),(2,2,249.0),(2,1,231.0), \
  (3,6,66.9),(3,5,215.0),(3,4,939.0),(3,3,974.0),(3,2,907.0),(3,1,846.0), \
  (4,7,131.0),(4,6,407.0),(4,5,4.25e3),(4,4,2.13e4),(4,3,2.41e4),(4,2,2.29e4),(4,1,1.97e4), \
  (5,7,3.59e3),(5,6,9.72e4),(5,5,3.45e5),(5,4,1.53e6),(5,3,1.58e6),(5,2,1.58e6),(5,1,1.73e6), \
  (6,7,1.61e5),(6,6,4.05e5),(6,5,1.48e6),(6,4,3.04e6),(6,3,2.59e6),(6,2,3.00e6), \
  (7,7,4.7e5),(7,6,1.11e6),(7,5,2.90e6),(7,4,3.57e6), \
  (8,7,1.94e6),(8,6,2.68e6),(8,5,3.14e6),(8,4,3.37e6)        ]

###### set this true to run the simulations
if do3way_exchange:
    k3wayx = [ (n,m, exchange(n,m)) for (n,m,v) in k3wayx_exp ]

    print "Toehold exchange rate constants, c.f. Zhang & Winfree 2009, figure 4B.  TE+12.5mM Mg++ at 25C.  n=incoming, m=incumbent. "
    print "invading,incumbent : detailed+algebra  :      condensed      :     experimental"
    for (model,exp) in zip(k3wayx,k3wayx_exp):
        (n,m,(k_eff,k_con))=model
        (n_exp,m_exp,k_exp)=exp
        assert n==n_exp and m==m_exp
        print "   n=%2d, m=%2d      :  k_eff=%10.4g :  k_con = %10.4g :  k_exp = %10.4g" % (n,m,k_eff,k_con,k_exp)
    raw_input("Press enter to continue...")

def fourway(n,m):
    """run enumerator for Dabby 2013, tables 5.1 and 5.2"""

    assert n<7 and (m<7 or m==16)
    M = 0 if (m==6 or m==16) else 6-m
    N = 6-n
    
    if n==0 and m==0:
        return (0,0)   # Leak not modeled, don't bother.  detailed and condensed semantics both give 0.
    if n==15:
        assert m==0
        sys = "length a = 9\nlength B = 9\nlength c = %d\nB c\na B( + c* )\n" % n
    if m==0:
        sys = "length a = 9\nlength B = 9\nlength c = %d\nlength d = 9\nB c\na B( + d* c* )\n" % n
    if m>0:
        sys = "length x = 21\nlength b = %d\nlength B = 9\nlength c = %d\nlength d = 9\nB c\na b(B( + d* c* ))\n" % (m,n)

    sys = "length x = 21\n"
    if m>0: 
        sys += "length m = %d\n" % m
    if M>0: 
        sys += "length M = %d\n" % M
    if n>0: 
        sys += "length n = %d\n" % n
    if N>0: 
        sys += "length N = %d\n" % N
    if m>0 and M>0:
        sys += "x*( m* M* + "
    elif m>0:
        sys += "x*( m* + "
    else:
        sys += "x*( M* + "
    if n>0 and N>0:
        sys += "N* n* )\n"
    elif n>0:
        sys += "n* )\n"
    else:
        sys += "N* )\n"
    if m>0:
        sys += "m x( + "
    else:
        sys += "x( + "
    if n>0:
        sys += ") n\n"
    else:
        sys += ")\n"

    pil_enum = CMI_enum(sys,8,'detailed')
    rates = [s for s in pil_enum if len(s)>0 and s[0]=='k']
    if watch:
        for s in pil_enum:
            print s

    # trust that the enumerator always lists reactions in a consistent order!
    if len(rates)==1:  # must be condensed, then, so must not actually ever happen here.
        k_eff = float(rates[0].split()[1][1:])
    elif len(rates)==3:  # must be reversible toehold, detailed model  (i.e. just one toehold)
        k0=float(rates[0].split()[1][1:])   # forward binding rate
        k1=float(rates[1].split()[1][1:])   # branch migration & strand displacement step
        k2=float(rates[2].split()[1][1:])   # toehold dissociation
        k_eff = k0*k1/(k1+k2)
    elif len(rates)==13: # both toeholds bind reversibly, and hilarity ensues...
        k = [ float(r.split()[1][1:]) for r in rates ]
        # reactions come out in one of two possible orders, due to mysterious reasons...
        if rates[0].find("18 -> 23") != -1:
            assert rates[1].find("6 -> 10") != -1
            assert rates[2].find("5 -> 10") != -1
            assert rates[3].find("44 -> 19") != -1
            assert rates[4].find("2 + 1 -> 6") != -1
            assert rates[5].find("2 + 1 -> 5") != -1
            assert rates[6].find("6 -> 18 + 19") != -1
            assert rates[7].find("5 -> 23 + 44") != -1
            assert rates[8].find("10 -> 23 + 19") != -1
            assert rates[9].find("6 -> 2 + 1") != -1
            assert rates[10].find("5 -> 2 + 1") != -1
            assert rates[11].find("10 -> 6") != -1
            assert rates[12].find("10 -> 5") != -1
            b10 = k[8]/(k[8]+k[11]+k[12])
            y6  = k[11]/(k[8]+k[11]+k[12])
            y5  = k[12]/(k[8]+k[11]+k[12])
            a5  = k[2]/(k[10]+k[7]+k[2])
            b5  = k[7]/(k[10]+k[7]+k[2])
            a6  = k[1]/(k[9]+k[6]+k[1])
            b6  = k[6]/(k[9]+k[6]+k[1])
            P10 = (b10+y6*b6+y5*b5)/(1-(y6*a6+y5*a5))
            k_eff = (k[4]*a6+k[5]*a5)*P10 + (k[4]*b6+k[5]*b5)
        elif rates[0].find("5 -> 10") != -1:
            assert rates[1].find("6 -> 10") != -1
            assert rates[2].find("43 -> 18") != -1
            assert rates[3].find("19 -> 23") != -1
            assert rates[4].find("1 + 2 -> 5") != -1
            assert rates[5].find("1 + 2 -> 6") != -1
            assert rates[6].find("5 -> 43 + 23") != -1
            assert rates[7].find("6 -> 18 + 19") != -1
            assert rates[8].find("10 -> 18 + 23") != -1
            assert rates[9].find("5 -> 1 + 2") != -1
            assert rates[10].find("6 -> 1 + 2") != -1
            assert rates[11].find("10 -> 5") != -1
            assert rates[12].find("10 -> 6") != -1
            b10 = k[8]/(k[8]+k[11]+k[12])
            y6  = k[12]/(k[8]+k[11]+k[12])
            y5  = k[11]/(k[8]+k[11]+k[12])
            a5  = k[0]/(k[9]+k[6]+k[0])
            b5  = k[6]/(k[9]+k[6]+k[0])
            a6  = k[1]/(k[10]+k[7]+k[1])
            b6  = k[7]/(k[10]+k[7]+k[1])
            P10 = (b10+y6*b6+y5*b5)/(1-(y6*a6+y5*a5))
            k_eff = (k[5]*a6+k[4]*a5)*P10 + (k[5]*b6+k[4]*b5)
        else:
            print "CRAP.  Reactions are coming out in an unexpected order."

    pil_enum = CMI_enum(sys,8,'condensed')
    rates = [s for s in pil_enum if len(s)>0 and s[0]=='k']

    if len(rates)==1:  # irreversibble toehold-mediated strand displacemen
        k_con = float(rates[0].split()[1][1:])
    elif len(rates)==2:  # reversible toehold exchange   #### check by hand to make sure first one is always forward
        k_con = float(rates[0].split()[1][1:])   # forward binding rate
    # now must modify stuff below to output & compare condensed rates

    if watch:
        for s in pil_enum:
            print s
        print "Calculated k_eff = %f /M/s from detailed reactions and k_con = %f /M/s from condensed reactions." % (k_eff,k_con)
        raw_input("Press enter to continue...")  # in python 3, just input()

    return (k_eff,k_con)


# (n,m,k1_fit) from and Nadine Dabby, Caltech PhD Thesis, Table 5.2  (note m,n have consistent meaning, but order in table is swapped.)
k4way_exp = [ (0,0,0.034),(0,2,0.047),(2,2,0.10),(2,0,0.033),(4,2,0.93),(4,0,0.039),(0,4,0.97), \
              (2,4,56),(6,2,490),(0,6,58),(4,4,770),(6,0,5.0),(2,6,9.4e3),(4,6,7.0e4),(6,4,2.8e5),(6,6,6.9e5) ]

###### set this true to run the simulations
if do4way:
    k4way = [ (n,m, fourway(n,m)) for (n,m,v) in k4way_exp ]

    print "Toehold-mediated 4-way rate constants, c.f. Dabby's PhD thesis, table 5.2.  TAE+12.5mM Mg++ at 25C."
    print "  toehold lengths  : detailed+algebra  :      condensed      :     experimental"
    for (model,exp) in zip(k4way,k4way_exp):
        (n,m,(k_eff,k_con))=model
        (n_exp,m_exp,k_exp)=exp
        assert n==n_exp and m==m_exp
        print "   n=%2d, m=%2d      :  k_eff=%10.4g :  k_con = %10.4g :  k_exp = %10.4g" % (n,m,k_eff,k_con,k_exp)
    raw_input("Press enter to continue...")

