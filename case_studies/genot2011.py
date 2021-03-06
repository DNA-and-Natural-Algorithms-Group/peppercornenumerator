
def genot2011_F3_pil(mn) :
    (m,n) = mn
    return """
    length d1 = 15
    length d2 = 22
    length d3 = 14
    
    length d05 = {}
    length d04 = {}

    S = d1 d2( + d3* d05* )
    I = d2 d04 d3

    T = d2( d04 d3( + ) d05* )
    W = d1 d2
    """.format(m, n)

def genot2011_F4A_pil(_) :
    return """
    length d1 = 15
    length d2 = 22
    length d3 = {}
    
    S = d1 d2( + d3* )
    I = d2 d3

    T = d2( d3( + ) )
    W = d1 d2
    """.format(_)

def genot2011_F4B_pil(_) :
    return """
    length d1 = 15
    length d2 = 22
    length d8 = 1
    length d3 = {}
    
    S = d1 d2( + d3* d8 )
    I = d2 d3

    T = d2( d3( + ) d8 )
    W = d1 d2
    """.format(_)

def genot2011_F4C_pil(_) :
    return """
    length d1 = 15
    length d2 = 22
    length d3 = 6 # original 
    #length d3 = 4 # try-to-fit
    
    S = d1 d2( + d3* )
    I = d2 d3

    T = d2( d3( + ) )
    W = d1 d2
    """

def genot2011_F4D_pil(_) :
    return """
    length d1 = 15
    length d2 = 22
    length d3 = 14 # original
    #length d3 = 7 # try-to-fit
    length d4 = 17
    length d5 = 17
    
    S = d1 d2( + d3* d5 )  # d3: 14
    I = d2 d4 d3

    T = d2( d4 d3( + ) d5 )
    W = d1 d2
    """

def genot2011_SF4A_pil(dl) :
    return """
    \rlength d1 = 15
    \rlength d2a = 5
    \rlength d2b = 7
    \rlength d2c = 10
    \rlength d3 = 14 # n
    \rlength d4 = {} # m
    \rlength d5 = {} # m

    \rsup-sequence d2 = d2a d2b d2c
    \r
    \rS = d1 d2( + d3* d5 )
    \rX = d2 d4 d3

    \rW = d2( d4 d3( + ) d5 )
    \rP = d1 d2

    \rR = d1( d2a( + d2b* ) )

    \rF = d1 d2a
    \rFW = d1( d2a( d2b( d2c + ) ) )
    """.format(dl, dl)

def genot2011_SF5_pil(dl) :
    return """
    \rlength d1 = 15
    \rlength d2a = 5
    \rlength d2b = 7
    \rlength d2c = 10
    \rlength d3 = {} # n
    \rlength d4 = 20 # m
    \rlength d5 = 20 # m

    \rsup-sequence d2 = d2a d2b d2c
    \r
    \rS = d1 d2( + d3* d5 )
    \rX = d2 d4 d3

    \rW = d2( d4 d3( + ) d5 ) 
    \rP = d1 d2

    \rR = d1( d2a( + d2b* ) )

    \rF = d1 d2a
    \rFW = d1( d2a( d2b( d2c + ) ) )
    """.format(dl)


def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData
    psim = "pilsimulator --no-jacobian --header --nxy --atol 1e-13 --rtol 1e-13 --mxstep 10000 "
    psim += "--t-lin 10000 --t8 1000"
    rates = {'k_slow': 1e-10, 'k_fast': 0.01}

    G11_F3 = FigureData('Fig. 3: |a|=14, |b|=22, varying: |n|, |m|')
    G11_F3.fname = 'Genot2011-F3'
    current = G11_F3
    template = genot2011_F3_pil

    pilp = [(8,7), (17,17), (20,20), (23,23)] # original
    # After 30% of the reaction has gone to completion, the signal dropped to 70%. 
    litr = [(42, 1.98), (48, 1.98), (100, 1.98), (91, 1.98)] # seconds

    for (pil, res) in zip(pilp, litr):
        pilstring  = template(pil)
        simulation = psim + " --pyplot-labels S I T W --p0 S=6.6 I=660"
        reporter = 'T'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '1.98', res)#, simargs='pilname')

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)


    G11_F4A = FigureData('Fig. 4A: |n|=|m|=0, either |a|=11 or |a|=9')
    G11_F4A.fname = 'Genot2011-F4A'
    current = G11_F4A
    template = genot2011_F4A_pil

    pilp = [11, 9]
    #litr = [(32, 0.58*6.6), (38, 0.58*6.6)]
    litr = [(2, 1.98), (3, 1.98)]

    for (pil, res) in zip(pilp, litr):
        pilstring  = template(pil)
        simulation = psim + " --pyplot-labels S I T W --p0 S=6.6 I=22"
        reporter = 'T'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '1.98', res)#, simargs='pilname')

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)


    G11_F4B = FigureData('Fig. 4B: |n|=|m|=1, either |a|=11 or |a|=9')
    G11_F4B.fname = 'Genot2011-F4B'
    current = G11_F4B
    template = genot2011_F4B_pil

    pilp = [11, 9]
    #litr = [(44, 3.3), (182, 3.3)]
    litr = [(15, 1.98), (75, 1.98)]

    for (pil, res) in zip(pilp, litr):
        pilstring  = template(pil)
        simulation = psim + " --pyplot-labels S I T W --p0 S=6.6 I=22"
        reporter = 'T'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '1.98', res)#, simargs='pilname')

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    G11_F4C = FigureData('Fig. 4C: |n|=|m|=0, |a|=6, varying: [invader]$_0$')
    G11_F4C.fname = 'Genot2011-F4C'
    current = G11_F4C
    template = genot2011_F4C_pil

    sims = [psim + " --pyplot-labels S I T W --p0 S=6.6 I=330",
            psim + " --pyplot-labels S I T W --p0 S=6.6 I=145",
            psim + " --pyplot-labels S I T W --p0 S=6.6 I=66"]
    #litr = [(69, 3.3), (177, 3.3), (347, 3.3)]
    litr = [(25, 1.98), (75, 1.98), (150, 1.98)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'T'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '1.98', res)#, simargs=sim[sim.find('I='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    G11_F4D = FigureData('Fig. 4D: |n|=|m|=17, |a|=14, varying: [invader]$_0$')
    G11_F4D.fname = 'Genot2011-F4D'
    current = G11_F4D
    template = genot2011_F4D_pil
    sims = [psim + " --pyplot-labels S I T W --p0 S=6.6 I=330",
            psim + " --pyplot-labels S I T W --p0 S=6.6 I=145",
            psim + " --pyplot-labels S I T W --p0 S=6.6 I=66"]
    #litr = [(176, 3.3), (213, 3.3), (248, 3.3)]
    litr =  [(51, 1.98), (55, 1.98), (60, 1.98)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'T'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '1.98', res)#, simargs=sim[sim.find('I='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    G11_SF4A = FigureData('Sup. Fig. 4A: |a|=14, |b|=22 varying: |n|,|m|')
    G11_SF4A.fname = 'Genot2011-SF4A'
    current = G11_SF4A
    template = genot2011_SF4A_pil

    pilp = [17, 20, 23]
    #litr = [(9088, 5), (13927, 5), (28681, 5)]
    litr = [(4249, 3), (6492, 3), (12157,3)]

    psim = "pilsimulator --nxy --header --atol 1e-12 --rtol 1e-12 --mxstep 10000 --t-lin 100000 --t8 1e5"
    for (pil, res) in zip(pilp, litr):
        pilstring  = template(pil)
        simulation = psim + " --pyplot-labels S X R F --p0 R=30 X=40 S=10"
        reporter = 'F'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '3', res)

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)
    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    G11_SF5 = FigureData('Sup. Fig. 5: |m|=|n|=20, varying: |a|')
    G11_SF5.fname = 'Genot2011-SF5'
    current = G11_SF5
    template = genot2011_SF5_pil

    pilp = [14, 10, 8, 6, 4, 2]
    litr = [(25382, 6.94), 
            (35273, 5.76), 
            (68655, 1.87), 
            (77018, 0.91), 
            (76218, 0.99), 
            (77818, 0.79)]

    psim = "pilsimulator --nxy --header --atol 1e-12 --rtol 1e-12 --mxstep 10000 --t-lin 100000 --t8 1e5"
    for (pil, res) in zip(pilp, litr):
        pilstring  = template(pil)
        simulation = psim + " --pyplot-labels S X R F --p0 R=30 X=40 S=10"
        reporter = 'F'
        metric = 'diagonal-crossing-time'
        tmax = '86400'
        cmax = '10'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric,
                ':'.join([tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['detailed'].copy()
    current.pepperargs['default'].update({'k_slow': 1e-10, 'k_fast': 1e-5})
    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)


    return [G11_F3, G11_F4A, G11_F4B, G11_F4C, G11_F4D, G11_SF4A]#, G11_SF5]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)

