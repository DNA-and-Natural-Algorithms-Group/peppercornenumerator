
def zhang2007_F1_pil(x):
    return """
    # Domains
    length d1  = 10
    length d2a =  6
    length d2b =  6
    length d2c = 12
    length d3  =  4
    length d4  = 16 
    length d5  =  6
    length d6  = 16

    # Species
    C = d4 d5
    OB = d1 d2a d2b d2c 
    ROX = d1 d2a       
    S = d1 d2a( d2b( d2c( + d6 d3( d4( + d5* ) ) ) ) ) 
    F = d2a d2b d2c d3 d4 # fuel
    OR = d1( d2a( + d2b* ) )
    SB = d6 d3 d4          

    I1 = d1 d2a( d2b( d2c( + d6 d3( d4( + d4 d5( + ) ) ) ) ) ) 
    I2 = d1 d2a( d2b( d2c( + d6 d3( d4 + d4( d5( + ) ) ) ) ) ) 
    I3 = d1 d2a( d2b( d2c( + d4( d5( + ) ) d3* ) ) ) 
    I4 = d1 d2a( d2b( d2c( + d2a d2b d2c d3( d4 + d4( d5( + ) ) ) ) ) ) 
    I4a = d1 d2a( d2b( d2c( + d2a d2b d2c d3( d4( + d4 d5( + ) ) ) ) ) ) 
    I5 = d2a( d2b( d2c( d3( d4 + d4( d5( + ) ) ) ) ) ) 
    W = d2a( d2b( d2c( d3( d4( + d5* ) ) ) ) )
    RQ = d1( d2a( d2b( d2c + ) ) )
    """.format(x)

def zhang2007_F3_pil(x):
    return """
    # Domains
    length d1  = 10
    length d2a =  6
    length d2b =  6
    length d2c = 12
    length d3  =  4
    length d4a = 10
    length d4b =  6 # = 2a
    length d5  =  6 # = 2b
    length d6  = 16
    length d7  = 12
    length d8  =  4
    length d9  = 16
    length d10 =  6 # = 2b

    # Species
    F0 = d4b d5 d7 d8 d9
    S0 = d4a d4b( d5( d7( + d8( d9( + d10* ) ) ) ) )
    OB0 = d4a d4b d5 d7
    C0 = d9 d10

    F1 = d2a d2b d2c d3 d4a d4b
    S1 = d1 d2a( d2b( d2c( + d6 d3( d4a( d4b( + d5* ) ) ) ) ) )
    OB1 = d1 d2a d2b d2c

    OR = d1( d2a( + d2b* ) )
    ROX = d1 d2a       
    RQ = d1( d2a( d2b( d2c + ) ) )
    """.format(x)

def zhang2007_F4_pil(x):
    return """
    # Domains
    length d2b =  6
    length d2c = 12
    length d3  =  4
    length d4t =  7
    length d4a =  3
    length d4b =  6 
    length d6  = 16

    # Species
    A = d4t d4a d4b d2b d2c
    SB = d6 d3 d4t d4a d4b
    TET = d6 d3
    I = d4t d4a d4b( d2b( d2c( + d4t( d4a( d4b( d2b( d2c + ) ) ) ) d3* ) ) )

    S = d4t d4a d4b( d2b( d2c( + d6 d3( d4t( d4a( d4b( + d2b* ) ) ) ) ) ) )
    SR = d6( d3( + d4t* ) )
    F = d4b d2b d2c d3 d4t d4a d4b

    W = d4b( d2b( d2c( d3( d4t( d4a( d4b( + d2b* ) ) ) ) ) ) )

    FQ = d6( d3( d4t( d4a d4b +  ) ) )
    """.format(x)

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    # Default pilsimulator call
    psim = "pilsimulator --nxy --header --atol 1e-13 --rtol 1e-13 --mxstep 10000"
    psim += " --t8 18000 --t-lin 18000"

    # Setup
    F1E = FigureData('Zhang2007-F1')
    current = F1E
    template = zhang2007_F1_pil
    sims = [psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=10',
            psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=5',
            psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=2',
            psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=1',
            psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=0.5',
            psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=0.2',
            psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=0.1',
            psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=0.05',
            psim + ' --pyplot-labels ROX S F OR C --p0 S=10 F=13 OR=30 C=0.02']
    litr = [(1206, 7.42), (1662, 6.83), (2450, 5.79), (3420, 4.47), (4350, 3.25), (5414, 1.84), (6070,1.0), (6383,0.58), (6602,0.27)] # seconds

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '7200'
        cmax = '10'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C='):])

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F3 = FigureData('Zhang2007-F3')
    current = F3
    template = zhang2007_F3_pil
    sims = [psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=10',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=5',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=2',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=1',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=0.5',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=0.2',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=0.1',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=0.05',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=0.02']
    litr = [(1235, 5.0), (1813, 5.0), (2708, 5.5), (3041, 5.0), (3837, 4.0), (4513, 3.1), (4901, 2.6), (5605, 1.7), (5901, 1.2), (6030, 1.1)] 

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '7200'
        cmax = '10'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C0='):])

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F4 = FigureData('Zhang2007-F4')
    current = F4
    template = zhang2007_F4_pil
    sims = [psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=10',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=7',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=5',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=3',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=1',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.7',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.5',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.3',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.2',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.1']
    litr = [(240, 5), (342, 5), (394, 5), (609, 5), (688, 5), (884, 5), (1080, 5), (1122, 5), (1273, 5), (1369, 5), (1473, 5)] # seconds

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'TET'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, res, simargs=sim[sim.find('A='):])

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)


    return [F1E, F3, F4]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)


