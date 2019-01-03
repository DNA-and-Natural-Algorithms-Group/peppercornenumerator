
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

def zhang2010_F10A_pil(x):
    return """
    # Domains
    length d41 = 16
    length d42a = 4
    length d42b = 7
    length d42c = 9
    length d43  = 5
    length d44 = 15
    length d45 =  8
    length d46 = 15

    F4 = d42a d42b d42c d43 d44
    S4 = d41 d42a( d42b( d42c( + d46 d43( d44( + d45* ) ) ) ) )
    C4 = d44 d45

    W4 = d42a( d42b( d42c( d43( d44( + d45* ) ) ) ) )
    OP4 = d41 d42a d42b d42c
    SP4 = d46 d43 d44

    RQ = d41( d42a( d42b( d42c + ) ) )
    ROX = d41 d42a
    OR4 = d41( d42a( + d42b* ) )
    """.format(x)

def zhang2010_F10D_pil(x):
    return """
    # Domains
    length d51  = 15
    length d52a =  5
    length d52b =  7
    length d52c = 10
    length d63 =  5
    length d64 = 11
    length d65 =  8
    length d66 = 10

    F43 = d52a d52b d52c d63 d64
    S43 = d51 d52a( d52b( d52c( + d66 d63( d64( + d65* ) ) ) ) ) 
    C43 = d64 d65
    W43 = d52a( d52b( d52c( d63( d64( + d65* ) ) ) ) )
    OP43 = d51 d52a d52b d52c
    SP43 = d66 d63 d64
    RQ = d51( d52a( d52b( d52c + ) ) )
    RG = d51 d52a
    OR43 = d51( d52a( + d52b* ) )
    """.format(x)

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    # Default pilsimulator call
    psim = "pilsimulator --no-jacobian --nxy --header"
    psim += " --atol 1e-10 --rtol 1e-10 --mxstep 1000"
    psim += " --t8 86500 --t-lin 86500"

    rates = {'k_slow': 1e-5, 'k_fast': 0.1}

    # Setup
    F3A = FigureData('Zhang2010-F3A')
    current = F3A
    template = zhang2007_F1_pil
    sims = [psim + ' --pyplot-labels S F OR C ROX --p0 S=100 F=200 OR=300 C=10',
            psim + ' --pyplot-labels S F OR C ROX --p0 S=100 F=200 OR=300 C=1']
    litr = [(2235, 86.25), (11389, 43.93)] 

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '21600'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C='):])

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F3B = FigureData('Zhang2010-F3B')
    current = F3B
    template = zhang2007_F1_pil
    sims = [psim + ' --pyplot-labels S F OR C ROX --p0 S=30 F=60 OR=90 C=3',
            psim + ' --pyplot-labels S F OR C ROX --p0 S=30 F=60 OR=90 C=0.9']
    litr = [(3887, 23.70), (8687, 17.04)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '21600'
        cmax = '30'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C='):])

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F3C = FigureData('Zhang2010-F3C')
    current = F3C
    template = zhang2007_F1_pil
    sims = [psim + ' --pyplot-labels S F OR C ROX --p0 S=3 F=6 OR=9 C=0.9',
            psim + ' --pyplot-labels S F OR C ROX --p0 S=3 F=6 OR=9 C=0.3']
    litr = [(9600, 2.24), (16066, 1.76)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '43200'
        cmax = '3'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F3D = FigureData('Zhang2010-F3D')
    current = F3D
    template = zhang2007_F1_pil
    sims = [psim + ' --pyplot-labels S F OR C ROX --p0 S=1 F=2 OR=3 C=1',
            psim + ' --pyplot-labels S F OR C ROX --p0 S=1 F=2 OR=3 C=0.1']
    litr = [(13332, 0.82), (37775, 0.51)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '86400'
        cmax = '1'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F10C = FigureData('Zhang2010-F10C')
    current = F10C
    template = zhang2010_F10A_pil
    sims = [psim + ' --pyplot-labels S4 F4 OR4 C4 ROX --p0 S4=30 F4=60 OR4=90 C4=3',
            psim + ' --pyplot-labels S4 F4 OR4 C4 ROX --p0 S4=30 F4=60 OR4=90 C4=0.09',
            psim + ' --pyplot-labels S4 F4 OR4 C4 ROX --p0 S4=30 F4=60 OR4=90 C4=0.03']
    litr = [(14475, 11.55), (24750, 5.88), (30425, 2.75)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '36000'
        cmax = '20'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C4='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F10F = FigureData('Zhang2010-F10F')
    current = F10F
    template = zhang2010_F10D_pil
    sims = [psim + ' --pyplot-labels S43 F43 OR43 C43 RG --p0 S43=10 F43=20 OR43=30 C43=3',
            psim + ' --pyplot-labels S43 F43 OR43 C43 RG --p0 S43=10 F43=20 OR43=30 C43=0.09',
            psim + ' --pyplot-labels S43 F43 OR43 C43 RG --p0 S43=10 F43=20 OR43=30 C43=0.03']
    litr = [(3358, 3.98), (5139, 2.04), (6075, 1.04)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'RG'
        metric = 'diagonal-crossing-time'
        tmax = '7200'
        cmax = '8'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C43='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    return [F3A, F3B, F3C, F3D, F10C, F10F]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)



