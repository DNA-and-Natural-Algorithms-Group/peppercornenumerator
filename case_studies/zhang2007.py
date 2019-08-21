
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
    F1E = FigureData('Zhang et al. (2007) Fig. 1 - Single-layer catalytic DSD system (varying catalyst)')
    F1E.fname = 'Zhang2007-F1'
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

    # DIAGONAL CROSSING
    diagX = [
        (1216.44, 7.44),
        (1660.27, 6.83),
        (2457.53, 5.78),
        (3435.62, 4.49),
        (4364.38, 3.26),
        (5424.66, 1.88),
        (6073.97, 1.03),
        (6386.30, 0.57),
        (6550.68, 0.32)]
        #(6715.07, 0.09) 0 input
    
    for (sim, res) in zip(sims, diagX):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '7200'
        cmax = '10'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, ':'.join([tmax, cmax]), res, simargs=sim[sim.find('C='):])

    # Half completion
    halfC = [
        (550.68, 5.00),
        (953.42, 5.00),
        (1906.85, 5.00),
        (4076.71, 5.00)]

    for (sim, res) in zip(sims, halfC):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '5', res, simargs=sim[sim.find('C='):])

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F3 = FigureData('Zhang et al. (2007) Fig. 3 - Two-layer feedforward DSD system (varying catalyst)')
    F3.fname = 'Zhang2007-F3'
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
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=0.02',
            psim + ' --pyplot-labels ROX S0 S1 F0 F1 OR C0 --p0 S1=10 F1=13 S0=10 F0=13 OR=30 C0=0.01']

    #Diagonal Dataset
    diagX = [
        (1860.53, 6.56),
        (2302.30, 5.98),
        (2744.07, 5.40),
        (3066.90, 4.98),
        (3856.99, 3.92),
        (4519.65, 3.07),
        (4910.44, 2.54),
        (5607.08, 1.63),
        (5904.42, 1.23),
        (6065.84, 1.02)]
        #(6320.71, 0.66), # 0x
        #(6720.00, 0.15)] # control

    for (sim, res) in zip(sims, diagX):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'diagonal-crossing-time'
        tmax = '7200'
        cmax = '10'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, ':'.join([tmax, cmax]), res, simargs=sim[sim.find('C0='):])


    # Half completion
    halfC = [
        (1257.35, 5.00), 
        (1852.04, 5.00),
        (2514.69, 5.00),
        (3066.90, 5.00),
        (4850.97, 5.00)]

    for (sim, res) in zip(sims, halfC):
        pilstring  = template(None)
        simulation = sim
        reporter = 'ROX'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '5', res, simargs=sim[sim.find('C0='):])

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F4 = FigureData('Zhang et al. (2007) Fig. 4 - Autocatalytic DSD system (varying autocatalyst)')
    F4.fname = 'Zhang2007-F4'
    current = F4
    template = zhang2007_F4_pil
    sims = [psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=10',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=7',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=5',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=3',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=2',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=1',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.7',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.5',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.3',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.2',
            psim + ' --pyplot-labels TET S F SR A --p0 S=10 F=13 SR=20 A=0.1']

    #Diagonal Dataset
    diagX = [
        (496.08, 7.62), 
        (600.76, 7.32),
        (655.37, 7.16),
        (855.63, 6.61),
        (910.24, 6.45),
        (1064.98, 6.04), 
        (1206.07, 5.64),
        (1237.93, 5.54),
        (1333.50, 5.26),
        (1392.67, 5.10),
        (1465.49, 4.91)]
        #(1756.76, 4.10)]

    for (sim, res) in zip(sims, diagX):
        pilstring  = template(None)
        simulation = sim
        reporter = 'TET'
        metric = 'diagonal-crossing-time'
        tmax = '3600'
        cmax = '10'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, ':'.join([tmax, cmax]), res, simargs=sim[sim.find('A='):])

    #Half completion
    halfC = [
        (241.21, 5.0),
        (350.44, 5.0), 
        (393.68, 5.0),
        (609.86, 5.0),
        (694.06, 5.0),
        (892.04, 5.0),
        (1087.74, 5.0), 
        (1121.87, 5.0),
        (1276.61, 5.0),
        (1376.74, 5.0),
        (1470.04, 5.0)]
        #(1966.12, 5.0)]
 
    for (sim, res) in zip(sims, halfC):
        pilstring  = template(None)
        simulation = sim
        reporter = 'TET'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '5', res, simargs=sim[sim.find('A='):])

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)


    return [F1E, F3, F4]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)


