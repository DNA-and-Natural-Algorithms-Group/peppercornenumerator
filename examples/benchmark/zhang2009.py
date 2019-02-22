def zhang2009_F1DF_exchange_pil((m,n)) :
    assert m != 0
    assert n != 0
    return """
    # Experiments Zhang 2009, 3-way strand displacement toehold exchange

    # For the reporter, assuming that beta^6 is actually beta^m, and then there
    # is a 7nt reporter toehold section (br) after bt...  not pretty, but
    # should work just fine...

    length a  = 16
    length bt = {:d} # m
    length br = 7    # 7
    length b  = {:d} # 13 - m
    length ct = {:d} # n
    length c  = {:d} # 15 - n
    
    X = br b ct
    S = a bt( br( b( + c* ct* ) ) )
    
    I = a bt( br( b( + br b ct( + c* ) ) ) )
    J = a bt( br b + br( b( ct( + c* ) ) ) )
    
    Y = a bt br b
    L = br( b( ct( + c* ) ) ) bt*


    R = a( bt( + br* ) )
    V = a( bt( br( b + ) ) )
    F = a bt
    """.format(m, 13-m, n, 15-n)

def zhang2009_F1DF_displacement_pil((n)) :
    assert n != 0
    return """
    # Toehold mediated strand displacement (no exchange)

    length a  = 16
    length b6 = 6
    length br = 7 
    length b  = 7
    length ct = {:d} # n
    length c  = {:d} # 15 - n
    
    X = b6 br b ct
    S = a b6( br( b( + c* ct* ) ) )
    
    I = a b6( br( b( + br b ct( + c* ) ) ) )
    J = a b6( br b + br( b( ct( + c* ) ) ) )
    
    Y = a b6 br b
    L = b6( br( b( ct( + c* ) ) ) )


    R = a( b6( + br* ) )
    V = a( b6( br( b + ) ) )
    F = a b6
    """.format(n, 15-n)

def zhang2009_F5_pil((n)):
    return """
    # Domains
    length a = 16
    length bm = 5  # m
    length br = 1  # m
    length bc = 7 
    length bt  = 7
    length n = {}  # n
    length c = {}  # 15-n
    
    # Species
    S = a bm( br( bc( bt( + c* n* ) ) ) )
    X = br bc bt n
    L = br( bc( bt( n( + c* ) ) ) ) bm*
    Y = a bm br bc bt
    Z = bm br bc bt
    W = bm( br( bc( bt( + c* n* ) ) ) )

    R = a( bm( br( + bc* ) ) )
    YW = a( bm( br( bc( bt + ) ) ) )
    F = a bm br
    """.format(n, 15-n)


def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    # Default pilsimulator call
    psim = "pilsimulator --no-jacobian --nxy --header"
    psim += " --atol 1e-12 --rtol 1e-12 --mxstep 1000"
    psim += " --t8 18000 --t-lin 18000"

    rates = {'k_slow': 1e-5, 'k_fast': 0.1}
    rc = {'release_cutoff': 15}

    # Setup
    F3 = FigureData('Zhang & Winfree (2009) Fig. 3 - Single strand displacement reactions (varying toehold length)')
    F3.fname = 'Zhang2009-F3'
    current = F3
    template = zhang2009_F1DF_displacement_pil
    sims = [psim + ' --pyplot-labels F R S X --p0 R=3 S=1 X=0.6',
            psim + ' --pyplot-labels F R S X --p0 R=3 S=1 X=0.4',
            psim + ' --pyplot-labels F R S X --p0 R=3 S=1 X=0.2']
    litr = [(1260, 0.31), (1867, 0.27), (3498, 0.17)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(5)
        simulation = sim
        reporter = 'F'
        metric = 'diagonal-crossing-time'
        tmax = '7200'
        cmax = '0.4'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('X='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F4 = FigureData('Zhang & Winfree (2009) Fig. 4 - Single toehold exchange reactions (varying toehold lengths)')
    F4.fname = 'Zhang2009-F4'
    current = F4
    template = zhang2009_F1DF_exchange_pil
    pilp = [(4,7), (5,7), (6,7), (7,7)]
    litr = [(785, 0.28), (823, 0.27), (1192, 0.23), (1622, 0.18)]

    for (pip, res) in zip(pilp, litr):
        pilstring  = template(pip)
        simulation = psim + ' --pyplot-labels F R S X --p0 R=3 S=1 X=0.4'
        reporter = 'F'
        metric = 'diagonal-crossing-time'
        tmax = '3600'
        cmax = '0.4'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('X='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F5 = FigureData('Zhang & Winfree (2009) Fig. 5 - Catalytic DSD system (varying toehold lengths)')
    F5.fname = 'Zhang2009-F5'
    current = F5
    template = zhang2009_F5_pil
    pilp = [6,7,5,8,4,9,3,2]
    litr = [(2822, 7.55), (3440, 7.11), (4716, 6.24), (9103, 3.23), 
            (10727, 2.13), (11249, 1.75), (13298, 0.34), (13626, 0.11)]

    for (pip, res) in zip(pilp, litr):
        pilstring  = template(pip)
        simulation = psim + ' --pyplot-labels F R S Z X --p0 R=30 S=10 Z=100 X=1'
        reporter = 'F'
        metric = 'diagonal-crossing-time'
        tmax = '14400'
        cmax = '10'
        current.add_system_simulation_setup(pilstring, simulation, reporter, 
                ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('X='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    return [F3, F4, F5]
    #return [F5]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)



