
def kotani2017_F2_pil(x):
    return """
    length a   = 22
    length b   = 22
    length c   = 22
    length t1  = 6
    length t2  = 6
    length t3  = 10
    length T2  = 2

    length d1s = 16
    length d2  = 6

    S1 = d1s T2 b( a( t2( + ) ) c*( t1* + ) )
    S2 = t1( c( a( + t2* ) b*( d2 t3 + ) ) )
    C1 = t1 c a

    P1 = t2* a*( c*( t1*( + ) ) )
    I1 = d1s T2 b( a t2 + c )
    I2 = d1s T2 b( a( t2( + ) ) b*( d2 t3 + ) c*( t1* + ) )

    P2 = d1s T2 b( a( t2( + ) ) ) d2 t3
    P3 = b( c*( t1* + ) )

    R = d1s( d2( + t3* ) )

    D = d1s d2
    RW = d1s( T2 b( a( t2( + ) ) ) d2( t3( + ) ) )

    """.format(x)

def kotani2017_F3_pil(x):
    return """
    length a   = 22
    length b   = 22
    length c1   = 11
    length c2   = 11
    length e   = 22
    length f   = 22
    length g   = 22
    length t1  = 6
    length t2  = 6
    length t3  = 10
    length t4 = 6
    length t5 = 6
    length T2  = 2

    length d1s = 16
    length d2  = 6

    S3 = t1 c1 c2( f( e( t5( + ) ) g*( t4* + ) ) )
    S4 = t4( g( e( + t5* ) f*( T2 a + c2 ) ) )
    C2 = t4 g e
    P4 = t5* e*( g*( t4*( + ) ) ) 
    P5 = t1 c1 c2 f( e( t5( + ) ) ) T2 a
    P6 = c2( f( g*( t4* + ) ) )


    S1 = d1s T2 b( a( t2( + ) ) c2*( c1*( t1* + ) ) )
    S2 = t1( c1( c2( a( + t2* ) b*( d2 t3 + ) ) ) )

    P1 = t2* a*( c2*( c1*( t1*( + ) ) ) )
    P2 = d1s T2 b( a( t2( + ) ) ) d2 t3
    P3 = b( c2*( c1*( t1* + ) ) )

    R = d1s( d2( + t3* ) )

    C1 = t1 c1 c2 a

    D = d1s d2
    RW = d1s( T2 b( a( t2( + ) ) ) d2( t3( + ) ) )
    """.format(x)

def kotani2017_F4_pil(x):
    return """
    \tlength a   = 22
    \tlength b   = 22
    \tlength o   = 22
    \tlength c1  = 11
    \tlength c2  = 11
    \tlength t1  = 6
    \tlength t2  = 6
    \tlength t3  = 10
    \tlength T2  = 2
    \tlength x   = 2
    \tlength y   = 2

    \tlength d1s = 16
    \tlength d1r = 2
    \tlength d2  = 6

    \tS5 = o( b*( T2 a + c2 ) a( t2( y + ) ) c2*( c1*( t1* x* + ) ) ) d2 t3
    \tS6 = y* t2* a*( b*( c2*( + x t1 c1 ) ) o*( + d1s T2 ) c2*( c1*( t1*( + x ) ) ) )
    \t
    \tC1 = x t1 c1 c2 a

    \tP1  = t2* a*( c2*( c1*( t1*( x*( + ) ) ) ) )           @i 0 M
    \tP2 = c2( b( a( t2( y( + ) ) ) ) )                     @i 0 M

    \tI5  = c1 c2 o*( d2 t3 + ) b*( T2 a + c2 ) a t2 y        @i 0 M
    \tI6 = c1( c2( o*( d2 t3 + ) b*( T2 a + c2 ) a( t2( y( + ) ) ) b*( c2*( + x t1 c1 ) ) o*( + d1s T2 ) ) ) t1* @i 0 M
    \tP10 = x( t1( c1( c2( b( o*( + ) ) T2 a( + t2* ) ) ) ) ) @i 0 M
 
    \tP8 = x t1 c1 c2 b( o*( + ) ) T2 a                     @i 0 M
    \tP9 = d1s T2 o( c2*( c1*( t1* + ) ) ) d2 t3            @i 0 M

    \tR = d1r( d1s( d2( + t3* ) ) )
    \tD = d1r d1s d2                                              @i 0 M
    \tRW = d1s( T2 o( c2*( c1*( t1* + ) ) ) d2( t3( + ) ) ) d1r*  @i 0 M
    """.format(x)

def kotani2017_F4_pil_simple(x):
    return """
    \tlength a   = 22
    \tlength b   = 22
    \tlength o   = 22
    \tlength c1  = 11
    \tlength c2  = 11
    \tlength t1  = 6
    \tlength t2  = 6
    \tlength t3a  = 6
    \tlength t3b  = 4
    \tlength T2  = 2
    \tlength x   = 2
    \tlength y   = 2

    \tlength d1s = 16
    \tlength d2  = 6

    \tS5 = o( b*( T2 a + c2 ) a( t2( y + ) ) c2*( c1*( t1* x* + ) ) ) d2 t3a t3b
    \tS6 = y* t2* a*( b*( c2*( + x t1 c1 ) ) o*( + d1s T2 ) c2*( c1*( t1*( + x ) ) ) )
    \t
    \tC1 = x t1 c1 c2 a

    \tP1  = t2* a*( c2*( c1*( t1*( x*( + ) ) ) ) )           @i 0 M
    \tP2 = c2( b( a( t2( y( + ) ) ) ) )                     @i 0 M

    \tI5  = c1 c2 o*( d2 t3a t3b + ) b*( T2 a + c2 ) a t2 y        @i 0 M
    \tI6 = c1( c2( o*( d2 t3a t3b + ) b*( T2 a + c2 ) a( t2( y( + ) ) ) b*( c2*( + x t1 c1 ) ) o*( + d1s T2 ) ) ) t1* @i 0 M
    \tP10 = x( t1( c1( c2( b( o*( + ) ) T2 a( + t2* ) ) ) ) ) @i 0 M
 
    \tP8 = x t1 c1 c2 b( o*( + ) ) T2 a                     @i 0 M
    \tP9 = d1s T2 o( c2*( c1*( t1* + ) ) ) d2 t3a t3b       @i 0 M

    \tR = d1s( d2( + t3a* ) )
    \tD = d1s d2                                              @i 0 M
    \tRW = d1s( T2 o( c2*( c1*( t1* + ) ) ) d2( t3a( t3b + ) ) ) @i 0 M
    """.format(x)

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    # Default pilsimulator call
    # If you use --no-jacobian, #7 doesn't show the 0 nM trajectory
    psim = "pilsimulator --nxy --header --atol 1e-12 --rtol 1e-12 --mxstep 1000"
    h10l = " --t0 0 --t8 36000 --t-lin 18000"
    h20l = " --t0 0 --t8 72000 --t-lin 18000"
    h30log = " --t0 0.1 --t8 108000 --t-log 18000"

    # Setup
    k17_F2 = FigureData('Kotani & Hughes (2017) Fig. 2 - Single-layer catalytic DSD with 4-way branch migration (varying catalyst)')
    k17_F2.fname = 'Kotani2017-F2'
    current = k17_F2
    template = kotani2017_F2_pil
    sims = [psim + h10l + " --pyplot-labels D S1 S2 R C1 --p0 S1=10 S2=10 R=20 C1=1",
            psim + h10l + " --pyplot-labels D S1 S2 R C1 --p0 S1=10 S2=10 R=20 C1=0.5",
            psim + h10l + " --pyplot-labels D S1 S2 R C1 --p0 S1=10 S2=10 R=20 C1=0.05"]
    litr = [(7733, 7.42), (11333, 6.18), (25533, 1.40)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'D'
        metric = 'diagonal-crossing-time'
        cmax = '10'
        tmax = '32400'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C1='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default']['release_cutoff'] = 7
    current.pepperargs['default']['max_complex_size'] = 8

    if evaluate:
        current.eval(cmpfig=True)

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    k17_F3 = FigureData('Kotani & Hughes (2017) Fig. 3 - Two-layer feedforward DSD system with 4-way branch migration (varying catalyst)')
    k17_F3.fname = 'Kotani2017-F3'
    current = k17_F3
    template = kotani2017_F3_pil
    sims = [psim + h30log + " --pyplot-labels D S1 S2 S3 S4 R C1 --p0 S1=10 S2=10 S3=10 S4=10 R=20 C1=0.1",
            psim + h30log + " --pyplot-labels D S1 S2 S3 S4 R C1 --p0 S1=10 S2=10 S3=10 S4=10 R=20 C1=0.01",
            psim + h30log + " --pyplot-labels D S1 S2 S3 S4 R C1 --p0 S1=10 S2=10 S3=10 S4=10 R=20 C1=0.001"]
    litr = [(21220, 7.72), (64203, 3.12), (86996, 0.69)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'D'
        metric = 'diagonal-crossing-time'
        cmax = '10'
        tmax = '97200'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res, simargs=sim[sim.find('C1='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default']['release_cutoff'] = 7
    current.pepperargs['default']['max_complex_size'] = 8

    if evaluate:
        current.eval(cmpfig=True)

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    k17_F4 = FigureData('Kotani & Hughes (2017) Fig. 4 - Autocatalytic DSD system with 4-way branch migration (varying autocatalyst)')
    k17_F4.fname = 'Kotani2017-F4'
    current = k17_F4
    template = kotani2017_F4_pil
    sims = [psim + h20l + " --pyplot-labels D S5 S6 R C1 --p0 S5=10 S6=10 R=20 C1=0.1",
            psim + h20l + " --pyplot-labels D S5 S6 R C1 --p0 S5=10 S6=10 R=20 C1=0.01",
            psim + h20l + " --pyplot-labels D S5 S6 R C1 --p0 S5=10 S6=10 R=20 C1=0.001",
            psim + h20l + " --pyplot-labels D S5 S6 R C1 --p0 S5=10 S6=10 R=20 C1=0"]
    litr = [(6136, 5), (9150, 5), (10776, 5), (11637, 5)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'D'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, res, simargs=sim[sim.find('C1='):])

    current.pepperargs['default'] = current.pepperargs['DETAILED'].copy()
    current.pepperargs['default']['release_cutoff'] = 8
    current.pepperargs['default']['max_complex_size'] = 24
    #current.pepperargs['default']['k_slow'] = 1e-10
    #current.pepperargs['default']['k_fast'] = 1e-2

    if evaluate:
        current.eval(cmpfig=True)

    # current.pepperargs['#1-detailed']  = {'condensed': False, 'conc': 'nM', 'release_cutoff': 8}
    # current.pepperargs['#1-condensed'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 8}
    # current.pepperargs['#2-condensed'] = {'condensed': True, 'conc': 'nM', 'k_slow': 1e-3, 'max_complex_size': 20} # doesn't work, no 4-way branch migration, no D
    # current.pepperargs['#3-condensed'] = {'condensed': True, 'conc': 'nM', 'k_slow': 1e-4, 'max_complex_size': 10} # works! but not as good as detailed...
    # current.pepperargs['#4-condensed'] = {'condensed': True, 'conc': 'nM', 'k_slow': 1e-4, 'k_fast': 1e-3, 'max_complex_size': 16}
    # current.pepperargs['#5-condensed'] = {'condensed': True, 'conc': 'nM', 'k_slow': 1e-4, 'k_fast': 1e-2, 'max_complex_size': 24} # solver cannot handle leak
    # current.pepperargs['#6-condensed'] = {'condensed': True, 'conc': 'nM', 'k_slow': 1e-5, 'k_fast': 1e-2, 'max_complex_size': 24}
    # current.pepperargs['#7-condensed'] = {'condensed': True, 'conc': 'nM', 'k_slow': 1e-10, 'k_fast': 1e-2, 'max_complex_size': 24}

    if verbose:
        for df in current.get_dataframes():
            print(df)

    return [k17_F2, k17_F3, k17_F4]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)

