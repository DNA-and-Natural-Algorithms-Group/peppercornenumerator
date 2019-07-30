
def yin2008_F3_pil(x):
    return """
    length a = 6
    length b = 6
    length c = 6
    length d = 6
    length u = 6
    length v = 6
    length x = 6
    length y = 6

    # Initially present
    I = a* x* v* b* y*
    A = x*( v*( b*( y*( u* c* a* x* ) ) ) ) a
    B = v* d* y*( u*( c*( a*( x*( v* b* y* ) ) ) ) ) b
    C = c u( y( d( v( u* c* a* x* ) ) ) ) 
    D = d v( x( a( c( u( v* d* y* ) ) ) ) ) b* y*
    
    # Not initially present
    IA   = a*( x*( v*( b*( y*( +                                                   x* v* b* y* u* c* a* x* ) ) ) ) ) @i 0 M
    AB   =                       v* d* y* u* c* a*( x*( v*( b*( y*( x( a( c( u( y( b( + x* v* ) ) ) ) ) ) ) ) ) ) ) @i 0 M
    ABC  = c( u( y( d( v( u* c* a* x* v* d* y* u* + ) ) ) ) ) a*( x*( v*( b*( y*( x( a( c( u( y( b( + x* v* ) ) ) ) ) ) ) ) ) ) )  @i 0 M
    CDB  = c( u( y( d( v( u*( c*( a*( x*( v*( d*( y* u* + ) ) ) ) ) ) ) ) ) ) ) a* x* v* b*( y*( + v* d* y* u*( c*( a*( x*( v* b* y* ) ) ) ) ) ) @i 0 M
    CD   = c( u( y( d( v( u*( c*( a*( x*( v*( d*( y* u* + ) ) ) ) ) ) ) ) ) ) ) a* x* v* b* y* @i 0 M
    CDA  = c( u( y( d( v( u*( c*( a*( x*( v*( d*( y* u* + ) ) ) ) ) ) ) ) ) ) ) a*( x*( v*( b*( y*( + x* v* b* y* u* c* a* x* ) ) ) ) ) @i 0 M

    """.format(x)

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    # Default pilsimulator call
    psim = "pilsimulator --no-jacobian --nxy --header --atol 1e-10 --rtol 1e-10 --mxstep 10000 --t8 18000 --t-lin 18000"
    rates = {'k_slow': 1e-5, 'k_fast': 0.1}

    # Setup
    F3 = FigureData('Yin et al. (2008) Fig. 3 - Autocatalytic hairpin system (varying initiator)')
    F3.fname = 'Yin2008-F3'
    current = F3
    template = yin2008_F3_pil
    sims = [psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=20',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=6',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=2',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=1',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=0.6',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=0.4',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=0.2',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=0.1',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=0.06',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=0.02',
            psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=0.01']
            #psim + ' --pyplot-labels A B C D I --p0 A=20 B=20 C=20 D=20 I=0']

    # Diagonal crossing
    diagX = [
    (297.38, 3.30),
    (2256.56, 5.23), 
    (4411.66, 7.37),
    (5811.08, 8.72),
    (6930.61, 9.83),
    (7546.36, 10.45), 
    (8302.04, 11.17),
    (8973.76, 11.76),
    (9281.63, 12.11),
    (9561.52, 12.38),
    (9701.46, 12.52),
    (9841.40, 12.66)]

    for (sim, res) in zip(sims, diagX):
        pilstring  = template(None)
        simulation = sim
        reporter = 'A'
        metric = 'diagonal-crossing-time'
        tmax = '18000'
        cmax = '2;-20'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, ':'.join([tmax, cmax]), res, simargs=sim[sim.find('I='):])


    # Half completion
    halfC = [
    (45.48, -10.00),
    (1067.06, -10.00), 
    (3488.05, -10.00),
    (5321.28, -10.00),
    (6874.64, -10.00),
    (7700.29, -10.00),
    (8847.81, -10.00),
    (9729.45, -10.00),
    (10247.23, -10.00), 
    (10681.05, -10.00),
    (10876.97, -10.00),
    (11086.88, -10.00)]

    for (sim, res) in zip(sims, halfC):
        pilstring  = template(None)
        simulation = sim
        reporter = 'A'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric, '-10', res, simargs=sim[sim.find('I='):])

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    return [F3]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)



