
def zhang2011_F1_pil(x):
    return """
    # Domains
    length d1 = 8
    length d2a = 7
    length d2b = 7
    length d2c = 4
    length d3a = 14
    length d3b = 4
    length d4 = 8

    sup-sequence d2 = d2a d2b d2c : 18
    sup-sequence d3 = d3a d3b : 18

    # Setup Figure 1
    D1 = d2( d3( + d4* ) ) d1*
    T2 = d3 d4
    II = d1( d2( + d2 d3( + d4* ) ) )
    IJ = d2( d3 + d3( d4( + ) ) ) d1*
    T1 = d1 d2
    H1 = d1( d2( + d3( d4( + ) ) ) ) @initial 0 M

    # Result
    P1 = d2 d3 @initial 0 M

    # Reporter System Figure 2
    R = d2c( d3a( + ) ) d2b*
    F = d2c d3a @initial 0 M
    FQ = d2a d2b( d2c( d3a( d3b + ) ) ) @initial 0 M

    """.format(x)

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    # Default pilsimulator call
    psim = "pilsimulator --no-jacobian --nxy --header"
    psim += " --atol 1e-10 --rtol 1e-10 --mxstep 1000"
    psim += " --t8 54000 --t-lin 54000"
    rates = {'k_slow': 0.001, 'k_fast': 0.01}

    # Setup
    F3A = FigureData('Zhang (2011) Fig. 3A - Cooperative strand displacement')
    F3A.fname = 'Zhang2011-F3A'
    current = F3A
    template = zhang2011_F1_pil
    sims = [psim + ' --pyplot-labels R D1 T1 T2 F --p0 R=60 D1=20 T1=18 T2=18',
            psim + ' --pyplot-labels R D1 T1 T2 F --p0 R=60 D1=20 T1=12 T2=12',
            psim + ' --pyplot-labels R D1 T1 T2 F --p0 R=60 D1=20 T1=6 T2=6']

    diagX = [(409.46, 9.15), (425.68, 6.14), (300.00, 3.25)]
    diagXe = [('1800','12'), ('1800','8'), ('1800','4')]

    for (sim, res, tc) in zip(sims, diagX, diagXe):
        pilstring  = template(None)
        simulation = sim
        reporter = 'F'
        metric = 'diagonal-crossing-time'
        tmax = tc[0]
        cmax = tc[1]
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric,
                ':'.join([tmax, cmax]), res, simargs=sim[sim.find('T1='):])

    halfC = [(381.08, 9.02), (385.14, 6.01), (263.51, 3.02)]

    for (sim, res) in zip(sims, halfC):
        pilstring  = template(None)
        simulation = sim
        reporter = 'F'
        metric = 'completion-time'
        current.add_system_simulation_setup(pilstring, simulation, reporter, metric,
                '50%', res, simargs=sim[sim.find('T1='):])

    current.pepperargs['default'] = current.pepperargs['Detailed'].copy()
    current.pepperargs['default'].update(rates)

    if evaluate:
        current.eval(verbose=verbose)

    if verbose:
        for df in current.get_dataframes():
            print(df)

    return [F3A]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)




