
def dabby2013_4way_displacement((n,m)):
    assert n<7 and (m<7 or m==16)
    M = 0 if (m==6 or m==16) else 6-m
    N = 6-n
    return """# Nadine Dabby experiments
    # (n,m,k1_fit) from and Nadine Dabby, 
    # Caltech PhD Thesis, Table 5.2  
    # (note m,n have consistent meaning, but order in table is swapped.)

    length x = 21
    {}length m = {} # m
    {}length M = {} # M
    {}length n = {} # n
    {}length N = {} # N

    # starting states

    # x*( m* M* + N* n* )
    rep = x*( {} {} + {} {} )

    # m x( + ) n
    clx = {} x( + ) {}

    # x*( m*( M* + ) )
    pr1 = x*( {} {} + {} ) 

    # x*( n( + N* ) )
    pr2 = x*( {} + {} {} )

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
                '' if m==0 else 'm*(',
                '' if M==0 else 'M*',
                '' if m==0 else ')',
                '' if n==0 else 'n(',
                '' if N==0 else 'N*',
                '' if n==0 else ')')

def setups():
    """Returns a list of hardcoded dictionaries for every experimental setup.

    Provide DNA strands in form of a kernel string. Parameters to
    describe variations in the setup and a target value.

    Provide options for enumeration, such as condensation of the CRN or a
    release cutoff.

    Provide options for simulation, such as the initial concentrations.

    Provide completion threshold for simulation, such as target concentration.
    """
    setups = []

    dabby2013 = dict()
    dabby2013['name'] = 'Dabby2013-4way'
    dabby2013['piltemplate'] = dabby2013_4way_displacement
    dabby2013['pilparams'] = [(0, 2), (2, 2), (2, 0), 
                              (4, 2), (4, 0), (0, 4), (2, 4), 
                              (6, 2), (0, 6), (4, 4), (6, 0), 
                              (2, 6), (4, 6), (6, 4), (6, 6)]
    dabby2013['pepperargs'] = [('', {'condensed': True})]
    dabby2013['rates'] = None
    dabby2013['exp_results'] = [(0.047), (0.10), (0.033), (0.93), (0.039), (0.97), (56), (490), (58), (770), (5.0), (9.4e3), (7.0e4), (2.8e5), (6.9e5)]
    setups.append(dabby2013)

    return setups

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    d13 = FigureData('Dabby2013-4way')
    current = d13

    pilp = [(0, 2), (2, 2), (2, 0), (4, 2), (4, 0), (0, 4), (2, 4), (6, 2), (0, 6), (4, 4), (6, 0), (2, 6), (4, 6), (6, 4), (6, 6)]
    litr = [(0.047), (0.10), (0.033), (0.93), (0.039), (0.97), (56), (490), (58), (770), (5.0), (9.4e3), (7.0e4), (2.8e5), (6.9e5)]

    for (pil, rate) in zip(pilp, litr):
        pilstring  = dabby2013_4way_displacement(pil)
        exp_result = 'rep + clx -> pr1 + pr2 [k={}]'.format(rate)
        current.add_reaction_rate_setup(pilstring, exp_result)

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    return [d13]

if __name__ == '__main__':
    data()

