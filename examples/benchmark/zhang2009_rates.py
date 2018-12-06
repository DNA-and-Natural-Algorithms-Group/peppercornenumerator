
def zhang2009_3way_displacement((c,d)):
    if c == 0 : # special leak case
        raise NotImplementedError()
    return """# Experiments Zhang 2009, 3-way strand displacement
        \rlength a = 16
        \rlength b1 = 6 
        \rlength b2 = 7 
        \rlength b3 = 7
        \rlength c = {:d}   # c
        \r{}length d = {:d} # d

        \r# Chopped sequence for reporter reaction.
        \rsup-sequence b = b1 b2 b3

        \r# Inputs
        \rX = b c
        \rS = a b( + {} c* )
        \rR = a( b1( + b2* ) )

        \r# Outputs
        \rL = b( c( + {} ) )        @initial 0 M
        \rY = a b                   @initial 0 M
        \rF = a b1                  @initial 0 M
        \rW = a( b1( b2( b3 + ) ) ) @initial 0 M

        \r# Intermediates
        \ri1 = a b( + b c( + {} ) ) @initial 0 M

        \r# Expected detailed reactions:
        \r# X + S <=> i1
        \r# i1 -> L + Y
        \r# Target condensed reactions:
        \r# X + S -> L + Y
        \r# R + Y -> F + W
        \r""".format(c,
                '#' if d==0 else '', d,
                '' if d==0 else 'd*',
                '' if d==0 else 'd*',
                '' if d==0 else 'd*')

def zhang2009_3way_exchange((n,m)):
    return """# Experiments Zhang 2009, 3-way strand displacement toehold exchange
        \rlength A = 16
        \rlength b = {:d} # m
        \rlength B = {:d} # 20 - m
        \rlength c = {:d} # n
        \rlength D = {:d} # 15 - n

        \rS = A b( B( + D* c* ) )
        \rX = B c

        \ri1 = A b( B( + B c( + D* ) ) ) @initial 0 M
        \ri2 = A b( B + B( c( + D* ) ) ) @initial 0 M
        \rL  = B( c( + D* ) ) b*         @initial 0 M
        \rY  = A b B                     @initial 0 M

        \r# potential reactions:
        \r# inv + clx -> i1
        \r# i1 -> inv + clx
        \r# i1 -> i2
        \r# i2 -> i1
        \r# i2 -> i3 + i4
        \r# i3 + i4 -> i2
        \r""".format(m, 20-m, n, 15-n)

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

    k3way_exp = [ 
        #((n,  m),  k)
        #((0,  15), 1.40), 
        ((1,  14), 8.17), 
        ((2,  13), 144), 
        ((3,  12), 1.08e3), 
        ((4,  11), 5.05e4), 
        ((5,  10), 9.64e5), 
        ((6,  9),  2.36e6), 
        ((7,  8),  3.22e6), 
        ((8,  9),  3.15e6), 
        ((9,  8),  2.77e6),
        ((10, 7),  2.83e6),
        ((15, 0),  4.78e6)]

    zhang2009_3way = dict()
    zhang2009_3way['name'] = 'Zhang2009-3way'
    zhang2009_3way['piltemplate'] = zhang2009_3way_displacement
    zhang2009_3way['pilparams'] = [x[0] for x in k3way_exp]
    zhang2009_3way['pepperargs'] = [('condensed', {'condensed': True}), ('param2', {'condensed': True, 'release_cutoff': 8})]
    zhang2009_3way['rates'] = None
    zhang2009_3way['exp_results'] = [x[1] for x in k3way_exp]
    setups.append(zhang2009_3way)

    k3wayX_exp = [
                                                              ((1, 4), 7.70),   ((1, 3),   5.48), ((1, 2),   23.5), ((1, 1),   18.9),
                                            ((2, 5), 43.6),   ((2, 4), 214.05), ((2, 3),  273.0), ((2, 2),  249.0), ((2, 1),  231.0), 
                          ((3, 6), 66.9),   ((3, 5), 215.0),  ((3, 4), 939.0),  ((3, 3),  974.0), ((3, 2),  907.0), ((3, 1),  846.0),
        ((4, 7), 131.0),  ((4, 6), 407.0),  ((4, 5), 4.25e3), ((4, 4), 2.13e4), ((4, 3), 2.41e4), ((4, 2), 2.29e4), ((4, 1), 1.97e4), 
        ((5, 7), 3.59e3), ((5, 6), 9.72e4), ((5, 5), 3.45e5), ((5, 4), 1.53e6), ((5, 3), 1.58e6), ((5, 2), 1.58e6), ((5, 1), 1.73e6),
        ((6, 7), 1.61e5), ((6, 6), 4.05e5), ((6, 5), 1.48e6), ((6, 4), 3.04e6), ((6, 3), 2.59e6), ((6, 2), 3.00e6), 
        ((7, 7), 4.7e5),  ((7, 6), 1.11e6), ((7, 5), 2.90e6), ((7, 4), 3.57e6), 
        ((8, 7), 1.94e6), ((8, 6), 2.68e6), ((8, 5), 3.14e6), ((8, 4), 3.37e6)]

    zhang2009_3wayX = dict()
    zhang2009_3wayX['name'] = 'Zhang2009-3wayX'
    zhang2009_3wayX['piltemplate'] = zhang2009_3way_exchange
    zhang2009_3wayX['pilparams'] = [x[0] for x in k3wayX_exp]
    zhang2009_3wayX['pepperargs'] = [('condensed', {'condensed': True})]
    zhang2009_3wayX['rates'] = {'reactants': ['inv', 'clx']}
    zhang2009_3wayX['exp_results'] = [x[1] for x in k3wayX_exp]
    setups.append(zhang2009_3wayX)

    return setups

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    # Setup Experiment 1
    z09_3way = FigureData('Zhang2009-3way-displacement')
    current = z09_3way

    k3way_exp = [ 
        #((n,  m),  k)
        #((0,  15), 1.40), 
        ((1,  14), 8.17), 
        ((2,  13), 144), 
        ((3,  12), 1.08e3), 
        ((4,  11), 5.05e4), 
        ((5,  10), 9.64e5), 
        ((6,  9),  2.36e6), 
        ((7,  8),  3.22e6), 
        ((8,  9),  3.15e6), 
        ((9,  8),  2.77e6),
        ((10, 7),  2.83e6),
        ((15, 0),  4.78e6)]

    for (pil, rate) in k3way_exp:
        pilstring  = zhang2009_3way_displacement(pil)
        exp_result = 'X + S -> Y + L [k={}]'.format(rate)
        current.add_reaction_rate_setup(pilstring, exp_result)

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup Experiment 2
    z09_3wayX = FigureData('Zhang2009-3way-exchange')
    current = z09_3wayX

    k3wayX_exp = [
                                                              ((1, 4), 7.70),   ((1, 3),   5.48), ((1, 2),   23.5), ((1, 1),   18.9),
                                            ((2, 5), 43.6),   ((2, 4), 214.05), ((2, 3),  273.0), ((2, 2),  249.0), ((2, 1),  231.0), 
                          ((3, 6), 66.9),   ((3, 5), 215.0),  ((3, 4), 939.0),  ((3, 3),  974.0), ((3, 2),  907.0), ((3, 1),  846.0),
        ((4, 7), 131.0),  ((4, 6), 407.0),  ((4, 5), 4.25e3), ((4, 4), 2.13e4), ((4, 3), 2.41e4), ((4, 2), 2.29e4), ((4, 1), 1.97e4), 
        ((5, 7), 3.59e3), ((5, 6), 9.72e4), ((5, 5), 3.45e5), ((5, 4), 1.53e6), ((5, 3), 1.58e6), ((5, 2), 1.58e6), ((5, 1), 1.73e6),
        ((6, 7), 1.61e5), ((6, 6), 4.05e5), ((6, 5), 1.48e6), ((6, 4), 3.04e6), ((6, 3), 2.59e6), ((6, 2), 3.00e6), 
        ((7, 7), 4.7e5),  ((7, 6), 1.11e6), ((7, 5), 2.90e6), ((7, 4), 3.57e6), 
        ((8, 7), 1.94e6), ((8, 6), 2.68e6), ((8, 5), 3.14e6), ((8, 4), 3.37e6)]

    for (pil, rate) in k3wayX_exp:
        pilstring  = zhang2009_3way_exchange(pil)
        exp_result = 'X + S -> Y + L [k={}]'.format(rate)
        current.add_reaction_rate_setup(pilstring, exp_result)

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    return [z09_3way, z09_3wayX]

if __name__ == '__main__':
    data()

