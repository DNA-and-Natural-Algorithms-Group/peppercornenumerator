
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

    zhang2011_F1 = dict()
    zhang2011_F1['name'] = 'Zhang2011c-F1'
    zhang2011_F1['piltemplate'] = zhang2011_F1_pil
    zhang2011_F1['pilparams'] = [None]
    zhang2011_F1['pepperargs'] = {'condensed': False, 'conc': 'nM', 
            'max_complex_count': 10000, 
            'max_reaction_count': 10000, 
            'k_fast': 0.01, 'k_slow': 0.001}
    zhang2011_F1['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '54000', '--p0', 'R=60', 'D1=20', 'T1=18', 'T2=18'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '54000', '--p0', 'R=60', 'D1=20', 'T1=12', 'T2=12'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '54000', '--p0', 'R=60', 'D1=20', 'T1=6', 'T2=6')]
    zhang2011_F1['reporter'] = 'F'
    zhang2011_F1['exp_results'] = [(463, 9.45), (751, 7.28), (1104, 4.73)] # NOTE: It is actually not clear to me if these are the correct values...
    setups.append(zhang2011_F1)

    return setups


