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

    zhang2009_F3 = dict()
    zhang2009_F3['name'] = 'zhang2009-F3'
    zhang2009_F3['piltemplate'] = zhang2009_F1DF_displacement_pil
    zhang2009_F3['pilparams'] = [5]
    zhang2009_F3['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2009_F3['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '36000', '--t8', '36000', '--p0', 'R=3', 'S=1', 'X=0.6'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '36000', '--t8', '36000', '--p0', 'R=3', 'S=1', 'X=0.4'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '36000', '--t8', '36000', '--p0', 'R=3', 'S=1', 'X=0.2')]
    zhang2009_F3['reporter'] = 'F'
    zhang2009_F3['exp_results'] = [(1260, 0.31), (1867, 0.27), (3498, 0.17)]
    setups.append(zhang2009_F3)

    zhang2009_F4 = dict()
    zhang2009_F4['name'] = 'zhang2009-F4'
    zhang2009_F4['piltemplate'] = zhang2009_F1DF_exchange_pil
    zhang2009_F4['pilparams'] = [(4,7), (5,7), (6,7), (7,7)]
    zhang2009_F4['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2009_F4['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '3600', '--t8', '3600', '--p0', 'R=3', 'S=1', 'X=0.4')]
    zhang2009_F4['reporter'] = 'F'
    zhang2009_F4['exp_results'] = [(785, 0.28), (823, 0.27), (1192, 0.23), (1622, 0.18)]
    setups.append(zhang2009_F4)

    zhang2009_F5 = dict()
    zhang2009_F5['name'] = 'zhang2009-F5'
    zhang2009_F5['piltemplate'] = zhang2009_F5_pil
    zhang2009_F5['pilparams'] = [6,7,5,8,4,9,3,2]
    zhang2009_F5['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2009_F5['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-12', '--rtol', '1e-12', '--mxstep', '1000', '--t-lin', '36000', '--t8', '360000', '--p0', 'R=30', 'S=10', 'Z=100', 'X=1')]
    zhang2009_F5['reporter'] = 'F'
    zhang2009_F5['exp_results'] = [(2822, 7.55), (3440, 7.11), (4716, 6.24), (9103, 3.23), (10727, 2.13), (11249, 1.75), (13298, 0.34), (13626, 0.11)]
    setups.append(zhang2009_F5)

    return setups

