
def genot2011_F3_pil((m,n)) :
    return """
    length d1 = 15
    length d2 = 22
    length d3 = 14
    
    length d05 = {}
    length d04 = {}

    S = d1 d2( + d3* d05* )
    I = d2 d04 d3

    T = d2( d04 d3( + ) d05* )
    W = d1 d2
    """.format(m, n)

def genot2011_F4A_pil(d3) :
    return """
    length d1 = 15
    length d2 = 22
    length d3 = {}
    
    S = d1 d2( + d3* )
    I = d2 d3

    T = d2( d3( + ) )
    W = d1 d2
    """.format(d3)

def genot2011_F4B_pil(d3) :
    return """
    length d1 = 15
    length d2 = 22
    length d8 = 1
    length d3 = {}
    
    S = d1 d2( + d3* d8 )
    I = d2 d3

    T = d2( d3( + ) d8 )
    W = d1 d2
    """.format(d3)

def genot2011_F4C_pil(d3) :
    return """
    length d1 = 15
    length d2 = 22
    length d3 = 6 # original 
    #length d3 = 4 # try-to-fit
    
    S = d1 d2( + d3* )
    I = d2 d3

    T = d2( d3( + ) )
    W = d1 d2
    """

def genot2011_F4D_pil(d3) :
    return """
    length d1 = 15
    length d2 = 22
    length d3 = 14 # original
    #length d3 = 7 # try-to-fit
    length d4 = 17
    length d5 = 17
    
    S = d1 d2( + d3* d5 )  # d3: 14
    I = d2 d4 d3

    T = d2( d4 d3( + ) d5 )
    W = d1 d2
    """

def genot2011_SF4A_pil(dl) :
    return """
    length d1 = 15
    length d2a = 5
    length d2b = 7
    length d2c = 10
    length d3 = 14 # original 
    #length d3 = 6 # try-to-fit
    length d4 = {}
    length d5 = {}

    #composite d2 = d2a d2b d2c
    
    S = d1 d2a( d2b( d2c( + d3* d5 ) ) )
    X = d2a d2b d2c d4 d3

    W = d2a( d2b( d2c( d4 d3( + ) d5 ) ) )
    P = d1 d2a d2b d2c

    R = d1( d2a( + d2b* ) )

    F = d1 d2a
    FW = d1( d2a( d2b( d2c + ) ) )


    """.format(dl, dl)


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

    genot2011_F3 = dict()
    genot2011_F3['name'] = 'Genot2011-F3'
    genot2011_F3['piltemplate'] = genot2011_F3_pil
    genot2011_F3['pilparams'] = [ (8,7), (17,17), (20,20), (23,23)] # original
    genot2011_F3['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_fast': 20, 'k_slow': 1e-10}
    #genot2011_F3['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_slow': 0.0001}
    genot2011_F3['simulation'] = [('pilsimulator', '--nxy', '--atol', '1e-12', '--rtol', '1e-12', '--t-lin', '10000', '--t8', '1000', '--p0', 'S=6.6', 'I=660')]
    genot2011_F3['reporter'] = 'T'
    genot2011_F3['exp_results'] = [(42, 2.2), (48, 2.2), (100, 2.2), (91, 2.2)] # seconds
    setups.append(genot2011_F3)

    genot2011_F4A = dict()
    genot2011_F4A['name'] = 'Genot2011-F4A'
    genot2011_F4A['piltemplate'] = genot2011_F4A_pil
    genot2011_F4A['pilparams'] = [ 11, 9] # original
    #genot2011_F4A['pilparams'] = [ 9, 7]
    genot2011_F4A['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_fast': 20, 'k_slow': 1e-10}
    #genot2011_F4A['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_slow': 0.0001}
    genot2011_F4A['simulation'] = [('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '1000', '--p0', 'S=6.6', 'I=22')]
    genot2011_F4A['reporter'] = 'T'
    #genot2011_F4A['exp_results'] = [(32, 0.58*6.6), (38, 0.58*6.6)]
    genot2011_F4A['exp_results'] = [(2, 2.2), (3, 2.2)]
    setups.append(genot2011_F4A)


    genot2011_F4B = dict()
    genot2011_F4B['name'] = 'Genot2011-F4B'
    genot2011_F4B['piltemplate'] = genot2011_F4B_pil
    genot2011_F4B['pilparams'] = [ 11, 9] # original
    #genot2011_F4B['pilparams'] = [ 9, 7]
    genot2011_F4B['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_fast': 20, 'k_slow': 1e-10}
    #genot2011_F4B['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_slow': 0.0001}
    genot2011_F4B['simulation'] = [('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '1000', '--p0', 'S=6.6', 'I=22')]
    genot2011_F4B['reporter'] = 'T'
    #genot2011_F4B['exp_results'] = [(44, 3.3), (182, 3.3)]
    genot2011_F4B['exp_results'] = [(15, 2.2), (75, 2.2)]
    setups.append(genot2011_F4B)


    genot2011_F4C = dict()
    genot2011_F4C['name'] = 'Genot2011-F4C'
    genot2011_F4C['piltemplate'] = genot2011_F4C_pil
    genot2011_F4C['pilparams'] = [None]
    genot2011_F4C['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_fast': 20, 'k_slow': 1e-10}
    #genot2011_F4C['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_slow': 0.0001}
    genot2011_F4C['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '10000',  '--p0', 'S=6.6', 'I=330'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '10000',  '--p0', 'S=6.6', 'I=145'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '10000',  '--p0', 'S=6.6', 'I=66')]
    genot2011_F4C['reporter'] = 'T'
    #genot2011_F4C['exp_results'] = [(69, 3.3), (177, 3.3), (347, 3.3)]
    genot2011_F4C['exp_results'] = [(25, 2.2), (75, 2.2), (150, 2.2)]
    setups.append(genot2011_F4C)

    genot2011_F4D = dict()
    genot2011_F4D['name'] = 'Genot2011-F4D'
    genot2011_F4D['piltemplate'] = genot2011_F4D_pil
    genot2011_F4D['pilparams'] = [None]
    genot2011_F4D['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_fast': 20, 'k_slow': 1e-10}
    #genot2011_F4D['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_slow': 0.0001}
    genot2011_F4D['simulation'] = [
            ('pilsimulator', '--force', '--nxy', '--atol', '1e-12', '--rtol', '1e-12', '--t0', '0.1', '--t-log', '10000', '--t8', '1000', '--p0', 'S=6.6', 'I=330'),
            ('pilsimulator', '--force', '--nxy', '--atol', '1e-12', '--rtol', '1e-12', '--t0', '0.1', '--t-log', '10000', '--t8', '1000', '--p0', 'S=6.6', 'I=145'),
            ('pilsimulator', '--force', '--nxy', '--atol', '1e-12', '--rtol', '1e-12', '--t0', '0.1', '--t-log', '10000', '--t8', '1000', '--p0', 'S=6.6', 'I=66')]
    genot2011_F4D['reporter'] = 'T'
    #genot2011_F4D['exp_results'] = [(176, 3.3), (213, 3.3), (248, 3.3)]
    genot2011_F4D['exp_results'] = [(51, 2.2), (55, 2.2), (60, 2.2)]
    setups.append(genot2011_F4D)

    genot2011_SF4A = dict()
    genot2011_SF4A['name'] = 'Genot2011-SF4A'
    genot2011_SF4A['piltemplate'] = genot2011_SF4A_pil
    genot2011_SF4A['pilparams'] = [17,20,23]
    genot2011_SF4A['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_fast': 20, 'k_slow': 1e-10}
    #genot2011_SF4A['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_slow': 0.0001}
    genot2011_SF4A['simulation'] = [
            ('pilsimulator', '--nxy', '--p0', 'R=30', 'X=40', 'S=10', '--t8', '1e5', '--t-lin', '100000', '--atol', '1e-10', '--rtol', '1e-10')]
    genot2011_SF4A['reporter'] = 'F'
    genot2011_SF4A['exp_results'] = [(9088, 5), (13927, 5), (28681, 5)]
    setups.append(genot2011_SF4A)


    return setups

