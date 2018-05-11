
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

    zhang2007_F1E = dict()
    zhang2007_F1E['name'] = 'zhang2007_F1'
    zhang2007_F1E['piltemplate'] = zhang2007_F1_pil
    zhang2007_F1E['pilparams'] = [None]
    zhang2007_F1E['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2007_F1E['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=5'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=2'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=0.5'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=0.2'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=0.1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=0.05'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '100000', '--t8', '360000', '--p0', 'S=10', 'F=13', 'OR=30', 'C=0.02')]
    zhang2007_F1E['reporter'] = 'ROX'
    zhang2007_F1E['exp_results'] = [(1206, 7.42), (1662, 6.83), (2450, 5.79), (3420, 4.47), (4350, 3.25), (5414, 1.84), (6070,1.0), (6383,0.58), (6602,0.27)] # seconds
    setups.append(zhang2007_F1E)

    zhang2007_F3 = dict()
    zhang2007_F3['name'] = 'zhang2007_F3'
    zhang2007_F3['piltemplate'] = zhang2007_F3_pil
    zhang2007_F3['pilparams'] = [None]
    zhang2007_F3['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2007_F3['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=5'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=2'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=0.5'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=0.2'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=0.1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=0.05'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'F1=13', 'S0=10', 'F0=13', 'OR=30', 'C0=0.02')]
    zhang2007_F3['reporter'] = 'ROX'
    zhang2007_F3['exp_results'] = [(1235, 5.0), (1813, 5.0), (2708, 5.5), (3041, 5.0), (3837, 4.0), (4513, 3.1), (4901, 2.6), (5605, 1.7), (5901, 1.2), (6030, 1.1)] # seconds
    setups.append(zhang2007_F3)

    zhang2007_F4 = dict()
    zhang2007_F4['name'] = 'zhang2007_F4'
    zhang2007_F4['piltemplate'] = zhang2007_F4_pil
    zhang2007_F4['pilparams'] = [None]
    zhang2007_F4['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2007_F4['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=7'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=5'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=3'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=0.7'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=0.5'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=0.3'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=0.2'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=10', 'F=13', 'SR=20', 'A=0.1')]
    zhang2007_F4['reporter'] = 'TET'
    zhang2007_F4['exp_results'] = [(240, 5), (342, 5), (394, 5), (609, 5), (688, 5), (884, 5), (1080, 5), (1122, 5), (1273, 5), (1369, 5), (1473, 5)] # seconds
    setups.append(zhang2007_F4)

    return setups


