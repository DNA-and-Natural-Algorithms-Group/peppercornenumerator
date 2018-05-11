
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

def zhang2010_F10A_pil(x):
    return """
    # Domains
    length d41 = 16
    length d42a = 4
    length d42b = 7
    length d42c = 9
    length d43  = 5
    length d44 = 15
    length d45 =  8
    length d46 = 15

    F4 = d42a d42b d42c d43 d44
    S4 = d41 d42a( d42b( d42c( + d46 d43( d44( + d45* ) ) ) ) )
    C4 = d44 d45

    W4 = d42a( d42b( d42c( d43( d44( + d45* ) ) ) ) )
    OP4 = d41 d42a d42b d42c
    SP4 = d46 d43 d44

    RQ = d41( d42a( d42b( d42c + ) ) )
    ROX = d41 d42a
    OR4 = d41( d42a( + d42b* ) )
    """.format(x)

def zhang2010_F10D_pil(x):
    return """
    # Domains
    length d51  = 15
    length d52a =  5
    length d52b =  7
    length d52c = 10
    length d63 =  5
    length d64 = 11
    length d65 =  8
    length d66 = 10

    F43 = d52a d52b d52c d63 d64
    S43 = d51 d52a( d52b( d52c( + d66 d63( d64( + d65* ) ) ) ) ) 
    C43 = d64 d65
    W43 = d52a( d52b( d52c( d63( d64( + d65* ) ) ) ) )
    OP43 = d51 d52a d52b d52c
    SP43 = d66 d63 d64
    RQ = d51( d52a( d52b( d52c + ) ) )
    RG = d51 d52a
    OR43 = d51( d52a( + d52b* ) )
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

    zhang2010_F3A = dict()
    zhang2010_F3A['name'] = 'zhang2010_F3A'
    zhang2010_F3A['piltemplate'] = zhang2007_F1_pil
    zhang2010_F3A['pilparams'] = [None]
    zhang2010_F3A['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2010_F3A['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=200', 'F=200', 'OR=300', 'C=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=200', 'F=200', 'OR=300', 'C=1')]
    zhang2010_F3A['reporter'] = 'ROX'
    zhang2010_F3A['exp_results'] = [(2235, 86.25), (11389, 43.93)] 
    setups.append(zhang2010_F3A)


    zhang2010_F3B = dict()
    zhang2010_F3B['name'] = 'zhang2010_F3B'
    zhang2010_F3B['piltemplate'] = zhang2007_F1_pil
    zhang2010_F3B['pilparams'] = [None]
    zhang2010_F3B['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2010_F3B['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=30', 'F=60', 'OR=90', 'C=3'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=30', 'F=60', 'OR=90', 'C=0.9')]
    zhang2010_F3B['reporter'] = 'ROX'
    zhang2010_F3B['exp_results'] = [(3887, 23.70), (8687, 17.04)]
    setups.append(zhang2010_F3B)

    zhang2010_F3C = dict()
    zhang2010_F3C['name'] = 'zhang2010_F3C'
    zhang2010_F3C['piltemplate'] = zhang2007_F1_pil
    zhang2010_F3C['pilparams'] = [None]
    zhang2010_F3C['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2010_F3C['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=3', 'F=6', 'OR=9', 'C=0.9'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'S=3', 'F=6', 'OR=9', 'C=0.3')]
    zhang2010_F3C['reporter'] = 'ROX'
    zhang2010_F3C['exp_results'] = [(9600, 2.24), (16066, 1.76)]
    setups.append(zhang2010_F3C)

    zhang2010_F3D = dict()
    zhang2010_F3D['name'] = 'zhang2010_F3D'
    zhang2010_F3D['piltemplate'] = zhang2007_F1_pil
    zhang2010_F3D['pilparams'] = [None]
    zhang2010_F3D['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2010_F3D['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '360000', '--p0', 'S=1', 'F=2', 'OR=3', 'C=1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '360000', '--p0', 'S=1', 'F=2', 'OR=3', 'C=0.1')]
    zhang2010_F3D['reporter'] = 'ROX'
    zhang2010_F3D['exp_results'] = [(13332, 0.82), (37775, 0.51)]
    setups.append(zhang2010_F3D)

    zhang2010_F10C = dict()
    zhang2010_F10C['name'] = 'zhang2010_F10C'
    zhang2010_F10C['piltemplate'] = zhang2010_F10A_pil
    zhang2010_F10C['pilparams'] = [None]
    zhang2010_F10C['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2010_F10C['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '36000', '--p0', 'S4=30', 'F4=60', 'OR4=90', 'C4=3'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '36000', '--p0', 'S4=30', 'F4=60', 'OR4=90', 'C4=0.09'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '36000', '--p0', 'S4=30', 'F4=60', 'OR4=90', 'C4=0.03')]
    zhang2010_F10C['reporter'] = 'ROX'
    zhang2010_F10C['exp_results'] = [(14475, 11.55), (24750, 5.88), (30425, 2.75)]
    setups.append(zhang2010_F10C)

    zhang2010_F10F = dict()
    zhang2010_F10F['name'] = 'zhang2010_F10F'
    zhang2010_F10F['piltemplate'] = zhang2010_F10D_pil
    zhang2010_F10F['pilparams'] = [None]
    zhang2010_F10F['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    zhang2010_F10F['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '360000', '--p0', 'S43=10', 'F43=20', 'OR43=30', 'C43=3'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '360000', '--p0', 'S43=10', 'F43=20', 'OR43=30', 'C43=0.09'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--t-lin', '10000', '--t8', '360000', '--p0', 'S43=10', 'F43=20', 'OR43=30', 'C43=0.03')]
    zhang2010_F10F['reporter'] = 'RG'
    zhang2010_F10F['exp_results'] = [(3358, 3.98), (5139, 2.04), (6075, 1.04)]
    setups.append(zhang2010_F10F)

    return setups

