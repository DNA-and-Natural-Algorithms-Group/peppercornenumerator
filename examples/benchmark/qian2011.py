
def qian2011_F1(x):
    return """
    sequence T   = TCT : 3
    sequence S2a  = CAAAACA : 7
    sequence S2b = AAACCTCA : 8
    sequence S5  = CACCACCAAACTTCA : 15
    sequence S6  = CATAACACAATCACA : 15
    sequence S7  = CAACATATCAATTCA : 15

    sup-sequence S2  = S2a S2b : 15

    I = S5 T S2
    G = S6 T( S5( + T* ) ) # @initial 100 nM
    T = S5( + S2a* T* )    # @initial  50 nM
    F = S7 T S5            # @initial 200 nM
    R = S6( + T* )         # @initial 150 nM

    Q = S6
    """.format(x)


def qian2011_F2(x):
    return """
    length T = 3
    length S1 = 15
    length S2a = 7
    length S2b = 8
    length S3 = 15
    length S5 = 15
    length S6 = 15
    length S7 = 15

    sup-sequence S2 = S2a S2b

    W12 = S2 T S1
    W32 = S2 T S3

    G2_25 = S5 T( S2( + T* ) )
    Th2_55 = S5( + S2a* T* )
    G5_56 = S6 T( S5( + T* ) )
    W57 = S7 T S5

    R6 = S6( + T* )
    Q = S6
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

    qian2011_F1E = dict()
    qian2011_F1E['name'] = 'qian2011_F1E'
    qian2011_F1E['piltemplate'] = qian2011_F1
    qian2011_F1E['pilparams'] = [None]
    qian2011_F1E['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    qian2011_F1E['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G=100', 'T=50', 'F=200', 'R=150', 'I=100'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G=100', 'T=50', 'F=200', 'R=150', 'I=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G=100', 'T=50', 'F=200', 'R=150', 'I=80'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G=100', 'T=50', 'F=200', 'R=150', 'I=70'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G=100', 'T=50', 'F=200', 'R=150', 'I=60'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G=100', 'T=50', 'F=200', 'R=150', 'I=50')]
    qian2011_F1E['reporter'] = 'Q'
    qian2011_F1E['exp_results'] = [(1779, 0.84), (2101, 0.81), (2609, 0.76), (3546, 0.67), (6624, 0.39), (9809, 0.08)]
    setups.append(qian2011_F1E)

    qian2011_F2C_OR = dict()
    qian2011_F2C_OR['name'] = 'qian2011_F2C_OR'
    qian2011_F2C_OR['piltemplate'] = qian2011_F2
    qian2011_F2C_OR['pilparams'] = [None]
    qian2011_F2C_OR['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    qian2011_F2C_OR['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G2_25=200', 'Th2_55=60', 'G5_56=100', 'W57=200', 'R6=150', 'W12=90', 'W32=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G2_25=200', 'Th2_55=60', 'G5_56=100', 'W57=200', 'R6=150', 'W12=10', 'W32=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G2_25=200', 'Th2_55=60', 'G5_56=100', 'W57=200', 'R6=150', 'W12=90', 'W32=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G2_25=200', 'Th2_55=60', 'G5_56=100', 'W57=200', 'R6=150', 'W12=10', 'W32=10')]
    qian2011_F2C_OR['reporter'] = 'Q'
    qian2011_F2C_OR['exp_results'] = [(3245, 0.85), (4733, 0.68), (4733, 0.68), (10530, 0.04)]
    setups.append(qian2011_F2C_OR)

    qian2011_F2C_AND = dict()
    qian2011_F2C_AND['name'] = 'qian2011_F2C_AND'
    qian2011_F2C_AND['piltemplate'] = qian2011_F2
    qian2011_F2C_AND['pilparams'] = [None]
    qian2011_F2C_AND['pepperargs'] = {'condensed': True, 'conc': 'nM'}
    qian2011_F2C_AND['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G2_25=200', 'Th2_55=120', 'G5_56=100', 'W57=200', 'R6=150', 'W12=90', 'W32=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G2_25=200', 'Th2_55=120', 'G5_56=100', 'W57=200', 'R6=150', 'W12=10', 'W32=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G2_25=200', 'Th2_55=120', 'G5_56=100', 'W57=200', 'R6=150', 'W12=90', 'W32=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'G2_25=200', 'Th2_55=120', 'G5_56=100', 'W57=200', 'R6=150', 'W12=10', 'W32=10')]
    qian2011_F2C_AND['reporter'] = 'Q'
    qian2011_F2C_AND['exp_results'] = [(13797, 0.68), (41022, 0.07), (41022, 0.07), (42239, 0.04)]
    setups.append(qian2011_F2C_AND)

    return setups


