
def qian2011_F1(x):
    return """
    INPUT(I) = w[2,5]
    OUTPUT(O) = Fluor[6]
    seesaw[5, {2}, {6,7}]
    reporter[6,5]

    conc[w[5,7], 2*c]
    conc[g[5, w[5,6]], 1*c]
    conc[th[w[2,5],5], 0.5*c]
    """

def qian2011_F2(x):
    return """
    INPUT(x1) = w[1,2]
    INPUT(x2) = w[3,2]
    OUTPUT(y) = Fluor[6]
    seesaw[2, {1,3}, {5}]
    seesaw[5, {2}, {6,7}]
    reporter[6,5]

    conc[w[5,7], 2*c]
    conc[g[5, w[5,6]], 1*c]
    conc[th[w[2,5],5], 0.6*c]
    #conc[th[w[2,5],5], 1.2*c]
    conc[g[2,w[2,5]], 2*c]
    """

def qian2011_SF26(x):
    return """
    INPUT(x1) = w[17,16]
    INPUT(x2) = w[18,16]
    INPUT(x3) = w[9,4]
    INPUT(x4) = w[3,2]
    OUTPUT(y) = Fluor[6]
    seesawOR[2,5,{1,3},{6}]
    seesawOR[4,1,{8,9},{2}]
    seesawOR[16,8,{17,18},{4}]
    reporter[6,5]
    """

def qian2011_SF27(x):
    return """
    INPUT(x1) = w[21,20]
    INPUT(x2) = w[22,20]
    INPUT(x3) = w[18,16]
    INPUT(x4) = w[9,4]
    INPUT(x5) = w[3,2]
    OUTPUT(y) = Fluor[6]

    seesawOR[2,5,{1,3},{6}]
    seesawOR[4,1,{8,9},{2}]
    seesawOR[16,8,{17,18},{4}]
    seesawOR[20,17,{21,22},{16}]
    reporter[6,5]
    """

def qian2011_SF28(x):
    return """
    INPUT(x1) = w[21,20]
    INPUT(x2) = w[22,20]
    INPUT(x3) = w[18,16]
    INPUT(x4) = w[9,4]
    INPUT(x5) = w[13,12]
    INPUT(x6) = w[14,12]
    OUTPUT(y) = Fluor[6]

    seesawOR[2,5,{1,3},{6}]
    seesawAND[12,3,{13,14},{2}]
    seesawOR[4,1,{8,9},{2}]
    seesawOR[16,8,{17,18},{4}]
    seesawOR[20,17,{21,22},{16}]
    reporter[6,5]
    """

def qian2011_SF29(x):
    return """
    INPUT(x1) = w[23,4]
    INPUT(x2) = w[9,4]
    INPUT(x3) = w[13,12]
    INPUT(x4) = w[14,12]
    INPUT(x5) = w[24,16]
    INPUT(x6) = w[18,16]
    INPUT(x7) = w[21,20]
    INPUT(x8) = w[22,20]
    OUTPUT(y) = Fluor[6]

    seesawOR[2,5,{1,3,8,17},{6}]
    seesawOR[4,1,{23,9},{2}]
    seesawOR[12,3,{13,14},{2}]
    seesawOR[16,8,{24,18},{2}]
    seesawOR[20,17,{21,22},{2}]
    reporter[6,5]
    """

def qian2011_SF29_AND(x):
    return """
    INPUT(x1) = w[23,4]
    INPUT(x2) = w[9,4]
    INPUT(x3) = w[13,12]
    INPUT(x4) = w[14,12]
    INPUT(x5) = w[24,16]
    INPUT(x6) = w[18,16]
    INPUT(x7) = w[21,20]
    INPUT(x8) = w[22,20]
    OUTPUT(y) = Fluor[6]

    seesawAND[2,5,{1,3,8,17},{6}]
    seesawOR[4,1,{23,9},{2}]
    seesawOR[12,3,{13,14},{2}]
    seesawOR[16,8,{24,18},{2}]
    seesawOR[20,17,{21,22},{2}]
    reporter[6,5]
    """

def qian2011_SF30_OR(x):
    return """
    INPUT(x1) = w[13,12]
    INPUT(x2) = w[14,12]
    OUTPUT(y1) = Fluor[6]
    OUTPUT(y2) = Fluor[23]
    OUTPUT(y3) = Fluor[24]
    OUTPUT(y4) = Fluor[25]

    seesawOR[12,3,{13,14},{2,4,16,20}]
    seesawOR[2,5,{3,22},{6}]
    seesawOR[4,1,{3,9},{23}]
    seesawOR[16,8,{3,18},{24}]
    seesawOR[20,17,{3,21},{25}]
    reporter[6,5]
    """


def setups():
    """Returns a list of hardcoded dictionaries for every experimental setup.

    Provide DNA strands in form of a kernel string. Parameters to
    describe variations in the setup and a target value.

    Provide options for enumeration, such as condensation of the CRN or a
    release cutoff.

    Provide options for simulation, such as the initial concentrations.

    Provide completion threshold for simulation, such as target concentration.
    """
    ddG_bind = 0
    setups = []

    qian2011_F1E = dict()
    qian2011_F1E['name'] = 'Qian2011-F1E'
    qian2011_F1E['piltemplate'] = qian2011_F1
    qian2011_F1E['pilparams'] = [None]
    qian2011_F1E['pepperargs'] = {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind, 'k_slow' : 0.01, 'k_fast': 1}
    qian2011_F1E['simulation'] = [
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=100'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=90'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=80'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=70'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=60'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=50'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=40'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=30'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=20'),
            ('pilsimulator', '--nxy', '--force', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'I=10')]
    qian2011_F1E['reporter'] = 'O'
    qian2011_F1E['exp_results'] = [(1779, 84), (2101, 81), (2609, 76), (3546, 67), (6624, 39)]
    qian2011_F1E['off_results'] = [(9809, 0.8), (9809, 0.7), (9809, 0.8), (9809, 0.7), (9809, 0.7)]
    setups.append(qian2011_F1E)

    qian2011_F2C_OR = dict()
    qian2011_F2C_OR['name'] = 'Qian2011-F2C-OR'
    qian2011_F2C_OR['piltemplate'] = qian2011_F2
    qian2011_F2C_OR['pilparams'] = [None]
    qian2011_F2C_OR['pepperargs'] = {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}
    qian2011_F2C_OR['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=10')]
    qian2011_F2C_OR['reporter'] = 'y'
    qian2011_F2C_OR['exp_results'] = [(3245, 85), (4733, 68), (4733, 68)]
    qian2011_F2C_OR['off_results'] = [(10530, 0.4)]
    setups.append(qian2011_F2C_OR)

    qian2011_F2C_AND = dict()
    qian2011_F2C_AND['name'] = 'Qian2011-F2C-AND'
    qian2011_F2C_AND['piltemplate'] = qian2011_F2
    qian2011_F2C_AND['pilparams'] = [None]
    qian2011_F2C_AND['pepperargs'] = {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}
    qian2011_F2C_AND['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '54000', '--p0', 'T2_5=120', 'x1=90', 'x2=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '54000', '--p0', 'T2_5=120', 'x1=10', 'x2=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '54000', '--p0', 'T2_5=120', 'x1=90', 'x2=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '54000', '--p0', 'T2_5=120', 'x1=10', 'x2=10')]
    qian2011_F2C_AND['reporter'] = 'y'
    qian2011_F2C_AND['exp_results'] = [(13797, 68)]
    qian2011_F2C_AND['off_results'] = [(41022, 0.7), (41022, 0.7), (42239, 0.4)]
    setups.append(qian2011_F2C_AND)

    qian2011_FS26 = dict()
    qian2011_FS26['name'] = 'Qian2011-FS26'
    qian2011_FS26['piltemplate'] = qian2011_SF26
    qian2011_FS26['pilparams'] = [None]
    qian2011_FS26['pepperargs'] = {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}
    qian2011_FS26['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=90', 'x3=90', 'x4=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=10', 'x3=90', 'x4=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=90', 'x3=10', 'x4=10')]
    qian2011_FS26['reporter'] = 'y'
    qian2011_FS26['exp_results'] = [(60*15, 50), (60*60, 50), (3600*3, 50)]
    setups.append(qian2011_FS26)

    qian2011_FS27 = dict()
    qian2011_FS27['name'] = 'Qian2011-FS27'
    qian2011_FS27['piltemplate'] = qian2011_SF27
    qian2011_FS27['pilparams'] = [None]
    qian2011_FS27['pepperargs'] = {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}
    qian2011_FS27['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=90', 'x3=90', 'x4=90', 'x5=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=10', 'x3=90', 'x4=10', 'x5=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=90', 'x3=10', 'x4=10', 'x5=10')]
    qian2011_FS27['reporter'] = 'y'
    qian2011_FS27['exp_results'] = [(60*15, 50), (3600*2.5, 50), (3600*3.5, 50)]
    setups.append(qian2011_FS27)

    qian2011_FS28 = dict()
    qian2011_FS28['name'] = 'Qian2011-FS28'
    qian2011_FS28['piltemplate'] = qian2011_SF28
    qian2011_FS28['pilparams'] = [None]
    qian2011_FS28['pepperargs'] = {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}
    qian2011_FS28['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=90', 'x3=90', 'x4=90', 'x5=90', 'x6=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=10', 'x3=90', 'x4=90', 'x5=10', 'x6=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=90', 'x3=10', 'x4=90', 'x5=10', 'x6=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=10', 'x3=10', 'x4=90', 'x5=10', 'x6=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=10', 'x3=10', 'x4=10', 'x5=90', 'x6=90'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=10', 'x3=90', 'x4=10', 'x5=90', 'x6=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=90', 'x3=90', 'x4=10', 'x5=10', 'x6=10')]
    qian2011_FS28['reporter'] = 'y'
    qian2011_FS28['exp_results'] = [(60*15, 50), (3600*1.2, 50), (3600*1.5, 50), (3600*1.7, 50), (3600*3, 50), (3600*3.5, 50), (3600*5, 50)]
    setups.append(qian2011_FS28)

    qian2011_FS29 = dict()
    qian2011_FS29['name'] = 'Qian2011-FS29-OR'
    qian2011_FS29['piltemplate'] = qian2011_SF29
    qian2011_FS29['pilparams'] = [None]
    qian2011_FS29['pepperargs'] = {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}
    qian2011_FS29['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=10', 'x3=90', 'x4=10', 'x5=90', 'x6=10', 'x7=90', 'x8=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=10', 'x3=90', 'x4=10', 'x5=90', 'x6=10', 'x7=10', 'x8=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=10', 'x3=90', 'x4=10', 'x5=10', 'x6=10', 'x7=10', 'x8=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=10', 'x3=10', 'x4=10', 'x5=10', 'x6=10', 'x7=90', 'x8=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=10', 'x3=10', 'x4=10', 'x5=90', 'x6=10', 'x7=10', 'x8=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=10', 'x2=10', 'x3=90', 'x4=10', 'x5=10', 'x6=10', 'x7=10', 'x8=10'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=10', 'x3=10', 'x4=10', 'x5=10', 'x6=10', 'x7=10', 'x8=10')]
    qian2011_FS29['reporter'] = 'y'
    qian2011_FS29['exp_results'] = [(60*15, 50), (60*20, 50), (60*30, 50), (60*45, 50), (60*50, 50), (60*60, 50), (60*65, 50), (60*70, 50)]
    setups.append(qian2011_FS29)

    qian2011_FS29AND = dict()
    qian2011_FS29AND['name'] = 'Qian2011-FS29-AND'
    qian2011_FS29AND['piltemplate'] = qian2011_SF29_AND
    qian2011_FS29AND['pilparams'] = [None]
    qian2011_FS29AND['pepperargs'] = {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}
    qian2011_FS29AND['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '36000', '--p0', 'x1=90', 'x2=10', 'x3=90', 'x4=10', 'x5=90', 'x6=10', 'x7=90', 'x8=10')]
    qian2011_FS29AND['reporter'] = 'y'
    qian2011_FS29AND['exp_results'] = [(3600*4, 50)]
    setups.append(qian2011_FS29AND)


    return setups


