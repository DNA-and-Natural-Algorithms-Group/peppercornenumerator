def qian2011_SF22(x):
    return """
    \r# Qian & Winfree, Science, 2011
    \r# 1x = 100 nM, T = 20 *C

    \rINPUT(I) = w[2,5]
    \rOUTPUT(O) = Fluor[6]
    \rseesaw[5, {2}, {6,7}]
    \rreporter[6,5]

    \rconc[w[5,7], 2*c]
    \rconc[g[5, w[5,6]], 1*c]
    """

def qian2011_SF23(x):
    return """
    \rINPUT(I) = w[2,5]
    \rOUTPUT(O) = Fluor[6]
    \rseesaw[5, {2}, {6,7}]
    \rreporter[6,5]

    \rconc[w[5,7], 2*c]
    \rconc[g[5, w[5,6]], 1*c]
    \rconc[th[w[2,5],5], 0.55*c] # simulated using * 1.1
    """

def qian2011_F2(x):
    return """
    \rINPUT(x1) = w[1,2]
    \rINPUT(x2) = w[3,2]
    \rOUTPUT(y) = Fluor[6]
    \rseesaw[2, {1,3}, {5}]
    \rseesaw[5, {2}, {6,7}]
    \rreporter[6,5]

    \rconc[w[5,7], 2*c]
    \rconc[g[5, w[5,6]], 1*c]
    \rconc[th[w[2,5],5], 0.66*c] # OR gate simulated using * 1.1
    \r#conc[th[w[2,5],5], 1.32*c] # AND gate simulated using * 1.1
    \rconc[g[2,w[2,5]], 2*c]
    """

def qian2011_SF26(x):
    return """
    \rINPUT(x1) = w[17,16]
    \rINPUT(x2) = w[18,16]
    \rINPUT(x3) = w[9,4]
    \rINPUT(x4) = w[3,2]
    \rOUTPUT(y) = Fluor[6]
    \rseesawOR[2,5,{1,3},{6}]
    \rseesawOR[4,1,{8,9},{2}]
    \rseesawOR[16,8,{17,18},{4}]
    \rreporter[6,5]
    \r"""

def qian2011_SF27(x):
    return """
    \rINPUT(x1) = w[21,20]
    \rINPUT(x2) = w[22,20]
    \rINPUT(x3) = w[18,16]
    \rINPUT(x4) = w[9,4]
    \rINPUT(x5) = w[3,2]
    \rOUTPUT(y) = Fluor[6]

    \rseesawOR[2,5,{1,3},{6}]
    \rseesawOR[4,1,{8,9},{2}]
    \rseesawOR[16,8,{17,18},{4}]
    \rseesawOR[20,17,{21,22},{16}]
    \rreporter[6,5]
    \r"""

def qian2011_SF28(x):
    return """
    \rINPUT(x1) = w[21,20]
    \rINPUT(x2) = w[22,20]
    \rINPUT(x3) = w[18,16]
    \rINPUT(x4) = w[9,4]
    \rINPUT(x5) = w[13,12]
    \rINPUT(x6) = w[14,12]
    \rOUTPUT(y) = Fluor[6]

    \rseesawOR[2,5,{1,3},{6}]
    \rseesawAND[12,3,{13,14},{2}]
    \rseesawOR[4,1,{8,9},{2}]
    \rseesawAND[16,8,{17,18},{4}]
    \rseesawOR[20,17,{21,22},{16}]
    \rreporter[6,5]
    \r"""

def qian2011_SF29_OR(x):
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

    INPUT(x3) = w[22,2]
    INPUT(x4) = w[9,4]
    INPUT(x5) = w[18,16]
    INPUT(x6) = w[21,20]

    seesawOR[12,3,{13,14},{2,4,16,20}]
    seesawOR[2,5,{3,22},{6}]
    seesawOR[4,1,{3,9},{23}]
    seesawOR[16,8,{3,18},{24}]
    seesawOR[20,17,{3,21},{25}]
    reporter[6,5]
    reporter[23,1]
    reporter[24,8]
    reporter[25,17]
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
    ddG_bind = 2.0
    setups = []

    # Default pilsimulator call
    psim = "pilsimulator --no-jacobian --nxy --atol 1e-13 --rtol 1e-13 --mxstep 10000 --t8 36000 --t-lin 50000".split()

    # Default peppercorn parameters?
    qian2011_SF22 = dict()
    qian2011_SF22['name'] = 'Qian2011-SF22'
    qian2011_SF22['piltemplate'] = qian2011_S22
    qian2011_SF22['pilparams'] = [None]
    qian2011_SF22['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind, 'k_slow' : 0.01, 'k_fast': 1}),
            ('condensed', {'condensed': True,  'conc': 'nM', 'ddG_bind': ddG_bind, 'k_slow' : 0.01, 'k_fast': 1}),
            #('seesaw-rxns-design', {'seesawrxns': 'seesaw-design'}),
            #('seesaw-rxns-leak', {'seesawrxns': 'seesaw-leak'}),
            #('seesaw-rxns-utbr', {'seesawrxns': 'seesaw-utbr'}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_SF22['simulation'] = [
            psim + "--pyplot-labels I O --p0 I=100".split(),
            psim + "--pyplot-labels I O --p0 I=90".split(),
            psim + "--pyplot-labels I O --p0 I=80".split(),
            psim + "--pyplot-labels I O --p0 I=70".split(),
            psim + "--pyplot-labels I O --p0 I=60".split(),
            psim + "--pyplot-labels I O --p0 I=50".split(),
            psim + "--pyplot-labels I O --p0 I=40".split(),
            psim + "--pyplot-labels I O --p0 I=30".split(),
            psim + "--pyplot-labels I O --p0 I=20".split(),
            psim + "--pyplot-labels I O --p0 I=10".split()]
    qian2011_SF22['reporter'] = 'O'
    qian2011_SF22['metric'] = 'diagonal-crossing-time'
    qian2011_SF22['cmax'] = 100
    qian2011_SF22['tmax'] = 2700
    qian2011_SF22['exp-results'] = [
            (738.55, 72.89), (766.03, 71.75), (796.95, 70.62), (834.73, 69.18), (886.26, 67.22), (951.53, 65.05), (1040.84, 61.75), (1157.63, 57.32), (1353.44, 50.31), (1700.38, 37.53)]
    #setups.append(qian2011_SF22)

    qian2011_F1E = dict()
    qian2011_F1E['name'] = 'Qian2011-F1E'
    qian2011_F1E['piltemplate'] = qian2011_F1
    qian2011_F1E['pilparams'] = [None]
    qian2011_F1E['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind, 'k_slow' : 0.01, 'k_fast': 1}),
            ('condensed', {'condensed': True,  'conc': 'nM', 'ddG_bind': ddG_bind, 'k_slow' : 0.01, 'k_fast': 1}),
            #('seesaw-rxns-design', {'seesawrxns': 'seesaw-design'}),
            #('seesaw-rxns-leak', {'seesawrxns': 'seesaw-leak'}),
            #('seesaw-rxns-utbr', {'seesawrxns': 'seesaw-utbr'}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_F1E['simulation'] = [
            psim + "--pyplot-labels I O --p0 I=100".split(),
            psim + "--pyplot-labels I O --p0 I=90".split(),
            psim + "--pyplot-labels I O --p0 I=80".split(),
            psim + "--pyplot-labels I O --p0 I=70".split(),
            psim + "--pyplot-labels I O --p0 I=60".split(),
            psim + "--pyplot-labels I O --p0 I=50".split(),
            psim + "--pyplot-labels I O --p0 I=40".split(),
            psim + "--pyplot-labels I O --p0 I=30".split(),
            psim + "--pyplot-labels I O --p0 I=20".split(),
            psim + "--pyplot-labels I O --p0 I=10".split()]
    qian2011_F1E['reporter'] = 'O'
    qian2011_F1E['metric'] = 'diagonal-crossing-time'
    qian2011_F1E['cmax'] = 100
    qian2011_F1E['tmax'] = 10800
    #qian2011_F1E['quadrant_data'] = [
    #        [(617.01, 50.01), 
    #         (797.01, 50.01), 
    #         (1118.45, 50.02), 
    #         (1979.89, 50.03)], 
    #        [(5402.07, 36.20), 
    #         (5380.55, 8.30), 
    #         (5378.41, 5.53), 
    #         (5377.01, 3.71), 
    #         (5362.76, 1.90), 
    #         (5361.82, 0.69)]]
    qian2011_F1E['exp-results'] = [
            (1808, 83.24), 
            (2081, 80.66), 
            (2587, 76.28), 
            (3516, 67.37), 
            (6590, 39.12), 
            (9828, 9.05), 
            (10115, 6.21), 
            (10265, 4.41), 
            (10497, 2.47), 
            (10648, 1.31)]
    #setups.append(qian2011_F1E)

    qian2011_F2C_OR = dict()
    qian2011_F2C_OR['name'] = 'Qian2011-F2C-OR'
    qian2011_F2C_OR['piltemplate'] = qian2011_F2
    qian2011_F2C_OR['pilparams'] = [None]
    qian2011_F2C_OR['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('condensed', {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_F2C_OR['simulation'] = [
            psim + "--pyplot-labels x1 x2 y --p0 x1=90 x2=90".split(),
            psim + "--pyplot-labels x1 x2 y --p0 x1=10 x2=90".split(),
            psim + "--pyplot-labels x1 x2 y --p0 x1=90 x2=10".split(),
            psim + "--pyplot-labels x1 x2 y --p0 x1=10 x2=10".split()]
    qian2011_F2C_OR['reporter'] = 'y'
    qian2011_F2C_OR['metric'] = 'diagonal-crossing-time'
    qian2011_F2C_OR['cmax'] = 100
    qian2011_F2C_OR['tmax'] = 18000
    qian2011_F2C_OR['exp-results'] = [(3245, 85), (4733, 68), (4733, 68), (10530, 0.4)]
    #setups.append(qian2011_F2C_OR)

    qian2011_F2C_AND = dict()
    qian2011_F2C_AND['name'] = 'Qian2011-F2C-AND'
    qian2011_F2C_AND['piltemplate'] = qian2011_F2
    qian2011_F2C_AND['pilparams'] = [None]
    qian2011_F2C_AND['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('condensed', {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_F2C_AND['simulation'] = [ #T2_5 = 120 * 1.1 = 132
            psim + "--t8 45000 --pyplot-labels x1 x2 y --p0 x1=90 x2=90 T2_5=132".split(),
            psim + "--t8 45000 --pyplot-labels x1 x2 y --p0 x1=10 x2=90 T2_5=132".split(),
            psim + "--t8 45000 --pyplot-labels x1 x2 y --p0 x1=90 x2=10 T2_5=132".split(),
            psim + "--t8 45000 --pyplot-labels x1 x2 y --p0 x1=10 x2=10 T2_5=132".split()]
    qian2011_F2C_AND['reporter'] = 'y'
    qian2011_F2C_AND['metric'] = 'diagonal-crossing-time'
    qian2011_F2C_AND['cmax'] = 100
    qian2011_F2C_AND['tmax'] = 43200
    qian2011_F2C_AND['exp-results'] = [(13797, 68), (41022, 0.7), (41022, 0.7), (42239, 0.4)]
    #setups.append(qian2011_F2C_AND)

    qian2011_FS26 = dict()
    qian2011_FS26['name'] = 'Qian2011-SF26'
    qian2011_FS26['piltemplate'] = qian2011_SF26
    qian2011_FS26['pilparams'] = [None]
    qian2011_FS26['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('condensed', {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_FS26['simulation'] = [
            psim + "--pyplot-labels x1 x2 x3 x4 y --p0 x1=90 x2=90 x3=90 x4=90".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 y --p0 x1=10 x2=10 x3=90 x4=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 y --p0 x1=10 x2=90 x3=10 x4=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 y --p0 x1=10 x2=10 x3=10 x4=10".split()]
    qian2011_FS26['reporter'] = 'y'
    qian2011_FS26['metric'] = 'diagonal-crossing-time'
    qian2011_FS26['cmax'] = 100
    qian2011_FS26['tmax'] = 21600
    qian2011_FS26['exp-results'] = [(3326.32, 84.49), (7368.42, 65.71), (9557.89, 55.51), (20842.11, 3.67)]
    #setups.append(qian2011_FS26)

    qian2011_FS27 = dict()
    qian2011_FS27['name'] = 'Qian2011-SF27'
    qian2011_FS27['piltemplate'] = qian2011_SF27
    qian2011_FS27['pilparams'] = [None]
    qian2011_FS27['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('condensed', {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_FS27['simulation'] = [
            psim + "--pyplot-labels x1 x2 x3 x4 x5 y --p0 x1=90 x2=90 x3=90 x4=90 x5=90".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 y --p0 x1=10 x2=10 x3=90 x4=10 x5=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 y --p0 x1=10 x2=90 x3=10 x4=10 x5=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 y --p0 x1=10 x2=10 x3=10 x4=10 x5=10".split()]
    qian2011_FS27['reporter'] = 'y'
    qian2011_FS27['metric'] = 'diagonal-crossing-time'
    qian2011_FS27['cmax'] = 100
    qian2011_FS27['tmax'] = 21600
    qian2011_FS27['exp-results'] = [(3558.62, 82.89), (9848.28, 53.83), (12082.76, 43.88), (20937.93, 3.28)]
    #setups.append(qian2011_FS27)

    qian2011_FS28 = dict()
    qian2011_FS28['name'] = 'Qian2011-SF28'
    qian2011_FS28['piltemplate'] = qian2011_SF28
    qian2011_FS28['pilparams'] = [None]
    qian2011_FS28['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('condensed', {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('seesaw-design', {'seesawrxns': 'seesaw-design', 'conc': 'nM'}),
            ('seesaw-utbr', {'seesawrxns': 'seesaw-utbr', 'conc': 'nM'}),
            ('seesaw-leak', {'seesawrxns': 'seesaw-leak', 'conc': 'nM'}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw', 'conc': 'nM'})]
    qian2011_FS28['simulation'] = [
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=90 x3=90 x4=90 x5=90 x6=90".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=10 x3=90 x4=90 x5=10 x6=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=90 x3=10 x4=90 x5=10 x6=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=10 x3=10 x4=90 x5=10 x6=90".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=10 x3=10 x4=10 x5=90 x6=90".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=90 x3=90 x4=10 x5=10 x6=10".split(),

            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=90 x3=10 x4=10 x5=10 x6=90".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=10 x3=10 x4=10 x5=90 x6=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=90 x3=10 x4=10 x5=10 x6=90".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=10 x3=90 x4=10 x5=90 x6=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=10 x3=10 x4=10 x5=10 x6=10".split()]
    qian2011_FS28['reporter'] = 'y'
    qian2011_FS28['metric'] = 'diagonal-crossing-time'
    qian2011_FS28['cmax'] = 100
    qian2011_FS28['tmax'] = 28800
    qian2011_FS28['exp-results'] = [
            (4467.98, 83.76), (7541.50, 73.47), (8224.51, 71.49), (8822.13, 69.50), (12151.78, 57.23), (13944.66, 51.29), 
            (16420.55, 42.57), (28031.62, 2.57), (28031.62, 2.57), (28031.62, 2.57), (28031.62, 2.57), (28031.62, 2.57)]
    setups.append(qian2011_FS28)

    qian2011_FS29 = dict()
    qian2011_FS29['name'] = 'Qian2011-SF29-OR'
    qian2011_FS29['piltemplate'] = qian2011_SF29
    qian2011_FS29['pilparams'] = [None]
    qian2011_FS29['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('condensed', {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_FS29['simulation'] = [ # Threshold = 60 * 1.1
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10 x7=90 x8=10 T2_5=66".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10 x7=10 x8=10 T2_5=66".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=66".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=10 x2=10 x3=10 x4=10 x5=10 x6=10 x7=90 x8=10 T2_5=66".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=10 x2=10 x3=10 x4=10 x5=90 x6=10 x7=10 x8=10 T2_5=66".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=10 x2=10 x3=90 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=66".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=10 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=66".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=10 x2=10 x3=10 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=66".split()]
    qian2011_FS29['reporter'] = 'y'
    qian2011_FS29['metric'] = 'diagonal-crossing-time'
    qian2011_FS29['cmax'] = 100
    qian2011_FS29['tmax'] = 10800
    qian2011_FS29['exp-results'] = [(1925.58, 81.20), (2302.33, 78.40), (2679.07, 74.80), (3809.30, 63.60), (4102.33, 61.20), (4311.63, 59.20), (4479.07, 58.00), (10548.84, 3.60)]
    #setups.append(qian2011_FS29)

    qian2011_FS29AND = dict()
    qian2011_FS29AND['name'] = 'Qian2011-SF29-AND'
    qian2011_FS29AND['piltemplate'] = qian2011_SF29_AND
    qian2011_FS29AND['pilparams'] = [None]
    qian2011_FS29AND['pepperargs'] = [ # Threshold = 320 * 1.1
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('condensed', {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_FS29AND['simulation'] = [
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10 x7=90 x8=10 T2_5=352".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10 x7=10 x8=10 T2_5=352".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=352".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=10 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=352".split()]
    qian2011_FS29AND['reporter'] = 'y'
    qian2011_FS29AND['metric'] = 'diagonal-crossing-time'
    qian2011_FS29AND['cmax'] = 100
    qian2011_FS29AND['tmax'] = 36000
    qian2011_FS29AND['exp-results'] = [(15947.37, 55.33), (31894.74, 11.07), (34263.16, 5.79), (35052.63, 2.94)]
    #setups.append(qian2011_FS29AND)

    # ???
    qian2011_FS30 = dict()
    qian2011_FS30['name'] = 'Qian2011-SF30-OR'
    qian2011_FS30['piltemplate'] = qian2011_SF30_OR
    qian2011_FS30['pilparams'] = [None]
    qian2011_FS30['pepperargs'] = [
            ('detailed',  {'condensed': False, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('condensed', {'condensed': True, 'conc': 'nM', 'ddG_bind': ddG_bind}),
            ('seesaw-rxns', {'seesawrxns': 'seesaw'})]
    qian2011_FS30['simulation'] = [
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 --p0 x1=90 x2=90 x3=10 x4=10 x5=10 x6=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 --p0 x1=90 x2=10 x3=10 x4=10 x5=10 x6=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 --p0 x1=10 x2=90 x3=10 x4=10 x5=10 x6=10".split(),
            psim + "--pyplot-labels x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 --p0 x1=10 x2=10 x3=10 x4=10 x5=10 x6=10".split()]
    qian2011_FS30['reporter'] = None
    qian2011_FS30['metric'] = 'diagonal-crossing-time'
    qian2011_FS30['cmax'] = 100
    qian2011_FS30['tmax'] = 36000
    qian2011_FS30['exp-results'] = [
            ('y1', [(15947.37, 55.33), (31894.74, 11.07), (34263.16, 5.79), (35052.63, 2.94)]),
            ('y2', [(15947.37, 55.33), (31894.74, 11.07), (34263.16, 5.79), (35052.63, 2.94)]),
            ('y3', [(15947.37, 55.33), (31894.74, 11.07), (34263.16, 5.79), (35052.63, 2.94)]),
            ('y4', [(15947.37, 55.33), (31894.74, 11.07), (34263.16, 5.79), (35052.63, 2.94)])]
    #setups.append(qian2011_FS30)


    return setups

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    ddG_bind = {'ddG_bind': 0.0}
    rates = {'k_slow': 0.01, 'k_fast': 1}
    seesaw = {'seesaw-rxns': 'seesaw-T20', 'seesaw-conc': 100e-9}

    # Default pilsimulator call
    psim = "pilsimulator --no-jacobian --nxy --atol 1e-13 --rtol 1e-13 --mxstep 10000 --t8 36000 --t-lin 50000"

    # Setup
    SF22 = FigureData('Qian2011-SF22')
    current = SF22
    template = qian2011_SF22
    sims = [psim + " --pyplot-labels I O --p0 I=100",
            psim + " --pyplot-labels I O --p0 I=90",
            psim + " --pyplot-labels I O --p0 I=80",
            psim + " --pyplot-labels I O --p0 I=70",
            psim + " --pyplot-labels I O --p0 I=60",
            psim + " --pyplot-labels I O --p0 I=50",
            psim + " --pyplot-labels I O --p0 I=40",
            psim + " --pyplot-labels I O --p0 I=30",
            psim + " --pyplot-labels I O --p0 I=20",
            psim + " --pyplot-labels I O --p0 I=10"]
    litr = [(738.55, 72.89), (766.03, 71.75), (796.95, 70.62), (834.73, 69.18), 
            (886.26, 67.22), (951.53, 65.05), (1040.84, 61.75), (1157.63, 57.32), 
            (1353.44, 50.31), (1700.38, 37.53)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'O'
        metric = 'diagonal-crossing-time'
        tmax = '2700'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, 
                reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    SF23 = FigureData('Qian2011-SF23')
    current = SF23
    template = qian2011_SF23
    sims = [psim + " --pyplot-labels I O --p0 I=100",
            psim + " --pyplot-labels I O --p0 I=90",
            psim + " --pyplot-labels I O --p0 I=80",
            psim + " --pyplot-labels I O --p0 I=70",
            psim + " --pyplot-labels I O --p0 I=60",
            psim + " --pyplot-labels I O --p0 I=50",
            psim + " --pyplot-labels I O --p0 I=40",
            psim + " --pyplot-labels I O --p0 I=30",
            psim + " --pyplot-labels I O --p0 I=20",
            psim + " --pyplot-labels I O --p0 I=10"]
    litr = [(1808, 83.24), (2081, 80.66), (2587, 76.28), (3516, 67.37), 
            (6590, 39.12), (9828, 9.05), (10115, 6.21), (10265, 4.41), 
            (10497, 2.47), (10648, 1.31)]


    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'O'
        metric = 'diagonal-crossing-time'
        tmax = '10800'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F2C_OR = FigureData('Qian2011-F2C-OR')
    current = F2C_OR
    template = qian2011_F2
    sims = [psim + " --pyplot-labels x1 x2 y --p0 x1=90 x2=90",
            psim + " --pyplot-labels x1 x2 y --p0 x1=10 x2=90",
            psim + " --pyplot-labels x1 x2 y --p0 x1=90 x2=10",
            psim + " --pyplot-labels x1 x2 y --p0 x1=10 x2=10"]
    litr = [(3245, 85), (4733, 68), (4733, 68), (10530, 0.4)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'y'
        metric = 'diagonal-crossing-time'
        tmax = '18000'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    F2C_AND = FigureData('Qian2011-F2C-AND')
    current = F2C_AND
    template = qian2011_F2
    sims = [ #T2_5 = 120 * 1.1 = 132
            psim + " --t8 45000 --pyplot-labels x1 x2 y --p0 x1=90 x2=90 T2_5=132",
            psim + " --t8 45000 --pyplot-labels x1 x2 y --p0 x1=10 x2=90 T2_5=132",
            psim + " --t8 45000 --pyplot-labels x1 x2 y --p0 x1=90 x2=10 T2_5=132",
            psim + " --t8 45000 --pyplot-labels x1 x2 y --p0 x1=10 x2=10 T2_5=132"]
    litr = [(13797, 68), (41022, 0.7), (41022, 0.7), (42239, 0.4)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'y'
        metric = 'diagonal-crossing-time'
        tmax = '43200'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    SF26 = FigureData('Qian2011-SF26')
    current = SF26
    template = qian2011_SF26
    sims = [psim + " --pyplot-labels x1 x2 x3 x4 y --p0 x1=90 x2=90 x3=90 x4=90",
            psim + " --pyplot-labels x1 x2 x3 x4 y --p0 x1=10 x2=10 x3=90 x4=10",
            psim + " --pyplot-labels x1 x2 x3 x4 y --p0 x1=10 x2=90 x3=10 x4=10",
            psim + " --pyplot-labels x1 x2 x3 x4 y --p0 x1=10 x2=10 x3=10 x4=10"]
    litr = [(3326.32, 84.49), (7368.42, 65.71), (9557.89, 55.51), (20842.11, 3.67)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'y'
        metric = 'diagonal-crossing-time'
        tmax = '21600'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)


    # Setup
    SF27 = FigureData('Qian2011-SF27')
    current = SF27
    template = qian2011_SF27
    sims = [psim + " --pyplot-labels x1 x2 x3 x4 x5 y --p0 x1=90 x2=90 x3=90 x4=90 x5=90",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 y --p0 x1=10 x2=10 x3=90 x4=10 x5=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 y --p0 x1=10 x2=90 x3=10 x4=10 x5=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 y --p0 x1=10 x2=10 x3=10 x4=10 x5=10"]
    litr = [(3558.62, 82.89), (9848.28, 53.83), (12082.76, 43.88), (20937.93, 3.28)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'y'
        metric = 'diagonal-crossing-time'
        tmax = '21600'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    SF28 = FigureData('Qian2011-SF28')
    current = SF28
    template = qian2011_SF28
    sims = [psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=90 x3=90 x4=90 x5=90 x6=90",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=10 x3=90 x4=90 x5=10 x6=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=90 x3=10 x4=90 x5=10 x6=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=10 x3=10 x4=90 x5=10 x6=90",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=10 x3=10 x4=10 x5=90 x6=90",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=90 x3=90 x4=10 x5=10 x6=10",

            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=90 x3=10 x4=10 x5=10 x6=90",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=90 x2=10 x3=10 x4=10 x5=90 x6=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=90 x3=10 x4=10 x5=10 x6=90",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=10 x3=90 x4=10 x5=90 x6=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y --p0 x1=10 x2=10 x3=10 x4=10 x5=10 x6=10"]
    litr = [(4467.98, 83.76), (7541.50, 73.47), (8224.51, 71.49), (8822.13, 69.50), (12151.78, 57.23), (13944.66, 51.29), 
            (16420.55, 42.57), (28031.62, 2.57), (28031.62, 2.57), (28031.62, 2.57), (28031.62, 2.57), (28031.62, 2.57)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'y'
        metric = 'diagonal-crossing-time'
        tmax = '28800'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate: 
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    SF29_OR = FigureData('Qian2011-SF29-OR')
    current = SF29_OR
    template = qian2011_SF29_OR
    sims = [ # Threshold = 60 * 1.1
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10 x7=90 x8=10 T2_5=66",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10 x7=10 x8=10 T2_5=66",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=66",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=10 x2=10 x3=10 x4=10 x5=10 x6=10 x7=90 x8=10 T2_5=66",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=10 x2=10 x3=10 x4=10 x5=90 x6=10 x7=10 x8=10 T2_5=66",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=10 x2=10 x3=90 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=66",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=10 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=66",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=10 x2=10 x3=10 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=66"]
    litr = [(1925.58, 81.20), (2302.33, 78.40), (2679.07, 74.80), (3809.30, 63.60), (4102.33, 61.20), (4311.63, 59.20), (4479.07, 58.00), (10548.84, 3.60)]

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'y'
        metric = 'diagonal-crossing-time'
        tmax = '10800'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    # Setup
    SF29_AND = FigureData('Qian2011-SF29-AND')
    current = SF29_AND
    template = qian2011_SF29_AND
    sims = [ # Threshold = 320 * 1.1
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10 x7=90 x8=10 T2_5=352",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=90 x6=10 x7=10 x8=10 T2_5=352",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=90 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=352",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 x7 x8 y --p0 x1=90 x2=10 x3=10 x4=10 x5=10 x6=10 x7=10 x8=10 T2_5=352"]
    litr = [(15947.37, 55.33), (31894.74, 11.07), (34263.16, 5.79), (35052.63, 2.94)] 

    for (sim, res) in zip(sims, litr):
        pilstring  = template(None)
        simulation = sim
        reporter = 'y'
        metric = 'diagonal-crossing-time'
        tmax = '36000'
        cmax = '100'
        current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['Condensed'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    SF30_OR = FigureData('Qian2011-SF30-OR')
    current = SF30_OR
    template = qian2011_SF30_OR
    sims = [ # Threshold = 320 * 1.1
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 --p0 x1=90 x2=90 x3=10 x4=10 x5=10 x6=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 --p0 x1=90 x2=10 x3=10 x4=10 x5=10 x6=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 --p0 x1=10 x2=90 x3=10 x4=10 x5=10 x6=10",
            psim + " --pyplot-labels x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 --p0 x1=10 x2=10 x3=10 x4=10 x5=10 x6=10"]
    litr = [('y1', [(4034.90, 0.61), (5387.92, 0.48), (5387.92, 0.48), (10413.42, 0.01)])]

    for (rep, results) in litr:
        for (sim, res) in zip(sims, results):
            pilstring  = template(None)
            simulation = sim
            reporter = rep
            metric = 'diagonal-crossing-time'
            tmax = '10800'
            cmax = '1'
            current.add_system_simulation_setup(pilstring, simulation, reporter, ':'.join([metric, tmax, cmax]), res)

    current.pepperargs['default'] = current.pepperargs['CONDENSED'].copy()
    current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(ddG_bind)

    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)


    return [SF22, SF23, F2C_OR, F2C_AND, SF26, SF27, SF28, SF29_OR, SF29_AND, SF30_OR]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)


