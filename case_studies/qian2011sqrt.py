
def qian2011_SF31(x):
    return """
INPUT(x1_0) = w[45,42]
INPUT(x1_1) = w[46,41]
INPUT(x2_0) = w[47,42]
INPUT(x2_1) = w[48,41]
INPUT(x3_0) = w[49,33]
INPUT(x3_1) = w[50,35]
INPUT(x4_0) = w[51,37]
INPUT(x4_1) = w[52,38]

OUTPUT(y1_0) = Fluor[6]
OUTPUT(y1_1) = Fluor[23]
OUTPUT(y2_0) = Fluor[24]
OUTPUT(y2_1) = Fluor[25]

seesawOR[10, 1, {21, 27}, {23}]
seesawAND[53, 5, {18, 22}, {6}]
seesawOR[34, 18, {28, 33, 37}, {53}]
seesawAND[36, 21, {29, 35, 38}, {10}]
seesawOR[39, 22, {29, 31}, {53}]
seesawAND[40, 27, {30, 28}, {10}]

seesawOR[41, 28, {46, 48}, {34, 40}]
seesawAND[42, 29, {45, 47}, {36, 39}]
seesawOR[43, 30, {33, 38}, {40}]
seesawAND[44, 31, {35, 37}, {39}]
seesawOR[20, 8, {35, 38}, {25}]
seesawAND[26, 13, {33, 37}, {24}]

inputfanout[33, 49, {34, 43, 26}]
inputfanout[35, 50, {36, 44, 20}]
inputfanout[37, 51, {34, 44, 26}]
inputfanout[38, 52, {36, 43, 20}]

reporter[6, 5]
reporter[23, 1]
reporter[24, 13]
reporter[25, 8]
"""

def data(evaluate=False, verbose = 0):
    from figure_analysis import FigureData

    rates = {'k_slow': 0.01, 'k_fast': 1}
    seesaw_on = {'ssw_rxns': 'seesaw-T25-utbr-leak-reduced',
            'ssw_expl': False,
            'ssw_conc': 50e-9,
            'dry_run': True}
    seesaw_off = {'ssw_rxns': None,
            'ssw_expl': False,
            'ssw_conc': 50e-9,
            'dry_run': False}

    # Default pilsimulator call
    psim = "pilsimulator --no-jacobian --nxy --header --atol 1e-13 --rtol 1e-13 --mxstep 1000 --t8 36000 --t-lin 18000"

    # Setup
    SF31 = FigureData('Sup. Fig. 31: Integer 4-bit square-root circuit')
    SF31.fname = 'Qian2011-SF31'
    current = SF31
    template = qian2011_SF31
    sims = [psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=45 x4_1=5  x3_0=45 x3_1=5  x2_0=45 x2_1=5  x1_0=45 x1_1=5',  # 0000
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=45 x4_1=5  x3_0=45 x3_1=5  x2_0=45 x2_1=5  x1_0=5 x1_1=45',  # 0001
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=45 x4_1=5  x3_0=45 x3_1=5  x2_0=5 x2_1=45  x1_0=45 x1_1=5',  # 0010
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=45 x4_1=5  x3_0=45 x3_1=5  x2_0=5 x2_1=45  x1_0=5 x1_1=45',  # 0011
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=45 x4_1=5  x3_0=5 x3_1=45  x2_0=45 x2_1=5  x1_0=45 x1_1=5',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=45 x4_1=5  x3_0=5 x3_1=45  x2_0=45 x2_1=5  x1_0=5 x1_1=45',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=45 x4_1=5  x3_0=5 x3_1=45  x2_0=5 x2_1=45  x1_0=45 x1_1=5',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=45 x4_1=5  x3_0=5 x3_1=45  x2_0=5 x2_1=45  x1_0=5 x1_1=45',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=5 x4_1=45  x3_0=45 x3_1=5  x2_0=45 x2_1=5  x1_0=45 x1_1=5',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=5 x4_1=45  x3_0=45 x3_1=5  x2_0=45 x2_1=5  x1_0=5 x1_1=45',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=5 x4_1=45  x3_0=45 x3_1=5  x2_0=5 x2_1=45  x1_0=45 x1_1=5',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=5 x4_1=45  x3_0=45 x3_1=5  x2_0=5 x2_1=45  x1_0=5 x1_1=45',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=5 x4_1=45  x3_0=5 x3_1=45  x2_0=45 x2_1=5  x1_0=45 x1_1=5',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=5 x4_1=45  x3_0=5 x3_1=45  x2_0=45 x2_1=5  x1_0=5 x1_1=45',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=5 x4_1=45  x3_0=5 x3_1=45  x2_0=5 x2_1=45  x1_0=45 x1_1=5',
            psim + ' --pyplot-labels y1_0 y1_1 y2_0 y2_1 --p0 x4_0=5 x4_1=45  x3_0=5 x3_1=45  x2_0=5 x2_1=45  x1_0=5 x1_1=45']

    simargs = ['0000', 
               '0001', 
               '0010', 
               '0011', 
               '0100', 
               '0101', 
               '0110', 
               '0111', 
               '1000', 
               '1001', 
               '1010', 
               '1011', 
               '1100', 
               '1101', 
               '1110', 
               '1111']

    litr = [('y1_0', [(4.2*3600, 50*0.42), 
                      (7.1*3600, 50*0.1), 
                      (7.1*3600, 50*0.1), 
                      (7.1*3600, 50*0.1), 
                      (3.9*3600, 50*0.55), 
                      (4.1*3600, 50*0.45), 
                      (4.1*3600, 50*0.45), 
                      (4.1*3600, 50*0.45), 
                      (4.5*3600, 50*0.42), 
                      (7.3*3600, 50*0.1), 
                      (7.3*3600, 50*0.1), 
                      (7.3*3600, 50*0.1), 
                      (7.3*3600, 50*0.1), 
                      (7.3*3600, 50*0.1), 
                      (7.3*3600, 50*0.1), 
                      (7.3*3600, 50*0.1)]),
            ('y1_1', [(7.1*3600, 50*0.1), 
                      (5.0*3600, 50*0.4), 
                      (5.0*3600, 50*0.4), 
                      (4.2*3600, 50*0.42), 
                      (7.3*3600, 50*0.08), 
                      (7.3*3600, 50*0.08), 
                      (7.3*3600, 50*0.08), 
                      (7.3*3600, 50*0.08), 
                      (6.9*3600, 50*0.1), 
                      (4.0*3600, 50*0.5), 
                      (4.0*3600, 50*0.5), 
                      (3.8*3600, 50*0.6), 
                      (6.0*3600, 50*0.25), 
                      (4.0*3600, 50*0.51), 
                      (4.0*3600, 50*0.51), 
                      (3.5*3600, 50*0.6)]),
            ('y2_0', [(3.1*3600, 50*0.55), 
                      (3.0*3600, 50*0.6), 
                      (3.0*3600, 50*0.6), 
                      (3.0*3600, 50*0.6), 
                      (6.8*3600, 50*0.1), 
                      (6.8*3600, 50*0.1), 
                      (6.8*3600, 50*0.1), 
                      (6.8*3600, 50*0.1), 
                      (7.0*3600, 50*0.11), 
                      (7.0*3600, 50*0.11), 
                      (7.0*3600, 50*0.11), 
                      (7.0*3600, 50*0.11), 
                      (7.2*3600, 50*0.08), 
                      (7.2*3600, 50*0.08), 
                      (7.2*3600, 50*0.08), 
                      (7.2*3600, 50*0.08)]),
            ('y2_1', [(7.0*3600, 50*0.10), 
                      (7.0*3600, 50*0.10), 
                      (7.0*3600, 50*0.10), 
                      (7.0*3600, 50*0.10), 
                      (2.0*3600, 50*0.75), 
                      (2.0*3600, 50*0.75), 
                      (2.0*3600, 50*0.75), 
                      (2.0*3600, 50*0.75), 
                      (1.8*3600, 50*0.8), 
                      (1.8*3600, 50*0.8), 
                      (1.8*3600, 50*0.8), 
                      (1.8*3600, 50*0.8), 
                      (2.0*3600, 50*0.8), 
                      (2.0*3600, 50*0.8), 
                      (2.0*3600, 50*0.8), 
                      (2.0*3600, 50*0.8), 
                      (1.5*3600, 50*0.85), 
                      (1.5*3600, 50*0.85), 
                      (1.5*3600, 50*0.85), 
                      (1.5*3600, 50*0.85)])]

    for (rep, results) in litr:
        for (sim, arg, res) in zip(sims, simargs, results):
            pilstring  = template(None)
            simulation = sim
            reporter = rep
            metric = 'diagonal-crossing-time'
            tmax = '36000'
            cmax = '50'
            current.add_system_simulation_setup(pilstring, simulation, reporter, metric, ':'.join([tmax, cmax]), res, simargs=rep + '-' + arg)

    current.pepperargs['default'] = current.pepperargs['CONDENSED'].copy()
    #current.pepperargs['default'].update(rates)
    current.pepperargs['default'].update(seesaw_off)
    if evaluate:
        current.eval()

    if verbose:
        for df in current.get_dataframes():
            print(df)

    return [SF31]

if __name__ == '__main__':
    data(evaluate=True, verbose=1)


