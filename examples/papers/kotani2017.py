
def kotani2017_F2_pil(x):
    return """
    length a   = 22
    length b   = 22
    length c   = 22
    length t1  = 6
    length t2  = 6
    length t3  = 10
    length T2  = 2

    length d1s = 16
    length d2  = 6

    S1 = d1s T2 b( a( t2( + ) ) c*( t1* + ) )
    S2 = t1( c( a( + t2* ) b*( d2 t3 + ) ) )
    C1 = t1 c a

    P1 = t2* a*( c*( t1*( + ) ) )
    I1 = d1s T2 b( a t2 + c )
    I2 = d1s T2 b( a( t2( + ) ) b*( d2 t3 + ) c*( t1* + ) )

    P2 = d1s T2 b( a( t2( + ) ) ) d2 t3
    P3 = b( c*( t1* + ) )

    R = d1s( d2( + t3* ) )

    D = d1s d2
    RW = d1s( T2 b( a( t2( + ) ) ) d2( t3( + ) ) )

    """.format(x)

def kotani2017_F3_pil(x):
    return """
    length a   = 22
    length b   = 22
    length c1   = 11
    length c2   = 11
    length e   = 22
    length f   = 22
    length g   = 22
    length t1  = 6
    length t2  = 6
    length t3  = 10
    length t4 = 6
    length t5 = 6
    length T2  = 2

    length d1s = 16
    length d2  = 6

    S3 = t1 c1 c2( f( e( t5( + ) ) g*( t4* + ) ) )
    S4 = t4( g( e( + t5* ) f*( T2 a + c2 ) ) )
    C2 = t4 g e
    P4 = t5* e*( g*( t4*( + ) ) ) 
    P5 = t1 c1 c2 f( e( t5( + ) ) ) T2 a
    P6 = c2( f( g*( t4* + ) ) )


    S1 = d1s T2 b( a( t2( + ) ) c2*( c1*( t1* + ) ) )
    S2 = t1( c1( c2( a( + t2* ) b*( d2 t3 + ) ) ) )

    P1 = t2* a*( c2*( c1*( t1*( + ) ) ) )
    P2 = d1s T2 b( a( t2( + ) ) ) d2 t3
    P3 = b( c2*( c1*( t1* + ) ) )

    R = d1s( d2( + t3* ) )

    C1 = t1 c1 c2 a

    D = d1s d2
    RW = d1s( T2 b( a( t2( + ) ) ) d2( t3( + ) ) )
    """.format(x)

def kotani2017_F4_pil(x):
    return """
    length a   = 22
    length b   = 22
    length o   = 22
    length c1  = 11
    length c2  = 11
    length t1  = 6
    length t2  = 6
    length t3  = 10
    length T2  = 2
    length x   = 2
    length y   = 2

    length d1s = 16
    length d2  = 6

    S5au = o( b*( T2 a + c2 ) a( t2( y + ) ) c2*( c1*( t1* x* + ) ) ) d2 t3
    S6au = y* t2* a*( b*( c2*( + x t1 c1 ) ) o*( + d1s T2 ) c2*( c1*( t1*( + x ) ) ) )
    
    C1x = x t1 c1 c2 a

    P1x  = t2* a*( c2*( c1*( t1*( x*( + ) ) ) ) )
    P2au = c2( b( a( t2( y( + ) ) ) ) )

    P8au = x t1 c1 c2 b( o*( + ) ) T2 a
    P9au = d1s T2 o( c2*( c1*( t1* + ) ) ) d2 t3

    R = d1s( d2( + t3* ) )
    D = d1s d2
    RWau = d1s( T2 o( c2*( c1*( t1* + ) ) ) d2( t3( + ) ) ) 

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

    # If you run this in detailed mode, you need to set --t8 to 1e8
    kotani2017_F2 = dict()
    kotani2017_F2['name'] = 'kotani2017_F2'
    kotani2017_F2['piltemplate'] = kotani2017_F2_pil
    kotani2017_F2['pilparams'] = [None]
    kotani2017_F2['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 10}
    kotani2017_F2['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-13', '--rtol', '1e-13', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'S2=10', 'R=20', 'C1=1'),
            ('pilsimulator', '--nxy', '--atol', '1e-13', '--rtol', '1e-13', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'S2=10', 'R=20', 'C1=0.5'),
            ('pilsimulator', '--nxy', '--atol', '1e-13', '--rtol', '1e-13', '--mxstep', '10000', '--t8', '36000', '--p0', 'S1=10', 'S2=10', 'R=20', 'C1=0.05')]
    kotani2017_F2['reporter'] = 'D'
    kotani2017_F2['exp_results'] = [(7733, 7.42), (11333, 6.18), (25533, 1.40)]
    setups.append(kotani2017_F2)



    # If you run this in detailed mode, you need to set --t8 to 1e8
    kotani2017_F3 = dict()
    kotani2017_F3['name'] = 'kotani2017_F3'
    kotani2017_F3['piltemplate'] = kotani2017_F3_pil
    kotani2017_F3['pilparams'] = [None]
    kotani2017_F3['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 10}
    kotani2017_F3['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '360000', '--p0', 'S1=10', 'S2=10', 'S3=10', 'S4=10', 'R=20', 'C1=0.1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '360000', '--p0', 'S1=10', 'S2=10', 'S3=10', 'S4=10', 'R=20', 'C1=0.01'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '360000', '--p0', 'S1=10', 'S2=10', 'S3=10', 'S4=10', 'R=20', 'C1=0.001')]
    kotani2017_F3['reporter'] = 'D'
    kotani2017_F3['exp_results'] = [(21220, 7.72), (64203, 3.12), (86996, 0.69)]
    setups.append(kotani2017_F3)

    # If you run this in detailed mode, you need to set --t8 to 1e8
    kotani2017_F4 = dict()
    kotani2017_F4['name'] = 'kotani2017_F4'
    kotani2017_F4['piltemplate'] = kotani2017_F4_pil
    kotani2017_F4['pilparams'] = [None]
    kotani2017_F4['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 10}
    kotani2017_F4['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '360000', '--p0', 'S5au=10', 'S6au=10', 'R=20', 'C1x=0.1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '360000', '--p0', 'S5au=10', 'S6au=10', 'R=20', 'C1x=0.01'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '360000', '--p0', 'S5au=10', 'S6au=10', 'R=20', 'C1x=0.001'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '360000', '--p0', 'S5au=10', 'S6au=10', 'R=20', 'C1x=0')]
    kotani2017_F4['reporter'] = 'D'
    kotani2017_F4['exp_results'] = [(6815, 6.06), (9004, 4.78), (10278, 4.03), (10795, 3.73)]
    setups.append(kotani2017_F4)

    return setups

