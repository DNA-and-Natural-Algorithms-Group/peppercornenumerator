
def yin2008_F3_pil(x):
    return """
    length a = 6
    length b = 6
    length c = 6
    length d = 6
    length u = 6
    length v = 6
    length x = 6
    length y = 6

    # Initially present
    I = a* x* v* b* y*
    A = x*( v*( b*( y*( u* c* a* x* ) ) ) ) a
    B = v* d* y*( u*( c*( a*( x*( v* b* y* ) ) ) ) ) b
    C = c u( y( d( v( u* c* a* x* ) ) ) ) 
    D = d v( x( a( c( u( v* d* y* ) ) ) ) ) b* y*
    
    # Not initially present
    IA   = a*( x*( v*( b*( y*( +                                                   x* v* b* y* u* c* a* x* ) ) ) ) ) @i 0 M
    AB   =                       v* d* y* u* c* a*( x*( v*( b*( y*( x( a( c( u( y( b( + x* v* ) ) ) ) ) ) ) ) ) ) ) @i 0 M
    ABC  = c( u( y( d( v( u* c* a* x* v* d* y* u* + ) ) ) ) ) a*( x*( v*( b*( y*( x( a( c( u( y( b( + x* v* ) ) ) ) ) ) ) ) ) ) )  @i 0 M
    CDB  = c( u( y( d( v( u*( c*( a*( x*( v*( d*( y* u* + ) ) ) ) ) ) ) ) ) ) ) a* x* v* b*( y*( + v* d* y* u*( c*( a*( x*( v* b* y* ) ) ) ) ) ) @i 0 M
    CD   = c( u( y( d( v( u*( c*( a*( x*( v*( d*( y* u* + ) ) ) ) ) ) ) ) ) ) ) a* x* v* b* y* @i 0 M
    CDA  = c( u( y( d( v( u*( c*( a*( x*( v*( d*( y* u* + ) ) ) ) ) ) ) ) ) ) ) a*( x*( v*( b*( y*( + x* v* b* y* u* c* a* x* ) ) ) ) ) @i 0 M

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

    ddG_bind = 0

    yin2008_F3 = dict()
    yin2008_F3['name'] = 'Yin2008-F3'
    yin2008_F3['piltemplate'] = yin2008_F3_pil
    yin2008_F3['pilparams'] = [None]
    yin2008_F3['pepperargs'] = {'condensed': True, 'conc': 'nM', 'release_cutoff': 15, 'k_slow' : 0.00001, 'k_fast': 0.1, 'ddG_bind': ddG_bind}
    yin2008_F3['simulation'] = [
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=20'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=6'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=2'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=0.6'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=0.4'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=0.2'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=0.1'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=0.06'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=0.02'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=0.01'),
            ('pilsimulator', '--nxy', '--atol', '1e-10', '--rtol', '1e-10', '--mxstep', '10000', '--t8', '18000', '--p0', 'A=20', 'B=20', 'C=20', 'D=20', 'I=0')]
    yin2008_F3['reporter'] = 'A'
    yin2008_F3['exp_results'] = [(43,  -10), (1064, -10), (3505, -10), (5341, -10), (6893, -10), (7706, -10), (8842, -10), (9731, -10), (10242, -10), (10696, -10), (10885, -10), (11112, -10)]
    setups.append(yin2008_F3)

    return setups

def main():
    print yin2008_F3_pil(None)


if __name__ == '__main__':
    main()


