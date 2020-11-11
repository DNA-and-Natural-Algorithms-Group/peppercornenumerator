#!/usr/bin/env python
#
#  test_enumerator.py
#  EnumeratorProject
#
import unittest

from peppercornenumerator.input import read_pil, read_seesaw
from peppercornenumerator.objects import clear_memory, PepperComplex, PepperReaction
from peppercornenumerator.enumerator import (PeppercornUsageError,
                                             Enumerator, enumerate_pil, enumerate_ssw)

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestEnumeratorInterface(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_interface_01(self):
        PepperComplex.PREFIX = 'enum'
        complexes, reactions = read_pil("""
        # Domain Specifications
        length a = 8
        length b = 8
        length c = 8
        length d2 = 8
        length d3 = 8
        length t = 4

        # Resting-set Complexes
        e0 = d2( d3( + ) a* + a*( b*( c ) ) ) t* 
        e12 = d2 d3( + ) a*                     
        e13 = t( d2( d3 + a*( b*( c ) ) ) )    
        e21 = d2 d3                           
        e22 = t( d2( d3( + ) a*( + a* b*( c ) ) ) )
        e27 = t( d2( d3( + ) a* + a*( b*( c ) ) ) )
        gate = d2( d3( + ) a*( + a* b*( c ) ) ) t* 
        t23 = t d2 d3   
                                                                                                     
        # Transient Complexes    
        e5 = t( d2 d3 + d2( d3( + ) a*( + a* b*( c ) ) ) )  
        e7 = t( d2 d3 + d2( d3( + ) a* + a*( b*( c ) ) ) ) 
        e18 = t( d2( d3 + d2 d3( + ) a*( + a* b*( c ) ) ) )

        # Detailed Reactions  
        reaction [bind21         =      1.2e+06 /M/s ] e0 + t23 -> e7   
        reaction [branch-3way    =     0.122307 /s   ] e0 -> gate       
        reaction [branch-3way    =      41.6667 /s   ] e5 -> e7         
        reaction [branch-3way    =      41.6667 /s   ] e5 -> e18        
        reaction [open           =      306.345 /s   ] e5 -> t23 + gate 
        reaction [open           =      306.345 /s   ] e7 -> e0 + t23   
        reaction [branch-3way    =     0.122307 /s   ] e7 -> e5         
        reaction [branch-3way    =      41.6667 /s   ] e7 -> e12 + e13  
        reaction [branch-3way    =     0.122307 /s   ] e18 -> e5        
        reaction [branch-3way    =      41.6667 /s   ] e18 -> e12 + e13 
        reaction [branch-3way    =     0.122307 /s   ] e18 -> e22 + e21 
        reaction [branch-3way    =      41.6667 /s   ] e22 -> e27       
        reaction [branch-3way    =     0.122307 /s   ] e27 -> e22       
        reaction [branch-3way    =      41.6667 /s   ] gate -> e0       
        reaction [bind21         =      1.2e+06 /M/s ] t23 + gate -> e5 
        """)

        enum = Enumerator(complexes.values(), reactions)
        enum.release_cutoff = 7
        enum.enumerate()
        self.assertEqual(sorted(enum.reactions), sorted(reactions))

    def test_interface_02(self):
        complexes, _ = read_pil("""
        length a = 3
        length n = 1
        length b = 1
        length c = 4
        length ab = 4

        X = a( b( c( + ) ) )
        Y = ab( c( + ) )

        Xf = a b( c( + ) ) a*
        Xff= a b c( + ) b* a*
        Xb = a( b( c + c* ) )
        Y1 = ab c( + ) ab*
        Y2 = ab( c + c* )
        """)
        X = complexes['X']
        Y = complexes['Y']

        Xf = complexes['Xf']
        Xff= complexes['Xff']
        Xb = complexes['Xb']
        Y1 = complexes['Y1']
        Y2 = complexes['Y2']

        enum = Enumerator([X, Y], named_complexes = [X, Y, Xf, Xff, Xb, Y1, Y2])
        enum.max_helix = False
        enum.dry_run()
        self.assertEqual(sorted(enum.complexes), sorted([X,Y]))
        self.assertEqual(sorted(enum.resting_complexes), sorted([X,Y]))
        self.assertEqual(sorted(r.representative for r in enum.resting_macrostates), 
                         sorted([X, Y]))
        enum.enumerate()
        assert len(list(enum.complexes)) == 16
        with self.assertRaises(PeppercornUsageError) as e:
            enum.dry_run()

    def test_named_complexes(self):
        complexes, reactions = read_pil("""
        sequence t = GGAGCC
        sequence s = ATATAT
        sequence r = GCGCGC
        sequence d2 = GGCAAACAAG
        sequence d3 = CGGCAGAATT
        sequence a = CGCATTTGCC
        sequence b = TACCTTTTCC
        sequence c = CAAAGCCCTT

        A = t d2 s* d3
        B = d2( d3( + ) s*( a* + a*( b*( c ) ) ) ) t*
        B2 = d2( d3( + ) s* a* + a*( b*( c ) ) s ) t* @ initial 0 M
        B3 = d2( d3( + ) s*( a*( + a* b*( c ) ) ) ) t* @ initial 0 M
        B4 = d2( d3( + ) s* a*( + a* b*( c ) ) s ) t* @ initial 0 M

        C = d2 d3( + ) s* a* @initial 0 M
        D = t( d2( s*( d3 + a*( b*( c ) ) ) ) ) @i 0 M
        E = d2 d3 @i 0 M
        F = t( d2( s*( d3( + ) s* a*( + a* b*( c ) ) ) ) ) @i 0 M
        G = t d2 s* d3( + ) s* a* @i 0 M
        """)
        A = complexes['A']
        B = complexes['B']
        F = complexes['F']

        enum = Enumerator([A, B], named_complexes = [A, B, F])
        enum.max_complex_count = 1000
        enum.max_reaction_count = 5000
        enum.enumerate()
        self.assertTrue(F in [rms.representative for rms in enum.resting_macrostates])

@unittest.skipIf(SKIP, "skipping tests")
class TestWrappers(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_enumerate_pil(self):
        Zhang2009_F5 = """ 
        # Domains
        length a = 16
        length bm = 5  # m
        length br = 1  # m
        length bc = 7
        length bt  = 7
        length n = 2  # n
        length c = 13  # 15-n

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
        """
        enum, outp = enumerate_pil(Zhang2009_F5, 
                is_file = False, enumfile = None, 
                detailed = True, condensed = True, 
                enumconc = 'nM', dG_bp = -1.3, 
                k_fast = 0.1, k_slow = 0.001)

        assert isinstance(enum,  Enumerator)
        assert isinstance(outp,  str)

    def test_enumerate_seesaw(self):
        Qian2011_F2C_OR = """
        INPUT(x1) = w[1,2]
        INPUT(x2) = w[3,2]
        OUTPUT(y) = Fluor[6]
        seesaw[2, {1,3}, {5}]
        seesaw[5, {2}, {6,7}]
        reporter[6,5]
        
        conc[w[5,7], 2*c]
        conc[g[5, w[5,6]], 1*c]
        conc[th[w[2,5],5], 0.66*c] # OR gate simulated using * 1.1
        conc[g[2,w[2,5]], 2*c]
        """
        enum, outp = enumerate_ssw(Qian2011_F2C_OR, is_file = False,
            ssw_expl = False,
            ssw_conc = 100e-9,
            ssw_rxns = 'T20-utbr-leak-reduced',
            dry_run = True,
            enumfile = None,
            detailed = True, 
            condensed = False, 
            enumconc = 'nM')

        assert isinstance(enum,  Enumerator)
        assert isinstance(outp,  str)

@unittest.skipIf(SKIP, "skipping tests")
class TestEnumeratorOutput(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_bind_and_displace3way(self):
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length t = 5

        I = a b c
        J = b c
        C = b( c( + ) ) a*
        B = a( b c + b( c( + ) ) )
        D = a( b( c( + ) ) )
        """)
        I = complexes['I']
        J = complexes['J']
        C = complexes['C']
        B = complexes['B']
        D = complexes['D']

        # DSD-pathway "bind21"
        path1 = PepperReaction([I, C], [B], 'bind21')
        # DSD-pathway "branch3way"
        path2 = PepperReaction([B], [D, J], 'branch-3way')

        enum = Enumerator(list(complexes.values()))
        enum.enumerate()
        self.assertEqual(sorted(enum.reactions), sorted([path1, path2]))

    def test_cooperative_binding(self):
        complexes, reactions = read_pil("""
        length a = 5
        length x = 10
        length y = 10
        length b = 5

        C = x( y( + b* ) ) a*
        L = a x
        R = y b
        T = x y

        LC = a( x + x( y( + b* ) ) )
         CR = x( y( + y b( + ) ) ) a*
        LCR = a( x + x( y( + y b( + ) ) ) )
        LCF = a( x( + x y( + b* ) ) )
         CRF = x( y + y( b( + ) ) ) a*

        LCRF1 = a( x( + x y( + y b( + ) ) ) )
        LCRF2 = a( x + x( y + y( b( + ) ) ) )
        LR = a( x( + y( b( + ) ) ) )
        """)

        C = complexes['C']
        L = complexes['L']
        R = complexes['R']
        T = complexes['T']
        LC = complexes['LC']
        LCF = complexes['LCF']
        CR = complexes['CR']
        CRF = complexes['CRF']
        LCRF1 = complexes['LCRF1']
        LR = complexes['LR']
        
        path1 = PepperReaction([L, C], [LC], 'bind21')
        path1r= PepperReaction([LC], [L, C], 'open')
        path2 = PepperReaction([LC], [LCF], 'branch-3way')
        path3 = PepperReaction([R, LCF], [LCRF1], 'bind21')
        path4 = PepperReaction([LCRF1], [LR, T], 'branch-3way')

        enum = Enumerator(list(complexes.values()))
        enum.k_fast = float('inf')
        enum.k_slow = 0
        enum.max_helix = True
        enum.enumerate()
        self.assertEqual(len(list(enum.reactions)), 22)

    def test_simple(self):
        complexes, reactions = read_pil("""
        length a = 6
        length a1 = 2
        length a2 = 2
        length a3 = 2
        length b  = 24
        length b1 = 8
        length b2 = 8
        length b3 = 8
        length c  = 24
        length c1 = 8
        length c2 = 8
        length c3 = 8

        I = a b c
        C = b( c( + ) ) a*
        J = a( b c + b( c( + ) ) )
        D = a( b( c( + ) ) )

        cI = a1 a2 a3 b1 b2 b3 c1 c2 c3
        cC = b1( b2( b3( c1( c2( c3( + ) ) ) ) ) ) a3* a2* a1*
        cJ = a1( a2( a3( b1 b2 b3 c1 c2 c3 + b1( b2( b3( c1( c2( c3( + ) ) ) ) ) ) ) ) )
        cD = a1( a2( a3( b1( b2( b3( c1( c2( c3( + ) ) ) ) ) ) ) ) )
        """)

        enum = Enumerator([complexes['I'], complexes['C'], complexes['J'], complexes['D']],
                named_complexes = list(complexes.values()))
        enum.k_fast = 0
        enum.k_slow = 0
        enum.max_helix = True
        enum.enumerate()

        enum2 = Enumerator([complexes['cI'], complexes['cC'], complexes['cJ'], complexes['cD']],
                named_complexes = list(complexes.values()))
        enum.k_fast = 0
        enum2.k_fast = 0
        enum2.k_slow = 0
        enum2.max_helix = True
        enum2.enumerate()

        self.assertEqual(len(list(enum2.reactions)), len(list(enum.reactions)))

    def test_max_helix_01(self):
        complexes, reactions = read_pil("""
        length a = 15
        length x = 15
        length x1 = 15
        length x2 = 15
        length y = 15
        length y1 = 15
        length y2 = 15
        length z = 15
        length z1 = 15
        length z2 = 15

        # should be one reaction, is one
        A1 = x( y z + y( z( + ) ) )
        A1_2 = x( y( z( + ) ) )
        YZ = y z

        # should be one reactions, is one
        B1 = x1( x2( y1 y2 z1 z2 + y1( y2( z1( z2( + ) ) ) ) ) ) 
        B1_2 = x1( x2( y1( y2( z1( z2( + ) ) ) ) ) ) 
        YZ2 = y1 y2 z1 z2

        # should be two reactions, is two
        A2   = x( y z + y( + z( + ) ) )
        A2_1 = x( y( z + z( + ) ) )
        #A2_2 = x( y( z( + ) ) ) # = A1_2
        Y1 = y
        Z1 = z

        # should be two reactions, is two
        B2 = x1( x2( y1 y2 z1 z2 + y1( y2( + z1( z2( + ) ) ) ) ) ) 
        B2_1 = x1( x2( y1( y2( z1 z2 + z1( z2( + ) ) ) ) ) ) 
        Y2 = y1 y2
        Z2 = z1 z2

        # should be two reactions, is two
        C = x( y z + y( + a( + ) z( + ) ) )
        C1 = x( y( z + a( + ) z( + ) ) )
        r1 = a( + ) z

        """)

        A1 = complexes['A1']
        A1_2 = complexes['A1_2']
        YZ = complexes['YZ']

        A2 = complexes['A2']
        A2_1 = complexes['A2_1']
        #A2_2 = complexes['A2_2']
        Y1 = complexes['Y1']
        Z1 = complexes['Z1']

        enum = Enumerator([A1, A2])
        enum.k_fast = 0
        enum.k_slow = 0
        enum.max_helix = True
        enum.enumerate()

        path1 = PepperReaction([A1], sorted([A1_2, YZ]), 'branch-3way')
        path2 = PepperReaction([A2], sorted([A2_1, Y1]), 'branch-3way')
        path3 = PepperReaction([A2_1], sorted([A1_2, Z1]), 'branch-3way')
        self.assertEqual(sorted(enum.reactions), sorted([path1, path2, path3]))
 
        B1 = complexes['B1']
        B1_2 = complexes['B1_2']
        YZ2 = complexes['YZ2']

        B2 = complexes['B2']
        B2_1 = complexes['B2_1']
        Y2 = complexes['Y2']
        Z2 = complexes['Z2']

        enum = Enumerator([B1, B2])
        enum.k_fast = 0
        enum.k_slow = 0
        enum.max_helix = True
        enum.enumerate()

        path1 = PepperReaction([B1], sorted([B1_2, YZ2]), 'branch-3way')
        path2 = PepperReaction([B2], sorted([B2_1, Y2]), 'branch-3way')
        path3 = PepperReaction([B2_1], sorted([B1_2, Z2]), 'branch-3way')
        self.assertEqual(sorted(enum.reactions), sorted([path1, path2, path3]))

        C = complexes['C']
        enum = Enumerator([C])
        enum.k_fast = 0
        enum.k_slow = 0
        enum.max_helix = True
        enum.enumerate()
        self.assertEqual(len(list(enum.reactions)), 2)

    def test_max_helix_02(self):
        complexes, reactions = read_pil("""
        length d1 = 15
        length d4 = 15
        length d6 = 15
        length d7 = 15
        length h8 = 15
        length t0 = 6
        length t2 = 6
        length t3 = 6
        length t5 = 6

        # Initial Complexes
        B2 = d7 t3 d4 t5 
        helper = t3 d7 t3
        PR_FL_B2  = d1 t2( d6( + d7( t3( d4 t5 + ) ) t3* ) )  @ initial 0 M
        
        # Intermediate Complexes
        PR_FLh1B2 = d1 t2( d6( + t3( d7 t3 + d7( t3( d4 t5 + ) ) ) ) )  @ initial 0 M
        PR_FLh2B2 = d1 t2( d6( + t3 d7 t3( + d7( t3( d4 t5 + ) ) ) ) )  @ initial 0 M
        PR_FL_h1w = d1 t2( d6( + t3( d7( t3( + ) ) ) ) )  @ initial 0 M
        
        # sidestuff
        PR_FLB2B2 = d1 t2( d6( + d7 t3( d4 t5 + d7( t3( d4 t5 + ) ) ) ) )  @ initial 0 M

        # casey-semantics
        PR_FLh2B2_v2 = d1 t2( d6( + t3( d7( t3 + d7 t3( d4 t5 + ) ) ) ) )  @ initial 0 M
        PR_FLh2w =     d1 t2( d6( + t3( d7( t3 + t3* ) ) ) )  @ initial 0 M

        """)
        B2 = complexes['B2']
        helper = complexes['helper']
        PR_FL_B2 = complexes['PR_FL_B2']
        PR_FLh1B2 = complexes['PR_FLh1B2']
        PR_FLh2B2 = complexes['PR_FLh2B2']
        PR_FL_h1w = complexes['PR_FL_h1w']
        PR_FLB2B2 = complexes['PR_FLB2B2']
        PR_FLh2B2_v2 = complexes['PR_FLh2B2_v2']
        PR_FLh2w = complexes['PR_FLh2w']
 
        path1  = PepperReaction([PR_FL_B2, helper], [PR_FLh1B2], 'bind21')
        path1r = PepperReaction([PR_FLh1B2], [PR_FL_B2, helper], 'open')
        path2  = PepperReaction([PR_FL_B2, helper], [PR_FLh2B2], 'bind21')
        path2r = PepperReaction([PR_FLh2B2], [PR_FL_B2, helper], 'open')
        path3  = PepperReaction([PR_FLh1B2], [PR_FL_h1w, B2], 'branch-3way')
        path4  = PepperReaction([PR_FL_B2, B2], [PR_FLB2B2], 'bind21')
        path4r = PepperReaction([PR_FLB2B2], [PR_FL_B2, B2], 'open')
        path5  = PepperReaction([PR_FLh1B2], [PR_FLh2B2], 'branch-3way')
        path6  = PepperReaction([PR_FLh2B2], [PR_FLh1B2], 'branch-3way')

        enum = Enumerator([B2, helper, PR_FL_B2])
        enum.max_helix = True
        enum.enumerate()
        assert sorted(enum.reactions) == sorted([path1, path1r,
                                                 path2, path2r,
                                                 path3, 
                                                 path4, path4r,
                                                 path5, path6])

    def test_self_displacement_01(self):
        complexes, reactions = read_pil("""
        length x = 10
        length y = 10
        B1 = x( y( x( y x + ) ) )
        B2 = x y x( y( x( + ) ) )
        B3 = x( y( x y x( + ) ) )
        B4 = x( y x y( x( + ) ) )
        """)
        B1 = complexes['B1']
        B2 = complexes['B2']
        B3 = complexes['B3']
        B4 = complexes['B4']
        path1  = PepperReaction([B1], [B4], 'branch-3way')
        path1r = PepperReaction([B4], [B1], 'branch-3way')
        path2  = PepperReaction([B4], [B2], 'branch-3way')
        path3  = PepperReaction([B2], [B3], 'branch-3way')
        path3r = PepperReaction([B3], [B2], 'branch-3way')
        path4  = PepperReaction([B3], [B1], 'branch-3way')

        enum = Enumerator([B1])
        enum.max_helix = True
        enum.enumerate()
        self.assertEqual(sorted(enum.reactions),
                         sorted([path1, path1r, path2, path3, path3r, path4]))

    def test_self_displacement_02(self):
        complexes, reactions = read_pil("""
        length x = 10
        length y = 10
        T = x( y x + ) y* x*
        T1 = x( y x + x* y* )
        T2 = x y x( + ) y* x*
        T3 = x( y( x( + ) ) )
        T4 = x y x( + x* y* )
        """)
        T = complexes['T']
        T1 = complexes['T1']
        T2 = complexes['T2']
        T3 = complexes['T3']
        T4 = complexes['T4']
        path1  = PepperReaction([T], [T1], 'branch-3way')
        path1r = PepperReaction([T1], [T], 'branch-3way')
        path2  = PepperReaction([T], [T2], 'branch-3way')
        path2r = PepperReaction([T2], [T], 'branch-3way')
        path3  = PepperReaction([T1], [T3], 'bind11')
        path4  = PepperReaction([T2], [T3], 'bind11')
        path5  = PepperReaction([T1], [T4], 'branch-3way')
        path5r = PepperReaction([T4], [T1], 'branch-3way')
        path6  = PepperReaction([T2], [T4], 'branch-3way')
        path6r = PepperReaction([T4], [T2], 'branch-3way')

        enum = Enumerator([T])
        enum.max_helix = True
        enum.enumerate()
        self.assertEqual(sorted(enum.reactions),
                         sorted([path1, path1r, path2, path2r, path3, path4, 
                                 path5, path5r, path6, path6r]))

    def test_self_displacement_03(self):
        complexes, reactions = read_pil("""
        length x1 = 10
        length x2 = 10
        length y1 = 10
        length y2 = 10
        B1 = x1( x2( y1( y2( x1( x2( y1 y2 x1 x2 + ) ) ) ) ) )
        B2 = x1 x2 y1 y2 x1( x2( y1( y2( x1( x2( + ) ) ) ) ) )
        i1 = x1( x2( y1 y2 x1 x2 y1( y2( x1( x2( + ) ) ) ) ) )
        i2 = x1( x2( y1( y2( x1 x2 y1 y2 x1( x2( + ) ) ) ) ) )
        """)
        B1 = complexes['B1']
        B2 = complexes['B2']
        i1 = complexes['i1']
        i2 = complexes['i2']
        path1  = PepperReaction([B1], [i1], 'branch-3way')
        path1r = PepperReaction([i1], [B1], 'branch-3way')
        path1f = PepperReaction([i1], [B2], 'branch-3way')
        path2  = PepperReaction([B2], [i2], 'branch-3way')
        path2r = PepperReaction([i2], [B2], 'branch-3way')
        path2f = PepperReaction([i2], [B1], 'branch-3way')

        enum = Enumerator([B1])
        enum.max_helix = True
        enum.enumerate()
        self.assertEqual(sorted(enum.reactions),
                         sorted([path1, path1r, path1f, path2, path2r, path2f]))


if __name__ == '__main__':
  unittest.main()

