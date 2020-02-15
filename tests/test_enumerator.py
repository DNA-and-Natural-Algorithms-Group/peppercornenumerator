#
#  test_enumerator.py
#
# Unittests for the Enumerator Object,
#   I/O using pil / seesaw
#

import unittest
import logging
logging.disable(logging.CRITICAL)

from peppercornenumerator.enumerator import Enumerator, enumerate_pil, enumerate_ssw
from peppercornenumerator.input import read_pil, read_seesaw
from peppercornenumerator.objects import PepperReaction, clear_memory, DSDObjectsError
import peppercornenumerator.reactions as reactlib

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestEnumeratorInterface(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_rxn_input(self):
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

        enum = Enumerator(list(complexes.values()), reactions)
        enum.release_cutoff = 7
        enum.enumerate()
        self.assertEqual(len(enum.reactions), len(reactions))


    def test_interface(self):
        #   - initialization
        #   - initial-complexes
        #   - resting-macrostates
        #   - transient-macrostates
        complexes, _ = read_pil("""
        length a = 3
        length n = 1
        length b = 1
        length c = 4
        length ab = 4

        X = a( b( n c( + ) ) )
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

        enum = Enumerator([X, Y])

        enum.dry_run()
        self.assertEqual(sorted(enum.complexes), sorted([X,Y]))
        self.assertEqual(sorted(enum.resting_complexes), sorted([X,Y]))
        self.assertEqual(sorted(r.canonical for r in enum.resting_macrostates), sorted([X,Y]))

        with self.assertRaises(DSDObjectsError) as e:
            enum.enumerate()

        #enum = Enumerator([X, Y])
        #enum.enumerate()
        #with self.assertRaises(DSDObjectsError) as e:
        #    enum.dry_run()

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

        enum = Enumerator([A,B])
        enum.max_complex_count = 1000
        enum.max_reaction_count = 5000
        enum.enumerate()

        self.assertTrue(F in [rms.canonical for rms in enum.resting_macrostates])

class TestWrappers(unittest.TestCase):

    def setUp(self):
        pass

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

        enum, outp = enumerate_pil(Zhang2009_F5, is_file = False, enumfile = None, 
                detailed = True, condensed = True, enumconc = 'nM',
                dG_bp = -1.3, k_fast = 0.1, k_slow = 0.001)

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

if __name__ == '__main__':
  unittest.main()

