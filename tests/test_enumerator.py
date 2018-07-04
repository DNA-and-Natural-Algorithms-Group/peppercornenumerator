#
#  test_enumerator.py
#
# Unittests for the Enumerator Object,
#   I/O using pil / seesaw
#   local-elevation
#

import unittest
import logging
logging.disable(logging.CRITICAL)

from peppercornenumerator import Enumerator
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
        enum.enumerate()

        self.assertTrue(F in [rms.canonical for rms in enum.resting_macrostates])

    def test_local_folding(self):
        pass

if __name__ == '__main__':
  unittest.main()

