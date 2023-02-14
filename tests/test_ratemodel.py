import unittest
import peppercornenumerator.reactions as rxn
from peppercornenumerator.input import read_pil, read_pil_line
from peppercornenumerator.objects import clear_memory

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class RateBranch3WayTests(unittest.TestCase):
    def setUp(self):
        self.a = read_pil_line("domain a = 15")
        self.b = read_pil_line("domain b = 15")
        self.c = read_pil_line("domain c = 15")
        self.d = read_pil_line("domain d = 15")
        self.e = read_pil_line("domain e = 15")
        self.x = read_pil_line("domain x = 15")

    def tearDown(self):
        clear_memory()

    def test_max_helix_initiation_pen(self):
        """ 
        A series of 3-way branch migration reactions.
        """
        # INPUT
        complexes, reactions = read_pil("""
        X  = a( b( c + c( d( + ) ) ) )
        Y  = a( b( c( + c d( + ) ) ) )
        """)
        X = complexes['X']
        Y = complexes['Y']
        output = rxn.branch_3way(X)
        reaction = output[0]
        self.assertEqual(int(reaction.const), 22)

        output = rxn.branch_3way(Y)
        reaction = output[0]
        self.assertEqual(int(reaction.const), 22)

    def test_max_helix_initiation_pen2(self):
        """ 
        A series of 3-way branch migration reactions.
        """
        # INPUT
        complexes, reactions = read_pil("""
        X  = a( b( c e + x c( d( + ) ) ) )
        Y  = a( b( c( e + x c d( + ) ) ) )
        """)

        X = complexes['X']
        Y = complexes['Y']
        output = rxn.branch_3way(X)
        reaction = output[0]
        self.assertEqual(int(reaction.const), 666)

        output = rxn.branch_3way(Y)
        reaction = output[0]
        self.assertEqual(int(reaction.const), 666)
