import unittest

SKIP = True

@unittest.skipIf(SKIP, "skipping tests")
class RateBranch3WayTests(unittest.TestCase):
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
        (displacing, bound, before, after) = reaction.meta
        rate  = rxn.branch_3way_remote_rate(len(displacing), before, after)
        self.assertEqual(int(rate), 22)

        output = rxn.branch_3way(Y)
        reaction = output[0]
        (displacing, bound, before, after) = reaction.meta
        rate = rxn.branch_3way_remote_rate(len(displacing), before, after)
        self.assertEqual(int(rate), 22)


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
        (displacing, bound, before, after) = reaction.meta
        rate = rxn.branch_3way_remote_rate(len(displacing), before, after)
        self.assertEqual(int(rate), 666)

        output = rxn.branch_3way(Y)
        reaction = output[0]
        (displacing, bound, before, after) = reaction.meta
        rate = rxn.branch_3way_remote_rate(len(displacing), before, after)
        self.assertEqual(int(rate), 666)

