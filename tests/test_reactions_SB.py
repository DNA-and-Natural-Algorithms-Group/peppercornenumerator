#
#  test_reactions.py
#

import copy
import unittest

import peppercornenumerator.reactions as rxn
from peppercornenumerator.input import from_kernel


class NewBranch3WayTests(unittest.TestCase):
    def setUp(self):
        pass

    def test_single_migration(self):
        """ 
        A single 3-way branch migration reaction.
        """
        # INPUT
        (domains, strands, complexes) = from_kernel([
        "X = a( b x + b( d( + ) ) )",
        "Y = a( b( x + b d( + ) ) )",
        ])
        reactant = complexes['X']
        product = complexes['Y']

        # OUTPUT
        forward = rxn.ReactionPathway('branch_3way', [reactant], [product])
        backward = rxn.ReactionPathway('branch_3way', [product], [reactant])

        output = rxn.branch_3way(reactant, greedy=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, greedy=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        output = rxn.branch_3way(reactant, greedy=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, greedy=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

    def test_greedy_migration(self):
        """ 
        A series of 3-way branch migration reactions.
        """
        # INPUT
        (domains, strands, complexes) = from_kernel([
        "X  = a( x y z + x( y( z( b( + ) ) ) ) )",
        "I1 = a( x( y z + x y( z( b( + ) ) ) ) )",
        "I2 = a( x( y( z + x y z( b( + ) ) ) ) )",
        "Y  = a( x( y( z( + x y z b( + ) ) ) ) )",
        ])
        reactant = complexes['X']
        inter1 = complexes['I1']
        inter2 = complexes['I2']
        product = complexes['Y']

        # ~~~~~~~~~~~~~ #
        # OUTPUT GREEDY #
        # ~~~~~~~~~~~~~ #
        forward = rxn.ReactionPathway('branch_3way', [reactant], [product])
        backward = rxn.ReactionPathway('branch_3way', [product], [reactant])

        output = rxn.branch_3way(reactant, greedy=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, greedy=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        forward = rxn.ReactionPathway('branch_3way', [inter1], [product])
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])

        output = rxn.branch_3way(inter1, greedy=True)
        #for o in output: print 'greedy', o.kernel_string()
        self.assertEqual(output, [backward, forward])
 
        # ~~~~~~~~~~~~~~~~~~~ #
        # OUTPUT NO-MAX-HELIX #
        # ~~~~~~~~~~~~~~~~~~~ #
        forward = rxn.ReactionPathway('branch_3way', [reactant], [inter1])
        backward = rxn.ReactionPathway('branch_3way', [product], [inter2])

        output = rxn.branch_3way(reactant, greedy=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, greedy=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        # OUTPUT NO-MAX-HELIX
        forward = rxn.ReactionPathway('branch_3way', [inter1], [inter2])
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])

        output = rxn.branch_3way(inter1, greedy=False)
        #for o in output: print 'nmheli', o.kernel_string()
        self.assertEqual(sorted(output), sorted([forward, backward]))
 
    def test_remote_migration(self):
        """ 
        A remote 3way branch migration reaction.
        """
        # INPUT
        (domains, strands, complexes) = from_kernel([
        "X  = a( b x y z + x( y( z( c( + ) ) ) ) )",
        "I1 = a( b x( y z + x y( z( c( + ) ) ) ) )",
        "I2 = a( b x( y( z + x y z( c( + ) ) ) ) )",
        "Y  = a( b x( y( z( + x y z c( + ) ) ) ) )",
        ])
        reactant = complexes['X']
        inter1 = complexes['I1']
        inter2 = complexes['I2']
        product = complexes['Y']

        # ~~~~~~~~~~~~~ #
        # OUTPUT REMOTE #
        # ~~~~~~~~~~~~~ #
        forward = rxn.ReactionPathway('branch_3way', [reactant], [product])
        backward = rxn.ReactionPathway('branch_3way', [product], [reactant])

        output = rxn.branch_3way(reactant, greedy=True, remote=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, greedy=True, remote=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        forward = rxn.ReactionPathway('branch_3way', [inter1], [product])
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])

        output = rxn.branch_3way(inter1, greedy=True, remote=True)
        #for o in output: print 'greedy', o.kernel_string()
        self.assertEqual(sorted(output), sorted([forward, backward]))
 
        # ~~~~~~~~~~~~~~~~~~~~ #
        # OUTPUT REJECT REMOTE #
        # ~~~~~~~~~~~~~~~~~~~~ #
        backward = rxn.ReactionPathway('branch_3way', [product], [reactant])

        output = rxn.branch_3way(reactant, greedy=True, remote=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [])

        output = rxn.branch_3way(product, greedy=True, remote=False)
        #for o in output: print 'greedy', o.kernel_string()
        self.assertEqual(output, [backward])

        forward = rxn.ReactionPathway('branch_3way', [inter1], [product])
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])

        output = rxn.branch_3way(inter1, greedy=True, remote=False)
        #for o in output: print 'greedy', o.kernel_string()
        self.assertEqual(sorted(output), sorted([forward, backward]))
 
    def test_multiple_choice(self):
        """ 
        A multiple choice input for 3way branch migration reactions.

        There are two important realizations: 
            1) max-helix move-set *only* extend immediately adjacent helices
            and *not* remote branches.
            2) max-helix move-set *always* extends in both directions.

        Example:

                y/
               x/
              z/
             y/
          a b/x y z c
          _x/________ 
          ___________ 
          a* x*y*z*c*

        """
        # INPUT
        (domains, strands, complexes) = from_kernel([
        "X  = a( x b y z x y + x( y( z( c( + ) ) ) ) )",
        "I1 = a( x( b y z x y + x y( z( c( + ) ) ) ) )",
        "I2 = a( x( b y( z x y + x y z( c( + ) ) ) ) )", # no-max-helix
        "Y1 = a( x( b y( z( x y + x y z c( + ) ) ) ) )",
        "I4 = a( x( b y z x y( + x y z( c( + ) ) ) ) )",
        "I3 = a( x b y z x( y + x y( z( c( + ) ) ) ) )", # no-max-helix
        "Y2 = a( x b y z x( y( + x y z( c( + ) ) ) ) )",
        ])
        reactant = complexes['X']
        inter1 = complexes['I1']
        inter2 = complexes['I2']
        inter3 = complexes['I3']
        inter4 = complexes['I4']
        product1 = complexes['Y1']
        product2 = complexes['Y2']

        # ~~~~~~ #
        # OUTPUT #
        # ~~~~~~ #
        forward = rxn.ReactionPathway('branch_3way', [reactant], [inter1])
        output = rxn.branch_3way(reactant, greedy=True, remote=False)
        self.assertEqual(output, [forward])

        # NOTE: NO GREEDY REMOTE TOEHOLDS!
        forward1b = rxn.ReactionPathway('branch_3way', [reactant], [inter1])
        forward2 = rxn.ReactionPathway('branch_3way', [reactant], [product2])
        output = rxn.branch_3way(reactant, greedy=True, remote=True)
        self.assertEqual(output, [forward2, forward1b])

        backward1 = rxn.ReactionPathway('branch_3way', [product1], [reactant])
        output = rxn.branch_3way(product1, greedy=True, remote=True)
        self.assertEqual(output, [backward1])

        # NOTE: INTERNAL REARRANGEMENT LEADS TO I4!
        backward2 = rxn.ReactionPathway('branch_3way', [product2], [reactant])
        backward2b = rxn.ReactionPathway('branch_3way', [product2], [inter4])
        output = rxn.branch_3way(product2, greedy=True, remote=True)
        #for o in output: print 'ow', o.kernel_string()
        self.assertEqual(output, [backward2, backward2b])

        # NOTE: GREEDY BRANCH MIGRATION MAXIMIZES NEW FORMED HELIX (whattf)!
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])
        forward2 = rxn.ReactionPathway('branch_3way', [inter1], [product2])
        forward  = rxn.ReactionPathway('branch_3way', [inter1], [product1])
        output = rxn.branch_3way(inter1, greedy=True, remote=True)
        #for o in output: print 'ow', o.kernel_string()
        self.assertEqual(output, [backward, forward2, forward])


if __name__ == '__main__':
  unittest.main()

