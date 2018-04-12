#
#  test_reactions.py
#

import copy
import unittest
import logging
logging.disable(logging.CRITICAL)

from peppercornenumerator import Enumerator
from peppercornenumerator.input import read_kernel
from peppercornenumerator.objects import PepperReaction, clear_memory

import peppercornenumerator.reactions as rxn

SKIP = False

# Helper functions
def filter_bind11((dom1, struct1, loc1), (dom2, struct2, loc2)):
    return struct1 is None and struct2 is None and dom2.can_pair(dom1)

def filter_3way((dom1, struct1, loc1), (dom2, struct2, loc2)):
    return (struct1 is None) and (struct2 is not None) and dom1.can_pair(dom2)

def filter_4way((dom1, struct1, loc1), (dom2, struct2, loc2)):
    return struct1 is not None and struct2 is not None and dom1 == dom2

def triple(reactant, loc):
    return (reactant.get_domain(loc), reactant.get_structure(loc), loc)

@unittest.skipIf(SKIP, "skipping tests")
class FindOnLoop(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_bind11(self):
        complexes, reaction = read_kernel("""
        A = a( b c d e( ) d* c* b* c* )
        """)
        A = complexes['A']

        folA = rxn.find_on_loop(A, (0,1), filter_bind11, direction=1)
        [s,x,t,y] = folA[0]
        self.assertEqual(s, [triple(A, (0,1))]) # b
        self.assertEqual(x, [triple(A, (0,2)),  # c
                             triple(A, (0,3)),  # d
                             triple(A, (0,4)),  # e
                             triple(A, (0,6)),  # d*
                             triple(A, (0,7))]) # c*
        self.assertEqual(t, [triple(A, (0,8))]) # b*
        self.assertEqual(y, [triple(A, (0,9)),  # c*
                             triple(A, (0,10))])# a*

        [s,x,t,y] = rxn.zipper(A, s[0], x, t[0], y, filter_bind11)
        self.assertEqual(s, [triple(A, (0,1)),  # b
                             triple(A, (0,2)),  # c
                             triple(A, (0,3))]) # d
        self.assertEqual(x, [triple(A, (0,4))]) # e
        self.assertEqual(t, [triple(A, (0,8)),  # b*
                             triple(A, (0,7)),  # c*
                             triple(A, (0,6))]) # d*
        self.assertEqual(y, [triple(A, (0,9)),  # c*
                             triple(A, (0,10))])# a*

        folA = rxn.find_on_loop(A, (0,2), filter_bind11, direction=1)
        [s,x,t,y] = folA[0]
        self.assertEqual(s, [triple(A, (0,2))]) # c
        self.assertEqual(x, [triple(A, (0,3)),  # d
                             triple(A, (0,4)),  # e
                             triple(A, (0,6))]) # d*
        self.assertEqual(t, [triple(A, (0,7))]) # c*
        self.assertEqual(y, [triple(A, (0,8)),  # b*
                             triple(A, (0,9)),  # c*
                             triple(A, (0,10)), # a*
                             triple(A, (0,1))]) # b

        [s,x,t,y] = rxn.zipper(A, s[0], x, t[0], y, filter_bind11)
        self.assertEqual(s, [triple(A, (0,1)),  # b
                             triple(A, (0,2)),  # c
                             triple(A, (0,3))]) # d
        self.assertEqual(x, [triple(A, (0,4))]) # e
        self.assertEqual(t, [triple(A, (0,8)),  # b*
                             triple(A, (0,7)),  # c*
                             triple(A, (0,6))]) # d*
        self.assertEqual(y, [triple(A, (0,9)),  # c*
                             triple(A, (0,10))])# a*

        [s,x,t,y] = folA[1]
        self.assertEqual(s, [triple(A, (0,2))]) # c
        self.assertEqual(x, [triple(A, (0,3)),  # d
                             triple(A, (0,4)),  # e
                             triple(A, (0,6)),  # d*
                             triple(A, (0,7)),  # c*
                             triple(A, (0,8))]) # b*
        self.assertEqual(t, [triple(A, (0,9))]) # c*
        self.assertEqual(y, [triple(A, (0,10)), # a*
                             triple(A, (0,1))]) # b

        [s,x,t,y] = rxn.zipper(A, s[0], x, t[0], y, filter_bind11)
        self.assertEqual(s, [triple(A, (0,2))]) # c
        self.assertEqual(x, [triple(A, (0,3)),  # d
                             triple(A, (0,4)),  # e
                             triple(A, (0,6)),  # d*
                             triple(A, (0,7)),  # c*
                             triple(A, (0,8))]) # b*
        self.assertEqual(t, [triple(A, (0,9))]) # c*
        self.assertEqual(y, [triple(A, (0,10)), # a*
                             triple(A, (0,1))]) # b

    def test_bind11_ms(self):
        complexes, reaction = read_kernel("""
        B = a( b c d e( + ) d* + c* b* c* )
        """)
        B = complexes['B']
        folB = rxn.find_on_loop(B, (0,1), filter_bind11, direction=1)

        # Find on Loop results
        [s,x,t,y] = folB[0]
        self.assertEqual(s, [triple(B, (0,1))]) # b
        self.assertEqual(x, [triple(B, (0,2)),  # c
                             triple(B, (0,3)),  # d
                             triple(B, (0,4)),  # e
                             triple(B, (1,1)),  # d*
                             None,              # +
                             triple(B, (2,0))]) # c*
        self.assertEqual(t, [triple(B, (2,1))]) # b*
        self.assertEqual(y, [triple(B, (2,2)),  # c*
                             triple(B, (2,3))]) # a*

        # Zipper results
        s,x,t,y = rxn.zipper(B, s[0], x, t[0], y, filter_bind11)
        self.assertEqual(s, [triple(B, (0,1)),  # b
                             triple(B, (0,2))]) # c
        self.assertEqual(x, [triple(B, (0,3)),  # d
                             triple(B, (0,4)),  # e
                             triple(B, (1,1)),  # d*
                             None])             # +
        self.assertEqual(t, [triple(B, (2,1)),  # b*
                             triple(B, (2,0))]) # c*
        self.assertEqual(y, [triple(B, (2,2)),  # c*
                             triple(B, (2,3))]) # a*

        folB = rxn.find_on_loop(B, (0,2), filter_bind11, direction=1)

        self.assertEqual(len(folB), 2)

        [s,x,t,y] = folB[0]
        self.assertEqual(s, [triple(B, (0,2))]) # c
        self.assertEqual(x, [triple(B, (0,3)),  # d
                             triple(B, (0,4)),  # e
                             triple(B, (1,1)),  # d*
                             None])             # +
        self.assertEqual(t, [triple(B, (2,0))]) # c*
        self.assertEqual(y, [triple(B, (2,1)),  # b*
                             triple(B, (2,2)),  # c*
                             triple(B, (2,3)),  # a*
                             triple(B, (0,1))]) # b

        s,x,t,y = rxn.zipper(B, s[0], x, t[0], y, filter_bind11)
        self.assertEqual(s, [triple(B, (0,1)),  # b
                             triple(B, (0,2))]) # c
        self.assertEqual(x, [triple(B, (0,3)),  # d
                             triple(B, (0,4)),  # e
                             triple(B, (1,1)),  # d*
                             None])             # +
        self.assertEqual(t, [triple(B, (2,1)),  # b*
                             triple(B, (2,0))]) # c*
        self.assertEqual(y, [triple(B, (2,2)),  # c*
                             triple(B, (2,3))]) # a*

        [s,x,t,y] = folB[1]
        self.assertEqual(s, [triple(B, (0,2))]) # c
        self.assertEqual(x, [triple(B, (0,3)),  # d
                             triple(B, (0,4)),  # e
                             triple(B, (1,1)),  # d*
                             None,              # +
                             triple(B, (2,0)),  # c*
                             triple(B, (2,1))]) # b*
        self.assertEqual(t, [triple(B, (2,2))]) # c*
        self.assertEqual(y, [triple(B, (2,3)),  # a*
                             triple(B, (0,1))]) # b

        s,x,t,y = rxn.zipper(B, s[0], x, t[0], y, filter_bind11)
        self.assertEqual(s, [triple(B, (0,2))]) # c
        self.assertEqual(x, [triple(B, (0,3)),  # d
                             triple(B, (0,4)),  # e
                             triple(B, (1,1)),  # d*
                             None,              # +
                             triple(B, (2,0)),  # c*
                             triple(B, (2,1))]) # b*
        self.assertEqual(t, [triple(B, (2,2))]) # c*
        self.assertEqual(y, [triple(B, (2,3)),  # a*
                             triple(B, (0,1))]) # b

    def test_3way_matching(self):
        complexes, reaction = read_kernel("""
        A = a( b( c( x a b c b y c( b( + ) ) + b c ) ) )
        #                  ^(0.6)                  ^(2,2)
        """)
        A = complexes['A']
        folA = rxn.find_on_loop(A, (0,6), filter_3way, direction=1)
        self.assertEqual(len(folA), 1)
        [s,x,t,y] = folA[0]
        self.assertEqual(s, [triple(A, (0,6))]) # c
        self.assertEqual(x, [triple(A, (0,7)),  # b
                             triple(A, (0,8)),  # y
                             triple(A, (0,9)),  # c
                             None,              # +
                             triple(A, (2,0)),  # b
                             triple(A, (2,1))]) # c
        self.assertEqual(t, [triple(A, (2,2))]) # c*
        self.assertEqual(y, [triple(A, (0,3)),  # x
                             triple(A, (0,4)),  # a
                             triple(A, (0,5))]) # b

        [s,x,t,y] = rxn.zipper(A, s[0], x, t[0], y, filter_3way)
        self.assertEqual(s, [triple(A, (0,4)),  # a
                             triple(A, (0,5)),  # b
                             triple(A, (0,6))]) # c

        self.assertEqual(x, [triple(A, (0,7)),  # b
                             triple(A, (0,8)),  # y
                             triple(A, (0,9)),  # c
                             None,              # +
                             triple(A, (2,0)),  # b
                             triple(A, (2,1))]) # c

        self.assertEqual(t, [triple(A, (2,4)),  # a*
                             triple(A, (2,3)),  # b*
                             triple(A, (2,2))]) # c*
        self.assertEqual(y, [triple(A, (0,3))]) # x



        folA = rxn.find_on_loop(A, (0,6), filter_3way, direction=-1)
        self.assertEqual(len(folA), 1)
        [s,x,t,y] = folA[0]
        self.assertEqual(s, [triple(A, (0,6))]) # c
        self.assertEqual(x, [triple(A, (0,5)),  # b
                             triple(A, (0,4)),  # a
                             triple(A, (0,3)),  # x
                             triple(A, (0,2)),  # c
                             triple(A, (2,1)),  # c
                             triple(A, (2,0)),  # b
                             None])             # +
        self.assertEqual(t, [triple(A, (1,1))]) # c*
        self.assertEqual(y, [triple(A, (0,8)),  # y
                             triple(A, (0,7))]) # b

        [s,x,t,y] = rxn.zipper(A, s[0], x, t[0], y, filter_3way)
        self.assertEqual(s, [triple(A, (0,6)),  # c
                             triple(A, (0,7))]) # b
        self.assertEqual(x, [triple(A, (0,5)),  # b
                             triple(A, (0,4)),  # a
                             triple(A, (0,3)),  # x
                             triple(A, (0,2)),  # c
                             triple(A, (2,1)),  # c
                             triple(A, (2,0)),  # b
                             None])             # +

        self.assertEqual(t, [triple(A, (1,1)),  # c*
                             triple(A, (1,0))]) # b*
        self.assertEqual(y, [triple(A, (0,8)),  # y
                             triple(A, (0,7))]) # b

    def test_4way_matching(self):
        complexes, reaction = read_kernel("""
        A = c* b*( X2*( X1*( + ) ) ) c y* b*( a* Y1* + a ) c( + ) b*( X2*( X1*( + ) ) ) c y* b*( a*( Y1* + ) ) c
        #      ^(0,1)                     ^(1,5)                  ^(3,3)
        """)
        A = complexes['A']
        folA = rxn.find_on_loop(A, (0,1), filter_4way)
        self.assertEqual(len(folA), 3)

        [s,x,t,y] = folA[0]
        self.assertEqual(s, [triple(A, (0, 1))])
        self.assertEqual(x, [triple(A, (1, 3)), triple(A, (1, 4))])
        self.assertEqual(t, [triple(A, (1, 5))])
        self.assertEqual(y, [triple(A, (2, 2)), 
                             triple(A, (3, 1)), 
                             triple(A, (4, 3)), 
                             triple(A, (4, 4)), 
                             triple(A, (4, 5)), 
                             triple(A, (5, 2)), 
                             None, 
                             triple(A, (0, 0))])

        [s,x,t,y] = folA[1]
        self.assertEqual(s, [triple(A, (0, 1))] )
        self.assertEqual(t, [triple(A, (3, 1))] )

        [s,x,t,y] = rxn.zipper(A, s[0], x, t[0], y, filter_4way)
        self.assertEqual(s, [triple(A, (0, 1)), 
                             triple(A, (0, 2)), 
                             triple(A, (0, 3))])
        self.assertEqual(x, [triple(A, (1, 3)), 
                             triple(A, (1, 4)), 
                             triple(A, (1, 5)), 
                             triple(A, (2, 2))])
        self.assertEqual(t, [triple(A, (3, 1)), 
                             triple(A, (3, 2)), 
                             triple(A, (3, 3))])
        self.assertEqual(y, [triple(A, (4, 3)), 
                             triple(A, (4, 4)), 
                             triple(A, (4, 5)), 
                             triple(A, (5, 2)), 
                             None, 
                             triple(A, (0, 0))])

@unittest.skipIf(SKIP, "skipping tests")
class NewOpenTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_basic_open(self):
        """ 
        A basic open reaction.

        Testing max-helix-semantics and release-cutoff 5, 8, 13
        """
        # INPUT
        complexes, reactions = read_kernel("""
        length a = 8
        length t = 5

        X = a( t( + ) )
        Y = a( t + t* )
        Z = a t( + ) a*
        S1 = a t 
        S2 = t* a*
        """)
        reactant = complexes['X']
        product1 = complexes['Y']
        product2 = complexes['Z']
        product_set = sorted([complexes['S1'], complexes['S2']])

        # max helix semantics ON -> no reactions vs dissociation
        output = rxn.open(reactant, max_helix=True, release_11=7, release_1N=7)
        self.assertEqual(output, [])
        output = rxn.open(reactant, max_helix=True, release_11=8, release_1N=8)
        self.assertEqual(output, [])

        forward = PepperReaction([reactant], product_set, 'open', memorycheck=False)
        forward.rate = rxn.opening_rate(13)

        forward = PepperReaction([reactant], product_set, 'open', memorycheck=False)
        output = rxn.open(reactant, max_helix=True, release_11=13, release_1N=13)
        #for o in output: print 'ow', o.kernel_string
        self.assertEqual(output, [forward])

        # max helix semantics OFF -> domains dissociate, but at most one at a time
        forward1 = PepperReaction([reactant], [product1], 'open', memorycheck=False)
        forward1.rate = rxn.opening_rate(5)
        forward2 = PepperReaction([reactant], [product2], 'open', memorycheck=False)
        forward2.rate = rxn.opening_rate(8)

        output = rxn.open(reactant, max_helix=False, release_11=7, release_1N=7)
        self.assertEqual(output, [forward1])

        output = rxn.open(reactant, max_helix=False, release_11=8, release_1N=8)
        self.assertEqual(output, sorted([forward2, forward1]))

        output = rxn.open(reactant, max_helix=False, release_11=13, release_1N=13)
        self.assertEqual(output, sorted([forward2, forward1]))

    def test_multiple_choice(self):
        # TODO: len(a) + len(b) = len(ab) - but the behavior is different!
        #
        ## think of hierarchical domain-level displacement.
        ## Upper  class: max-helix
        ## Middle class: no-max-helix
        ## Lower  class: sequence-level
        #
        # Should release-cutoff be dependent on domain-level representation?
        # That means, should release-cutoff of 4 mean "max-helix-semantics up
        # to length 4"? What would that mean if you start in 
        # "a b c( + ) b* a*" vs "ab c( + ) ab*"
        complexes, reactions = read_kernel("""
        length a = 3
        length b = 1
        length c = 3
        length ab = 4

        X = a( b( c( + ) ) )
        Y = ab( c( + ) )

        """)
        pass

@unittest.skipIf(SKIP, "skipping tests")
class NewBindTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_binding(self):
        complexes, reactions = read_kernel("""
        length a = 10
        length t = 5

        S1 = a t 
        S2 = t* a*
        X = a( t + t* )
        Y = a t( + ) a*
        Z = a( t( + ) )

        SB = a a t t* a* a*
        SG1 = a( a t t* ) a* 
        SG2 = a a( t t* a* ) 
        SG3 = a( a( t(  ) ) )

        SI1 = a( a t t* a* )
        SI2 = a a( t t* ) a*
        SI3 = a a t( ) a* a*
        """)
        S1 = complexes['S1']
        S2 = complexes['S2']
        X = complexes['X']
        Y = complexes['Y']
        Z = complexes['Z']
        singles = sorted([S1, S2])
        SB = complexes['SB']
        SG1 = complexes['SG1']
        SG2 = complexes['SG2']
        SG3 = complexes['SG3']
        SI1 = complexes['SI1']
        SI2 = complexes['SI2']
        SI3 = complexes['SI3']

        path = PepperReaction([X], [Z], 'bind11', memorycheck=False)
        output = rxn.bind11(X, max_helix=True)
        #for o in output: print 'mh', o.kernel_string()
        self.assertEqual(output, [path])

        path = PepperReaction([Y], [Z], 'bind11', memorycheck=False)
        output = rxn.bind11(Y, max_helix=True)
        #for o in output: print 'mh', o.kernel_string()
        self.assertEqual(output, [path])

        output = rxn.bind11(Z, max_helix=True)
        #for o in output: print 'mh', o.kernel_string()
        self.assertEqual(output, [])

        path1 = PepperReaction([S2, S1], [X], 'bind21', memorycheck=False)
        path2 = PepperReaction([S2, S1], [Y], 'bind21', memorycheck=False)
        output = rxn.bind21(S1, S2, max_helix=False)
        #for o in output: print 'bind21', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2]))

        path = PepperReaction([S2, S1], [Z], 'bind21', memorycheck=False)
        output = rxn.bind21(S1, S2, max_helix=True)
        #for o in output: print 'bind21g', o.kernel_string()
        self.assertEqual(output, [path])

        path1 = PepperReaction([SB], [SG1], 'bind11', memorycheck=False)
        path2 = PepperReaction([SB], [SG2], 'bind11', memorycheck=False)
        path3 = PepperReaction([SB], [SG3], 'bind11', memorycheck=False)
        output = rxn.bind11(SB, max_helix=True)
        #for o in output: print 'open11g', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3]))

        #path1 = PepperReaction([SB], [SG1], 'bind11', memorycheck=False)
        #path2 = PepperReaction([SB], [SG2], 'bind11', memorycheck=False)
        path3 = PepperReaction([SB], [SI1], 'bind11', memorycheck=False)
        path4 = PepperReaction([SB], [SI2], 'bind11', memorycheck=False)
        path5 = PepperReaction([SB], [SI3], 'bind11', memorycheck=False)
        output = rxn.bind11(SB, max_helix=False)
        #for o in output: print 'open11f', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3, path4, path5]))

    def test_multiple_choice(self):
        #TODO
        pass
 
@unittest.skipIf(SKIP, "skipping tests")
class NewBranch3WayTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_break_3way(self):
        complexes, reactions = read_kernel("""
        length b = 7
        A = b( b*( b b*( + ) b* ) )
        """)

        A1 = complexes['A']
        output = rxn.branch_3way(A1, max_helix=False, remote=True)
        #for o in output: print 'branch_3way', o.kernel_string, o.rate

    def test_single_migration(self):
        """ 
        A single 3-way branch migration reaction.
        """
        # INPUT
        complexes, reactions = read_kernel("""
        X = a( b x + b( d( + ) ) )
        Y = a( b( x + b d( + ) ) )
        """)
        reactant = complexes['X']
        product = complexes['Y']

        # OUTPUT
        forward  = PepperReaction([reactant], [product], 'branch-3way', memorycheck=False)
        backward = PepperReaction([product], [reactant], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(reactant, max_helix=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, max_helix=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        output = rxn.branch_3way(reactant, max_helix=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, max_helix=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

    def test_max_helix_migration(self):
        """ 
        A series of 3-way branch migration reactions.
        """
        # INPUT
        complexes, reactions = read_kernel("""
        X  = a( x y z + x( y( z( b( + ) ) ) ) )
        I1 = a( x( y z + x y( z( b( + ) ) ) ) )
        I2 = a( x( y( z + x y z( b( + ) ) ) ) )
        Y  = a( x( y( z( + x y z b( + ) ) ) ) )
        """)
        reactant = complexes['X']
        inter1 = complexes['I1']
        inter2 = complexes['I2']
        product = complexes['Y']

        # ~~~~~~~~~~~~~~~~ #
        # OUTPUT max-helix #
        # ~~~~~~~~~~~~~~~~ #
        forward =  PepperReaction([reactant], [product], 'branch-3way', memorycheck=False)
        backward = PepperReaction([product], [reactant], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(reactant, max_helix=True)
        #print output[0].kernel_string
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, max_helix=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        forward = PepperReaction([inter1], [product], 'branch-3way', memorycheck=False)
        backward = PepperReaction([inter1], [reactant], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(inter1, max_helix=True)
        #for o in output: print 'max_helix', o.kernel_string()
        self.assertEqual(output, sorted([backward, forward]))
 
        # ~~~~~~~~~~~~~~~~~~~ #
        # OUTPUT NO-MAX-HELIX #
        # ~~~~~~~~~~~~~~~~~~~ #
        forward = PepperReaction([reactant], [inter1], 'branch-3way', memorycheck=False)
        backward = PepperReaction([product], [inter2], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(reactant, max_helix=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, max_helix=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        # OUTPUT NO-MAX-HELIX
        forward = PepperReaction([inter1], [inter2], 'branch-3way', memorycheck=False)
        backward = PepperReaction([inter1], [reactant], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(inter1, max_helix=False)
        #for o in output: print 'nmheli', o.kernel_string()
        self.assertEqual(sorted(output), sorted([forward, backward]))
 
    def test_remote_migration(self):
        """ 
        A remote 3way branch migration reaction.
        """
        # INPUT
        complexes, reactions = read_kernel("""
        X  = a( b x y z + x( y( z( c( + ) ) ) ) )
        I1 = a( b x( y z + x y( z( c( + ) ) ) ) )
        I2 = a( b x( y( z + x y z( c( + ) ) ) ) )
        Y  = a( b x( y( z( + x y z c( + ) ) ) ) )
        """)
        reactant = complexes['X']
        inter1 = complexes['I1']
        inter2 = complexes['I2']
        product = complexes['Y']

        # ~~~~~~~~~~~~~ #
        # OUTPUT REMOTE #
        # ~~~~~~~~~~~~~ #
        forward = PepperReaction([reactant], [product], 'branch-3way', memorycheck=False)
        backward = PepperReaction([product], [reactant], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(reactant, max_helix=True, remote=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, max_helix=True, remote=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        forward = PepperReaction([inter1], [product], 'branch-3way', memorycheck=False)
        backward = PepperReaction([inter1], [reactant], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(inter1, max_helix=True, remote=True)
        #for o in output: print 'max_helix', o.kernel_string()
        self.assertEqual(sorted(output), sorted([backward, forward]))
 
        # ~~~~~~~~~~~~~~~~~~~~ #
        # OUTPUT REJECT REMOTE #
        # ~~~~~~~~~~~~~~~~~~~~ #
        backward = PepperReaction([product], [reactant], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(reactant, max_helix=True, remote=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [])

        output = rxn.branch_3way(product, max_helix=True, remote=False)
        #for o in output: print 'max_helix', o.kernel_string()
        self.assertEqual(output, [backward])

        forward = PepperReaction([inter1], [product], 'branch-3way', memorycheck=False)
        backward = PepperReaction([inter1], [reactant], 'branch-3way', memorycheck=False)

        output = rxn.branch_3way(inter1, max_helix=True, remote=False)
        #for o in output: print 'max_helix', o.kernel_string()
        self.assertEqual(sorted(map(str,output)), sorted(map(str,[forward, backward])))
 
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
        complexes, reactions = read_kernel("""
        X  = a( x b y z x y + x( y( z( c( + ) ) ) ) )
        I1 = a( x( b y z x y + x y( z( c( + ) ) ) ) )
        I2 = a( x( b y( z x y + x y z( c( + ) ) ) ) ) # no-max-helix
        Y1 = a( x( b y( z( x y + x y z c( + ) ) ) ) )
        Y2 = a( x b y z x( y( + x y z( c( + ) ) ) ) )
        Y3 = a( x( b y z x y( + x y z( c( + ) ) ) ) )
        Y4 = a( x b y z x( y + x y( z( c( + ) ) ) ) )
        """)
        reactant = complexes['X']
        inter1 = complexes['I1']
        inter2 = complexes['I2']
        product1 = complexes['Y1']
        product2 = complexes['Y2']
        product3 = complexes['Y3']
        product4 = complexes['Y4']

        # ~~~~~~ #
        # OUTPUT #
        # ~~~~~~ #
        forward = PepperReaction([reactant], [inter1], 'branch-3way', memorycheck=False)
        output = rxn.branch_3way(reactant, max_helix=True, remote=False)
        self.assertEqual(output, [forward])

        # NOTE: NO max_helix REMOTE TOEHOLDS!
        forward1b = PepperReaction([reactant], [inter1], 'branch-3way', memorycheck=False)
        forward2 = PepperReaction([reactant], [product2], 'branch-3way', memorycheck=False)
        output = rxn.branch_3way(reactant, max_helix=True, remote=True)
        self.assertEqual(sorted(map(str,output)), sorted(map(str,[forward2, forward1b])))

        # NOTE: THIS behavior changed now!!!
        backward1 = PepperReaction([product1], [inter1], 'branch-3way', memorycheck=False)
        output = rxn.branch_3way(product1, max_helix=True, remote=True)
        #for o in output: print 'ow', o, o.kernel_string
        self.assertEqual(output, [backward1])

        # NOTE: THIS behavior changed!
        backward2 = PepperReaction([product2], [reactant], 'branch-3way', memorycheck=False)
        #backward2b = PepperReaction([product2], [product4], 'branch-3way', memorycheck=False)
        backward2b = PepperReaction([product2], [product3], 'branch-3way', memorycheck=False)
        output = rxn.branch_3way(product2, max_helix=True, remote=True)
        #for o in output: print 'ow', o.kernel_string()
        self.assertEqual(output, sorted([backward2, backward2b]))

        # NOTE: max_helix zippering cannot involve different strands than the initial step
        backward = PepperReaction([inter1], [reactant], 'branch-3way', memorycheck=False) #1
        forward  = PepperReaction([inter1], [product1], 'branch-3way', memorycheck=False) #4
        forward2 = PepperReaction([inter1], [product3], 'branch-3way', memorycheck=False)
        forward3 = PepperReaction([inter1], [product4], 'branch-3way', memorycheck=False)
        output = rxn.branch_3way(inter1, max_helix=True, remote=True)
        #for o in output: print 'fixme', o.kernel_string()
        self.assertEqual(output, sorted([backward, forward, forward2, forward3]))

    def test_multiple_choice_2(self):
        # INPUT
        complexes, reactions = read_kernel("""
        N = a( x( b x y + y( z( + ) ) ) )
        P1 = a( x b x( y + y( z( + ) ) ) )
        P2 = a( x( b x y( + y z( + ) ) ) )
        """)
        N = complexes['N']
        P1 = complexes['P1']
        P2 = complexes['P2']

        forward1 = PepperReaction([N], [P1], 'branch-3way', memorycheck=False)
        forward2 = PepperReaction([N], [P2], 'branch-3way', memorycheck=False)
        output = rxn.branch_3way(N, max_helix=True, remote=True)
        #for o in output: print 'new-max', o.kernel_string
        self.assertEqual(sorted(output), sorted([forward1, forward2]))


@unittest.skipIf(SKIP, "skipping tests")
class NewBranch4WayTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_break_casey_4way(self):
        # Note the new max-helix semantics also fixes the old 4way error,
        # so this test can be safely removed...
        complexes, reactions = read_kernel("""
        A1 = t0*( d3*( d4*( + ) ) + d3*( t0* d3*( d4*( + ) ) + ) )
        """)

        A1 = complexes['A1']
        output = rxn.branch_4way(A1, max_helix=True, remote=True)
        #for o in output: print 'branch_4way_bug', o.kernel_string()

    def test_break_reverse_4way(self):
        complexes, reactions = read_kernel("""
        B = y*( b*( a*( Y1* + ) b( c( + ) ) Y2* Y1*( + ) a ) )
        A = c* b* X2* X1*( + ) X2 b c y*( b*( a*( Y1* + ) b( c( + ) ) Y2* Y1*( + ) a ) )
        R = c* b* X2* X1*( + ) X2 b c y*( b*( a*( Y1* + ) ) c( + ) b*( Y2* Y1*( + ) a ) )
        """)

        A1 = complexes['A']
        R1 = complexes['R']
        output = rxn.branch_4way(A1, max_helix=True, remote=True)
        oR = rxn.branch_4way(R1, max_helix=True, remote=True)
        #for o in output: print 'branch_4way_reverse', o.kernel_string, o.rate
        #for o in oR: print 'branch_4way_reverse', o.kernel_string, o.rate


    def test_4wayfilter_bugfix(self):
        # a test to ensure the 4way filter includes struct1
        complexes, reactions = read_kernel("""
        # Domain Specifications
        length d1 = 15
        length d2 = 15
        length d3 = 15
        length d4 = 15
        length d5 = 15
        length d6 = 15
        length t0 = 6

        BUG = d4*( t0* + ) t0( + ) d4*( t0*( + ) ) 
        FIX = d4*( t0* + d4( t0( + ) ) t0*( + ) )
        """)
        BUG = complexes['BUG']
        FIX = complexes['FIX']

        path = PepperReaction([BUG], [FIX], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(BUG, max_helix=True, remote=True)
        #for o in output: print 'branch_4way_bug', o.kernel_string()
        self.assertEqual(output, [path])

    def test_branch4_way(self):
        # Standard 3state 4way junction, no end-dissociation
        complexes, reactions = read_kernel("""
        A1 = a( b( c( + ) ) x*( + ) b( c( y( + ) ) ) )
        A2 = a( b( c( + ) b*( x*( + ) ) c( y( + ) ) ) )
        A3 = a( b( c( + c*( b*( x*( + ) ) ) y( + ) ) ) )
        """)
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']

        path = PepperReaction([A1], [A2], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A1, max_helix=False)
        #for o in output: print 'branch_4way one-step', o.kernel_string, o.rate
        self.assertEqual(output, [path])

        path = PepperReaction([A3], [A2], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A3, max_helix=False)
        #for o in output: print 'branch_4way one-step', o.kernel_string(), o.rate
        self.assertEqual(output, [path])

        path1 = PepperReaction([A2], [A1], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A2], [A3], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A2, max_helix=False)
        #for o in output: print 'branch_4way one-step', o.kernel_string(), o.rate
        self.assertEqual(output, sorted([path1, path2]))

        path = PepperReaction([A1], [A3], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A1, max_helix=True)
        #for o in output: print 'branch_4way_two-step', o.kernel_string(), o.rate
        self.assertEqual(output, [path])

        path = PepperReaction([A3], [A1], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A3, max_helix=True)
        #for o in output: print 'branch_4way_two-step', o.kernel_string(), o.rate
        self.assertEqual(output, [path])

    def test_branch4_way_2(self):
        # Unconventional multi-state 4way junction, no end-dissociation
        complexes, reactions = read_kernel("""
        A1 = a( a*( b( + ) ) a*( x*( + ) ) b( c( + ) ) )
        A2 = a( ) b( + ) a( a*( x*( + ) ) b( c( + ) ) )
        A3 = a( a*( b( + ) a( ) x*( + ) ) b( c( + ) ) )
        A4 = a( a*( b( + ) ) ) x*( + ) a( b( c( + ) ) )
        """)
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']
        A4 = complexes['A4']

        path1 = PepperReaction([A1], [A2], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A1], [A3], 'branch-4way', memorycheck=False)
        path3 = PepperReaction([A1], [A4], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A1, max_helix=False)
        #for o in output: print 'branch_4way_next', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3]))

        path1 = PepperReaction([A1], [A2], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A1], [A3], 'branch-4way', memorycheck=False)
        path3 = PepperReaction([A1], [A4], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A1, max_helix=True)
        #for o in output: print 'branch_4way_next_mh', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3]))

    def test_branch4_way_long(self):
        # Unconventional multi-state 4way junction, no end-dissociation
        complexes, reactions = read_kernel("""
        A1 = a( x*( y*( z*( b z( y( x( c( + ) ) ) ) d z( y( x( e( + ) ) ) ) f z( y( x( g( + ) ) ) ) h ) ) ) )

        A2 = a( x*( y*( z*( b ) y( x( c( + ) ) ) z*( d z( y( x( e( + ) ) ) ) f z( y( x( g( + ) ) ) ) h ) ) ) )
        A3 = a( x*( y*( z*( b z( y( x( c( + ) ) ) ) d ) y( x( e( + ) ) ) z*( f z( y( x( g( + ) ) ) ) h ) ) ) )
        A4 = a( x*( y*( z*( b z( y( x( c( + ) ) ) ) d z( y( x( e( + ) ) ) ) f ) y( x( g( + ) ) ) z*( h ) ) ) )
        A5 = a( x*( y*( z*( b z( y( x( c( + ) ) ) ) d z( y( x( e( + ) ) ) z*( f ) y( x( g( + ) ) ) ) h ) ) ) )
        A6 = a( x*( y*( z*( b z( y( x( c( + ) ) ) z*( d ) y( x( e( + ) ) ) ) f z( y( x( g( + ) ) ) ) h ) ) ) )
        A7 = a( x*( y*( z*( b z( y( x( c( + ) ) ) z*( d z( y( x( e( + ) ) ) ) f ) y( x( g( + ) ) ) ) h ) ) ) )

        A2m = a( x*( y*( z*( b ) ) ) c( + ) x*( y*( z*( d z( y( x( e( + ) ) ) ) f z( y( x( g( + ) ) ) ) h ) ) ) )
        A3m = a( x*( y*( z*( b z( y( x( c( + ) ) ) ) d ) ) ) e( + ) x*( y*( z*( f z( y( x( g( + ) ) ) ) h ) ) ) )
        A4m = a( x*( y*( z*( b z( y( x( c( + ) ) ) ) d z( y( x( e( + ) ) ) ) f ) ) ) g( + ) x*( y*( z*( h ) ) ) )
        A5m = a( x*( y*( z*( b z( y( x( c( + ) ) ) ) d z( y( x( e( + ) x*( y*( z*( f ) ) ) g( + ) ) ) ) h ) ) ) )
        A6m = a( x*( y*( z*( b z( y( x( c( + ) x*( y*( z*( d ) ) ) e( + ) ) ) ) f z( y( x( g( + ) ) ) ) h ) ) ) )
        A7m = a( x*( y*( z*( b z( y( x( c( + ) x*( y*( z*( d z( y( x( e( + ) ) ) ) f ) ) ) g( + ) ) ) ) h ) ) ) )
        """)
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']
        A4 = complexes['A4']
        A5 = complexes['A5']
        A6 = complexes['A6']
        A7 = complexes['A7']
        A2m = complexes['A2m']
        A3m = complexes['A3m']
        A4m = complexes['A4m']
        A5m = complexes['A5m']
        A6m = complexes['A6m']
        A7m = complexes['A7m']

        path1 = PepperReaction([A1], [A2], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A1], [A3], 'branch-4way', memorycheck=False)
        path3 = PepperReaction([A1], [A4], 'branch-4way', memorycheck=False)
        path4 = PepperReaction([A1], [A5], 'branch-4way', memorycheck=False)
        path5 = PepperReaction([A1], [A6], 'branch-4way', memorycheck=False)
        path6 = PepperReaction([A1], [A7], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A1, max_helix=False)
        #for o in output: print 'b_4way', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3, path4, path5, path6]))

        path1 = PepperReaction([A1], [A2m], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A1], [A3m], 'branch-4way', memorycheck=False)
        path3 = PepperReaction([A1], [A4m], 'branch-4way', memorycheck=False)
        path4 = PepperReaction([A1], [A5m], 'branch-4way', memorycheck=False)
        path5 = PepperReaction([A1], [A6m], 'branch-4way', memorycheck=False)
        path6 = PepperReaction([A1], [A7m], 'branch-4way', memorycheck=False)

        output = rxn.branch_4way(A1, max_helix=True)
        #for o in output: print 'b_4way_mh', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3, path4, path5, path6]))

    def test_branch4_way_no_remote(self):
        # Standard 4state 4way junction, no end-dissociation
        complexes, reactions = read_kernel("""
        A0 = a( x( y( z( b(  + ) ) ) ) c*( + ) x( y( z( d( + ) ) ) ) )
        A1 = a( x( y( z( b(  + ) ) ) x*( c*( + ) ) y( z( d( + ) ) ) ) )
        A2 = a( x( y( z( b(  + ) ) y*( x*( c*( + ) ) ) z( d( + ) ) ) ) )
        A3 = a( x( y( z( b(  + ) z*( y*( x*( c*( + ) ) ) ) d( + ) ) ) ) )
        """)
        A0 = complexes['A0']
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']

        path1 = PepperReaction([A1], [A3], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A1], [A0], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A1, max_helix=True, remote = False)
        #for o in output: print 'branch-4way', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2]))
        path1 = PepperReaction([A2], [A0], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A2], [A3], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A2, max_helix=True, remote = False)
        #for o in output: print 'branch-4way', o.kernel_string(), o.rate
        self.assertEqual(output, sorted([path1, path2]))

        path = PepperReaction([A3], [A0], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A3, max_helix=True, remote = False)
        #for o in output: print 'branch-4way', o, o.kernel_string, o.rate
        self.assertEqual(output, [path])

        path1 = PepperReaction([A2], [A0], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A2], [A3], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A2, max_helix=True, remote = False)
        #for o in output: print 'branch-4way', o.kernel_string(), o.rate
        self.assertEqual(output, sorted([path1, path2]))

        path1 = PepperReaction([A1], [A3], 'branch-4way', memorycheck=False)
        path2 = PepperReaction([A1], [A0], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A1, max_helix=True, remote = True)
        #for o in output: print 'branch-4way', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2]))

        path = PepperReaction([A3], [A0], 'branch-4way', memorycheck=False)
        output = rxn.branch_4way(A3, max_helix=True, remote = False)
        #for o in output: print 'branch-4way', o.kernel_string(), o.rate
        self.assertEqual(output, [path])


@unittest.skipIf(SKIP, "skipping tests")
class DSD_PathwayTests(unittest.TestCase):
    def setUp(self):
        self.fast_reactions = [rxn.bind11, 
                    rxn.open, 
                    rxn.branch_3way, 
                    rxn.branch_4way]
        self.slow_reactions = {
                1: [],
                2: [rxn.bind21]
                }

        # Default options
        self._release_11 = 6
        self._release_1N = 6
        self._max_helix = True
        self._remote = True
        self._k_fast = 0.0
        self._k_slow = 0.0

    def tearDown(self):
        clear_memory()

    def test_bind_and_displace3way(self):
        # Skip the outer loop of the enumerator...
        complexes, reactions = read_kernel("""
        length a = 10
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
        path1 = PepperReaction(sorted([I, C]), [B], 'bind21', memorycheck=False)
        output = rxn.bind21(I, C, max_helix=True)
        self.assertEqual(output, [path1])

        # DSD-pathway "branch3way"
        path2 = PepperReaction([B], sorted([D, J]), 'branch-3way', memorycheck=False)
        output = rxn.branch_3way(B, max_helix=True)
        self.assertEqual(output, [path2])

        #enum = Enumerator(domains.values(), strands.values(), complexes.values())
        enum = Enumerator(complexes.values())
        enum.enumerate()
        self.assertEqual(sorted(enum.reactions), sorted([path1, path2]))

    def test_cooperative_binding(self):
        complexes, reactions = read_kernel("""
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
        
        self.k_fast = float('inf')
        self.k_slow = 0

        path1 = PepperReaction(sorted([L, C]), [LC], 'bind21', memorycheck=False)
        path1r= PepperReaction([LC], sorted([L, C]), 'open', memorycheck=False)
        path2 = PepperReaction([LC], [LCF], 'branch-3way', memorycheck=False)
        path3 = PepperReaction(sorted([R, LCF]), [LCRF1], 'bind21', memorycheck=False)
        path4 = PepperReaction([LCRF1], sorted([LR, T]), 'branch-3way', memorycheck=False)

        enum = Enumerator(complexes.values())
        enum.k_fast = self.k_fast
        enum.k_slow = self.k_slow
        enum.max_helix_migration = True
        enum.enumerate()
        #for r in enum.reactions:
        #    if r == path1:
        #        print 'bind21: L+C->LC', r.rate
        #    elif r == path1r:
        #        print 'open: LC->L+C', r.rate
        #    elif r == path2:
        #        print 'branch_3way: LC->LCF', r.rate
        #    elif r == path3:
        #        print 'bind21: R + LCF->LCRF1', r.rate
        #    elif r == path4:
        #        print 'brnach_3way: LCRF1 -> LR', r.rate
        #    else:
        #        print r.kernel_string(), r.rate
        self.assertEqual(len(enum.reactions), 22)


@unittest.skipIf(SKIP, "skipping tests")
class ReactionMatchingTests (unittest.TestCase):
    # A selection of Casey's reaction tests.
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def testBind11A(self):

        # bind11: a ? a* ? -> a( ? ) ?
        complexes, _ = read_kernel("""
        A1 = x( ) a  x( ) a* x( )
        A2 = x( ) a( x( ) )  x( )

        A3 = x( ) a  b  x( ) b* a* x( )
        A4 = x( ) a( b( x( ) )  )  x( )
        """)

        # No zipping possible
        rxns = rxn.bind11(complexes['A1'])
        self.assertEqual(rxns,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'bind11', memorycheck=False)])

        # Zipping possible
        rxns = rxn.bind11(complexes['A3'])
        self.assertEqual(rxns,
        [PepperReaction([complexes['A3']], [complexes['A4']], 'bind11', memorycheck=False)])

    def testBind21A(self):
        # bind11: a ? a* ? -> a( ? ) ?
        complexes, _ = read_kernel("""
        A1 = w( ) a  x( )
        A2 = y( ) a* z( )
        A3 = w( ) a( x( ) + y( ) ) z( )
        """)

        # No zipping possible
        rxns = rxn.bind21(complexes['A1'], complexes['A2'])
        self.assertEqual(rxns,
            [PepperReaction([complexes['A1'], complexes['A2']], 
            [complexes['A3']], 'bind21', memorycheck=False)])

    def testOpenA(self):
        # open:  ? a( ? ) ? -> ? a ? a* ?
        complexes, _ = read_kernel("""
        length a = 5
        length b = 5

        A1 = x( ) a( x )  x( )
        A2 = x( ) a  x a* x( )

        A3 = x a( b( y )  )  z
        A4 = x a  b  y b* a* z
        """)

        # No zipping possible
        self.assertEqual(rxn.open(complexes['A1']),
            [PepperReaction([complexes['A1']], [complexes['A2']], 'open', memorycheck=False)])

        # Zipping possible
        self.assertEqual(rxn.open(complexes['A3'], release_11 = 13, release_1N=13),
            [PepperReaction([complexes['A3']], [complexes['A4']], 'open', memorycheck=False)])

    def testOpenB(self):
        # open:  ? a( ? ) ? -> ? a ? a* ?
        complexes, _ = read_kernel("""
        length a = 5
        length b = 5

        A1 = x a( y )   z
        A2 = x a  y a* z

        A3 = x a( b( y )  )  z
        A4 = x a  b( y )  a* z
        A5 = x a( b  y b* )  z
        """)

        # No zipping possible
        rxns = rxn.open(complexes['A1'], max_helix=False, release_11 = 7, release_1N=7)
        self.assertEqual(rxns,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'open', memorycheck=False)])

        # Zipping possible
        rxns = rxn.open(complexes['A3'], max_helix=False, release_11 = 7, release_1N=7)
        self.assertEqual(rxns, sorted([
            PepperReaction([complexes['A3']], [complexes['A4']], 'open', memorycheck=False),
            PepperReaction([complexes['A3']], [complexes['A5']], 'open', memorycheck=False)]))

    def testOpenNoMaxHelix(self):
        # open:  ? a( ? ) ? -> ? a ? a* ?
        complexes, _ = read_kernel("""
        length a = 5
        length b = 5

        # A1 = a( b( ) )
        A1 = a( b( ) )
        A2 = a( b b* )
        A3 = a b( ) a*
        """)

        # Zipping possible
        rxns = rxn.open(complexes['A1'], max_helix=False, release_11 = 10, release_1N=10)
        self.assertEqual(rxns, sorted([
            PepperReaction([complexes['A1']], [complexes['A2']], 'open', memorycheck=False),
            PepperReaction([complexes['A1']], [complexes['A3']], 'open', memorycheck=False)]))

    def testBranch3wayA(self):
        # 3wayA: ? b ? b(?) ? <-> ? b(? b ?) ?
        complexes, _ = read_kernel("""
        A1 = d1( ) b  d2( ) b( d3( ) ) d4( )
        A2 = d1( ) b( d2( ) b  d3( ) ) d4( )
        """)
        forward = rxn.branch_3way(complexes['A1'])

        self.assertEqual(forward,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'branch-3way', memorycheck=False)])

        reverse = rxn.branch_3way(complexes['A2'])
        self.assertEqual(reverse,
        [PepperReaction([complexes['A2']], [complexes['A1']], 'branch-3way', memorycheck=False)])

    def testBranch3wayB(self):
        # 3wayB: ? b(?) ? b ? <-> ? b ? b*(?) ?
        complexes, _ = read_kernel("""
        A1 = d1( ) b( d2( )  )  d3( ) b d4( )
        A2 = d1( ) b  d2( ) b*( d3( ) ) d4( )
        """)
        forward = rxn.branch_3way(complexes['A1'])
        self.assertEqual(forward,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'branch-3way', memorycheck=False)])

        reverse = rxn.branch_3way(complexes['A2'])
        self.assertEqual(reverse,
        [PepperReaction([complexes['A2']], [complexes['A1']], 'branch-3way', memorycheck=False)])

    def testBranch3wayC(self):
        # 3wayC: ? b*(?) ? b ? <-> ? b*(? b ?) ?
        complexes, _ = read_kernel("""
        A1 = d1( ) b*( d2( ) ) d3( ) b d4( )
        A2 = d1( ) b*( d2( ) b d3( ) ) d4( )
        """)
        forward = rxn.branch_3way(complexes['A1'])
        self.assertEqual(forward,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'branch-3way', memorycheck=False)])

        reverse = rxn.branch_3way(complexes['A2'])
        self.assertEqual(reverse,
        [PepperReaction([complexes['A2']], [complexes['A1']], 'branch-3way', memorycheck=False)])

    def testBranch4wayA(self):
        # 4way: b( ? ) ? b ( ? ) --> b( ? b*( ? ) ? )
        complexes, _ = read_kernel("""
        A1 = d1( ) b( d2( )  )  d3( ) b( d4( ) ) d5( )
        A2 = d1( ) b( d2( ) b*( d3( ) )  d4( ) ) d5( )
        """)

        forward = rxn.branch_4way(complexes['A1'])
        self.assertEqual(forward,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'branch-4way', memorycheck=False)])

        reverse = rxn.branch_4way(complexes['A2'])
        self.assertEqual(reverse,
        [PepperReaction([complexes['A2']], [complexes['A1']], 'branch-4way', memorycheck=False)])

@unittest.skipIf(SKIP, "skipping tests")
class IsomorphicSets(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_simple(self):
        # works just fine
        complexes, reactions = read_kernel("""
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

        self.k_fast = 0
        self.k_slow = 0


        enum = Enumerator([complexes['I'], complexes['C'], complexes['J'], complexes['D']])
        enum.k_fast = self.k_fast
        enum.k_slow = self.k_slow
        enum.max_helix_migration = True
        enum.enumerate()

        enum2 = Enumerator([complexes['cI'], complexes['cC'], complexes['cJ'], complexes['cD']])
        enum2.k_fast = self.k_fast
        enum2.k_slow = self.k_slow
        enum2.max_helix_migration = True
        enum2.enumerate()

        self.assertEqual(len(enum2.reactions), len(enum.reactions))
        #for r in enum.reactions:
        #    print 'test_isomorph', r.kernel_string(), r.rate

    def test_old_vs_new_max_helix(self):
        pass

    def test_erik_max_helix_examples_3way(self):
        complexes, reactions = read_kernel("""

        # should be one reaction, is one
        A1 = x( y z + y( z( + ) ) )
        A1_2 = x( y( z( + ) ) )
        YZ = y z

        # should be one reactions, is one
        B1 = x1( x2( y1 y2 z1 z2 + y1( y2( z1( z2( + ) ) ) ) ) ) 
        B1_2 = x1( x2( y1( y2( z1( z2( + ) ) ) ) ) ) 
        YZ2 = y1 y2 z1 z2

        # should be two reactions, is one
        A2 = x( y z + y( + z( + ) ) )
        A2_1 = x( y( z + z( + ) ) )
        #A2_2 = x( y( z( + ) ) )
        Y1 = y
        Z1 = z

        # should be two reactions, is one
        B2 = x1( x2( y1 y2 z1 z2 + y1( y2( + z1( z2( + ) ) ) ) ) ) 
        B2_1 = x1( x2( y1( y2( z1 z2 + z1( z2( + ) ) ) ) ) ) 
        #B2_2 = x1( x2( y1( y2( z1( z2( + ) ) ) ) ) ) 
        Y2 = y1 y2
        Z2 = z1 z2

        # should be two reactions, is one
        C = x( y z + y( + a( + ) z( + ) ) )

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
        enum.max_helix_migration = True
        enum.enumerate()

        path1 = PepperReaction([A1], sorted([A1_2, YZ]), 'branch-3way', memorycheck=False)
        path2 = PepperReaction([A2], sorted([A2_1, Y1]), 'branch-3way', memorycheck=False)
        path3 = PepperReaction([A2_1], sorted([A1_2, Z1]), 'branch-3way', memorycheck=False)

        self.assertEqual(sorted(enum.reactions), sorted([path1, path2, path3]))
        #for r in enum.reactions:
        #    print 'invade', r, r.kernel_string(), r.rate
 
        B1 = complexes['B1']
        B1_2 = complexes['B1_2']
        YZ2 = complexes['YZ2']

        B2 = complexes['B2']
        B2_1 = complexes['B2_1']
        #B2_2 = complexes['B2_2']
        Y2 = complexes['Y2']
        Z2 = complexes['Z2']

        enum = Enumerator([B1, B2])
        enum.k_fast = 0
        enum.k_slow = 0
        enum.max_helix_migration = True
        enum.enumerate()

        path1 = PepperReaction([B1], sorted([B1_2, YZ2]), 'branch-3way', memorycheck=False)
        path2 = PepperReaction([B2], sorted([B2_1, Y2]), 'branch-3way', memorycheck=False)
        path3 = PepperReaction([B2_1], sorted([B1_2, Z2]), 'branch-3way', memorycheck=False)
        #for r in enum.reactions:
        #    print 'invade', r, r.kernel_string(), r.rate
        self.assertEqual(sorted(enum.reactions), sorted([path1, path2, path3]))

@unittest.skipIf(SKIP, "skipping tests")
class Compare_MaxHelix(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_self_displacement_bug(self):
        complexes, reactions = read_kernel("""
        B1 = x( y( x( y x + ) ) )
        B2 = x y x( y( x( + ) ) )

        B3 = x( y( x y x( + ) ) )
        B4 = x( y x y( x( + ) ) )
        """)
        B1 = complexes['B1']
        B2 = complexes['B2']
        B3 = complexes['B3']
        B4 = complexes['B4']

        forward = PepperReaction([B1], [B2], 'branch-3way', memorycheck=False)
        backward= PepperReaction([B2], [B1], 'branch-3way', memorycheck=False)

        path1  = PepperReaction([B1], [B4], 'branch-3way', memorycheck=False)
        path1r = PepperReaction([B4], [B1], 'branch-3way', memorycheck=False)
        path2  = PepperReaction([B4], [B2], 'branch-3way', memorycheck=False)

        path3  = PepperReaction([B2], [B3], 'branch-3way', memorycheck=False)
        path3r = PepperReaction([B3], [B2], 'branch-3way', memorycheck=False)
        path4  = PepperReaction([B3], [B1], 'branch-3way', memorycheck=False)

        enum = Enumerator([B1])
        enum.max_helix_migration = True
        enum.enumerate()
        #for r in sorted(enum.reactions): print 'max-helix', r, r.kernel_string(), r.rate
        self.assertEqual(sorted([path1, path1r, path2, path3, path3r, path4]), sorted(enum.reactions))
        #self.assertEqual(sorted([forward, backward]), sorted(enum.reactions))

    def test_self_displacement_bug_iso(self):
        complexes, reactions = read_kernel("""
        B1 = x1( x2( y1( y2( x1( x2( y1 y2 x1 x2 + ) ) ) ) ) )
        B2 = x1 x2 y1 y2 x1( x2( y1( y2( x1( x2( + ) ) ) ) ) )
        i1 = x1( x2( y1 y2 x1 x2 y1( y2( x1( x2( + ) ) ) ) ) )
        i2 = x1( x2( y1( y2( x1 x2 y1 y2 x1( x2( + ) ) ) ) ) )

        """)
        B1 = complexes['B1']
        B2 = complexes['B2']
        i1 = complexes['i1']
        i2 = complexes['i2']

        forward = PepperReaction([B1], [B2], 'branch-3way', memorycheck=False)
        backward= PepperReaction([B2], [B1], 'branch-3way', memorycheck=False)

        path1 = PepperReaction([B1], [i1], 'branch-3way', memorycheck=False)
        path1r = PepperReaction([i1], [B1], 'branch-3way', memorycheck=False)
        path1f = PepperReaction([i1], [B2], 'branch-3way', memorycheck=False)
        path2 = PepperReaction([B2], [i2], 'branch-3way', memorycheck=False)
        path2r = PepperReaction([i2], [B2], 'branch-3way', memorycheck=False)
        path2f = PepperReaction([i2], [B1], 'branch-3way', memorycheck=False)

        backward= PepperReaction([B2], [B1], 'branch-3way', memorycheck=False)

        enum = Enumerator([B1])
        enum.max_helix_migration = True
        enum.enumerate()
        #for r in sorted(enum.reactions): print 'max-helix', r, r.kernel_string, r.rtype, r.rate
        self.assertEqual(sorted([path1, path1r, path1f, path2, path2r, path2f]), sorted(enum.reactions))
        #self.assertEqual(sorted([forward, backward]), sorted(enum.reactions))


    def test_self_displacement(self):
        complexes, reactions = read_kernel("""
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
 
        enum = Enumerator([T])
        enum.max_helix_migration = True
        enum.enumerate()

        path1  = PepperReaction([T], [T1], 'branch-3way', memorycheck=False)
        path1r = PepperReaction([T1], [T], 'branch-3way', memorycheck=False)
        path2  = PepperReaction([T], [T2], 'branch-3way', memorycheck=False)
        path2r = PepperReaction([T2], [T], 'branch-3way', memorycheck=False)

        path3 = PepperReaction([T1], [T3], 'bind11', memorycheck=False)
        path4 = PepperReaction([T2], [T3], 'bind11', memorycheck=False)

        path5  = PepperReaction([T1], [T4], 'branch-3way', memorycheck=False)
        path5r = PepperReaction([T4], [T1], 'branch-3way', memorycheck=False)
        path6  = PepperReaction([T2], [T4], 'branch-3way', memorycheck=False)
        path6r = PepperReaction([T4], [T2], 'branch-3way', memorycheck=False)

        self.assertEqual(sorted([path1, path1r, path2, path2r, path3, path4, 
                                 path5, path5r, path6, path6r]), sorted(enum.reactions))


    def test_compare_semantics(self):
        complexes, reactions = read_kernel("""
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
 
        enum = Enumerator([B2, helper, PR_FL_B2])
        enum.max_helix_migration = True

        path1  = PepperReaction(sorted([PR_FL_B2, helper]), [PR_FLh1B2], 'bind21', memorycheck=False)
        path1r = PepperReaction([PR_FLh1B2], sorted([PR_FL_B2, helper]), 'open', memorycheck=False)
        path2  = PepperReaction(sorted([PR_FL_B2, helper]), [PR_FLh2B2], 'bind21', memorycheck=False)
        path2r = PepperReaction([PR_FLh2B2], sorted([PR_FL_B2, helper]), 'open', memorycheck=False)
        path3  = PepperReaction([PR_FLh1B2], sorted([PR_FL_h1w, B2]), 'branch-3way', memorycheck=False)
        path4  = PepperReaction(sorted([PR_FL_B2, B2]), [PR_FLB2B2], 'bind21', memorycheck=False)
        path4r = PepperReaction([PR_FLB2B2], sorted([PR_FL_B2, B2]), 'open', memorycheck=False)

        path5  = PepperReaction([PR_FLh1B2], [PR_FLh2B2], 'branch-3way', memorycheck=False)
        path6  = PepperReaction([PR_FLh2B2], [PR_FLh1B2], 'branch-3way', memorycheck=False)

        enum.enumerate()

        self.assertTrue(path1 in enum.reactions)
        self.assertTrue(path1r in enum.reactions)
        self.assertTrue(path2 in enum.reactions)
        self.assertTrue(path2r in enum.reactions)
        self.assertTrue(path3 in enum.reactions)
        self.assertTrue(path4 in enum.reactions)
        self.assertTrue(path4r in enum.reactions)
        self.assertTrue(path5 in enum.reactions)
        self.assertTrue(path6 in enum.reactions)

        # CASEY Semantics after fix
        #path6  = PepperReaction([PR_FLh2B2], sorted([PR_FL_h1w, B2]), 'branch-3way', memorycheck=False)

        #path10  = PepperReaction([PR_FLh2B2], [PR_FLh2B2_v2], 'branch-3way', memorycheck=False)
        #path11  = PepperReaction([PR_FLh2B2_v2], sorted([PR_FLh2w, B2]), 'open', memorycheck=False)
        #path12  = PepperReaction([PR_FLh2w], [PR_FL_h1w], 'bind11', memorycheck=False)
        #path13  = PepperReaction([PR_FLh2B2_v2], [PR_FLh1B2], 'branch-3way', memorycheck=False)
        #path14  = PepperReaction([PR_FLh2B2_v2], sorted([PR_FL_h1w, B2]), 'branch-3way', memorycheck=False)
        #self.assertTrue(path10 not in enum.reactions)
        #self.assertTrue(path11 not in enum.reactions)
        #self.assertTrue(path12 not in enum.reactions)
        #self.assertTrue(path13 not in enum.reactions)
        #self.assertTrue(path14 not in enum.reactions)

        #for r in enum.reactions:
        #    #if r not in [path1, path1r, path2, path2r, path3, path4, path4r, path5, path6, path10, path11, path12, path13, path14]:
        #    #if r not in [path1, path1r, path2, path2r, path3, path4, path4r, path5, path6]:
        #    if r not in [path1, path1r, path2, path2r, path3, path4, path4r, path5, path6]:
        #        print 'new-max-helix', r, r.kernel_string(), r.rate
        #    else :
        #        print 'both-max-helix', r, r.kernel_string(), r.rate


if __name__ == '__main__':
  unittest.main()

