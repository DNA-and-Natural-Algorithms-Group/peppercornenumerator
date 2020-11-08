#!/usr/bin/env python
#
#  test_reactions.py
#  EnumeratorProject
#
import unittest

from peppercornenumerator import Enumerator
from peppercornenumerator.input import read_pil_line, read_pil
from peppercornenumerator.objects import PepperReaction, clear_memory
from peppercornenumerator.reactions import (find_on_loop,
                                            filter_bind11,
                                            filter_3way,
                                            filter_4way,
                                            zipper,
                                            bind11,
                                            bind21,
                                            open1N,
                                            branch_3way,
                                            branch_4way,
                                            join_complexes_21)

SKIP = False

def triple(reactant, loc):
    return (reactant.get_domain(loc), reactant.get_paired_loc(loc), loc)

@unittest.skipIf(SKIP, "skipping tests")
class TestFindOnLoop(unittest.TestCase):
    def setUp(self):
        self.a = read_pil_line("domain a = 15")
        self.b = read_pil_line("domain b = 15")
        self.c = read_pil_line("domain c = 15")
        self.d = read_pil_line("domain d = 15")
        self.e = read_pil_line("domain e = 15")
        self.x = read_pil_line("domain x = 15")
        self.y = read_pil_line("domain y = 15")

    def tearDown(self):
        clear_memory()

    def test_bind11(self):
        # need to define domains
        A = read_pil_line("A = a( b c d e( ) d* c* b* c* )")

        folA = find_on_loop(A, (0,1), filter_bind11, direction=1)
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

        [s,x,t,y] = zipper(A, s[0], x, t[0], y, filter_bind11)
        self.assertEqual(s, [triple(A, (0,1)),  # b
                             triple(A, (0,2)),  # c
                             triple(A, (0,3))]) # d
        self.assertEqual(x, [triple(A, (0,4))]) # e
        self.assertEqual(t, [triple(A, (0,8)),  # b*
                             triple(A, (0,7)),  # c*
                             triple(A, (0,6))]) # d*
        self.assertEqual(y, [triple(A, (0,9)),  # c*
                             triple(A, (0,10))])# a*

        folA = find_on_loop(A, (0,2), filter_bind11, direction=1)
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

        [s,x,t,y] = zipper(A, s[0], x, t[0], y, filter_bind11)
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

        [s,x,t,y] = zipper(A, s[0], x, t[0], y, filter_bind11)
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
        B = read_pil_line("B = a( b c d e( + ) d* + c* b* c* )")
        folB = find_on_loop(B, (0,1), filter_bind11, direction=1)

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
        s,x,t,y = zipper(B, s[0], x, t[0], y, filter_bind11)
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

        folB = find_on_loop(B, (0,2), filter_bind11, direction=1)

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

        s,x,t,y = zipper(B, s[0], x, t[0], y, filter_bind11)
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

        s,x,t,y = zipper(B, s[0], x, t[0], y, filter_bind11)
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
        A = read_pil_line("""
        A = a( b( c( x a b c b y c( b( + ) ) + b c ) ) )
        #                  ^(0.6)                  ^(2,2) 
        """)
        folA = find_on_loop(A, (0,6), filter_3way, direction=1)
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

        [s,x,t,y] = zipper(A, s[0], x, t[0], y, filter_3way)
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



        folA = find_on_loop(A, (0,6), filter_3way, direction=-1)
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

        [s,x,t,y] = zipper(A, s[0], x, t[0], y, filter_3way)
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
        x2 = read_pil_line("domain X2 = 15")
        x1 = read_pil_line("domain X1 = 15")
        y1 = read_pil_line("domain Y1 = 15")
        A = read_pil_line("""
        A = c* b*( X2*( X1*( + ) ) ) c y* b*( a* Y1* + a ) c( + ) b*( X2*( X1*( + ) ) ) c y* b*( a*( Y1* + ) ) c
        #      ^(0,1)                     ^(1,5)                  ^(3,3)
        """)
        folA = find_on_loop(A, (0,1), filter_4way)
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

        [s,x,t,y] = zipper(A, s[0], x, t[0], y, filter_4way)
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
class TestBind11(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_bind11_01(self):
        complexes, reactions = read_pil("""
        length a = 10
        length t = 5
        X = a( t + t* )
        Y = a t( + ) a*
        Z = a( t( + ) )
        """)
        X = complexes['X']
        Y = complexes['Y']
        Z = complexes['Z']

        rxnXZ = PepperReaction([X], [Z], 'bind11')
        rxnYZ = PepperReaction([Y], [Z], 'bind11')
        assert bind11(X, max_helix = True) == [rxnXZ]
        assert bind11(Y, max_helix = True) == [rxnYZ]
        assert bind11(Z, max_helix = True) == []
        assert bind11(X, max_helix = True) == [rxnXZ]
        assert bind11(Y, max_helix = True) == [rxnYZ]
        assert bind11(Z, max_helix = True) == []

    def test_bind11_02(self):
        complexes, reactions = read_pil("""
        length a = 10
        length t = 5
        S0  = a a t t* a* a*
        SG1 = a( a t t* ) a* 
        SG2 = a a( t t* a* ) 
        SG3 = a( a( t( ) ) )
        SI1 = a( a t t* a* )
        SI2 = a a( t t* ) a*
        SI3 = a a t( ) a* a*
        """)
        S0  = complexes['S0']
        SG1 = complexes['SG1']
        SG2 = complexes['SG2']
        SG3 = complexes['SG3']
        SI1 = complexes['SI1']
        SI2 = complexes['SI2']
        SI3 = complexes['SI3']

        rg1 = PepperReaction([S0], [SG1], 'bind11')
        rg2 = PepperReaction([S0], [SG2], 'bind11')
        rg3 = PepperReaction([S0], [SG3], 'bind11')
        out = bind11(S0, max_helix = True)
        #for o in out: print(o.kernel_string)
        assert sorted(out) == sorted([rg1, rg2, rg3])

        ri1 = PepperReaction([S0], [SI1], 'bind11')
        ri2 = PepperReaction([S0], [SI2], 'bind11')
        ri3 = PepperReaction([S0], [SI3], 'bind11')
        out = bind11(S0, max_helix = False)
        #for o in out: print(o.kernel_string)
        assert sorted(out) == sorted([rg1, rg2, ri1, ri2, ri3])

@unittest.skipIf(SKIP, "skipping tests")
class TestBind21(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_bind21_01(self):
        complexes, reactions = read_pil("""
        length a = 10
        length t = 5
        S1 = a t 
        S2 = t* a*
        X = a( t + t* )
        Y = a t( + ) a*
        Z = a( t( + ) )
        """)
        S1 = complexes['S1']
        S2 = complexes['S2']
        X = complexes['X']
        Y = complexes['Y']
        Z = complexes['Z']
        singles = sorted([S1, S2])

        rxn1 = PepperReaction([S2, S1], [X], 'bind21')
        rxn2 = PepperReaction([S1, S2], [Y], 'bind21')
        out = bind21(S1, S2, max_helix = False)
        #for o in out: print(o.kernel_string)
        self.assertEqual(out, sorted([rxn1, rxn2]))

        rxn3 = PepperReaction([S2, S1], [Z], 'bind21')
        out = bind21(S1, S2, max_helix = True)
        self.assertEqual(out, [rxn3])

    def test_join_complexes_00(self):
        complexes, reactions = read_pil("""
        length a = 10
        length t = 5
        S1 = a t 
        S2 = t* a*
        """)
        S1 = complexes['S1']
        S2 = complexes['S2']

        c, l1, l2 = join_complexes_21(S1, (0, 1), S2, (0, 0))
        assert c.kernel_string == 'a t + t* a*'
        assert l1 == (0, 1)
        assert l2 == (1, 0)
        n, n1, n2 = join_complexes_21(S2, (0, 0), S1, (0, 1))
        assert n is c
        assert n1 == l1
        assert n2 == l2
        c, l1, l2 = join_complexes_21(S1, (0, 0), S2, (0, 1))
        assert c.kernel_string == 'a t + t* a*'
        assert l1 == (0, 0)
        assert l2 == (1, 1)
        n, n1, n2 = join_complexes_21(S2, (0, 1), S1, (0, 0))
        assert n is c
        assert n1 == l1
        assert n2 == l2

    def test_join_complexes_01(self):
        complexes, reactions = read_pil("""
        length a = 10
        length t = 5
        S1 = a( + a t )
        S2 = t* a*
        S3 = a( a( + t* a*( + ) ) + a( + ) )
        X = a t a*( + ) + t* a*
        XB = a t( a*( + ) + ) a*
        Y = a( a( + a t a*( + ) + t* a*( + ) ) + a( + ) )
        """)
        S1 = complexes['S1']
        S2 = complexes['S2']
        S3 = complexes['S3']
        X = complexes['X']
        Y = complexes['Y']

        c, l1, l2 = join_complexes_21(S1, (1, 1), S2, (0, 0))
        assert c is X
        assert l1 == (0, 1)
        assert l2 == (2, 0)
        n, n1, n2 = join_complexes_21(S2, (0, 0), S1, (1, 1))
        assert n is c
        assert n1 == l1
        assert n2 == l2
        c, l1, l2 = join_complexes_21(S1, (1, 1), S3, (1, 0))
        assert c is Y
        assert l1 == (1, 1)
        assert l2 == (3, 0)
        c, l1, l2 = join_complexes_21(S3, (1, 0), S1, (1, 1))
        assert c is Y
        assert l1 == (1, 1)
        assert l2 == (3, 0)

@unittest.skipIf(SKIP, "skipping tests")
class TestOpen1N(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_open1N_01(self):
        complexes, reactions = read_pil("""
        length a = 8
        length t = 5

        X = a( t( + ) )
        Y = a( t + t* )
        Z = a t( + ) a*
        S1 = a t 
        S2 = t* a*
        """)
        X = complexes['X']
        Y = complexes['Y']
        Z = complexes['Z']
        S1 = complexes['S1']
        S2 = complexes['S2']
        rxnYS = PepperReaction([Y], [S1, S2], 'open')
        rxnXS = PepperReaction([X], [S1, S2], 'open')
        # max helix semantics ON -> no reactions vs dissociation
        out = open1N(Y, max_helix = True, release_11 = 7, release_1N = 7)
        assert out == []
        out = open1N(Y, max_helix = True, release_11 = 8, release_1N = 8)
        assert out == [rxnYS]
        out = open1N(X, max_helix = True, release_11 = 13, release_1N = 13)
        assert out == [rxnXS]

        rxnXY = PepperReaction([X], [Y], 'open')
        rxnXZ = PepperReaction([X], [Z], 'open')
        # max helix semantics OFF -> domains dissociate, but at most one at a time
        out = open1N(X, max_helix = False, release_11 = 7, release_1N = 7)
        self.assertEqual(out, [rxnXY])
        out = open1N(X, max_helix = False, release_11 = 8, release_1N = 8)
        self.assertEqual(out, sorted([rxnXY, rxnXZ]))
        out = open1N(X, max_helix = False, release_11 = 13, release_1N = 13)
        self.assertEqual(out, sorted([rxnXY, rxnXZ]))

    def test_open1N_02(self):
        # NOTE this used to be a test for breathing ...
        # Should release-cutoff of 4 mean: 
        #   "max-helix-semantics up to length 4"? 
        complexes, reactions = read_pil("""
        length a = 3
        length b = 1
        length c = 4
        length ab = 4

        X = a( b( c( + ) ) )
        Y = ab( c( + ) )

        S1 = a b c 
        S2 = c* b* a*
        X1 = a b( c( + ) )  a*
        X2 = a( b  c( + ) b* )
        X3 = a( b( c + c* ) )
        Y1 = ab c( + ) ab*
        Y2 = ab( c + c* )
        """)
        X = complexes['X']
        Y = complexes['Y']
        S1 = complexes['S1']
        S2 = complexes['S2']
        X1 = complexes['X1']
        X2 = complexes['X2']
        X3 = complexes['X3']
        Y1 = complexes['Y1']
        Y2 = complexes['Y2']

        rxn1 = PepperReaction([X], [X1], 'open')
        rxn2 = PepperReaction([X], [X2], 'open')
        rxn3 = PepperReaction([X], [X3], 'open')
        out = open1N(X, max_helix = False, release_11 = 6, release_1N = 6)
        assert sorted(out) == sorted([rxn1, rxn2, rxn3])

        rxn4 = PepperReaction([X], [S1, S2], 'open')
        out = open1N(X, max_helix = True, release_11 = 8, release_1N = 8)
        assert out == [rxn4]

        out = open1N(Y, max_helix = True, release_11 = 6, release_1N = 6)
        assert out == []

@unittest.skipIf(SKIP, "skipping tests")
class TestBranch3Way(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_3way_01(self):
        complexes, reaction = read_pil("""
        length d11 = long
        length d12 = long
        length d7 = long
        length t6 = short
        length t8 = short
        A = d11  t6( d7 t8 + d12( t6( d7 t8 + ) ) ) d11*( t8*( + d7 ) )
        B = d11( t6( d7 t8 + d12( t6( d7 t8 + ) ) ) )     t8*( + d7 ) d11
        """)
        A = complexes['A']
        B = complexes['B']
        rxn = PepperReaction([A], [B], 'branch-3way')
        out = branch_3way(A, max_helix = True)
        assert out == [rxn]

    def test_3way_02(self):
        complexes, reactions = read_pil("""
        length b = 7
        A = b( b*( b b*( + ) b* ) )
        B = b( b*( b( ) + b b* ) )
        C = b( b*( ) b*( + ) b* b )
        D = b( b*( b b* + b( ) ) )
        """)
        A = complexes['A']
        B = complexes['B']
        C = complexes['C']
        D = complexes['D']
        rxn1 = PepperReaction([A], [B], 'branch-3way')
        rxn2 = PepperReaction([A], [C], 'branch-3way')
        rxn3 = PepperReaction([A], [D], 'branch-3way')
        out = branch_3way(A, max_helix = False)
        assert sorted(out) == sorted([rxn1, rxn2, rxn3])

    def test_3way_single(self):
        # A single 3-way branch migration reaction.
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length d = 10
        length x = 10
        X = a( b x + b( d( + ) ) )
        Y = a( b( x + b d( + ) ) )
        """)
        X = complexes['X']
        Y = complexes['Y']
        fw = PepperReaction([X], [Y], 'branch-3way')
        bw = PepperReaction([Y], [X], 'branch-3way')

        out = branch_3way(X, max_helix = True)
        self.assertEqual(out, [fw])
        out = branch_3way(Y, max_helix = True)
        self.assertEqual(out, [bw])
        out = branch_3way(X, max_helix = False)
        self.assertEqual(out, [fw])
        out = branch_3way(Y, max_helix = False)
        self.assertEqual(out, [bw])

    def test_3way_max_helix(self):
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length x = 10
        length y = 10
        length z = 10
        X  = a( x y z + x( y( z( b( + ) ) ) ) )
        I1 = a( x( y z + x y( z( b( + ) ) ) ) )
        I2 = a( x( y( z + x y z( b( + ) ) ) ) )
        Y  = a( x( y( z( + x y z b( + ) ) ) ) )
        """)
        X  = complexes['X']
        I1 = complexes['I1']
        I2 = complexes['I2']
        Y  = complexes['Y']

        fw = PepperReaction([X], [Y], 'branch-3way')
        bw = PepperReaction([Y], [X], 'branch-3way')
        out = branch_3way(X, max_helix = True)
        self.assertEqual(out, [fw])
        out = branch_3way(Y, max_helix = True)
        self.assertEqual(out, [bw])

        fw = PepperReaction([I1], [Y], 'branch-3way')
        bw = PepperReaction([I1], [X], 'branch-3way')
        out = branch_3way(I1, max_helix = True)
        self.assertEqual(out, sorted([bw, fw]))
 
        fw = PepperReaction([X], [I1], 'branch-3way')
        bw = PepperReaction([Y], [I2], 'branch-3way')
        out = branch_3way(X, max_helix = False)
        self.assertEqual(out, [fw])
        out = branch_3way(Y, max_helix = False)
        self.assertEqual(out, [bw])

        fw = PepperReaction([I1], [I2], 'branch-3way')
        bw = PepperReaction([I1], [X], 'branch-3way')
        out = branch_3way(I1, max_helix = False)
        self.assertEqual(sorted(out), sorted([fw, bw]))
 
    def test_3way_remote(self):
        """ 
        A remote 3way branch migration reaction.
        """
        # INPUT
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length x = 10
        length y = 10
        length z = 10
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
        forward = PepperReaction([reactant], [product], 'branch-3way')
        backward = PepperReaction([product], [reactant], 'branch-3way')

        output = branch_3way(reactant, max_helix = True, remote = True)
        self.assertEqual(output, [forward])

        output = branch_3way(product, max_helix = True, remote = True)
        self.assertEqual(output, [backward])

        forward = PepperReaction([inter1], [product], 'branch-3way')
        backward = PepperReaction([inter1], [reactant], 'branch-3way')

        output = branch_3way(inter1, max_helix = True, remote = True)
        self.assertEqual(sorted(output), sorted([backward, forward]))
 
        # ~~~~~~~~~~~~~~~~~~~~ #
        # OUTPUT REJECT REMOTE #
        # ~~~~~~~~~~~~~~~~~~~~ #
        backward = PepperReaction([product], [reactant], 'branch-3way')

        output = branch_3way(reactant, max_helix = True, remote = False)
        self.assertEqual(output, [])
        output = branch_3way(product, max_helix = True, remote = False)
        self.assertEqual(output, [backward])

        forward = PepperReaction([inter1], [product], 'branch-3way')
        backward = PepperReaction([inter1], [reactant], 'branch-3way')

        output = branch_3way(inter1, max_helix = True, remote = False)
        self.assertEqual(sorted(output), sorted([forward, backward]))
 
    def test_3way_multiple_choice_01(self):
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
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length x = 10
        length y = 10
        length z = 10
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
        forward = PepperReaction([reactant], [inter1], 'branch-3way')
        output = branch_3way(reactant, max_helix = True, remote = False)
        self.assertEqual(output, [forward])

        forward1b = PepperReaction([reactant], [inter1], 'branch-3way')
        forward2 = PepperReaction([reactant], [product2], 'branch-3way')
        output = branch_3way(reactant, max_helix = True, remote = True)
        self.assertEqual(sorted(map(str,output)), sorted(map(str,[forward2, forward1b])))

        backward1 = PepperReaction([product1], [inter1], 'branch-3way')
        output = branch_3way(product1, max_helix = True, remote = True)
        self.assertEqual(output, [backward1])

        backward2 = PepperReaction([product2], [reactant], 'branch-3way')
        backward2b = PepperReaction([product2], [product3], 'branch-3way')
        output = branch_3way(product2, max_helix = True, remote = True)
        self.assertEqual(output, sorted([backward2, backward2b]))

        # max_helix zipping cannot involve different strands than the initial step
        backward = PepperReaction([inter1], [reactant], 'branch-3way')
        forward  = PepperReaction([inter1], [product1], 'branch-3way')
        forward2 = PepperReaction([inter1], [product3], 'branch-3way')
        forward3 = PepperReaction([inter1], [product4], 'branch-3way')
        output = branch_3way(inter1, max_helix = True, remote = True)
        self.assertEqual(output, sorted([backward, forward, forward2, forward3]))

    def test_3way_multiple_choice_02(self):
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length x = 10
        length y = 10
        length z = 10
        N = a( x( b x y + y( z( + ) ) ) )
        P1 = a( x b x( y + y( z( + ) ) ) )
        P2 = a( x( b x y( + y z( + ) ) ) )
        """)
        N = complexes['N']
        P1 = complexes['P1']
        P2 = complexes['P2']
        forward1 = PepperReaction([N], [P1], 'branch-3way')
        forward2 = PepperReaction([N], [P2], 'branch-3way')
        output = branch_3way(N, max_helix = True, remote = True)
        self.assertEqual(sorted(output), sorted([forward1, forward2]))

    def test_seesaw_leak(self):
        complexes, reactions = read_pil("""
        # Domain specifications 
        length c = 2
        length s1 = 11
        length s2 = 11
        length s3 = 11
        length s4 = 11
        length t = 3

        # Resting complexes 
        G_g2_w2_3 = c s3 c( t( c( s2( c( + c* t* ) ) ) ) )
        G_w4_2_g2 = c( s2( c( t( c( s4 c + ) ) ) ) ) t* c*

        I12_223 = c s2 c t( c( s1 c + ) ) c*( s2*( c*( t*( c*( + c s3 ) ) ) ) )
        B12_223 = c( s2( c( t( c( s1 c + ) ) ) ) ) t*( c*( + c s3 ) ) c s2 c
        I12_422 = c s2 c( t( c s1 c + c( s2( c( t( c( s4 c + ) ) ) ) ) ) )

        w1_2 = c s2 c t c s1 c
        w2_3 = c s3 c t c s2 c
        w4_2 = c s2 c t c s4 c
        """)
        G223 = complexes['G_g2_w2_3']
        G422 = complexes['G_w4_2_g2']
        w12 = complexes['w1_2']
        w23 = complexes['w2_3']
        w42 = complexes['w4_2']
        I12_223 = complexes['I12_223']
        I12_422 = complexes['I12_422']

        out = bind21(w12, G223, max_helix = True)
        self.assertEqual(len(out), 4)
        self.assertTrue(I12_223 in r.products for r in out)

        out = bind21(w12, G422, max_helix = True)
        self.assertEqual(len(out), 4)
        self.assertTrue(I12_422 in r.products for r in out)

        out = branch_3way(I12_223, max_helix = True)
        self.assertEqual(len(out), 4)
        self.assertTrue(B12_422 in r.products for r in out)

        out = branch_3way(I12_422, max_helix = True)
        assert len(out) == 4
        #for o in out: print(o, o.kernel_string)

@unittest.skipIf(SKIP, "skipping tests")
class TestBranch4Way(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_4way_01(self):
        complexes, reactions = read_pil("""
        length t0 = short
        length d3 = long
        length d4 = long
        A1 = t0*( d3*( d4*( + ) ) + d3*( t0* d3*( d4*( + ) ) + ) )
        A2 = t0*( d3*( d4*( + ) d3( + ) t0* d3*( d4*( + ) ) + ) )
        """)
        A1 = complexes['A1']
        A2 = complexes['A2']
        rxn = PepperReaction([A1], [A2], 'branch-4way')
        out = branch_4way(A1, max_helix = True, remote = True)
        assert out == [rxn]
        #for o in out: print(o, o.kernel_string)

    def test_4way_02(self):
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length x = 10
        length y = 10
        length Y1 = 10
        length Y2 = 10
        length X1 = 10
        length X2 = 10
        B = y*( b*( a*( Y1* + ) b( c( + ) ) Y2* Y1*( + ) a ) )
        A = c* b* X2* X1*( + ) X2 b c y*( b*( a*( Y1* + ) b( c( + ) ) Y2* Y1*( + ) a ) )
        R = c* b* X2* X1*( + ) X2 b c y*( b*( a*( Y1* + ) ) c( + ) b*( Y2* Y1*( + ) a ) )
        """)
        A = complexes['A']
        R = complexes['R']
        fw = PepperReaction([A], [R], 'branch-4way')
        bw = PepperReaction([R], [A], 'branch-4way')
        out = branch_4way(A, max_helix = True, remote = True)
        assert out == [fw]
        #for o in out: print(o, o.kernel_string)
        out = branch_4way(R, max_helix = True, remote = True)
        assert out == [bw]
        #for o in out: print(o, o.kernel_string)

    def test_4way_03(self):
        complexes, reactions = read_pil("""
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
        rxn = PepperReaction([BUG], [FIX], 'branch-4way')
        output = branch_4way(BUG, max_helix = True, remote = True)
        self.assertEqual(output, [rxn])

    def test_4way_04(self):
        # Standard 4way junction, no end-dissociation
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length x = 10
        length y = 10
        A1 = a( b( c( + ) ) x*( + ) b( c( y( + ) ) ) )
        A2 = a( b( c( + ) b*( x*( + ) ) c( y( + ) ) ) )
        A3 = a( b( c( + c*( b*( x*( + ) ) ) y( + ) ) ) )
        """)
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']

        path = PepperReaction([A1], [A2], 'branch-4way')
        output = branch_4way(A1, max_helix=False)
        self.assertEqual(output, [path])

        path = PepperReaction([A3], [A2], 'branch-4way')
        output = branch_4way(A3, max_helix=False)
        self.assertEqual(output, [path])

        path1 = PepperReaction([A2], [A1], 'branch-4way')
        path2 = PepperReaction([A2], [A3], 'branch-4way')
        output = branch_4way(A2, max_helix=False)
        self.assertEqual(output, sorted([path1, path2]))

        path = PepperReaction([A1], [A3], 'branch-4way')
        output = branch_4way(A1, max_helix=True)
        self.assertEqual(output, [path])

        path = PepperReaction([A3], [A1], 'branch-4way')
        output = branch_4way(A3, max_helix=True)
        self.assertEqual(output, [path])

    def test_4way_05(self):
        # Unconventional multi-state 4way junction, no end-dissociation
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length x = 10
        A1 = a( a*( b( + ) ) a*( x*( + ) ) b( c( + ) ) )
        A2 = a( ) b( + ) a( a*( x*( + ) ) b( c( + ) ) )
        A3 = a( a*( b( + ) a( ) x*( + ) ) b( c( + ) ) )
        A4 = a( a*( b( + ) ) ) x*( + ) a( b( c( + ) ) )
        """)
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']
        A4 = complexes['A4']

        path1 = PepperReaction([A1], [A2], 'branch-4way')
        path2 = PepperReaction([A1], [A3], 'branch-4way')
        path3 = PepperReaction([A1], [A4], 'branch-4way')
        output = branch_4way(A1, max_helix = False)
        self.assertEqual(output, sorted([path1, path2, path3]))

        path1 = PepperReaction([A1], [A2], 'branch-4way')
        path2 = PepperReaction([A1], [A3], 'branch-4way')
        path3 = PepperReaction([A1], [A4], 'branch-4way')
        output = branch_4way(A1, max_helix = True)
        self.assertEqual(output, sorted([path1, path2, path3]))

    def dtest_branch4_way_long(self):
        # Unconventional multi-state 4way junction, no end-dissociation
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length d = 10
        length e = 10
        length f = 10
        length g = 10
        length h = 10
        length x = 10
        length y = 10
        length z = 10
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

        path1 = PepperReaction([A1], [A2], 'branch-4way')
        path2 = PepperReaction([A1], [A3], 'branch-4way')
        path3 = PepperReaction([A1], [A4], 'branch-4way')
        path4 = PepperReaction([A1], [A5], 'branch-4way')
        path5 = PepperReaction([A1], [A6], 'branch-4way')
        path6 = PepperReaction([A1], [A7], 'branch-4way')
        output = branch_4way(A1, max_helix = False)
        self.assertEqual(sorted(output), sorted([path1, path2, path3, path4, path5, path6]))

        path1 = PepperReaction([A1], [A2m], 'branch-4way')
        path2 = PepperReaction([A1], [A3m], 'branch-4way')
        path3 = PepperReaction([A1], [A4m], 'branch-4way')
        path4 = PepperReaction([A1], [A5m], 'branch-4way')
        path5 = PepperReaction([A1], [A6m], 'branch-4way')
        path6 = PepperReaction([A1], [A7m], 'branch-4way')
        output = branch_4way(A1, max_helix = True)
        self.assertEqual(sorted(output), sorted([path1, path2, path3, path4, path5, path6]))

    def test_4way_reject_remote(self):
        # Standard 4way junction, no end-dissociation
        complexes, reactions = read_pil("""
        length a = 10
        length b = 10
        length c = 10
        length d = 10
        length x = 10
        length y = 10
        length z = 10
        A0 = a( x( y( z( b(  + ) ) ) ) c*( + ) x( y( z( d( + ) ) ) ) )
        A1 = a( x( y( z( b(  + ) ) ) x*( c*( + ) ) y( z( d( + ) ) ) ) )
        A2 = a( x( y( z( b(  + ) ) y*( x*( c*( + ) ) ) z( d( + ) ) ) ) )
        A3 = a( x( y( z( b(  + ) z*( y*( x*( c*( + ) ) ) ) d( + ) ) ) ) )
        """)
        A0 = complexes['A0']
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']

        path1 = PepperReaction([A1], [A3], 'branch-4way')
        path2 = PepperReaction([A1], [A0], 'branch-4way')
        output = branch_4way(A1, max_helix=True, remote = False)
        self.assertEqual(output, sorted([path1, path2]))

        path1 = PepperReaction([A2], [A0], 'branch-4way')
        path2 = PepperReaction([A2], [A3], 'branch-4way')
        output = branch_4way(A2, max_helix=True, remote = False)
        self.assertEqual(output, sorted([path1, path2]))

        path = PepperReaction([A3], [A0], 'branch-4way')
        output = branch_4way(A3, max_helix=True, remote = False)
        self.assertEqual(output, [path])

        path1 = PepperReaction([A2], [A0], 'branch-4way')
        path2 = PepperReaction([A2], [A3], 'branch-4way')
        output = branch_4way(A2, max_helix=True, remote = False)
        self.assertEqual(output, sorted([path1, path2]))

        path1 = PepperReaction([A1], [A3], 'branch-4way')
        path2 = PepperReaction([A1], [A0], 'branch-4way')
        output = branch_4way(A1, max_helix=True, remote = True)
        self.assertEqual(output, sorted([path1, path2]))

        path = PepperReaction([A3], [A0], 'branch-4way')
        output = branch_4way(A3, max_helix=True, remote = False)
        self.assertEqual(output, [path])

@unittest.skipIf(SKIP, "skipping tests")
class ReactionMatchingTests (unittest.TestCase):
    # A selection of Casey's reaction tests.
    def tearDown(self):
        clear_memory()

    def testBind11A(self):
        # bind11: a ? a* ? -> a( ? ) ?
        complexes, _ = read_pil("""
        length a = 10
        length b = 10
        length x = 10
        A1 = x( ) a  x( ) a* x( )
        A2 = x( ) a( x( ) )  x( )
        A3 = x( ) a  b  x( ) b* a* x( )
        A4 = x( ) a( b( x( ) )  )  x( )
        """)

        # No zipping possible
        rxns = bind11(complexes['A1'])
        self.assertEqual(rxns, [PepperReaction([complexes['A1']], [complexes['A2']], 'bind11')])

        # Zipping possible
        rxns = bind11(complexes['A3'])
        self.assertEqual(rxns, [PepperReaction([complexes['A3']], [complexes['A4']], 'bind11')])

    def testBind21A(self):
        # bind11: a ? a* ? -> a( ? ) ?
        complexes, _ = read_pil("""
        length a = 10
        length w = 10
        length x = 10
        length y = 10
        length z = 10
        A1 = w( ) a  x( )
        A2 = y( ) a* z( )
        A3 = w( ) a( x( ) + y( ) ) z( )
        """)

        # No zipping possible
        rxns = bind21(complexes['A1'], complexes['A2'])
        self.assertEqual(rxns,
            [PepperReaction([complexes['A1'], complexes['A2']], 
            [complexes['A3']], 'bind21')])

    def testOpenA(self):
        # open:  ? a( ? ) ? -> ? a ? a* ?
        complexes, _ = read_pil("""
        length a = 5
        length b = 5
        length x = 15
        length y = 15
        length z = 15

        A1 = x( ) a( x )  x( )
        A2 = x( ) a  x a* x( )
        A3 = x a( b( y )  )  z
        A4 = x a  b  y b* a* z
        """)

        # No zipping possible
        self.assertEqual(open1N(complexes['A1']),
            [PepperReaction([complexes['A1']], [complexes['A2']], 'open')])

        # Zipping possible
        self.assertEqual(open1N(complexes['A3'], release_11 = 13, release_1N=13),
            [PepperReaction([complexes['A3']], [complexes['A4']], 'open')])

    def testOpenB(self):
        # open:  ? a( ? ) ? -> ? a ? a* ?
        complexes, _ = read_pil("""
        length a = 5
        length b = 5
        length x = 15
        length y = 15
        length z = 15

        A1 = x a( y )   z
        A2 = x a  y a* z

        A3 = x a( b( y )  )  z
        A4 = x a  b( y )  a* z
        A5 = x a( b  y b* )  z
        """)

        # No zipping possible
        rxns = open1N(complexes['A1'], max_helix=False, release_11 = 7, release_1N=7)
        self.assertEqual(rxns,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'open')])

        # Zipping possible
        rxns = open1N(complexes['A3'], max_helix=False, release_11 = 7, release_1N=7)
        self.assertEqual(rxns, sorted([
            PepperReaction([complexes['A3']], [complexes['A4']], 'open'),
            PepperReaction([complexes['A3']], [complexes['A5']], 'open')]))

    def testOpenNoMaxHelix(self):
        # open:  ? a( ? ) ? -> ? a ? a* ?
        complexes, _ = read_pil("""
        length a = 5
        length b = 5

        # A1 = a( b( ) )
        A1 = a( b( ) )
        A2 = a( b b* )
        A3 = a b( ) a*
        """)

        # Zipping possible
        rxns = open1N(complexes['A1'], max_helix=False, release_11 = 10, release_1N=10)
        self.assertEqual(rxns, sorted([
            PepperReaction([complexes['A1']], [complexes['A2']], 'open'),
            PepperReaction([complexes['A1']], [complexes['A3']], 'open')]))

    def testBranch3wayA(self):
        # 3wayA: ? b ? b(?) ? <-> ? b(? b ?) ?
        complexes, _ = read_pil("""
        length b = 15
        length d1 = 15
        length d2 = 15
        length d3 = 15
        length d4 = 15

        A1 = d1( ) b  d2( ) b( d3( ) ) d4( )
        A2 = d1( ) b( d2( ) b  d3( ) ) d4( )
        """)
        forward = branch_3way(complexes['A1'])

        self.assertEqual(forward,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'branch-3way')])

        reverse = branch_3way(complexes['A2'])
        self.assertEqual(reverse,
        [PepperReaction([complexes['A2']], [complexes['A1']], 'branch-3way')])

    def testBranch3wayB(self):
        # 3wayB: ? b(?) ? b ? <-> ? b ? b*(?) ?
        complexes, _ = read_pil("""
        length b = 15
        length d1 = 15
        length d2 = 15
        length d3 = 15
        length d4 = 15
        A1 = d1( ) b( d2( )  )  d3( ) b d4( )
        A2 = d1( ) b  d2( ) b*( d3( ) ) d4( )
        """)
        forward = branch_3way(complexes['A1'])
        self.assertEqual(forward,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'branch-3way')])

        reverse = branch_3way(complexes['A2'])
        self.assertEqual(reverse,
        [PepperReaction([complexes['A2']], [complexes['A1']], 'branch-3way')])

    def testBranch3wayC(self):
        # 3wayC: ? b*(?) ? b ? <-> ? b*(? b ?) ?
        complexes, _ = read_pil("""
        length b = 15
        length d1 = 15
        length d2 = 15
        length d3 = 15
        length d4 = 15
        A1 = d1( ) b*( d2( ) ) d3( ) b d4( )
        A2 = d1( ) b*( d2( ) b d3( ) ) d4( )
        """)
        forward = branch_3way(complexes['A1'])
        self.assertEqual(forward,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'branch-3way')])

        reverse = branch_3way(complexes['A2'])
        self.assertEqual(reverse,
        [PepperReaction([complexes['A2']], [complexes['A1']], 'branch-3way')])

    def testBranch4wayA(self):
        # 4way: b( ? ) ? b ( ? ) --> b( ? b*( ? ) ? )
        complexes, _ = read_pil("""
        length b = 15
        length d1 = 15
        length d2 = 15
        length d3 = 15
        length d4 = 15
        length d5 = 15
        A1 = d1( ) b( d2( )  )  d3( ) b( d4( ) ) d5( )
        A2 = d1( ) b( d2( ) b*( d3( ) )  d4( ) ) d5( )
        """)

        forward = branch_4way(complexes['A1'])
        self.assertEqual(forward,
        [PepperReaction([complexes['A1']], [complexes['A2']], 'branch-4way')])

        reverse = branch_4way(complexes['A2'])
        self.assertEqual(reverse,
        [PepperReaction([complexes['A2']], [complexes['A1']], 'branch-4way')])

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
        output = bind21(I, C)
        self.assertEqual(output, [path1])

        # DSD-pathway "branch3way"
        path2 = PepperReaction([B], [D, J], 'branch-3way')
        output = branch_3way(B)
        self.assertEqual(output, [path2])

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
        #print()
        #for r in enum.reactions:
        #    print(r, r.kernel_string)
        self.assertEqual(len(enum.reactions), 22)

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
        self.assertEqual(len(enum2.reactions), len(enum.reactions))

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
        self.assertEqual(len(enum.reactions), 2)

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

