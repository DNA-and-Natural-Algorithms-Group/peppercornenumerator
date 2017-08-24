#
#  test_reactions.py
#

import copy
import unittest


from peppercornenumerator import Enumerator
from peppercornenumerator.pil_parser import parse_pil_string
import peppercornenumerator.reactions as rxn

# Input parsing stuff
from peppercornenumerator.utils import PepperDomain, PepperComplex
from peppercornenumerator.dsdobjects import reset_names
from peppercornenumerator.input import from_kernel

SKIP = False

def nuskell_parser(pil_string, ddlen=15):
    # snatched from nuskell.objects.TestTube.load_pil_kernel
    ppil = parse_pil_string(pil_string)

    #for p in ppil: print p

    def resolve_loops(loop):
      """ Return a sequence, structure pair from kernel format with parenthesis. """
      sequen = []
      struct = []
      for dom in loop :
        if isinstance(dom, str):
          sequen.append(dom)
          if dom == '+' :
            struct.append('+')
          else :
            struct.append('.')
        elif isinstance(dom, list):
          struct[-1] = '('
          old = sequen[-1]
          se, ss = resolve_loops(dom)
          sequen.extend(se)
          struct.extend(ss)
          sequen.append(old + '*' if old[-1] != '*' else old[:-1])
          struct.append(')')
      return sequen, struct

    domains = {}
    strands = {}
    complexes = {}
    for line in ppil :
      if line[0] == 'domain':
          #domains[line[1]] = PepperDomain(line[1], int(line[2]), sequence = 'N'* int(line[2]))
          domains[line[1]] = PepperDomain(list('N'* int(line[2])), name=line[1])
      elif line[0] == 'complex':
        name = line[1]
        sequence, structure = resolve_loops(line[2])

        strand = [] # construct a strand
        cplx_strands = [] # store strands for Complex
        for e in range(len(sequence)):
          d = sequence[e]
          if d == '+': 
              sid = '_'.join(tuple(map(str,strand)))
              if sid not in strands :
                  strands[sid] = None
              cplx_strands.append(strands[sid])
              strand = [] # construct a strand
              continue
          if d[-1] == '*' : 
            dname = d[:-1]
            if dname in domains :
                cdom = domains[dname]
                #dom = PepperDomain(cdom.name, len(cdom), , is_complement=True)
                dom = cdom.get_ComplementDomain(list(sequence = 'N'* len(cdom)))
            else :
                #cdom = PepperDomain(dname, ddlen, sequence = 'N'* ddlen, is_complement=False)
                cdom = PepperDomain(list('N'* ddlen), name=dname)
                domains[dname] = cdom
                #dom = Domain(dname, ddlen, sequence = 'N'* ddlen, is_complement=True)
                dom = cdom.get_ComplementDomain(list(sequence = 'N'* len(cdom)))
          else :
            dname = d
            if dname in domains :
                dom = domains[dname]
            else :
                #dom = Domain(dname, ddlen, sequence = 'N'* ddlen)
                dom = PepperDomain(list('N'* ddlen), name=dname)
                domains[dname] = dom
          strand.append(dom)

        sid = '_'.join(tuple(map(str,strand)))
        if sid not in strands :
            strands[sid] = None
        cplx_strands.append(strands[sid])
        strand = [] # construct a strand

        for e, dom in enumerate(sequence):
            if dom == '+':
                domain = '+'
            elif dom[-1] == '*' : 
                domain = ~domains[dom[:-1]]
            else :
                domain = domains[dom]
            sequence[e] = domain

        complex = PepperComplex(sequence, structure, name=name)
        complexes[name] = complex
      else :
        raise NotImplementedError('Weird expression returned from pil_parser!')

    return (domains, strands, complexes)

@unittest.skipIf(SKIP, "skipping tests")
class NewOpenTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        reset_names()

    def test_basic_open(self):
        """ 
        A basic open reaction.

        Testing max-helix-semantics and release-cutoff 5, 8, 13
        """
        # INPUT
        (domains, strands, complexes) = nuskell_parser("""
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

        forward = rxn.ReactionPathway('open', [reactant], product_set)
        output = rxn.open(reactant, max_helix=True, release_11=13, release_1N=13)
        #for o in output: print 'ow', o.kernel_string()
        self.assertEqual(output, [forward])

        # max helix semantics OFF -> domains dissociate, but at most one at a time
        forward1 = rxn.ReactionPathway('open', [reactant], [product1])
        forward2 = rxn.ReactionPathway('open', [reactant], [product2])
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
        (domains, strands, complexes) = nuskell_parser("""
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
        reset_names()

    def test_binding(self):
        (domains, strands, complexes) = nuskell_parser("""
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

        path = rxn.ReactionPathway('bind11', [X], [Z])
        output = rxn.bind11(X, max_helix=True)
        #for o in output: print 'mh', o.kernel_string()
        self.assertEqual(output, [path])

        path = rxn.ReactionPathway('bind11', [Y], [Z])
        output = rxn.bind11(Y, max_helix=True)
        #for o in output: print 'mh', o.kernel_string()
        self.assertEqual(output, [path])

        output = rxn.bind11(Z, max_helix=True)
        #for o in output: print 'mh', o.kernel_string()
        self.assertEqual(output, [])

        path1 = rxn.ReactionPathway('bind21', [S2, S1], [X])
        path2 = rxn.ReactionPathway('bind21', [S2, S1], [Y])
        output = rxn.bind21(S1, S2, max_helix=False)
        #for o in output: print 'bind21', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2]))

        path = rxn.ReactionPathway('bind21', [S2, S1], [Z])
        output = rxn.bind21(S1, S2, max_helix=True)
        #for o in output: print 'bind21g', o.kernel_string()
        self.assertEqual(output, [path])

        path1 = rxn.ReactionPathway('bind11', [SB], [SG1])
        path2 = rxn.ReactionPathway('bind11', [SB], [SG2])
        path3 = rxn.ReactionPathway('bind11', [SB], [SG3])
        output = rxn.bind11(SB, max_helix=True)
        #for o in output: print 'open11g', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3]))

        path1 = rxn.ReactionPathway('bind11', [SB], [SG1])
        path2 = rxn.ReactionPathway('bind11', [SB], [SG2])
        path3 = rxn.ReactionPathway('bind11', [SB], [SI1])
        path4 = rxn.ReactionPathway('bind11', [SB], [SI2])
        path5 = rxn.ReactionPathway('bind11', [SB], [SI3])
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
        reset_names()

    def test_single_migration(self):
        """ 
        A single 3-way branch migration reaction.
        """
        # INPUT
        (domains, strands, complexes) = nuskell_parser("""
        X = a( b x + b( d( + ) ) )
        Y = a( b( x + b d( + ) ) )
        """)
        reactant = complexes['X']
        product = complexes['Y']

        # OUTPUT
        forward = rxn.ReactionPathway('branch_3way', [reactant], [product])
        backward = rxn.ReactionPathway('branch_3way', [product], [reactant])

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
        (domains, strands, complexes) = nuskell_parser("""
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
        forward = rxn.ReactionPathway('branch_3way', [reactant], [product])
        backward = rxn.ReactionPathway('branch_3way', [product], [reactant])

        output = rxn.branch_3way(reactant, max_helix=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, max_helix=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        forward = rxn.ReactionPathway('branch_3way', [inter1], [product])
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])

        output = rxn.branch_3way(inter1, max_helix=True)
        #for o in output: print 'max_helix', o.kernel_string()
        self.assertEqual(output, sorted([backward, forward]))
 
        # ~~~~~~~~~~~~~~~~~~~ #
        # OUTPUT NO-MAX-HELIX #
        # ~~~~~~~~~~~~~~~~~~~ #
        forward = rxn.ReactionPathway('branch_3way', [reactant], [inter1])
        backward = rxn.ReactionPathway('branch_3way', [product], [inter2])

        output = rxn.branch_3way(reactant, max_helix=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, max_helix=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        # OUTPUT NO-MAX-HELIX
        forward = rxn.ReactionPathway('branch_3way', [inter1], [inter2])
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])

        output = rxn.branch_3way(inter1, max_helix=False)
        #for o in output: print 'nmheli', o.kernel_string()
        self.assertEqual(sorted(output), sorted([forward, backward]))
 
    def test_remote_migration(self):
        """ 
        A remote 3way branch migration reaction.
        """
        # INPUT
        (domains, strands, complexes) = nuskell_parser("""
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
        forward = rxn.ReactionPathway('branch_3way', [reactant], [product])
        backward = rxn.ReactionPathway('branch_3way', [product], [reactant])

        output = rxn.branch_3way(reactant, max_helix=True, remote=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [forward])

        output = rxn.branch_3way(product, max_helix=True, remote=True)
        #print output[0].kernel_string()
        self.assertEqual(output, [backward])

        forward = rxn.ReactionPathway('branch_3way', [inter1], [product])
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])

        output = rxn.branch_3way(inter1, max_helix=True, remote=True)
        #for o in output: print 'max_helix', o.kernel_string()
        self.assertEqual(sorted(output), sorted([forward, backward]))
 
        # ~~~~~~~~~~~~~~~~~~~~ #
        # OUTPUT REJECT REMOTE #
        # ~~~~~~~~~~~~~~~~~~~~ #
        backward = rxn.ReactionPathway('branch_3way', [product], [reactant])

        output = rxn.branch_3way(reactant, max_helix=True, remote=False)
        #print output[0].kernel_string()
        self.assertEqual(output, [])

        output = rxn.branch_3way(product, max_helix=True, remote=False)
        #for o in output: print 'max_helix', o.kernel_string()
        self.assertEqual(output, [backward])

        forward = rxn.ReactionPathway('branch_3way', [inter1], [product])
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant])

        output = rxn.branch_3way(inter1, max_helix=True, remote=False)
        #for o in output: print 'max_helix', o.kernel_string()
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
        (domains, strands, complexes) = nuskell_parser("""
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
        forward = rxn.ReactionPathway('branch_3way', [reactant], [inter1])
        output = rxn.branch_3way(reactant, max_helix=True, remote=False)
        self.assertEqual(output, [forward])

        # NOTE: NO max_helix REMOTE TOEHOLDS!
        forward1b = rxn.ReactionPathway('branch_3way', [reactant], [inter1])
        forward2 = rxn.ReactionPathway('branch_3way', [reactant], [product2])
        output = rxn.branch_3way(reactant, max_helix=True, remote=True)
        self.assertEqual(output, sorted([forward2, forward1b]))

        backward1 = rxn.ReactionPathway('branch_3way', [product1], [reactant])
        output = rxn.branch_3way(product1, max_helix=True, remote=True)
        self.assertEqual(output, [backward1])

        # NOTE: THIS behavior changed!
        backward2 = rxn.ReactionPathway('branch_3way', [product2], [reactant])
        #backward2b = rxn.ReactionPathway('branch_3way', [product2], [product4])
        backward2b = rxn.ReactionPathway('branch_3way', [product2], [product3])
        output = rxn.branch_3way(product2, max_helix=True, remote=True)
        #for o in output: print 'ow', o.kernel_string()
        self.assertEqual(output, sorted([backward2, backward2b]))

        # NOTE: max_helix zippering cannot involve different strands than the initial step
        backward = rxn.ReactionPathway('branch_3way', [inter1], [reactant]) #1
        forward  = rxn.ReactionPathway('branch_3way', [inter1], [product1]) #4
        forward2 = rxn.ReactionPathway('branch_3way', [inter1], [product3])
        forward3 = rxn.ReactionPathway('branch_3way', [inter1], [product4])
        output = rxn.branch_3way(inter1, max_helix=True, remote=True)
        #for o in output: print 'fixme', o.kernel_string()
        self.assertEqual(output, sorted([backward, forward, forward2, forward3]))

    def test_multiple_choice_2(self):
        # INPUT
        (domains, strands, complexes) = nuskell_parser("""
        N = a( x( b x y + y( z( + ) ) ) )
        P1 = a( x b x( y + y( z( + ) ) ) )
        P2 = a( x( b x y( + y z( + ) ) ) )
        """)
        N = complexes['N']
        P1 = complexes['P1']
        P2 = complexes['P2']

        forward1 = rxn.ReactionPathway('branch_3way', [N], [P1])
        forward2 = rxn.ReactionPathway('branch_3way', [N], [P2])
        output = rxn.branch_3way(N, max_helix=True, remote=True)
        #for o in output: print 'new-max', o.kernel_string()
        self.assertEqual(output, sorted([forward1, forward2]))


@unittest.skipIf(SKIP, "skipping tests")
class NewBranch4WayTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        reset_names()

    def test_break_casey_4way(self):
        # Note the new max-helix semantics also fixes the old 4way error,
        # so this test can be safely removed...
        (domains, strands, complexes) = nuskell_parser("""
        A1 = t0*( d3*( d4*( + ) ) + d3*( t0* d3*( d4*( + ) ) + ) )
        """)

        A1 = complexes['A1']
        output = rxn.branch_4way(A1, max_helix=True, remote=True)
        #for o in output: print 'branch_4way_bug', o.kernel_string()

    def test_4wayfilter_bugfix(self):
        # a test to ensure the 4way filter includes struct1
        (domains, strands, complexes) = nuskell_parser("""

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

        path = rxn.ReactionPathway('branch_4way', [BUG], [FIX])
        output = rxn.branch_4way(BUG, max_helix=True, remote=True)
        #for o in output: print 'branch_4way_bug', o.kernel_string()
        self.assertEqual(output, [path])

    def test_branch4_way(self):
        # Standard 3state 4way junction, no end-dissociation
        (domains, strands, complexes) = nuskell_parser("""
        A1 = a( b( c( + ) ) x*( + ) b( c( y( + ) ) ) )
        A2 = a( b( c( + ) b*( x*( + ) ) c( y( + ) ) ) )
        A3 = a( b( c( + c*( b*( x*( + ) ) ) y( + ) ) ) )
        """)
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']

        path = rxn.ReactionPathway('branch_4way', [A1], [A2])
        output = rxn.branch_4way(A1, max_helix=False)
        #for o in output: print 'branch_4way one-step', o.kernel_string(), o.rate
        self.assertEqual(output, [path])

        path = rxn.ReactionPathway('branch_4way', [A3], [A2])
        output = rxn.branch_4way(A3, max_helix=False)
        #for o in output: print 'branch_4way one-step', o.kernel_string(), o.rate
        self.assertEqual(output, [path])

        path1 = rxn.ReactionPathway('branch_4way', [A2], [A1])
        path2 = rxn.ReactionPathway('branch_4way', [A2], [A3])
        output = rxn.branch_4way(A2, max_helix=False)
        #for o in output: print 'branch_4way one-step', o.kernel_string(), o.rate
        self.assertEqual(output, sorted([path1, path2]))

        path = rxn.ReactionPathway('branch_4way', [A1], [A3])
        output = rxn.branch_4way(A1, max_helix=True)
        #for o in output: print 'branch_4way_two-step', o.kernel_string(), o.rate
        self.assertEqual(output, [path])

        path = rxn.ReactionPathway('branch_4way', [A3], [A1])
        output = rxn.branch_4way(A3, max_helix=True)
        #for o in output: print 'branch_4way_two-step', o.kernel_string(), o.rate
        self.assertEqual(output, [path])

    def test_branch4_way_2(self):
        # Unconventional multi-state 4way junction, no end-dissociation
        (domains, strands, complexes) = nuskell_parser("""
        A1 = a( a*( b( + ) ) a*( x*( + ) ) b( c( + ) ) )
        A2 = a( ) b( + ) a( a*( x*( + ) ) b( c( + ) ) )
        A3 = a( a*( b( + ) a( ) x*( + ) ) b( c( + ) ) )
        A4 = a( a*( b( + ) ) ) x*( + ) a( b( c( + ) ) )
        """)
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']
        A4 = complexes['A4']

        path1 = rxn.ReactionPathway('branch_4way', [A1], [A2])
        path2 = rxn.ReactionPathway('branch_4way', [A1], [A3])
        path3 = rxn.ReactionPathway('branch_4way', [A1], [A4])
        output = rxn.branch_4way(A1, max_helix=False)
        #for o in output: print 'branch_4way_next', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3]))

        path1 = rxn.ReactionPathway('branch_4way', [A1], [A2])
        path2 = rxn.ReactionPathway('branch_4way', [A1], [A3])
        path3 = rxn.ReactionPathway('branch_4way', [A1], [A4])
        output = rxn.branch_4way(A1, max_helix=True)
        #for o in output: print 'branch_4way_next_mh', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3]))

    def test_branch4_way_long(self):
        # Unconventional multi-state 4way junction, no end-dissociation
        (domains, strands, complexes) = nuskell_parser("""
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

        path1 = rxn.ReactionPathway('branch_4way', [A1], [A2])
        path2 = rxn.ReactionPathway('branch_4way', [A1], [A3])
        path3 = rxn.ReactionPathway('branch_4way', [A1], [A4])
        path4 = rxn.ReactionPathway('branch_4way', [A1], [A5])
        path5 = rxn.ReactionPathway('branch_4way', [A1], [A6])
        path6 = rxn.ReactionPathway('branch_4way', [A1], [A7])
        output = rxn.branch_4way(A1, max_helix=False)
        #for o in output: print 'b_4way', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3, path4, path5, path6]))

        path1 = rxn.ReactionPathway('branch_4way', [A1], [A2m])
        path2 = rxn.ReactionPathway('branch_4way', [A1], [A3m])
        path3 = rxn.ReactionPathway('branch_4way', [A1], [A4m])
        path4 = rxn.ReactionPathway('branch_4way', [A1], [A5m])
        path5 = rxn.ReactionPathway('branch_4way', [A1], [A6m])
        path6 = rxn.ReactionPathway('branch_4way', [A1], [A7m])

        output = rxn.branch_4way(A1, max_helix=True)
        #for o in output: print 'b_4way_mh', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2, path3, path4, path5, path6]))

    def test_branch4_way_no_remote(self):
        # Standard 4state 4way junction, no end-dissociation
        (domains, strands, complexes) = nuskell_parser("""
        A0 = a( x( y( z( b(  + ) ) ) ) c*( + ) x( y( z( d( + ) ) ) ) )
        A1 = a( x( y( z( b(  + ) ) ) x*( c*( + ) ) y( z( d( + ) ) ) ) )
        A2 = a( x( y( z( b(  + ) ) y*( x*( c*( + ) ) ) z( d( + ) ) ) ) )
        A3 = a( x( y( z( b(  + ) z*( y*( x*( c*( + ) ) ) ) d( + ) ) ) ) )
        """)
        A0 = complexes['A0']
        A1 = complexes['A1']
        A2 = complexes['A2']
        A3 = complexes['A3']

        path1 = rxn.ReactionPathway('branch_4way', [A1], [A3])
        path2 = rxn.ReactionPathway('branch_4way', [A1], [A0])
        output = rxn.branch_4way(A1, max_helix=True, remote = False)
        #for o in output: print 'branch_4way', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2]))
        path1 = rxn.ReactionPathway('branch_4way', [A2], [A0])
        path2 = rxn.ReactionPathway('branch_4way', [A2], [A3])
        output = rxn.branch_4way(A2, max_helix=True, remote = False)
        #for o in output: print 'branch_4way', o.kernel_string(), o.rate
        self.assertEqual(output, sorted([path1, path2]))

        path = rxn.ReactionPathway('branch_4way', [A3], [A0])
        output = rxn.branch_4way(A3, max_helix=True, remote = False)
        #for o in output: print 'branch_4way', o.kernel_string(), o.rate
        self.assertEqual(output, [path])

        path1 = rxn.ReactionPathway('branch_4way', [A2], [A0])
        path2 = rxn.ReactionPathway('branch_4way', [A2], [A3])
        output = rxn.branch_4way(A2, max_helix=True, remote = False)
        #for o in output: print 'branch_4way', o.kernel_string(), o.rate
        self.assertEqual(output, sorted([path1, path2]))

        path1 = rxn.ReactionPathway('branch_4way', [A1], [A3])
        path2 = rxn.ReactionPathway('branch_4way', [A1], [A0])
        output = rxn.branch_4way(A1, max_helix=True, remote = True)
        #for o in output: print 'branch_4way', o.kernel_string()
        self.assertEqual(output, sorted([path1, path2]))

        path = rxn.ReactionPathway('branch_4way', [A3], [A0])
        output = rxn.branch_4way(A3, max_helix=True, remote = False)
        #for o in output: print 'branch_4way', o.kernel_string(), o.rate
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
        reset_names()

    def test_bind_and_displace3way(self):
        # Skip the outer loop of the enumerator...
        (domains, strands, complexes) = nuskell_parser("""
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
        path1 = rxn.ReactionPathway('bind21', sorted([I, C]), [B])
        output = rxn.bind21(I, C, max_helix=True)
        self.assertEqual(output, [path1])

        # DSD-pathway "branch3way"
        path2 = rxn.ReactionPathway('branch_3way', [B], sorted([D, J]))
        output = rxn.branch_3way(B, max_helix=True)
        #for o in output: print 'test', o.kernel_string()
        self.assertEqual(output, [path2])

        #enum = Enumerator(domains.values(), strands.values(), complexes.values())
        enum = Enumerator(complexes.values())
        enum.enumerate()
        self.assertEqual(sorted(enum.reactions), sorted([path1, path2]))

    def test_cooperative_binding(self):
        (domains, strands, complexes) = nuskell_parser("""
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

        path1 = rxn.ReactionPathway('bind21', sorted([L, C]), [LC])
        path1r = rxn.ReactionPathway('open', [LC], sorted([L, C]))
        path2 = rxn.ReactionPathway('branch_3way', [LC], [LCF])
        path3 = rxn.ReactionPathway('bind21', sorted([R, LCF]), [LCRF1])
        path4 = rxn.ReactionPathway('branch_3way', [LCRF1], sorted([LR, T]))

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

        # NOTE: condensation has no effect for cooperative binding with k_fast 
        from peppercornenumerator.condense import condense_resting_states
        condensed = condense_resting_states(enum, compute_rates = True, k_fast=self.k_fast)
        #for r in condensed['reactions']:
        #    print r.kernel_string(), r.rate
        self.assertEqual(len(enum.reactions), 22)


class NeighborhoodSearch(unittest.TestCase):
    # Test a basic move set and enumerate using k-fast/k-slow
    pass

@unittest.skipIf(SKIP, "skipping tests")
class IsomorphicSets(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        reset_names()

    def test_simple(self):
        # works just fine
        (domains, strands, complexes) = nuskell_parser("""
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
        (domains, strands, complexes) = nuskell_parser("""

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

        path1 = rxn.ReactionPathway('branch_3way', [A1], sorted([A1_2, YZ]))
        path2 = rxn.ReactionPathway('branch_3way', [A2], sorted([A2_1, Y1]))
        path3 = rxn.ReactionPathway('branch_3way', [A2_1], sorted([A1_2, Z1]))

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

        path1 = rxn.ReactionPathway('branch_3way', [B1], sorted([B1_2, YZ2]))
        path2 = rxn.ReactionPathway('branch_3way', [B2], sorted([B2_1, Y2]))
        path3 = rxn.ReactionPathway('branch_3way', [B2_1], sorted([B1_2, Z2]))
        #for r in enum.reactions:
        #    print 'invade', r, r.kernel_string(), r.rate
        self.assertEqual(sorted(enum.reactions), sorted([path1, path2, path3]))

@unittest.skipIf(SKIP, "skipping tests")
class Compare_MaxHelix(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        reset_names()

    def test_self_displacement_bug(self):
        (domains, strands, complexes) = nuskell_parser("""
        B1 = x( y( x( y x + ) ) )
        B2 = x y x( y( x( + ) ) )

        B3 = x( y( x y x( + ) ) )
        B4 = x( y x y( x( + ) ) )
        """)
        B1 = complexes['B1']
        B2 = complexes['B2']
        B3 = complexes['B3']
        B4 = complexes['B4']

        forward = rxn.ReactionPathway('branch_3way', [B1], [B2])
        backward= rxn.ReactionPathway('branch_3way', [B2], [B1])

        path1  = rxn.ReactionPathway('branch_3way', [B1], [B4])
        path1r  = rxn.ReactionPathway('branch_3way', [B4], [B1])
        path2 = rxn.ReactionPathway('branch_3way', [B4], [B2])

        path3  = rxn.ReactionPathway('branch_3way', [B2], [B3])
        path3r = rxn.ReactionPathway('branch_3way', [B3], [B2])
        path4 = rxn.ReactionPathway('branch_3way', [B3], [B1])

        enum = Enumerator([B1])
        enum.max_helix_migration = True
        enum.enumerate()
        #for r in sorted(enum.reactions): print 'max-helix', r, r.kernel_string(), r.rate
        #self.assertEqual(sorted([path1, path1r, path2, path3, path3r, path4]), sorted(enum.reactions))
        self.assertEqual(sorted([forward, backward]), sorted(enum.reactions))

    def test_self_displacement_bug_iso(self):
        (domains, strands, complexes) = nuskell_parser("""
        B1 = x1( x2( y1( y2( x1( x2( y1 y2 x1 x2 + ) ) ) ) ) )
        B2 = x1 x2 y1 y2 x1( x2( y1( y2( x1( x2( + ) ) ) ) ) )

        """)
        B1 = complexes['B1']
        B2 = complexes['B2']

        forward = rxn.ReactionPathway('branch_3way', [B1], [B2])
        backward= rxn.ReactionPathway('branch_3way', [B2], [B1])

        enum = Enumerator([B1])
        enum.max_helix_migration = True
        enum.enumerate()
        #for r in sorted(enum.reactions): print 'max-helix', r, r.kernel_string(), r.rate
        self.assertEqual(sorted([forward, backward]), sorted(enum.reactions))


    def test_self_displacement(self):
        (domains, strands, complexes) = nuskell_parser("""
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

        path1  = rxn.ReactionPathway('branch_3way', [T], [T1])
        path1r = rxn.ReactionPathway('branch_3way', [T1], [T])
        path2  = rxn.ReactionPathway('branch_3way', [T], [T2])
        path2r = rxn.ReactionPathway('branch_3way', [T2], [T])

        path3 = rxn.ReactionPathway('bind11', [T1], [T3])
        path4 = rxn.ReactionPathway('bind11', [T2], [T3])

        path5  = rxn.ReactionPathway('branch_3way', [T1], [T4])
        path5r = rxn.ReactionPathway('branch_3way', [T4], [T1])
        path6  = rxn.ReactionPathway('branch_3way', [T2], [T4])
        path6r = rxn.ReactionPathway('branch_3way', [T4], [T2])

        self.assertEqual(sorted([path1, path1r, path2, path2r, path3, path4, 
                                 path5, path5r, path6, path6r]), sorted(enum.reactions))


    def test_compare_semantics(self):
        (domains, strands, complexes) = nuskell_parser("""
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
        enum.enumerate()

        path1  = rxn.ReactionPathway('bind21', sorted([PR_FL_B2, helper]), [PR_FLh1B2])
        path1r = rxn.ReactionPathway('open', [PR_FLh1B2], sorted([PR_FL_B2, helper]))
        path2  = rxn.ReactionPathway('bind21', sorted([PR_FL_B2, helper]), [PR_FLh2B2])
        path2r = rxn.ReactionPathway('open', [PR_FLh2B2], sorted([PR_FL_B2, helper]))
        path3  = rxn.ReactionPathway('branch_3way', [PR_FLh1B2], sorted([PR_FL_h1w, B2]))
        path4  = rxn.ReactionPathway('bind21', sorted([PR_FL_B2, B2]), [PR_FLB2B2])
        path4r = rxn.ReactionPathway('open', [PR_FLB2B2], sorted([PR_FL_B2, B2]))

        path5  = rxn.ReactionPathway('branch_3way', [PR_FLh1B2], [PR_FLh2B2])
        path6  = rxn.ReactionPathway('branch_3way', [PR_FLh2B2], [PR_FLh1B2])

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
        #path6  = rxn.ReactionPathway('branch_3way', [PR_FLh2B2], sorted([PR_FL_h1w, B2]))

        #path10  = rxn.ReactionPathway('branch_3way', [PR_FLh2B2], [PR_FLh2B2_v2])
        #path11  = rxn.ReactionPathway('open', [PR_FLh2B2_v2], sorted([PR_FLh2w, B2]))
        #path12  = rxn.ReactionPathway('bind11', [PR_FLh2w], [PR_FL_h1w])
        #path13  = rxn.ReactionPathway('branch_3way', [PR_FLh2B2_v2], [PR_FLh1B2])
        #path14  = rxn.ReactionPathway('branch_3way', [PR_FLh2B2_v2], sorted([PR_FL_h1w, B2]))
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

