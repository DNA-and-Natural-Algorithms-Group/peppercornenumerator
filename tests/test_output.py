import unittest

from peppercornenumerator import enumerate_pil
from peppercornenumerator.output import write_sbml
from peppercornenumerator.objects import clear_memory

class Test_SBML_output(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def test_SCL_system(self):
        SCL_input = """
        # This file describes the suppressed-leak catalyst system
        # by DY Zhang and K Sarma
        
        # Domains (13)
        length d1 = 5
        length d2 = 5
        length d3 = 5
        length d4 = 15
        length d5 = 5
        length d6 = 15
        length d7 = 5
        
        # Strands or composite domains (10)
        sup-sequence PS = d3* d2* d1* d5 d6 : 35
        sup-sequence SP = d5 d6 : 20
        sup-sequence Cat = d6 d7 : 20
        sup-sequence BS = d7* d6* d5* d1 d2 d3 : 40
        sup-sequence OP = d1 d2 d3 d4 : 30
        
        # Resting complexes (7)
        C1 = d3*( d2*( d1*( d5 d6 + ) ) ) d4
        C2 = d5( d6( + d7* ) ) d1 d2 d3
        Cat = d6 d7
        I3 = d7*( d6*( d5* d1 d2 d3 + ) )
        OP = d1 d2 d3 d4
        SP = d5 d6
        W = d7* d6*( d5*( d1( d2( d3( + ) ) ) ) )
        
        # Transient complexes (7)
        new1 = d7*( d6*( d5*( d1 d2 d3 + d1( d2( d3( d4 + ) ) ) ) ) + d6 )
        new2 = d7* d6*( d5*( d1 d2 d3 + d1( d2( d3( d4 + ) ) ) ) )
        I1 = d5( d6( + d6 d7( + ) ) ) d1 d2 d3
        I2 = d5( d6 + d6( d7( + ) ) ) d1 d2 d3
        I4 = d7*( d6*( d5*( d1 d2 d3 + d1( d2( d3( d4 + ) ) ) ) d6 + ) )
        I5 = d7*( d6*( d5*( d1( d2( d3( + ) ) ) ) d6 + ) )
        I6 = d7*( d6*( d5*( d1( d2( d3( + ) ) ) ) ) + d6 )
        """

        enum, out = enumerate_pil(SCL_input, is_file = False)
        #print(out)
        print(write_sbml(enum))
        print()
        print(write_sbml(enum, condensed = True))

    def test_SCL_system(self):
        ARM3J = """
        # Domains (12)
        length a = 6
        length b = 6
        length c = 6
        length x = 6
        length y = 6
        length z = 6
        
        # Strands or composite domains (8)
        sup-sequence A = a x b y z* c* y* b* x* : 54
        sup-sequence C = c z a x y* b* x* a* z* : 54
        sup-sequence B = b y c z x* a* z* c* y* : 54
        sup-sequence I = y* b* x* a* : 24
        
        # Complexes (8)
        A = a x( b( y( z* c* ) ) ) @i 1e-7M
        B = b y( c( z( x* a* ) ) ) @i 1e-7M
        C = c z( a( x( y* b* ) ) ) @i 1e-7M
        I = y* b* x* a* @i 1e-9M
        IA = y*( b*( x*( a*( + ) ) ) ) z* c* y* b* x* @i 0 M
        IAB = y*( b*( x*( a*( + ) ) ) ) z*( c*( y*( b*( x* + ) ) ) ) x* a* z* c* y* @i 0 M
        IABC = y*( b*( x*( a*( + ) ) ) ) z*( c*( y*( b*( x* + ) ) ) ) x*( a*( z*( c*( y* + ) ) ) ) y* b* x* a* z* @i 0 M
        ABC = a( x( b( y( z*( c*( y*( b*( x* + ) ) ) ) x*( a*( z*( c*( y* + ) ) ) ) ) ) ) ) z* @i 0 M
        """

        enum, out = enumerate_pil(ARM3J, is_file = False, k_fast = 100, k_slow = 20)
        print(out)
        print(write_sbml(enum))
        print()
        print(write_sbml(enum, condensed = True))


 

if __name__ == '__main__':
  unittest.main()


