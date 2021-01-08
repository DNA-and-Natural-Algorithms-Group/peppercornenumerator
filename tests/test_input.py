#
# tests/test_input.py
#
import logging
logger = logging.getLogger('peppercornenumerator')
logger.setLevel(logging.INFO)
import unittest

from dsdobjects import SingletonError, clear_singletons
from peppercornenumerator.input import read_pil
from peppercornenumerator.objects import PepperDomain, PepperComplex, PepperMacrostate, PepperReaction

SKIP = False

@unittest.skipIf(SKIP, "skipping tests.")
class TestPilIO(unittest.TestCase):
    def tearDown(self):
        clear_singletons(PepperDomain)
        clear_singletons(PepperComplex)
        clear_singletons(PepperMacrostate)
        clear_singletons(PepperReaction)
 
    def test_read_pil_01(self):
        complexes, reactions = read_pil(
        """
        # Domains (12) 
        sequence a = NNNNNN
        sequence b = NNNNNN 
        sequence c = NNNNNN
        sequence x = NNNNNN
        sequence y = NNNNNN 
        sequence z = NNNNNN 

        # Resting complexes (8) 
        A = a x( b( y( z* c* ) ) ) @i 1e-08 M
        B = b y( c( z( x* a* ) ) ) @i 1e-08 M
        C = c z( a( x( y* b* ) ) ) @i 1e-08 M
        I = y* b* x* a* @i 1e-08 M

        IA = a( x( b( y( z* c* y* b* x* + ) ) ) ) @i 0.0 M
        IAB = y*( b*( x*( a*( + ) ) ) ) z*( c*( y*( b*( x* + ) ) ) ) x* a* z* c* y* @i 0.0 M
        IABC = y*( b*( x*( a*( + ) ) ) ) z*( c*( y*( b*( x* + ) ) ) ) x*( a*( z*( c*( y* + ) ) ) ) y* b* x* a* z* @i 0.0 M
        ABC = a( x( b( y( z*( c*( y*( b*( x* + ) ) ) ) x*( a*( z*( c*( y* + ) ) ) ) ) ) ) ) z* @i 0.0 M

        # Resting macrostates (8) 
        macrostate A = [A]
        macrostate B = [B]
        macrostate C = [C]
        macrostate I = [I]
        macrostate IA = [IA]
        macrostate IAB = [IAB]
        macrostate IABC = [IABC]
        macrostate ABC = [ABC]

        # Condensed reactions (10) 
        reaction [condensed      =  1.66666e+06 /M/s ] A + I -> IA
        reaction [condensed      =  1.66666e+06 /M/s ] IA + B -> IAB
        reaction [condensed      =  1.66646e+06 /M/s ] IAB + C -> IABC
        reaction [condensed      =    0.0261637 /s   ] IABC -> ABC + I
        """)

        # A preliminary interface to start testing prototype functions.
        #a, b, c = dom['a'], dom['b'], dom['c']
        #assert a is PepperDomain('a')
        #assert b is PepperDomain('b')
        #assert c is PepperDomain('c')
        assert PepperComplex(None, None, 'IA').kernel_string == 'a( x( b( y( z* c* y* b* x* + ) ) ) )'
        #assert PepperMacrostate(None, 'A').kernel_string == 'a x( b( y( z* c* ) ) )'
        #PepperReaction(None, None, None, name = '[condensed] A + I -> IA')
