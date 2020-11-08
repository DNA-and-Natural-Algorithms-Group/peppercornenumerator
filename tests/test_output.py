import unittest
import warnings

try:
    from roadrunner import RoadRunner
    SKIP_ROADRUNNER = False
except ImportError as err:
    warnings.warn('Testing the SBML output requires libroadrunner')
    SKIP_ROADRUNNER = True

from peppercornenumerator import enumerate_pil
from peppercornenumerator.output import write_sbml, write_crn
from peppercornenumerator.objects import clear_memory

SKIP = False

@unittest.skipIf(SKIP or SKIP_ROADRUNNER, "skipping tests that require libroadrunner")
class Test_SBML_roadrunner(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_delayed_choice(self):
        inpil = """
        length a = 5
        length b = 15
        length c = 15
        length d = 15
        length e = 15
        length z = 15
        
        RC1 = d( e( z ) ) b*( a* + ) c( + ) d @i 1e-8 M
        RC1b = d e( z ) d*( b*( a* + ) c( + ) ) @i 1e-8 M
        R0 = a b c @i 1e-7 M
        
        R21 = b c( + ) d @i 0 M
        R22 = d( e( z ) ) b*( a*( + ) ) c @i 0 M
        
        R31 = b c @i 0 M
        RC32b = d e( z ) d*( b*( a*( + ) ) c( + ) ) @i 0 M
        RC32 = d( e( z ) ) b*( a*( + ) ) c( + ) d @i 0 M
        
        TC1 = a( b c + b( c( + ) d + d( e( z ) ) ) ) @i 0 M
        TC2 = a( b c + b( c( + ) d( + d e( z ) ) ) ) @i 0 M
        TC3 = a( b( c + b c( + ) d( + d e( z ) ) ) ) @i 0 M
        """

        enum, out = enumerate_pil(inpil, is_file = False, condensed = True)
        xml = write_sbml(enum, condensed = True) 

        rr = RoadRunner()
        rr.load(xml)
        ExecutableModel = rr.model

        RC32 = 'RC32' if 'RC32' in rr.model.getFloatingSpeciesIds() else 'RC32b'
        RC1 = 'RC1' if 'RC1' in rr.model.getFloatingSpeciesIds() else 'RC1b'

        # Compartment
        assert ExecutableModel.getCompartmentVolumes() == 1.
        assert ExecutableModel.getNumCompartments() == 1
        assert ExecutableModel.getCompartmentIds() == ['TestTube']

        # Species
        speciesIDs = ['R0', RC32, 'R22', 'R31', 'R21', RC1]
        ini_spIDs = sorted('init({})'.format(x) for x in speciesIDs)
        concIDs = ['[R0]', '[' + RC32 + ']', '[R22]', '[R31]', '[R21]', '[' + RC1 + ']']
        ini_coIDs = sorted('init({})'.format(x) for x in concIDs)
        assert ExecutableModel.getNumFloatingSpecies() == len(speciesIDs)
        assert sorted(ExecutableModel.getFloatingSpeciesIds()) == sorted(speciesIDs)
        assert sorted(ExecutableModel.getFloatingSpeciesInitAmountIds()) == sorted(ini_spIDs)
        assert sorted(ExecutableModel.getFloatingSpeciesInitConcentrationIds()) == sorted(ini_coIDs)
        #print(sorted(ini_coIDs))
        #print(ExecutableModel.getFloatingSpeciesInitAmounts())
        #print(ExecutableModel.getFloatingSpeciesInitConcentrations())

        # Reactions
        rxnIDs = ['R0_' + RC1 + '__' + RC32 + '_R31', 'R0_' + RC1 + '__R22_R21']
        assert ExecutableModel.getNumReactions() == len(rxnIDs)
        assert sorted(ExecutableModel.getReactionIds()) == sorted(rxnIDs)
        assert rr.model[RC32] == 0
        assert rr.model[RC1] == 2e-8

        # simulate deterministic
        rr.model['init([R0])'] = 1e-8
        assert rr.model['init(R0)'] == 1e-8

        Vol = 1.66e-15
        rr.model.setCompartmentVolumes([Vol])
        rr.model['init([R0])'] = 1e-8 * Vol
        rr.model[f'init({RC1})'] = 2e-8 * Vol
        rr.integrator.absolute_tolerance = 1e-12 * Vol
        rr.integrator.relative_tolerance = 1e-12 * Vol
        rr.integrator.initial_time_step = 0.00001
        result = rr.simulate(0, 500, steps=100)
        #print(result)
        #rr.plot() # look at it if you like!

        # NOTE: something is off with the units with stochastic simulations, weird...
        #print(rr.model.getCompartmentVolumes())
        #print(rr.model.getFloatingSpeciesInitConcentrations())
        #print(rr.model.getFloatingSpeciesInitAmounts())
        #print(rr.model.getFloatingSpeciesConcentrations())
        #print(rr.model.getFloatingSpeciesAmounts())
        rr.reset()

@unittest.skipIf(SKIP, "skipping tests")
class Test_SBML_output(unittest.TestCase):
    def tearDown(self):
        clear_memory()

    def test_SCL_system_01(self):
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
        sup-sequence sPS = d3* d2* d1* d5 d6 : 35
        sup-sequence sSP = d5 d6 : 20
        sup-sequence sCat = d6 d7 : 20
        sup-sequence sBS = d7* d6* d5* d1 d2 d3 : 40
        sup-sequence sOP = d1 d2 d3 d4 : 30
        
        # Resting complexes (7)
        U = d1 d2 d3 @initial 0 M
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
        assert 'length' in out
        assert 'sup-sequence' in out
        assert write_crn(enum)
        assert write_crn(enum, condensed = True)
        assert write_sbml(enum)
        assert write_sbml(enum, condensed = True)

    def test_SCL_system_02(self):
        ARM3J = """
        # Domains (12)
        length a = 6
        length b = 6
        length c = 6
        length x = 6
        length y = 6
        length z = 6
        
        # Strands or composite domains (8)
        sup-sequence sA = a x b y z* c* y* b* x* : 54
        sup-sequence sC = c z a x y* b* x* a* z* : 54
        sup-sequence sB = b y c z x* a* z* c* y* : 54
        sup-sequence sI = y* b* x* a* : 24
        
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
        assert 'length' in out
        assert 'sup-sequence' in out
        assert write_crn(enum)
        assert write_crn(enum, condensed = True)
        assert write_sbml(enum)
        assert write_sbml(enum, condensed = True)

if __name__ == '__main__':
  unittest.main()

