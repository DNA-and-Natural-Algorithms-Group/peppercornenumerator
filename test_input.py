#
#  test_input.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 1/8/11.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#

import unittest
from utils import *
from input import *
from nose.tools import *
import copy

class InputStandardTests(unittest.TestCase):
	def testStandard_SLC(self):
		enumerator = input_standard('test_files/test_input_standard_SLC.in')
		
		d1 = Domain('1', 'short')
		d1c = Domain('1', 'short', True)
		d2 = Domain('2', 'short')
		d2c = Domain('2', 'short', True)
		d3 = Domain('3', 'short')
		d3c = Domain('3', 'short', True)
		d4 = Domain('4', 'long')
		d4c = Domain('4', 'long', True)
		d5 = Domain('5', 'short')
		d5c = Domain('5', 'short', True)
		d6 = Domain('6', 'long')
		d6c = Domain('6', 'long', True)
		d7 = Domain('7', 'short')
		d7c = Domain('7', 'short', True)
		
		enum_doms = enumerator.domains[:]
		doms = [d1, d1c, d2, d2c, d3, d3c, d4, d4c, d5, d5c, d6, d6c, d7, d7c]
		
		enum_doms.sort()
		doms.sort()
		
		assert (enum_doms == doms)
		
		PS = Strand('PS', [d3c, d2c, d1c, d5, d6])
		OP = Strand('OP', [d1, d2, d3, d4])
		SP = Strand('SP', [d5, d6])
		BS = Strand('BS', [d7c, d6c, d5c, d1, d2, d3])
		Cat = Strand('Cat', [d6, d7])
		
		enum_strands = enumerator.strands[:]
		strands = [PS, OP, SP, BS, Cat]
		
		enum_strands.sort()
		strands.sort()
		
		assert (enum_strands == strands)
		
		C1 = Complex('C1', [PS, OP], [[(1, 2), (1, 1), (1, 0), None, None], [(0, 2), (0, 1), (0, 0), None]])
		C2 = Complex('C2', [SP, BS], [[(1, 2), (1, 1)], [None, (0, 1), (0, 0), None, None, None]])

		CatC = Complex('Cat', [Cat], [[None, None]])
		OPC = Complex('OP', [OP], [[None, None, None, None]])
		SPC = Complex('SP', [SP], [[None, None]])
		
		I1 = Complex('I1', [SP, Cat, BS], [[(2, 2), (2, 1)], [None, (2, 0)], [(1, 1), (0, 1), (0, 0), None, None, None]])
		I2 = Complex('I2', [SP, Cat, BS], [[(2, 2), None], [(2, 1), (2, 0)], [(1, 1), (1, 0), (0, 0), None, None, None]])
		I3 = Complex('I3', [BS, Cat], [[(1, 1), (1, 0), None, None, None, None], [(0, 1), (0, 0)]])
		I4 = Complex('I4', [BS, OP, PS, Cat], [[(3,1), (3, 0), (2, 3), None, None, None],[(2, 2), (2, 1), (2, 0), None],[(1, 2), (1, 1), (1, 0), (0, 2), None],[(0, 1), (0, 0)]])
		I5 = Complex('I5', [BS, PS, Cat], [[(2, 1), (2, 0), (1, 3), (1, 2), (1, 1), (1, 0)], [(0, 5), (0, 4), (0, 3), (0, 2), None], [(0, 1), (0, 0)]])
		I6 = Complex('I6', [BS, PS, Cat], [[(2, 1), (1, 4), (1, 3), (1, 2), (1, 1), (1, 0)], [(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)], [None, (0, 0)]])
		W = Complex('W', [BS, PS], [[None, (1, 4), (1, 3), (1, 2), (1, 1), (1, 0)], [(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)]])
		
		enum_complexes = enumerator.initial_complexes[:]
		complexes = [C1, C2, CatC, OPC, SPC, I1, I2, I3, I4, I5, I6, W]
		
		
		enum_complexes.sort()
		complexes.sort()
		
		assert (enum_complexes == complexes)

	def testStandard_3arm_junction(self):
		enumerator = input_standard('test_files/test_input_standard_3arm_junction.in')
		
		da = Domain('a', 6)		
		db = Domain('b', 6)
		dc = Domain('c', 6)
		dx = Domain('x', 6)
		dy = Domain('y', 6)
		dz = Domain('z', 6)
		dac = Domain('a', 6, True)
		dbc = Domain('b', 6, True)
		dcc = Domain('c', 6, True)
		dxc = Domain('x', 6, True)
		dyc = Domain('y', 6, True)
		dzc = Domain('z', 6, True)
		
		enum_doms = enumerator.domains[:]
		doms = [da, dac, db, dbc, dc, dcc, dx, dxc, dy, dyc, dz, dzc]
				
		enum_doms.sort()
		doms.sort()
		
		assert (enum_doms == doms)
		
		I = Strand('I', [dyc, dbc, dxc, dac])
		A = Strand('A', [da, dx, db, dy, dzc, dcc, dyc, dbc, dxc])
		B = Strand('B', [db, dy, dc, dz, dxc, dac, dzc, dcc, dyc])
		C = Strand('C', [dc, dz, da, dx, dyc, dbc, dxc, dac, dzc])
		
		enum_strands = enumerator.strands[:]
		strands = [I, A, B, C]
		
		enum_strands.sort()
		strands.sort()
		
		assert (enum_strands == strands)
		
		IC = Complex('I', [I], [[None, None, None, None]])
		AC = Complex('A', [A], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		BC = Complex('B', [B], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		CC = Complex('C', [C], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		
		IA = Complex('IA', [I, A], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), None, None, None, None, None]])
		IAB = Complex('IAB', [I, A, B], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), None, None, None, None, None]])
		IABC = Complex('IABC', [I, A, B, C], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (3, 3), (3, 2), (3, 1), (3, 0), None], [(2, 7), (2, 6), (2, 5), (2, 4), None, None, None, None, None]])
		ABC = Complex('ABC', [A, B, C], [[(2, 7), (2, 6), (2, 5), (2, 4), (1, 3), (1, 2), (1, 1), (1, 0), None], [(0, 7), (0, 6), (0, 5), (0, 4), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (0, 3), (0, 2), (0, 1), (0, 0), None]])
		
		enum_complexes = enumerator.initial_complexes[:]
		complexes = [IC, AC, BC, CC, IA, IAB, IABC, ABC]
		
		
		enum_complexes.sort()
		complexes.sort()
		
		
		assert (enum_complexes == complexes)
	
	def testPil_3arm_junction(self):
		enumerator = input_pil('test_files/test_input_pil_3arm_junction.pil')
		
		da = Domain('a', 6, sequence="CTACTC")		
		db = Domain('b', 6, sequence="TCCTCA")
		dc = Domain('c', 6, sequence="TTTCCA")
		dx = Domain('x', 6, sequence="CTACTC")
		dy = Domain('y', 6, sequence="TCCTCA")
		dz = Domain('z', 6, sequence="TTTCCA")
		dac = Domain('a', 6, True, sequence="CTACTC")
		dbc = Domain('b', 6, True, sequence="TCCTCA")
		dcc = Domain('c', 6, True, sequence="TTTCCA")
		dxc = Domain('x', 6, True, sequence="CTACTC")
		dyc = Domain('y', 6, True, sequence="TCCTCA")
		dzc = Domain('z', 6, True, sequence="TTTCCA")
		
		enum_doms = enumerator.domains[:]
		doms = [da, dac, db, dbc, dc, dcc, dx, dxc, dy, dyc, dz, dzc]
				
		enum_doms.sort()
		doms.sort()
		
		assert (enum_doms == doms)
		
		I = Strand('I', [dyc, dbc, dxc, dac])
		A = Strand('A', [da, dx, db, dy, dzc, dcc, dyc, dbc, dxc])
		B = Strand('B', [db, dy, dc, dz, dxc, dac, dzc, dcc, dyc])
		C = Strand('C', [dc, dz, da, dx, dyc, dbc, dxc, dac, dzc])
		
		enum_strands = enumerator.strands[:]
		strands = [I, A, B, C]
		
		enum_strands.sort()
		strands.sort()
		
		assert (enum_strands == strands)
		
		IC = Complex('I', [I], [[None, None, None, None]])
		AC = Complex('A', [A], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		BC = Complex('B', [B], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		CC = Complex('C', [C], [[None, (0, 8), (0, 7), (0, 6), None, None, (0, 3), (0, 2), (0, 1)]])
		
		IA = Complex('IA', [I, A], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), None, None, None, None, None]])
		IAB = Complex('IAB', [I, A, B], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), None, None, None, None, None]])
		IABC = Complex('IABC', [I, A, B, C], [[(1, 3), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (0, 0), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (3, 3), (3, 2), (3, 1), (3, 0), None], [(2, 7), (2, 6), (2, 5), (2, 4), None, None, None, None, None]])
		ABC = Complex('ABC', [A, B, C], [[(2, 7), (2, 6), (2, 5), (2, 4), (1, 3), (1, 2), (1, 1), (1, 0), None], [(0, 7), (0, 6), (0, 5), (0, 4), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (0, 3), (0, 2), (0, 1), (0, 0), None]])
		
		enum_complexes = enumerator.initial_complexes[:]
		complexes = [IC, AC, BC, CC, IA, IAB, IABC, ABC]
		
		
		enum_complexes.sort()
		complexes.sort()
		
		
		assert (enum_complexes == complexes)
	
	
	def testErrors(self):
		
		def testDuplicateDomain():
			enum = input_standard('test_files/test_input_errors/test_input_duplicate_domain.in')
		assert_raises(Exception, testDuplicateDomain)
		
		def testDuplicateStrand():
			enum = input_standard('test_files/test_input_errors/test_input_duplicate_strand.in')	
		assert_raises(Exception, testDuplicateStrand)
		
		def testDuplicateComplex():
			enum = input_standard('test_files/test_input_errors/test_input_duplicate_complex.in')	
		assert_raises(Exception, testDuplicateComplex)

		def testMissingDomain():
			enum = input_standard('test_files/test_input_errors/test_input_missing_domain.in')	
		assert_raises(Exception, testMissingDomain)

		def testMissingStrand():
			enum = input_standard('test_files/test_input_errors/test_input_missing_strand.in')	
		assert_raises(Exception, testMissingStrand)
		
		def testComplexError():
			enum = input_standard('test_files/test_input_errors/test_input_size_mismatch.in')
		assert_raises(Exception, testComplexError)
		
		enum = input_standard('test_files/test_input_errors/test_input_warnings.in')	
		assert (enum != None)
			
		
			