#
#  test_enumerator.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 2/6/11.
#

import unittest
from utils import *
from enumerator import *
from nose.tools import *
from input import input_standard
import copy

class EnumeratorTests(unittest.TestCase):
	def setUp(self):
		self.SLC_enumerator = input_standard('test_files/test_input_standard_SLC.in')
		self.domains = {}
		self.strands = {}
		self.complexes = {}		
		
		for domain in self.SLC_enumerator.domains:
			self.domains[domain.name] = domain
		
		for strand in self.SLC_enumerator.strands:
			self.strands[strand.name] = strand
		
		for complex in self.SLC_enumerator.initial_complexes:
			self.complexes[complex.name] = complex
			
	def testDomains(self):
		exp_domains = [Domain('1', 'short'), Domain('1', 'short', True), Domain('2', 'short'), Domain('2', 'short', True), Domain('3', 'short'), Domain('3', 'short', True), Domain('4', 'long'), Domain('4', 'long', True), Domain('5', 'short'), Domain('5', 'short', True), Domain('6', 'long'), Domain('6', 'long', True), Domain('7', 'short'), Domain('7', 'short', True)]
		act_domains = self.SLC_enumerator.domains
		
		exp_domains.sort()
		act_domains.sort()
		
		
		assert exp_domains == act_domains
		act_domains[0] = Domain('99', 'short')
		
		assert exp_domains != act_domains
		
		def assnDomains(self):
			self.SLC_enumerator.domains = []
			
		assert_raises(AttributeError, assnDomains, self)
		
	def testStrands(self):
		exp_strands = [Strand('PS', [self.domains['3*'], self.domains['2*'], self.domains['1*'], self.domains['5'], self.domains['6']]), Strand('OP', [self.domains['1'], self.domains['2'], self.domains['3'], self.domains['4']]), Strand('SP', [self.domains['5'], self.domains['6']]), Strand('BS', [self.domains['7*'], self.domains['6*'], self.domains['5*'], self.domains['1'], self.domains['2'], self.domains['3']]), Strand('Cat', [self.domains['6'], self.domains['7']])]
		act_strands = self.SLC_enumerator.strands
		
		exp_strands.sort()
		act_strands.sort()
		
		assert exp_strands == act_strands
		
		act_strands[0] = Strand('A', [self.domains['1']])
		
		assert exp_strands != act_strands
		
		def assnStrands(self):
			self.SLC_enumerator.domains = []
			
		assert_raises(AttributeError, assnStrands, self)
		
	def testInitialComplexes(self):
		PS = self.strands['PS']
		OP = self.strands['OP']
		Cat = self.strands['Cat']
		SP = self.strands['SP']
		BS = self.strands['BS']
	
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

		exp_complexes = [C1, C2, CatC, OPC, SPC, I1, I2, I3, I4, I5, I6, W]
		act_complexes = self.SLC_enumerator.initial_complexes
		
		exp_complexes.sort()
		act_complexes.sort()
		
		assert exp_complexes == act_complexes
		
		act_complexes[0] = Complex('A', [Cat], [[None, None]])
		
		assert exp_complexes != act_complexes
		
		def assnComplexes(self):
			self.SLC_enumerator.initial_complexes = []
			
		assert_raises(AttributeError, assnComplexes, self)
		

	def testProcessNeighborhood1(self):
		enum = Enumerator(self.SLC_enumerator._domains, self.SLC_enumerator._strands, [self.complexes['Cat']])
		enum._N = []
		enum._S = []
		enum._T = []
		enum._resting_states = []
		enum._reactions = []
		enum._complexes = []
		enum.process_neighborhood(self.complexes['Cat'])
		assert enum._S == [self.complexes['Cat']]

	def testEnumeration1(self):
		enum = input_standard('test_files/test_input_standard_simple.in')
		strands = {}
		for strand in enum.strands:
			strands[strand.name] = strand
		
		expected_complexes = []
		for complex in enum.initial_complexes:
			expected_complexes.append(complex)
			
		new_complex = Complex('C3', [strands['S1'], strands['S2']], [[(1, 0), None],[(0, 1)]])
		expected_complexes.append(new_complex)
		enum.enumerate()
		resting_complexes = enum.resting_complexes
		
		expected_complexes.sort()
		resting_complexes.sort()
		assert expected_complexes == resting_complexes
		
