#
#  test_output.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 2/6/11.
#

import unittest
from utils import *
from enumerator import *
from nose.tools import *
from input import *
from output import *
import copy

class OutputTests(unittest.TestCase):
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

		self.SLC_enumerator_reduced = Enumerator(self.SLC_enumerator.domains, self.SLC_enumerator.strands, [self.complexes['Cat'], self.complexes['C1'], self.complexes['C2']])
			
		self.three_arm_enumerator = input_standard('test_files/test_input_standard_3arm_junction.in')
		
				
		for complex in self.three_arm_enumerator.initial_complexes:
			self.complexes[complex.name] = complex
				
		self.three_arm_enumerator_reduced = Enumerator(self.three_arm_enumerator.domains, self.three_arm_enumerator.strands, [self.complexes['I'], self.complexes['A'], self.complexes['B'], self.complexes['C']])

	def testOutputLegacy(self):
		self.SLC_enumerator_reduced.enumerate()
		output_legacy(self.SLC_enumerator_reduced, 'test_files/testOutputLegacy.out', output_condensed=False)

	def testOutputLegacy2(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_legacy(self.three_arm_enumerator_reduced, 'test_files/testOutputLegacy2.out', output_condensed=False)
			
	def testOutputFullGraph(self):
		self.SLC_enumerator_reduced.enumerate()
		output_full_graph(self.SLC_enumerator_reduced, 'test_files/testOutputFullGraph.out')
	
	def testOutputFullGraph2(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_full_graph(self.three_arm_enumerator_reduced, 'test_files/testOutputFullGraph2.out')
		
	def testOutputJSON(self):
		self.SLC_enumerator_reduced.enumerate()
		output_json(self.SLC_enumerator_reduced, 'test_files/testOutputJSON.out')
			
	def testJSONInputOutputLoop(self):
		self.SLC_enumerator_reduced.enumerate()
		output_json(self.SLC_enumerator_reduced, 'test_files/testJSONInputOutputLoop.out')
		enumerator = load_json('test_files/testJSONInputOutputLoop.out')
		assert enumerator == self.SLC_enumerator_reduced
			
	def testOutputLegacyCondensed(self):
		return True
		self.SLC_enumerator_reduced.enumerate()
		output_legacy(self.SLC_enumerator_reduced, 'test_files/testOutputLegacyCondensed.out', output_condensed=True)
		assert False