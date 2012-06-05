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
	# Disable output tests until all other functionality is working!
	
	# __test__ = False
	
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

		self.simple_enumerator = input_standard('test_files/test_input_standard_simple.in')
				
		self.simple2_enumerator = input_standard('test_files/test_input_standard_simple2.in')

	def testOutputLegacy(self):
		self.SLC_enumerator_reduced.enumerate()
		output_legacy(self.SLC_enumerator_reduced, 'test_files/testOutputLegacy.out', output_condensed=False)

	def testOutputLegacy2(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_legacy(self.three_arm_enumerator_reduced, 'test_files/testOutputLegacy2.out', output_condensed=False)
		
	def testOutputLegacy3(self):
		self.simple2_enumerator.enumerate()
		output_legacy(self.simple2_enumerator, 'test_files/testOutputLegacy3.out', output_condensed=False)
	
	def testOutputFullGraph(self):
		self.SLC_enumerator_reduced.enumerate()
		output_full_graph(self.SLC_enumerator_reduced, 'test_files/testOutputFullGraph.out')
	
	def testOutputFullGraph2(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_full_graph(self.three_arm_enumerator_reduced, 'test_files/testOutputFullGraph2.out')
	
	def testOutputFullGraph3(self):
		self.simple2_enumerator.enumerate()
		output_full_graph(self.simple2_enumerator, 'test_files/testOutputFullGraph3.out')
		
	def testOutputCondensedGraph(self):
		self.SLC_enumerator_reduced.enumerate()
		output_condensed_graph(self.SLC_enumerator_reduced, 'test_files/testOutputCondensedGraph.out')
	
	def testOutputCondensedGraph2(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_condensed_graph(self.three_arm_enumerator_reduced, 'test_files/testOutputCondensedGraph2.out')
	
	def testOutputCondensedGraph3(self):
		self.simple2_enumerator.enumerate()
		output_condensed_graph(self.simple2_enumerator, 'test_files/testOutputCondensedGraph3.out')
		
	def testOutputGraph(self):
		self.simple2_enumerator.enumerate()
		output_graph(self.simple2_enumerator,'test_files/testOutputGraph.out',output_condensed=True)
		output_graph(self.simple2_enumerator,'test_files/testOutputGraphCondensed.out',output_condensed=False)
		
		
	def testOutputJSON(self):
		self.SLC_enumerator_reduced.enumerate()
		output_json(self.SLC_enumerator_reduced, 'test_files/testOutputJSON.out')
		enumerator = load_json('test_files/testOutputJSON.out')
		assert enumerator == self.SLC_enumerator_reduced
		
	def testOutputJSON2(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_json(self.three_arm_enumerator_reduced, 'test_files/testOutputJSON2.out')
		enumerator = load_json('test_files/testOutputJSON2.out')
		assert enumerator == self.three_arm_enumerator_reduced
		
	def testOutputJSONCondensed(self):
		self.SLC_enumerator_reduced.enumerate()
		output_json(self.SLC_enumerator_reduced, \
				'test_files/testOutputJSONCondensed.out',output_condensed=True)
		enumerator = load_json('test_files/testOutputJSON.out')
		assert enumerator == self.SLC_enumerator_reduced
		
	def testOutputJSONCondensed2(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_json(self.three_arm_enumerator_reduced, \
				'test_files/testOutputJSONCondensed2.out',output_condensed=True)
		enumerator = load_json('test_files/testOutputJSON2.out')
		assert enumerator == self.three_arm_enumerator_reduced
			
	def testJSONInputOutputLoop(self):
		self.SLC_enumerator_reduced.enumerate()
		output_json(self.SLC_enumerator_reduced, 'test_files/testJSONInputOutputLoop.out')
		enumerator = load_json('test_files/testJSONInputOutputLoop.out')
		assert enumerator == self.SLC_enumerator_reduced

	def testOutputSBML(self):
		self.SLC_enumerator_reduced.enumerate()
		output_sbml(self.SLC_enumerator_reduced, 'test_files/testOutputSBML')
		
	def testOutputSBML2(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_sbml(self.three_arm_enumerator_reduced, 'test_files/testOutputSBML2')
		
	def testOutputSBMLCondensed(self):
		self.SLC_enumerator_reduced.enumerate()
		output_sbml(self.SLC_enumerator_reduced, 'test_files/testOutputSBMLCondensed',output_condensed=True)
		
	def testOutputSBML2Condensed(self):
		self.three_arm_enumerator_reduced.enumerate()
		output_sbml(self.three_arm_enumerator_reduced, 'test_files/testOutputSBML2Condensed',output_condensed=True)
			
	def testOutputLegacyCondensed(self):
		self.simple_enumerator.enumerate()
		output_legacy(self.simple_enumerator, 'test_files/testOutputLegacyCondensed.out', output_condensed=True)
		
		
	def testCondenseRestingStates(self):
		for enum in [self.SLC_enumerator, self.three_arm_enumerator]:
			enum.enumerate()
			condensed = condense_resting_states(enum)
				
			assert (sorted(condensed['resting_states']) == \
					sorted(enum.resting_states))
			
			# Test that transient complexes have been eliminated from complexes
			resting_complexes = set([complex for state in condensed['resting_states'] for complex in state.complexes])
			predicted_resting_complexes = set(enum.complexes) - set(enum.transient_complexes)
			
			print resting_complexes
			print predicted_resting_complexes
			
			assert (resting_complexes == predicted_resting_complexes)
			
			# Test transient complexes have been removed from reactions
			transients = set(enum.transient_complexes)
			for reaction in condensed['reactions']:
				for reactant in reaction.reactants:
					assert reactant not in transients
					
				for product in reaction.products:
					assert product not in transients
					
			# Test that no reactions are duplicated
			assert (sorted(condensed['reactions']) == sorted(list(set(condensed['reactions']))))
					
			
				
				
		
