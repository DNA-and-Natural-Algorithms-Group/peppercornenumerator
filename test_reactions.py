#
#  test_reactions.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 5/19/10.
#

from utils import *
import reactions
from reactions import *
from input import input_standard

import unittest
from nose.tools import *
import copy

class BindTests(unittest.TestCase):
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
		
		
	def testFindExternalStrandBreak(self):		
		I4 = self.complexes['I4']
		
		for (i, strand) in enumerate(I4.strands):
			if (strand.name == 'BS'):
				BS_index = i
			elif (strand.name == 'PS'):
				PS_index = i
		
		assert PS_index == find_external_strand_break(I4, (PS_index, 4))
		
		assert BS_index == find_external_strand_break(I4, (BS_index, 3))
		
		I1 = self.complexes['I1']
		
		
		for (i, strand) in enumerate(I1.strands):
			if (strand.name == 'SP'):
				SP_index = i
			elif (strand.name == 'Cat'):
				Cat_index = i
			elif (strand.name == 'BS'):
				BS_index = i
		
		assert SP_index == find_external_strand_break(I1, (Cat_index, 0))
		
		assert BS_index == find_external_strand_break(I1, (BS_index, 4))
		
	def testBind11(self):
		strand = Strand('A', [self.domains['1'], self.domains['2'], self.domains['1*']])
		complex = Complex('C', [strand], [[None, None, None]])
		
		out_list = bind11(complex)
		
		exp_complex = Complex('C-New', [strand], [[(0, 2), None, (0, 0)]])
		exp_out = [ReactionPathway('bind11', [complex], [exp_complex])]
		
		out_list.sort()
		exp_out.sort()
		
		assert out_list == exp_out
		
	def test_combine_complexes_21(self):
		PS_index = self.complexes['C1'].strand_index('PS')
		BS_index = self.complexes['I3'].strand_index('BS')
		
		out_complex = combine_complexes_21(self.complexes['C1'], (PS_index, 3),
										   self.complexes['I3'], (BS_index, 2))
										   
		print out_complex
		print out_complex.strands
		print out_complex.structure
		
		exp_complex = Complex('I4', [self.strands['PS'], self.strands['Cat'], self.strands['BS'], self.strands['OP']], [[(3, 2), (3, 1), (3, 0), (2, 2), None], [(2, 1), (2, 0)], [(1, 1), (1, 0), (0, 3), None, None, None], [(0, 2), (0, 1), (0, 0), None]])

		print exp_complex
		print exp_complex.strands
		print exp_complex.structure

		assert exp_complex == out_complex
		
	def testBind21(self):
		out_list = bind21(self.complexes['C1'], self.complexes['I3'])
		
		exp_complex = Complex('I4', [self.strands['PS'], self.strands['Cat'], self.strands['BS'], self.strands['OP']], [[(3, 2), (3, 1), (3, 0), (2, 2), None], [(2, 1), (2, 0)], [(1, 1), (1, 0), (0, 3), None, None, None], [(0, 2), (0, 1), (0, 0), None]])

		exp_pathway = ReactionPathway('bind21', [self.complexes['C1'], self.complexes['I3']], [exp_complex])
		
		print out_list
		
		assert out_list == [exp_pathway]
		
class OpenTests(unittest.TestCase):
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
	
	def testSplitComplex1(self):
		complex = self.complexes['I3'].clone()
		cat_index = complex.strand_index('Cat')
		bs_index = complex.strand_index('BS')
		
		complex.structure[cat_index][0] = None
		complex.structure[cat_index][1] = None
		complex.structure[bs_index][0] = None
		complex.structure[bs_index][1] = None
		
		res_list = split_complex(complex, (bs_index, -1), (bs_index, -1))
		bs_complex = Complex('BS', [self.strands['BS']], [[None, None, None, None, None, None]])
		
		exp_list = [bs_complex, self.complexes['Cat']]
				
					
		res_list.sort()
		exp_list.sort()
				
		assert exp_list == res_list
		
	def testSplitComplex2(self):
		complex = self.complexes['I4'].clone()
				
		BS_index = complex.strand_index('BS')
		PS_index = complex.strand_index('PS')
		
		complex.structure[BS_index][2] = None
		complex.structure[PS_index][3] = None
		
		res_list = split_complex(complex, (1, -1), (2, 5))
		
		exp_list = [self.complexes['C1'], self.complexes['I3']]
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list

		
	def testFindReleases1(self):
		# first test no releases case
		complex = self.complexes['I3'].clone()
		res_list = find_releases(complex)
		assert res_list == [complex]
			

	def testFindReleases2(self):			
		# now test 1 release case
		complex = self.complexes['I3'].clone()
		cat_index = complex.strand_index('Cat')
		bs_index = complex.strand_index('BS')
		
		complex.structure[cat_index][0] = None
		complex.structure[cat_index][1] = None
		complex.structure[bs_index][0] = None
		complex.structure[bs_index][1] = None
		
		res_list = find_releases(complex)
		bs_complex = Complex('BS', [self.strands['BS']], [[None, None, None, None, None, None]])
		
		exp_list = [bs_complex, self.complexes['Cat']]
		res_list.sort()
		exp_list.sort()
		
		assert exp_list == res_list

	def testFindReleases3(self):		
		complex = self.complexes['I4'].clone()
				
		BS_index = complex.strand_index('BS')
		PS_index = complex.strand_index('PS')
		
		complex.structure[BS_index][2] = None
		complex.structure[PS_index][3] = None
		
		res_list = find_releases(complex)
		
		exp_complex1 = self.complexes['C1']
		exp_complex2 = self.complexes['I3']
		
		exp_list = [exp_complex1, exp_complex2]
				
		res_list.sort()
		exp_list.sort()
		
		assert exp_list == res_list

	def testFindReleases4(self):
		complex = self.complexes['I3'].clone()
		cat_index = complex.strand_index('Cat')
		bs_index = complex.strand_index('BS')
		
		complex.structure[cat_index][0] = None
		complex.structure[cat_index][1] = None
		complex.structure[bs_index][0] = None
		complex.structure[bs_index][1] = None
		
		res_list = find_releases(complex)
		bs_complex = Complex('BS', [self.strands['BS']], [[None, None, None, None, None, None]])
		
		exp_list = [bs_complex, self.complexes['Cat']]
									
		res_list.sort()
		exp_list.sort()
				
		assert exp_list == res_list

	def testFindReleases5(self):
		# Test multiple releases case
		
		complex = self.complexes['I4'].clone()
		BS_index = complex.strand_index('BS')
		PS_index = complex.strand_index('PS')
		cat_index = complex.strand_index('Cat')
		
		complex.structure[BS_index][2] = None
		complex.structure[PS_index][3] = None
		complex.structure[cat_index][0] = None
		complex.structure[cat_index][1] = None
		complex.structure[BS_index][0] = None
		complex.structure[BS_index][1] = None
		
		res_list = find_releases(complex)
		
		exp_list = [self.complexes['C1'], self.complexes['Cat'], Complex('BS', [self.strands['BS']], [[None, None, None, None, None, None]])]
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list
		
	def testFindReleases6(self):
		# Test another no releases case
		
		complex = Complex('T1', [self.strands['BS'], self.strands['OP'], self.strands['PS'], self.strands['Cat']], [[(3, 1), (3, 0), (2, 3), None, None, None], [(2, 2), (2, 1), None, None], [None, (1, 1), (1, 0), (0, 2), None], [(0, 1), (0, 0)]])
		
		res_list = find_releases(complex)
		
		exp_list = [complex]
		
		res_list.sort()
		exp_list.sort()
		
		
		assert exp_list == res_list
		
	def testOpen1(self):
		# first test no reaction case
		res_list = open(self.complexes['I3'])
		
		assert res_list == []
		
	def testOpen2(self):
		# next test a simple case
		
		complex = Complex('t1', [self.strands['PS'], self.strands['BS']], [[None, None, None, (1, 2), None], [None, None, (0, 3), None, None, None]])
		res_list = open(complex)
		
		exp_rp = ReactionPathway('open', [complex], [Complex('PS', [self.strands['PS']], [[None, None, None, None, None]]), Complex('BS', [self.strands['BS']], [[None, None, None, None, None, None]])])
		
		assert res_list == [exp_rp]
		
	def testOpen3(self):
		# test a real case
		
		res_list = open(self.complexes['I4'])
		exp_list = [ReactionPathway('open', [self.complexes['I4']], [self.complexes['I3'], self.complexes['C1']])]
		
		res_list.sort()
		exp_list.sort()
		
		print res_list
		print exp_list
		
		assert res_list == exp_list
		