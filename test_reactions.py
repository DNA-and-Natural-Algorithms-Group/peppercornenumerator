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
from enumerator import Enumerator

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
	
		self.three_arm_enumerator = input_standard('test_files/test_input_standard_3arm_junction.in')
	
		
		for complex in self.three_arm_enumerator.initial_complexes:
			self.complexes[complex.name] = complex
	
		self.three_arm_enumerator_reduced = Enumerator(self.three_arm_enumerator.domains, self.three_arm_enumerator.strands, [self.complexes['I'], self.complexes['A'], self.complexes['B'], self.complexes['C']])
		
		
		self.index_parts = index_parts
		
	def testFindExternalStrandBreak(self):

		#	complex I4 :
		#	BS OP PS Cat
		#	(((... + (((. + )))). + ))	
		I4 = self.complexes['I4']
		
		# find index of these strands
		for (i, strand) in enumerate(I4.strands):
			if (strand.name == 'BS'):
				BS_index = i
			elif (strand.name == 'PS'):
				PS_index = i
		
		assert PS_index == find_external_strand_break(I4, (PS_index, 4))
		
		assert BS_index == find_external_strand_break(I4, (BS_index, 3))
		
		
		#	complex I1 :
		#	SP Cat BS
		#	(( + .( + )))...
		I1 = self.complexes['I1']
		
		# find index of these strands
		for (i, strand) in enumerate(I1.strands):
			if (strand.name == 'SP'):
				SP_index = i
			elif (strand.name == 'Cat'):
				Cat_index = i
			elif (strand.name == 'BS'):
				BS_index = i
		
		assert SP_index == find_external_strand_break(I1, (Cat_index, 0))
		
		assert BS_index == find_external_strand_break(I1, (BS_index, 4))
		
		
		# complex C :
		# s1 s2
		# .( + ).
		
		s1 = Strand('s1', [self.domains['1'], self.domains['2']])
		s2 = Strand('s2', [self.domains['2*'], self.domains['3']])
		
		c = Complex('C', [s1, s2], [[None, (1, 0)], [(0, 1), None]])
		
		assert find_external_strand_break(c, (c.strand_index('s1'), 0)) == 1
		
		
		# complex C : 
		# s1 s2 s3
		# (. + ( + ).)
		
		s1 = Strand('s1', [self.domains['4'], self.domains['2']])
		s2 = Strand('s2', [self.domains['1']])
		s3 = Strand('s3', [self.domains['1*'], self.domains['3*'], self.domains['4*']])
		
		c = Complex('C', [s1, s2, s3], [[(2, 2), None], [(2, 0)], [(1, 0), None, (0, 0)]])
		
		assert find_external_strand_break(c, (c.strand_index('s1'), 1)) == 0
		
		
		# complex C :
		# s1 s2 s3
		# ( + .(. + )..(.))
		
		s1 = Strand('s1', [self.domains['1']])
		s2 = Strand('s2', [self.domains['5'], self.domains['2'], self.domains['6'], self.domains['3'], self.domains['4'], self.domains['3*'], self.domains['1*']])
		s3 = Strand('s3', [self.domains['4*'], self.domains['5*'], self.domains['6']])
		
		c = Complex('C', [s1, s3, s2], [[(2, 6)], [None, (2, 0), None], [(1, 1), None, None, (2, 5), None, (2, 3), (0, 0)]])
		
		assert find_external_strand_break(c, (c.strand_index('s2'), 2)) == 0
		
		def try_pseudoknot(self):
			c.structure[c.strand_index('s3')][0] = (c.strand_index('s2'), 4)
			c.structure[c.strand_index('s2')][4] = (c.strand_index('s3'), 0)
			
			find_external_strand_break(c, (c.strand_index('s2'), 2))
			
		assert_raises(Exception, try_pseudoknot, self)
		
		
				
	def testBind11(self):
		# Test with no reaction
		
		out_list = bind11(self.complexes['I4'])
		exp_list = []
		
		assert out_list == exp_list
	
		# Test simple
		strand = Strand('A', [self.domains['1'], self.domains['2'], self.domains['1*']])
		complex = Complex('C', [strand], [[None, None, None]])
		
		out_list = bind11(complex)
		
		exp_complex = Complex('C-New', [strand], [[(0, 2), None, (0, 0)]])
		exp_out = [ReactionPathway('bind11', [complex], [exp_complex])]
		
		out_list.sort()
		exp_out.sort()
		
		assert out_list == exp_out
		
		# Test complicated
		
		c = self.complexes['I4'].clone()
		BS_index = c.strand_index('BS')
		OP_index = c.strand_index('OP')
		PS_index = c.strand_index('PS')
		Cat_index = c.strand_index('Cat')
		
		c.structure[OP_index][0] = None
		c.structure[PS_index][2] = None
		
		c1 = c.clone()
		c1.structure[PS_index][2] = (BS_index, 3)
		c1.structure[BS_index][3] = (PS_index, 2)
		
		out_list = bind11(c)
		exp_out = [ReactionPathway('bind11', [c], [self.complexes['I4']]), ReactionPathway('bind11', [c], [c1])]
		
		out_list.sort()
		exp_out.sort()
		
		assert out_list == exp_out
		
	def test_combine_complexes_21(self):
		PS_index = self.complexes['C1'].strand_index('PS')
		BS_index = self.complexes['I3'].strand_index('BS')
		
		out_complex = combine_complexes_21(self.complexes['C1'], (PS_index, 3),
										   self.complexes['I3'], (BS_index, 2))
		
		exp_complex = Complex('I4', [self.strands['PS'], self.strands['Cat'], self.strands['BS'], self.strands['OP']], [[(3, 2), (3, 1), (3, 0), (2, 2), None], [(2, 1), (2, 0)], [(1, 1), (1, 0), (0, 3), None, None, None], [(0, 2), (0, 1), (0, 0), None]])
		
		
		assert exp_complex == out_complex
		
		BS_index = self.complexes['W'].strand_index('BS')
		
		out_complex = combine_complexes_21(self.complexes['W'], (BS_index, 0),
										   self.complexes['Cat'], (0, 1))
										
		exp_complex = Complex('C', [self.strands['BS'], self.strands['PS'], self.strands['Cat']], [[(2, 1), (1, 4), (1, 3), (1, 2), (1, 1), (1, 0)], [(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)], [None, (0, 0)]])
		
		assert exp_complex == out_complex
		
		s1 = Strand('A', [self.domains['1'], self.domains['2']])
		s2 = Strand('B', [self.domains['2*'], self.domains['3']])
		
		c1 = Complex('A', [s1], [[None, None]])
		c2 = Complex('B', [s2], [[None, None]])
		
		out_complex = combine_complexes_21(c2, (0, 0), c1, (0, 1))
		
		exp_complex = Complex('AB', [s1, s2], [[None, (1, 0)], [(0, 1), None]])
		
		print "Output: "
		print out_complex
		print out_complex.strands
		print out_complex.structure
		print "Expected: "
		print exp_complex
		print exp_complex.strands
		print exp_complex.structure		
		
		assert out_complex == exp_complex
		
	def test_combine_complexes_21_seesaw(self):
		self.seesaw_enum = input_standard('test_files/examples/seesaw/seesaw.enum')
		(domains,strands,complexes) = self.index_parts(self.seesaw_enum)
		
		exp_complex = Complex('complex',[strands['S2_T_S3'],strands['S2_T_S3'],strands['T_S3_T']],[[None,(2,2),(2,1)],[None,(2,0),None],[(1,1),(0,2),(0,1)]])
		out_complex = combine_complexes_21(complexes['Waste'], (1,0), complexes['Fuel'], (0,1))
		
		print "Waste: "
		print
		print complexes['Waste']
		print complexes['Waste'].strands
		print complexes['Waste'].structure
		print
		print "Fuel: "
		print complexes['Fuel']
		print complexes['Fuel'].strands
		print complexes['Fuel'].structure
		print
		
		print "Output: "
		print out_complex
		print out_complex.strands
		print out_complex.structure
		print
		print "Expected: "
		print exp_complex
		print exp_complex.strands
		print exp_complex.structure		
		
		assert out_complex == exp_complex
		
	def test_combine_complexes_21_seesaw_2(self):
		self.seesaw_enum = input_standard('test_files/examples/seesaw/seesaw2.enum')
		(domains,strands,complexes) = self.index_parts(self.seesaw_enum)
		out_complex = reactions.combine_complexes_21(complexes["C1"], (0,3), complexes["C2"], (0,1))
		exp_complex = complexes['C3']
		
		print "Output: "
		print out_complex
		print out_complex.strands
		print out_complex.structure
		print
		print "Expected: "
		print exp_complex
		print exp_complex.strands
		print exp_complex.structure	
		
		assert out_complex == exp_complex
		
		
		
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
		# Test another multiple releases case
		
		complex = self.complexes['I4'].clone()
		BS_index = complex.strand_index('BS')
		PS_index = complex.strand_index('PS')
		OP_index = complex.strand_index('OP')
		cat_index = complex.strand_index('Cat')
		
		complex.structure[BS_index][2] = None
		complex.structure[PS_index][3] = None
		
		complex.structure[OP_index][0] = None
		complex.structure[OP_index][1] = None
		complex.structure[OP_index][2] = None
		complex.structure[PS_index][0] = None
		complex.structure[PS_index][1] = None
		complex.structure[PS_index][2] = None
		
		res_list = find_releases(complex)
		
		exp_list = [self.complexes['I3'], Complex('OP', [self.strands['OP']], [[None, None, None, None]]), \
				Complex('PS', [self.strands['PS']], [[None, None, None, None, None]])]
		
		res_list.sort()
		exp_list.sort()
		
		
		print "Computed Result (res_list):"	
		print res_list
		for complex in res_list:
			print complex.full_string()
		
		print "Expected result (exp_list):"
		print exp_list
		for complex in exp_list:
			print complex.full_string()
		
		assert res_list == exp_list

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
		
		complex = Complex('T1', [self.strands['BS'], self.strands['OP'], self.strands['PS'], self.strands['Cat']], \
						 [[(3, 1), (3, 0), (2, 3), None, None, None], [(2, 2), (2, 1), None, None], [None, (1, 1), (1, 0), (0, 2), None], [(0, 1), (0, 0)]])
		
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
		
		assert res_list == exp_list
	
	def testOpen4(self):
		# another real case
		
		res_list = open(self.complexes['I1'])
		exp_list = [ReactionPathway('open', [self.complexes['I1']], [self.complexes['C2'], self.complexes['Cat']])]
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list
		
	def testOpen5(self):
		
		res_list = open(self.complexes['I2'])
		exp_list = [ReactionPathway('open', [self.complexes['I2']], [self.complexes['I3'], self.complexes['SP']])]
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list		
		
	def testOpen6(self):
		# test an Open1-1 reaction
		strand = Strand('A', [self.domains['1'], self.domains['2'], self.domains['1*']])
		complex = Complex('A', [strand], [[(0, 2), None, (0, 0)]])
		
		res_list = open(complex)
		exp_list = [ReactionPathway('open', [complex], [Complex('B', [strand], [[None, None, None]])])]
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list

	def testOpen7(self):
		S1 = Strand('S1', [self.domains['1'], self.domains['4']])
		S2 = Strand('S2', [self.domains['4*']])
		S3 = Strand('S3', [self.domains['1*'], self.domains['4*']])
		complex = Complex('C', [S1, S2, S3], [[(2, 0), (1, 0)], [(0, 1)], [(0, 0), None]])

		print "Input Complex"
		print complex.full_string()
		print "Strand:", complex.strands[0], " ", complex.strands[0].domains
		for domain in complex.strands[0].domains:
			print "Domain:", domain, " ", domain.length
		
		print
		
		
		res_list = open(complex)
		exp_list = [ReactionPathway('open', [complex], sorted([Complex('C1', [S3], [[None, None]]), Complex('C2', [S1, S2], [[None, (1, 0)], [(0, 1)]])]))]

		print "Computed Result (res_list):"	
		print res_list
		print "Outputs:"
		for complex in res_list[0].products:
			print complex.full_string()
		
		print
		print "Expected result (exp_list):"
		print exp_list
		for complex in exp_list[0].products:
			print complex.full_string()
		assert res_list == exp_list
		
class BranchMigrationTests(unittest.TestCase):
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



	def testDo3wayMigration1(self):
		# test a 3-way producing 1 outcomplex
		res_list = do_3way_migration(self.complexes['I1'], (2, 0), (0, 1))
		exp_list = [self.complexes['I2']]
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list
		
	def testDo3wayMigration2(self):
		# test a 3-way producing 2 outcomplexes
		complex = self.complexes['I5'].clone()
				
		complex.structure[complex.strand_index('BS')][0] = None
		complex.structure[complex.strand_index('Cat')][1] = None
		
		res_list = do_3way_migration(complex, (complex.strand_index('PS'), 4), (complex.strand_index('BS'), 1))
		exp_list = [self.complexes['W'], self.complexes['Cat']]
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list
		
	def testBranch3way1(self):
		s1 = Strand('S1', [self.domains['1'], self.domains['2']])
		s2 = Strand('S2', [self.domains['1*']])
		s3 = Strand('S3', [self.domains['2*'], self.domains['1*']])
		
		c1 = Complex('c1', [s1, s3, s2], [[(2, 0), (1, 0)], [(0, 1), None], [(0, 0)]])
				
		res_list = branch_3way(c1)
		
		c2 = Complex('c2', [s1, s3], [[(1, 1), (1, 0)], [(0, 1), (0, 0)]])
		c3 = Complex('c3', [s2], [[None]])
		
		exp_list = [ReactionPathway('branch_3way', [c1], [c2, c3])]
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list
	
	def testBranch3way2(self):
		res_list = branch_3way(self.complexes['I1'])
		exp_list = [ReactionPathway('branch_3way', [self.complexes['I1']], [self.complexes['I2']])]
		
		res_list.sort()
		exp_list.sort()
		
		print res_list
		print exp_list
		
		assert res_list == exp_list
		
	def testBranch3way3(self):
		c = self.complexes['I4']
		res_list = branch_3way(c)
		
		BS_index = c.strand_index('BS')
		PS_index = c.strand_index('PS')
		Cat_index = c.strand_index('Cat')
		OP_index = c.strand_index('OP')
		
		o1 = c.clone()
		o1.structure[BS_index][1] = (PS_index, 4)
		o1.structure[PS_index][4] = (BS_index, 1)
		o1.structure[Cat_index][0] = None
		
		exp_list = []
		exp_list.append(ReactionPathway('branch_3way', [c], [o1]))
		
		# These results are expected with UNZIP=False
		#o2 = c.clone()
		#o2.structure[BS_index][3] = (PS_index, 2)
		#o2.structure[PS_index][2] = (BS_index, 3)
		#o2.structure[OP_index][0] = None
		
		#exp_list.append(ReactionPathway('branch_3way', [c], [o2]))
		
		exp_list.append(ReactionPathway('branch_3way', [c], [self.complexes['OP'], self.complexes['I5']]))
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list
		
	# Test with remote toehold
	
	def testBranch3way4(self):
		enumerator = input_standard('test_files/test_input_standard_3arm_junction.in')
		
		self.domains = {}
		self.strands = {}
		self.complexes = {}		
		
		for domain in enumerator.domains:
			self.domains[domain.name] = domain
		
		for strand in enumerator.strands:
			self.strands[strand.name] = strand
		
		for complex in enumerator.initial_complexes:
			self.complexes[complex.name] = complex
			
		res_list = branch_3way(self.complexes['IABC'])
		
		# These results are expected with UNZIP=False
		#out_complex = Complex('IABC-new', [self.strands['I'], self.strands['A'], self.strands['B'], self.strands['C']], [[None, (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), (3, 4), (2, 3), (2, 2), (2, 1), (2, 0), None], [(1, 7), (1, 6), (1, 5), (1, 4), (3, 3), (3, 2), (3, 1), (3, 0), None], [(2, 7), (2, 6), (2, 5), (2, 4), (1, 3), None, None, None, None]])
		
		
		#exp_list = [ReactionPathway('branch_3way', [self.complexes['IABC']], [out_complex])]
		
				
		exp_list = [ReactionPathway('branch_3way', [self.complexes['IABC']], [self.complexes['I'], self.complexes['ABC']])]
		res_list.sort()
		exp_list.sort()
				
		print res_list
		print exp_list
		
		assert res_list == exp_list

	def testBranch3way5(self):
		
		enumerator = input_standard('test_files/test_input_standard_3arm_junction.in')
		
		self.domains = {}
		self.strands = {}
		self.complexes = {}		
		
		for domain in enumerator.domains:
			self.domains[domain.name] = domain
		
		for strand in enumerator.strands:
			self.strands[strand.name] = strand
		
		for complex in enumerator.initial_complexes:
			self.complexes[complex.name] = complex
		
		IAbind = Complex('IAbind', [self.strands['I'], self.strands['A']], [[None, None, None, (1, 0)], [(0, 3), (1, 8), (1, 7), (1, 6), None, None, (1, 3), (1, 2), (1, 1)]])
		
		res_list = branch_3way(IAbind)
		exp_list = [ReactionPathway('branch_3way', [IAbind], [self.complexes['IA']])]
		
		print res_list
		print exp_list
		assert res_list == exp_list

	def testBranch3way6(self):
		enumerator = input_standard('test_files/test_input_standard_remote.in')
		
		self.domains = {}
		self.strands = {}
		self.complexes = {}		
		
		for domain in enumerator.domains:
			self.domains[domain.name] = domain
		
		for strand in enumerator.strands:
			self.strands[strand.name] = strand
		
		for complex in enumerator.initial_complexes:
			self.complexes[complex.name] = complex
		
		res_list = branch_3way(self.complexes['C1'])
		exp_list = [ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2'],self.complexes['C3']] )]
		
		print res_list
		print exp_list
		assert res_list == exp_list
		
	def testDo4wayMigration1(self):	
		s1 = Strand('s1', [self.domains['1*'], self.domains['2*'], self.domains['3']])
		s2 = Strand('s2', [self.domains['3*'], self.domains['2'], self.domains['4']])
		s3 = Strand('s3', [self.domains['4*'], self.domains['2*'], self.domains['5*']])
		s4 = Strand('s4', [self.domains['5'], self.domains['2'], self.domains['1']])
		
		c1 = Complex('c1', [s1, s2, s3, s4], [[(3, 2), (3, 1), (1, 0)], [(0, 2), (2, 1), (2, 0)], [(1, 2), (1, 1), (3, 0)], [(2, 2), (0, 1), (0, 0)]])
		
		c2 = Complex('c2', [s1, s2, s3, s4], [[(3, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (2, 0)], [(1, 2), (3, 1), (3, 0)], [(2, 2), (2, 1), (0, 0)]])
		
		res_list = do_4way_migration(c1, (0, 1), (3, 1), (1, 1), (2, 1))
		exp_list = [c2]
		
		res_list.sort()
		exp_list.sort()
		
		assert exp_list == res_list
		
	def testDo4wayMigration2(self):		
		s1 = Strand('s1', [self.domains['1*'], self.domains['2*'], self.domains['3']])
		s2 = Strand('s2', [self.domains['3*'], self.domains['2'], self.domains['1']])
		s3 = Strand('s3', [self.domains['1*'], self.domains['2*'], self.domains['5*']])
		s4 = Strand('s4', [self.domains['5'], self.domains['2'], self.domains['1']])
		
		c2 = Complex('c2', [s1, s2, s3, s4], [[(3, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (2, 0)], [(1, 2), (3, 1), (3, 0)], [(2, 2), (2, 1), (0, 0)]])
		
		c3 = Complex('c3', [s1, s2], [[(1, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (0, 0)]])
		c4 = Complex('c4', [s3, s4], [[(1, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (0, 0)]])
		
		res_list = do_4way_migration(c2, (0, 0), (3, 2), (1, 2), (2, 0))
		exp_list = [c3, c4]
		
		res_list.sort()
		exp_list.sort()
		
		assert exp_list == res_list
		
	def testBranch4way1(self):
		# test with no reaction
		
		res_list = branch_4way(self.complexes['C1'])
		exp_list = []
		
		res_list.sort()
		exp_list.sort()
		
		assert res_list == exp_list
		
	def testBranch4way2(self):
		s1 = Strand('s1', [self.domains['1*'], self.domains['2*'], self.domains['3']])
		s2 = Strand('s2', [self.domains['3*'], self.domains['2'], self.domains['4']])
		s3 = Strand('s3', [self.domains['4*'], self.domains['2*'], self.domains['5*']])
		s4 = Strand('s4', [self.domains['5'], self.domains['2'], self.domains['1']])
		
		c1 = Complex('c1', [s1, s2, s3, s4], [[(3, 2), (3, 1), (1, 0)], [(0, 2), (2, 1), (2, 0)], [(1, 2), (1, 1), (3, 0)], [(2, 2), (0, 1), (0, 0)]])
		
		c2 = Complex('c2', [s1, s2, s3, s4], [[(3, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (2, 0)], [(1, 2), (3, 1), (3, 0)], [(2, 2), (2, 1), (0, 0)]])
		
		res_list = branch_4way(c1)
		exp_list = [ReactionPathway('branch_4way', [c1], [c2])]
		
		res_list.sort()
		exp_list.sort()
				
		assert exp_list == res_list		
		
	def testBranch4way3(self):
		s1 = Strand('s1', [self.domains['1*'], self.domains['2*'], self.domains['3']])
		s2 = Strand('s2', [self.domains['3*'], self.domains['2'], self.domains['4']])
		s3 = Strand('s3', [self.domains['4*'], self.domains['2*'], self.domains['5*']])
		s4 = Strand('s4', [self.domains['5'], self.domains['2'], self.domains['1']])
		
		c1 = Complex('c1', [s1, s2, s3, s4], [[(3, 2), (3, 1), (1, 0)], [(0, 2), (2, 1), (2, 0)], [(1, 2), (1, 1), (3, 0)], [(2, 2), (0, 1), (0, 0)]])
		
		c2 = Complex('c2', [s1, s2, s3, s4], [[(3, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (2, 0)], [(1, 2), (3, 1), (3, 0)], [(2, 2), (2, 1), (0, 0)]])
		
		res_list = branch_4way(c2)
		exp_list = [ReactionPathway('branch_4way', [c2], [c1])]
		
		res_list.sort()
		exp_list.sort()
				
		assert exp_list == res_list			
		
	def testBranch4way4(self):
		s1 = Strand('s1', [self.domains['1*'], self.domains['2*'], self.domains['3']])
		s2 = Strand('s2', [self.domains['3*'], self.domains['2'], self.domains['1']])
		s3 = Strand('s3', [self.domains['1*'], self.domains['2*'], self.domains['5*']])
		s4 = Strand('s4', [self.domains['5'], self.domains['2'], self.domains['1']])
		
		c1 = Complex('c1', [s1, s2, s3, s4], [[(3, 2), (3, 1), (1, 0)], [(0, 2), (2, 1), (2, 0)], [(1, 2), (1, 1), (3, 0)], [(2, 2), (0, 1), (0, 0)]])
		
		c2 = Complex('c2', [s1, s2, s3, s4], [[(3, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (2, 0)], [(1, 2), (3, 1), (3, 0)], [(2, 2), (2, 1), (0, 0)]])
		
		c3 = Complex('c3', [s1, s2], [[(1, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (0, 0)]])
		c4 = Complex('c4', [s3, s4], [[(1, 2), (1, 1), (1, 0)], [(0, 2), (0, 1), (0, 0)]])
		
		res_list = branch_4way(c2)
		exp_list = [ReactionPathway('branch_4way', [c2], [c3, c4]), ReactionPathway('branch_4way', [c2], [c1])]
		
		print res_list
		print exp_list
		
		res_list.sort()
		exp_list.sort()
		
		assert exp_list == res_list	

class ReactionPathwayTests(unittest.TestCase):
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
	
	def testHash(self):
		assert hash(ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']])) == hash(ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']]))
		assert hash(ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']])) != hash(ReactionPathway('branch_4way', [self.complexes['C2']], [self.complexes['Cat']]))
	
	def testName(self):
		self.rp = ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']])
		
		assert self.rp.name == 'branch_3way'
		def assnName(self):
			self.rp.name = 'branch_1way'
			
		assert_raises(AttributeError, assnName, self)
		
	def testStr(self):
		self.rp = ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']])
		assert str(self.rp) == self.rp.name
		
	def testEq(self):
		rp1 = ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']])
		
		rp2 = ReactionPathway('branch_4way', [self.complexes['C2']], [self.complexes['Cat']])
		
		assert not rp1 == rp2
		
		rp3 = ReactionPathway('branch_3way', [Complex('C1', [self.strands['PS'], self.strands['OP']], [[(1, 2), (1, 1), (1, 0), None, None], [(0, 2), (0, 1), (0, 0), None]])
], [self.complexes['C2']])

		assert rp1 == rp3
		
		rp4 = ReactionPathway('branch_4way', [self.complexes['C1']], [self.complexes['C2']])
		
		assert not rp1 == rp4
		
		rp5 = ReactionPathway('branch_3way', [self.complexes['Cat']], [self.complexes['C2']])
		
		assert not rp1 == rp5
		
		rp6 = ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['Cat']])
		
		assert not rp1 == rp6
					
	def testCmp(self):
		rp1 = ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']])		
		rp2 = ReactionPathway('branch_4way', [self.complexes['C2']], [self.complexes['Cat']])		
		rp3 = ReactionPathway('branch_3way', [Complex('C1', [self.strands['PS'], self.strands['OP']], [[(1, 2), (1, 1), (1, 0), None, None], [(0, 2), (0, 1), (0, 0), None]])
], [self.complexes['C2']])
		rp4 = ReactionPathway('branch_4way', [self.complexes['C1']], [self.complexes['C2']])
		rp5 = ReactionPathway('branch_3way', [self.complexes['Cat']], [self.complexes['C2']])
		rp6 = ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['Cat']])
		
		assert cmp(rp1, rp1) == 0
		assert cmp(rp1, rp3) == 0
		assert cmp(rp1, rp2) == cmp('branch_3way', 'branch_4way')
		assert cmp(rp1, rp5) == cmp([self.complexes['C1']], [self.complexes['Cat']])
		assert cmp(rp1, rp6) == cmp([self.complexes['C2']], [self.complexes['Cat']])
		
	def testNormalize(self):
		rp1 = ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']])		
		rp2 = ReactionPathway('branch_3way', [self.complexes['C1']], [self.complexes['C2']])
		rp2.normalize()
		
		assert rp1 == rp2		
		
		rp3 = ReactionPathway('branch_3way', [self.complexes['C1'], self.complexes['Cat'], self.complexes['Cat']], [self.complexes['C2'], self.complexes['Cat'], self.complexes['Cat']])
		rp3.normalize()
		
		assert rp1 == rp3
		
		