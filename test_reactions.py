#
#  test_reactions.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 5/19/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#

from utils import *
from reactions import *
from test_utils import setUpSLC

import unittest
from nose.tools import *
import copy

class BindTests(unittest.TestCase):
	def setUp(self):
		setUpSLC(self)
		
	def testFindExternalStrandBreak(self):
		A = Strand('A', [self.domains['1'], self.domains['2']])
		B = Strand('B', [self.domains['3'], self.domains['4']])
		C = Strand('C', [self.domains['4*'], self.domains['2*']])
		
		D = Strand('D', [self.domains['5'], self.domains['3*']])
		E = Strand('E', [self.domains['5*']])
		
		C1 = Complex('C1', [A, B, C], [[None, (2, 1)], [None, (2, 0)], [(1, 1), (0, 1)]])
		C2 = Complex('C2', [D, E], [[(1, 0), None], [(0, 0)]])
		
		b1 = find_external_strand_break(C1, (1, 0))
		b2 = find_external_strand_break(C2, (0, 1))
		
		print (b1, b2)
		assert False