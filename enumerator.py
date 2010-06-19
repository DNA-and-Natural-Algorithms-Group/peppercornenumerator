#
#  enumerator.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/18/10.
#  Copyright (c) 2010 __MyCompanyName__. All rights reserved.
#

import utils
import reactions

# Fast reactions cannot be bimolecular!
fast_reactions = {
	1 = [bind11, open, branch_3way, branch_4way]
}

slow_reactions = {
	1 = [],
	2 = [bind21]
}


class Enumerator(object):
	"""
	Represents a single enumerator instance, consisting of all the information
	required for a reaction graph.
	"""

	def __init__(self, domains, strands, initial_complexes):
		self._domains = domains
		self._strands = strands
		self._initial_complexes = initial_complexes
		
	@property
	def domains(self):
		return _domains[:]
		
	@property
	def strands(self):
		return _strands[:]
	
	@property
	def initial_complexes(self):
		return _initial_complexes[:]