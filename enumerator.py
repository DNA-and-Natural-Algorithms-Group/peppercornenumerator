#
#  enumerator.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/18/10.
#

import utils
import reactions
import logging
import itertools

# These are sanity checks to prevent infinite looping
MAX_COMPLEX_SIZE = 6
MAX_REACTION_COUNT = 1000
MAX_COMPLEX_COUNT = 200


# Fast reactions cannot be bimolecular!
fast_reactions = {
	1 = [bind11, open, branch_3way, branch_4way]
}

# Slow reactions can only be unimolecular or bimolecular, though
# make_slow_reactions below could be changed in order to lift this
# restriction
slow_reactions = {
	1 = [],
	2 = [bind21]
}

class Enumerator(object):
	"""
	Represents a single enumerator instance, consisting of all the information
	required for a reaction graph. This class is the coordinator for the state
	enumerator. Enumerators have immutable starting conditions.
	"""

	def __init__(self, domains, strands, initial_complexes):
		"""
		Initializes the enumerator. Takes a list of domains, a list of strands
		made up of those domains, and a list of the initial complexes, made
		of the strands.		
		"""
		self._domains = domains
		self._strands = strands
		self._initial_complexes = initial_complexes
		self._reactions = None
		self._resting_states = None
		self._complexes = None
		self._transient_complexes = None
		self._resting_complexes = None
		
	@property
	def domains(self):
		return self._domains[:]
		
	@property
	def strands(self):
		return self._strands[:]
	
	@property
	def initial_complexes(self):
		return self._initial_complexes[:]
		
	@property
	def reactions(self):
		if self._reactions == None:
			raise Exception("enumerate not yet called!")
		return self._reactions[:]
		
	@property
	def resting_states(self):
		if self._resting_states == None:
			raise Exception("enumerate not yet called!")
		return self._resting_states[:]
		
	@property
	def complexes(self):
		if self._complexes == None:
			raise Exception("enumerate not yet called!")
		return self._complexes[:]
			
	@property
	def resting_complexes(self):
		if self._resting_complexes == None:
			raise Exception("enumerate not yet called!")
		return self._resting_complexes[:]
		
	@property
	def transient_complexes(self):
		if self._transient_complexes == None:
			raise Exception("enumerate not yet called!")
		return self._transient_complexes[:]
		
			
	def enumerate(self):
		"""
		Generates the reaction graph consisting of all complexes reachable from
		the initial set of complexes. Produces a full list of complexes, resting
		states, and reactions.
		"""
		
		# List E contains enumerated resting states. Only cross-reactions with
		# other end states need to be considered for these complexes. These
		# complexes will remain in this list throughout function execution.
		self._E = []
		
		# List S contains resting states which have not yet had cross-reactions
		# with set E enumerated yet. All self-interactions for these complexes
		# have been enumerated
		self._S = []
		
		# List T contains transient states which have had their self-reactions
		# enumerated. These complexes will remain in this list throughout
		# function execution.
		self._T = []
		
		# List N contains self-enumerated components of the current 
		# 'neighborhood' consisting of states which are connected via fast 
		# reactions to the current complex of interest, but have not yet been 
		# characterized as transient or resting states.
		self._N = []
		
		# List F contains components of the current 'neighborhood' which have
		# not yet had their self-reactions enumerated. They will be moved to N
		# when they are enumerated.
		self._F = []
		
		# List B contains products of bimolecular reactions that have had no
		# reactions enumerated yet. They will be moved to F when their
		# 'neighborhood' is to be considered.
		self._B = self.initial_complexes
		
		self._reactions = []
		self._complexes = []
		self._resting_states = []
		
		# We first generate the states reachable by fast reactions from the
		# initial complexes
		while (len(self._B) > 0):
			# source is the complex from which we will generate a neighborhood
			source = self._B.pop()			
			process_neighborhood(source)
		
		# We now consider slow reactions
		while len(self._S) > 0:
			# element is the complex for which we will consider slow reactions
			element = self._S.pop()
			slow_reactions = self.get_slow_reactions(element)
			self._E.append(element)
			
			# Find the new complexes which were generated
			self._B = self.get_new_products(slow_reactions)
			
			self._reactions.extend(slow_reactions)
			
			while len(self._B) > 0:
				if (len(self._E) + len(self._T) + len(self._S) > self.MAX_COMPLEX_COUNT):
					logging.error("Too many complexes enumerated!")
					raise Exception("Too many complexes generated, aborting...")
				
				if (len(reactions) > self.MAX_REACTION_COUNT):
					logging.error("Too many reactions enumerated!")
					raise Exception("Too many reactions generated, aborting...")
			
				source = self._B.pop()
				process_neighborhood(source)
				
		self._complexes.extend(self._E)
		self._complexes.extend(self._T)
		self._transient_complexes = self._T
		self._resting_complexes = self._E
				
	def process_neighborhood(self, source):
		"""
		Takes a single complex, generates the 'neighborhood' of complexes
		reachable from that complex through fast reactions, classifies these
		complexes as transient or resting state, and modifies the lists and
		list of reactions accordingly.
		"""
		
		# N_reactions holds reactions which are part of the current
		# neighborhood
		N_reactions = []
		
		self._F = [source]
		
		# First we find all of the complexes accessible through fast
		# reactions starting with the source
		while (len(self._F) > 0):
			curr_element = self._F.pop()
			curr_reactions = self.get_fast_reactions(fast_element)								
			new_products = self.get_new_products(curr_reactions)
			self._F.extend(new_products)
			N_reactions.extend(curr_reactions)
			self._N.append(curr_element)
		
		# Now we segment the neighborhood into transient and resting states
		# by finding the strongly connected components
		segmented_neighborhood = self.segment_neighborhood(self._N, N_reactions)
		
		# Resting state complexes are added to S
		self._S.extend(segmented_neighborhood['resting_state_complexes'])
		
		# Transient state complexes are added to T
		self._T.extend(segmented_neighborhood['transient_state_complexes'])
		
		# Resting states are added to the list
		self._resting_states.extend(segmented_neighborhood['resting_states'])
		
		# Reactions are added to the list
		self._reactions.extend(N_reactions)
		
		# Reset neighborhood
		self._N = []
					
	def get_slow_reactions(self, complex):
		"""
		Returns a list of slow reactions possible using complex and other
		complexes in list E as reagents.
		
		This only supports unimolecular and bimolecular reactions. Could be
		extended to support arbitrary reactions.
		"""

		reactions = []
		
		# Do unimolecular reactions
		for function in slow_reactions[1]:
			reactions.extend(function(complex))
			
		# Do bimolecular reactions
		for function in slow_reactions[2]:
			reactions.extend(function(complex, complex))
			for complex2 in self._E:
				reactions.extend(function(complex, complex2))
				
		return reactions			

	def get_fast_reactions(self, complex):
		"""
		Returns a list of fast reactions possible using complex as a reagent.
		"""
	
		reactions = []
		for reaction in fast_reactions[1]:
			reactions.extend(reaction(fast_element))
		return reactions
				
	def get_new_products(self, reactions):
		"""
		Checks the products in the list of reactions. Updates the pointers in
		those reactions to point to pre-existing complexes if necessary. Else,
		returns the new complexes in a list.
		"""
		new_products = []
		
		# Loop over every reaction
		for reaction in reactions:
			# Check every product of the reaction to see if it is new
			for (i, product) in enumerate(reaction.products):
			
				if (len(product.strands) > MAX_COMPLEX_SIZE):
					logging.warning("Complex %s too large, ignoring!" % len(product.strands))
					continue
			
				# This will be set to True if we've already seen this complex
				enumerated = False
				
				# If the product is in any of these lists, we don't need to
				# deal with it, so just update the reaction to point correctly
				for complex in self._E + self._S + self._T + self._N + self._F:
					if (product == complex):
						enumerated = True
						reaction.products[i] = complex
						break
												
				if not enumerated:
					# If the product is in list B, then we need to remove it from
					# that list so that it can be enumerated for self-interactions
					# as part of this neighborhood
					for complex in self._B:
						if (product == complex):
							reaction.products[i] = complex
							self._B.remove(complex)
							break
							
					# If the product has already been seen in this loop, update
					# the pointer appropriately
					for complex in new_products:
						if (products == complex):
							enumerated = True
							reaction.products[i] = complex
							break
				
				if not enumerated:
					new_products.append(product)
					
		return new_products
				
	def segment_neighborhood(self, complexes, reactions):
		"""
		Segments a set of complexes and reactions between them representing a
		neighborhood into resting states and transient states. Returns the set
		of complexes which are transient states, complexes which are in resting
		states, and the set of resting states, all in a dictionary.
		"""
		
		# First we initialize the graph variables that will be used for
		# Tarjan's algorithm
		
		self._tarjan_index = 0
		self._tarjan_stack = []
		self._SCC_stack = []
		
		for node in complexes:
			node._outward_edges = []
			node._full_outward_edges = []
			node._index = -1
			node._lowlink = -1
			
		# Detect which products are actually in the neighborhood	
		for reaction in reactions:
			for product in reaction.products:
				product_in_N = False
				for complex in complexes:
					if (complex == product):
						product_in_N = True
						break
				
				# If this product is in the neighborhood, we have an edge
				if product_in_N:
					# We know all these reactions are unimolecular
					reaction.reactants[0]._outward_edges.append(product)
					
				reaction.reactants[0]._full_outward_edges
									 .extend(reaction.products)
									 
		# We now perform Tarjan's algorithm, marking nodes as appropriate
		for node in complexes:
			if node._index == -1:
				self.tarjans(node)
		
		# Now check to see which of the SCCs are resting states
		resting_states = []
		resting_state_complexes = []
		transient_state_complexes = []
		for scc in self._SCC_stack:
			scc_products = []
			is_resting_state = True
			
			for node in scc:
				for product in node._full_outward_edges:
					scc_products.append(product)
			
			for product in scc_products:
				product_in_scc = False
				for complex in scc:
					if product == complex:
						product_in_scc = True
						break
				
				# If the product is not in the SCC, then there is a fast edge
				# leading out of the SCC, so this is not a resting state
				if not product_in_scc:
					is_resting_state = False
					break
			
			if is_resting_state:
				resting_state_complexes.extend(scc)
				resting_state = RestingState(str(auto_name), scc[:])
				auto_name += 1
				resting_states.append(resting_state)
				
			else:
				transient_state_complexes.extend(scc)
		
		return {
				'resting_states': resting_states, 
			    'resting_state_complexes': resting_state_complexes,
				'transient_states': transient_state_complexes
				}
		
	def tarjans(self, node):
		"""
		Executes an iteration of Tarjan's algorithm (a modified DFS) starting
		at the given node.
		"""
		
		# Set this node's tarjan numbers
		node._index = self._tarjan_index
		node._lowlink = self._tarjan_index
		self._tarjan_index += 1
		self._tarjan_stack.append(node)
		node._onStack = True
		
		# Process all connected nodes, setting lowlink as needed
		for next in node._outward_edges:
			if next._index == -1:
				self.tarjans(next)
				node._lowlink = min(node._lowlink, next._lowlink)
			elif next.onStack:
				node._lowlink = min(node._lowlink, next._lowlink)
				
		# This indicates that this node is a 'root' node, and children are
		# part of an SCC
		if node._lowlink == node._index:
			stop_flag = False
			scc = []
			while stop_flag == False:
				nextNode = self._tarjan_stack.pop()
				nextNode._onStack = False
				scc.append(nextNode)
				if nextNode == node:
					stop_flag = True
			
			# Add the SCC to the list of SCCs
			self._SCC_stack.append(scc)