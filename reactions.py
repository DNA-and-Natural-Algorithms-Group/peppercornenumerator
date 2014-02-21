#
#  reactions.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/18/10.
#

import copy
import utils
from utils import *

# TODO: fix naming!
auto_name = 0

def get_auto_name():
	global auto_name
	auto_name += 1
	return str(auto_name)

# TODO: refactor this
RELEASE_CUTOFF = 6

# If true, 3 way branch migrations are always greedy
UNZIP = True

class ReactionPathway(object):
	"""
	Represents a reaction node on a reaction graph. Has a list of reactants
	and a list of products, along with a name field to denote the type of
	reaction involved.
	"""
	# TODO: This needs to be actually used!
	
	def __init__(self, name, reactants, products):
		"""
		Default constructor. Takes a name, a list of reactants (Complexes) and
		a list of products (Complexes).
		"""
		self._name = name
		self._reactants = reactants
		self._reactants.sort()
		self._products = products
		self._products.sort()
		
	def __repr__(self):
		return self.full_string()
#		return "ReactionPathway(%s)" % (self.name)
	
	def __str__(self):
		return self.name

	def full_string(self):
		return "ReactionPathway(%s): %s -> %s" % (self.name, str(self.reactants), str(self.products))

	@property
	def name(self):
		return self._name
		
	@property
	def reactants(self):
		return self._reactants
		
	@property
	def products(self):
		return self._products
		
	@property
	def arity(self):
		return (len(self._reactants),len(self._products))	
	
	def __eq__(self, other):
		return (self.name == other.name) and \
			   (self.reactants == other.reactants) and \
			   (self.products == other.products)

	def __hash__(self):
		return hash(frozenset(self.reactants)) + hash(frozenset(self.products))

	def __cmp__(self, other):
		out = cmp(self.name, other.name)
		if (out != 0):
			return out
			
		out = cmp(self.reactants, other.reactants)
		if (out != 0):
			return out
			
		return cmp(self.products, other.products)
			   
	def normalize(self):
		"""
		Ensures that complexes appear on only one side of the reaction by
		removing them evenly from both sides until only one side has any.
		"""
		for reactant in self.reactants:
			while (reactant in self.reactants and \
				   reactant in self.products):
				self.reactants.remove(reactant)
				self.products.remove(reactant)
				
				


def bind11(reactant):
	"""
	Returns a list of reaction pathways which can be produced by 1-1 binding 
	reactions of the argument complex. The 1-1 binding reaction is the 
	hybridization of two complementary unpaired domains within a single complex 
	to produce a single unpseudoknotted product complex.
	"""
	# tuple: (strand, domain)
	outer_index = (0, 0)
	inner_index = (0, 0)
	struct = reactant.structure
	strands = reactant.strands
	reactions = []
		
	# We will loop over all unpaired domains
	# TODO: we can just use the available_domains field of complexes to
	# simplify this!
	for (strand_num, outer_strand) in enumerate(strands):
		for (domain_num, outer_domain) in enumerate(outer_strand.domains):
			
			if (struct[strand_num][domain_num] != None):
				# This is not an unpaired domain
				continue
			
			# For each unpaired domain, we loop over all higher-number unpaired
			# domains
			inner_strand = strand_num
			inner_domain = domain_num + 1
			
			while (inner_strand != len(strands)):
				# If we're back where we started, then we've traversed 
				# a whole loop
				if (inner_strand == strand_num) and \
				   (inner_domain == domain_num):
				   break
				
				# If we're at the end of a strand, move to the next
				if (inner_domain == len(strands[inner_strand].domains)):
					inner_domain = 0
					inner_strand = inner_strand + 1
				elif (struct[inner_strand][inner_domain] != None):
					# This is not an unpaired domain; skip around the loop
					# to return to an external loop by following the structure
					struct_element = struct[inner_strand][inner_domain]
					inner_strand = struct_element[0]
					inner_domain = struct_element[1] + 1
				elif not outer_domain.can_pair(strands[inner_strand].domains[inner_domain]):
					# These domains aren't complementary
					inner_domain = inner_domain + 1
				else:
					# These domains can pair, add to the list
					reactions.append(((strand_num, domain_num), 
									 (inner_strand, inner_domain)))
					# Continue to loop
					inner_domain = inner_domain + 1
	
	products = []
	for (d1, d2) in reactions:
		new_structure = copy.deepcopy(reactant.structure)
		new_structure[d1[0]][d1[1]] = d2
		new_structure[d2[0]][d2[1]] = d1
		product = Complex(get_auto_name(),  #reactant.name + "(" + str(d1) + "+" + str(d2) + ")",
						  reactant.strands, new_structure)		
		products.append(product)
	
	output = []
	for product in products:
		output.append(ReactionPathway('bind11', [reactant], [product]))

	# remove any duplicate reactions	
	if (len(output) == 0):
		return output
		
	output.sort()
	last = output[-1]
	for i in range(len(output) - 2, -1, -1):
		if last == output[i]:
			del output[i]
		else:
			last = output[i]

	return output
	
def bind21(reactant1, reactant2):
	"""
	Returns a list of reaction pathways which can be produced by 2-1 binding 
	reactions of the argument complexes. The 2-1 binding reaction is the 
	hybridization of two complementary unpaired domains, each in a different 
	complex, to produce a single, unpseudoknotted product complex containing 
	all of the strands contained in either of the original complexes.
	"""
	r1_doms = reactant1.available_domains
	r2_doms = reactant2.available_domains
	
	interaction_list = []
	
	# TODO: significant optimization likely possible here because the lists
	# are sorted
	
	new_complexes = []
	
	# Iterate through all the free domains in reactant1
	for (i, (dom1, strand_num1, dom_num1)) in enumerate(r1_doms):
		# For each, find any domains in reactant2 that could bind
		for (j, (dom2, strand_num2, dom_num2)) in enumerate(r2_doms):
			# If it can pair, this is one possible reaction (this kind of
			# reaction cannot possibly produce a pseudoknotted structure)
			if (dom1.can_pair(dom2)):
				new_complexes.append(combine_complexes_21(
									 reactant1, (strand_num1, dom_num1), 
								     reactant2, (strand_num2, dom_num2)))
									 
	output = []
	for complex in new_complexes:
		output.append(ReactionPathway('bind21', [reactant1, reactant2], 
								      [complex]))
	
	return output
	
	
def find_external_strand_break(complex, location):
	"""
	Takes a complex and a location (strand index, domain index). This location
	refers to a domain in the `complex` which is on an external loop. Returns the 
	index of the last strand _before_ the strand break which would put `location` 
	on an external loop.
	
	Finds the location of a strand break on an external loop.  Used to 
	determine where to split a complex when merging it with another complex at
	`location`.
	"""
	
	strand_num = location[0]
	domain_num = location[1]
	
	stop_flag = False
	search_dom_index = domain_num + 1
	search_strand_index = strand_num
	insertion_strand_index = None
	
	# We will first search to the 'right', looking for the first external
	# strand break
	while insertion_strand_index == None:
		# First check to see if we've run off the end of a strand
		# in which case the external break is between this strand and the
		# next one
		if (search_dom_index == len(complex.strands[search_strand_index])):
			insertion_strand_index = search_strand_index
		else:
			paired_dom = complex.structure[search_strand_index]\
										  [search_dom_index]
			if (paired_dom == None):			
				# If the current domain is unpaired, move to the next domain
				# to the right
				search_dom_index += 1
			
			# elif ((paired_dom[0] > search_strand_index) or \
			# 		 ((paired_dom[0] == search_strand_index) and \
			# 		  (paired_dom[1] > search_dom_index))):
			elif (paired_dom > (search_strand_index,search_dom_index)):
				# If the current domain is paired to something which is to the
				# right of the current search domain, then we jump past that
				# loop
				search_strand_index = paired_dom[0]
				search_dom_index = paired_dom[1] + 1
			else:
				# Otherwise the current domain is paired to something to the
				# left of the current search domain, which means that the
				# external strand break must be to the left
				break
				
	search_dom_index = domain_num - 1
	search_strand_index = strand_num
	
	# We now search to the 'left' if we haven't found the break yet
	while insertion_strand_index == None:
		# First check to see if we've run off the end of a strand, in which
		# case the external break is between the previous strand and this one
		if (search_dom_index == -1):
			insertion_strand_index = search_strand_index - 1
		else:
			paired_dom = complex.structure[search_strand_index] \
										  [search_dom_index]
			if (paired_dom == None):			
				# If the current domain is unpaired, move to the next domain
				# to the left
				search_dom_index -= 1
			# elif ((paired_dom[0] < search_strand_index) or 
			# 		 ((paired_dom[0] == search_strand_index) and
			# 		  (paired_dom[1] < search_dom_index))):
			elif (paired_dom < (search_strand_index,search_dom_index)):
				# If the current domain is paired to something which is to the
				# left of the current search domain, then we jump past that
				# loop
				search_strand_index = paired_dom[0]
				search_dom_index = paired_dom[1] - 1
			else:
				# Otherwise the current domain is paired to something to the
				# left of the current search domain, which means that the
				# external strand break must be to the right; however we
				# checked that, so there must have been an error...
				raise \
				    Exception("Unexpected error in find_external_strand_break")
	
	return insertion_strand_index


def combine_complexes_21(complex1, location1, complex2, location2):
	"""
	Combines two complexes to form one complex, binding the domain in
	complex1 at location1 to the domain in complex2 at location2.
	
	Returns the new complex.
	"""
	
	# First we need to find the external strand breaks where we will be
	# splitting the complexes
	
	insertion_index_1 = find_external_strand_break(complex1, location1)
	insertion_index_2 = find_external_strand_break(complex2, location2)
		
	# We then find the four parts, 1-4, which we will stick together
	# in order to create our final complex
	
	# If this condition is true, we will need to cut up complex1
	if insertion_index_1 >= 0:
		d1 = complex1.strands[0:(insertion_index_1+1)]
		d4 = complex1.strands[(insertion_index_1+1):]
		s1 = copy.deepcopy(complex1.structure[0:(insertion_index_1+1)])
		s4 = copy.deepcopy(complex1.structure[(insertion_index_1+1):])	
	else:
		# It's not clear to me whether this  condition will ever be executed, 
		# because find_external_strand_break should never return < 0 -CG 5/31
		d1 = complex1.strands[:]
		d4 = []
		s1 = copy.deepcopy(complex1.structure)
		s4 = []
	
	# If this condition is true, we will need to cut up complex2	
	if insertion_index_2 >= 0:
		d2 = complex2.strands[(insertion_index_2+1):]
		d3 = complex2.strands[0:(insertion_index_2+1)]
		s2 = copy.deepcopy(complex2.structure[(insertion_index_2+1):])
		s3 = copy.deepcopy(complex2.structure[0:(insertion_index_2+1)])
	else:
		# Likewise with this one; find_external_strand_break should always
		# produce an index >= 0
		d2 = []
		d3 = complex2.strands[:]
		s2 = []
		s3 = copy.deepcopy(complex2.structure)
	
	# We now need to update the structures based on the way things were shifted
	# We calculate the offsets for strand locations in each chunk
	s1_strand_offset = 0
	s2_strand_offset = len(d1) - len(d3)
	s3_strand_offset = len(d1) + len(d2)
	s4_strand_offset = len(d2) + len(d3)
	
	
	# We then iterate through the structure and update as needed
	if insertion_index_1 >= 0:
		for (strand_index, strand_list) in enumerate(s1):
			for (dom_index, pair) in enumerate(strand_list):
				# If the domain is paired, check if it is paired to something 
				# in s4
				if (pair and (pair[0] > insertion_index_1)):
					new_pair = (pair[0] + s4_strand_offset, pair[1])
					s1[strand_index][dom_index] = new_pair
	
	for (strand_index, strand_list) in enumerate(s2):
		for (dom_index, pair) in enumerate(strand_list):
			# If the domain is paired, check to see if it is paired to something
			# in s2 or s3
			if (pair):
				# Check if paired to s2
				if ((insertion_index_2 >= 0) and (pair[0] > insertion_index_2)):
					new_pair = (pair[0] + s2_strand_offset, pair[1])
					s2[strand_index][dom_index] = new_pair
				# Otherwise paired to s3
				else:
					new_pair = (pair[0] + s3_strand_offset, pair[1])
					s2[strand_index][dom_index] = new_pair
	
		
	for (strand_index, strand_list) in enumerate(s3):
		for (dom_index, pair) in enumerate(strand_list):
			# If the domain is paired, check to see if it is paired to something
			# in s2 or s3
			if (pair):
				# Check if paired to s2
				if ((insertion_index_2 >= 0) and (pair[0] > insertion_index_2)):
					new_pair = (pair[0] + s2_strand_offset, pair[1])
					s3[strand_index][dom_index] = new_pair
				# Otherwise paired to s3
				else:
					new_pair = (pair[0] + s3_strand_offset, pair[1])
					s3[strand_index][dom_index] = new_pair
	
	for (strand_index, strand_list) in enumerate(s4):
		for (dom_index, pair) in enumerate(strand_list):
			# Check if the domain is paired to something in s4
			if (pair and (pair[0] > insertion_index_1)):
				new_pair = (pair[0] + s4_strand_offset, pair[1])
				s4[strand_index][dom_index] = new_pair
	
	new_strands = d1 + d2 + d3 + d4
	new_structure = s1 + s2 + s3 + s4
	
	# Finally, we need to update given reaction indices to new indices
	# and then add the new pairing to the structure
	# if ((insertion_index_1 > 0) and (location1[0] > insertion_index_1)):
	if (location1[0] > insertion_index_1):
		location1 = (location1[0] + s4_strand_offset, location1[1])
	
	
	# NOTE: This case is unexercised by all tests written by Karthik, regardless
	# whether the first clause is included. Since this tests whether the 
	# location2 is within d2 or d3, I think the first clause is unnecessary.

	# if ((insertion_index_2 > 0) and (location2[0] > insertion_index_2)):
	if (location2[0] > insertion_index_2):
		location2 = (location2[0] + s2_strand_offset, location2[1])
	else:
		location2 = (location2[0] + s3_strand_offset, location2[1])
	
	
	new_structure[location1[0]][location1[1]] = location2
	new_structure[location2[0]][location2[1]] = location1
#	try:
#		new_structure[location1[0]][location1[1]] = location2
#		new_structure[location2[0]][location2[1]] = location1
#	except IndexError:
#		print "Agh!"
#		raise Exception()
		
	new_complex = Complex(get_auto_name(), new_strands, new_structure)
	
	return new_complex

def open(reactant):
	"""
	Returns a list of reaction product sets that can be produced by the
	'open' reaction, in which a short helix dissociates. Each product
	set are the results of one particular dissociation; each strand in the
	reactant occurs exactly once in one of the complexes in the product set.
	
	A dissociation can happen to any helix under the threshold length		
	"""

	product_sets = []
	
	
	structure = reactant.structure
	strands = reactant.strands
	
	
	# We loop through all stands, domains
	for (strand_index, strand) in enumerate(strands):
		for (domain_index, domain) in enumerate(strand.domains):
			# If the domain is unpaired, skip it
			if (structure[strand_index][domain_index] == None):
				continue
			
			# A: Strand/domain position on "top" strand - CG 5/21 
			helix_startA = [strand_index, domain_index]
			
			# B: Strand/domain position on "bottom" strand - CG 5/21 
			helix_startB = list(structure[strand_index][domain_index])
		
			# If the domain is bound to an earlier domain, then we have
			# already considered it, so skip it
			if (helix_startB < helix_startA):
				continue
			
			helix_endA = helix_startA[:]
			helix_endB = helix_startB[:]
						
			helix_length = domain.length
			
			# Now iterate through the whole helix to find the other end
			# of this one
			# (The helix ends at the first strand break from either direction)
			while True:
				
				# Strands run in opposite directions, so A must be incremented 
				# and B decremented in order that both pointers move "right" 
				# along the helix- CG 5/21
				helix_endA[1] += 1
				helix_endB[1] -= 1
				
				# If one of the strands has broken, the helix has ended
				if (helix_endA[1] >= strands[helix_endA[0]].length):
					break
				elif (helix_endB[1] < 0):
					break
				
				
				# If these domains aren't bound to each other, the helix
				# has ended
				if (tuple(helix_endA) != structure[helix_endB[0]][helix_endB[1]]):
					break
				
				# Add the current domain to the current helix
				helix_length += strands[helix_endA[0]].domains[helix_endA[1]]\
									   .length
				
			# We must also iterate in the other direction
			while True:
				
				# Now we move 
				helix_startA[1] -= 1
				helix_startB[1] += 1

				# If one of the strands has broken, the helix has ended
				if (helix_startA[1] < 0):
					break
				elif (helix_startB[1] >= strands[helix_startB[0]].length):
					break
				
				# If these domains aren't bound to each other, the helix
				# has ended
				if (tuple(helix_startA) != structure[helix_startB[0]][helix_startB[1]]):
					break
				
				# Add the current domain to the current helix
				helix_length += strands[helix_startA[0]].domains[helix_startA[1]]\
									   .length
				
			
			# Move start location to the first domain in the helix
			helix_startA[1] += 1
			helix_startB[1] -= 1
			
			# If the helix is short enough, we have a reaction	
			if (helix_length <= RELEASE_CUTOFF):


				release_reactant = Complex(get_auto_name(), reactant.strands[:], 
					copy.deepcopy(reactant.structure))
				
				# Delete all the pairs in the released helix
				for dom in range(helix_startA[1], helix_endA[1]):
					bound_loc = reactant.structure[helix_startA[0]][dom]
					release_reactant.structure[helix_startA[0]][dom] = None
					release_reactant.structure[bound_loc[0]][bound_loc[1]] = None
		
		
				product_set = find_releases(release_reactant)
				
				
				product_sets.append(product_set)
		
	output = []
	for product_set in product_sets:
		output.append(ReactionPathway('open', [reactant], sorted(product_set)))
	
	return output
	
	
def find_releases(reactant):
	"""
	Determines if reactant actually represents multiple unattached complexes.
	If so, returns a list of these complexes. Otherwise, returns [reactant].
	"""
	
	structure = reactant.structure
	strands = reactant.strands
	
	output_list = []
	
	# Iterate through the strands and determine if a split is necessary
	# Note that we don't iterate to the last strand because if the last
	# strand were free from the rest, we would catch that as the rest
	# of the strands being free from the last strand. Hence, we also
	# terminate loops if we ever get to the last strand
	for (strand_index, strand) in enumerate(strands[:-1]):
		
		domain_index = len(strand.domains)
		
		inner_index = (strand_index, domain_index - 1)
		
		# We now iterate through lower domains and see if we can find
		# a release point (Iterate while:
		#	0 > inner_index strand > index of last strand, and
		# 		inner_index strand < strand_index, or
		#		inner_index strand = strand_index and inner_index domain < domain_index)
		while (inner_index[0] >= 0) and (inner_index[0] < (len(strands) - 1)) and \
			  ( inner_index < (strand_index, domain_index)):


			# If we have run off of the end of a strand,
			# then we have found a release point  
			if (inner_index[1] == -1):
				split_start = inner_index
				split_end = (strand_index, domain_index)
				split_list = split_complex(reactant, split_start, split_end)
				# We check the two resulting complexes to see if they can
				# be split further
				for complex in split_list:
					output_list += (find_releases(complex))
				return output_list					
			
			# Otherwise decide where to go next
			curr_structure = structure[inner_index[0]][inner_index[1]]
						
			# If this domain is unpaired, move to the next
			if (curr_structure == None):
				inner_index = (inner_index[0], inner_index[1] - 1)
			
			# Check if the structure points to a higher domain 
			elif ((curr_structure[0] > inner_index[0]) or \
			     ((curr_structure[0] == inner_index[0]) and \
				  (curr_structure[1] > inner_index[1]))):

				# If the structure points to a domain above the start, this section
				# is connected to something higher, so abort this loop
				#if (curr_structure[0] > strand_index) and (curr_structure[1] > 0):
				if (curr_structure[0] > strand_index) and (curr_structure[1] >= 0):
					break
				  				  
				# Otherwise it points to a domain between this one and the start
				# -- we've already been there so move to the next
				inner_index = (inner_index[0], inner_index[1] - 1)			
				
			# Otherwise, follow the structure
			else:
				inner_index = curr_structure


		# Alright; the following statement precludes the while loop below from 
		# ever running, since the second clause will always fail. 
		# Then again, without ever running this while loop, all other tests pass.
		# OK, I'm fairly convinced this isn't needed, since the above loop 
		# effectively works by iterating (forwards) across each strand, then looking 
		# backwards -CG 5/29
		
		
		##### TODO: FIND OUT IF THIS IS NEEDED...
#		inner_index = (strand_index, 0)
#				
#		# If we didn't find a release point in the lower domains,
#		# we now try iterating through the higher domains
#		while (inner_index[0] < len(strands) - 1) and \
#			  ((inner_index[0] > strand_index) or \
#			  ((inner_index[0] == strand_index) and \
#			  (inner_index[1] > domain_index))):
#			# If we have run off of the end of a strand,
#			# then we have found a release point  
#			if (inner_index[1] == len(strands[inner_index[0]].domains)):
#				split_start = (strand_index, domain_index)
#				split_end = inner_index
#				split_list = split_complex(reactant, split_start, split_end)
#				# We check the two resulting complexes to see if they can
#				# be split further
#				for complex in split_list:
#					output_list.extend(find_releases(complex))
#				return output_list					
#							
#			# Otherwise decide where to go next
#			curr_structure = structure[inner_index[0]][inner_index[1]]	
#			
#			# If this domain is unpaired, move to the next
#			if (curr_structure == None):
#				inner_index = (inner_index[0], inner_index[1] + 1)
#
#			# Check if the structure points to a lower domain			
#			elif (curr_structure[0] < inner_index[0]) or \
#				 ((curr_structure[0] == inner_index[0]) and \
#				  (curr_structure[1] < inner_index[1])):
#				  
#				# If the structure points to a domain before the start,
#				# this section is connected to something lower, so abort this
#				# loop
#				if (curr_structure[0] < strand_index):
#					break
#				
#				# Otherwise it points to a domain between this one and the start
#				# -- we've already been there to move to the next
#				inner_index = (inner_index[0], inner_index[1] + 1)			
#			
#			
#			# Otherwise, follow the structure
#			else:
#				inner_index = curr_structure			
				
	# If we still haven't found any splits, then this complex cannot be split
	
	return [reactant]
	
def split_complex(reactant, split_start, split_end):
	"""
	Splits a disconnected complex between split_start and split_end (which
	should be at ends of strands), and returns the two resulting complexes.
	"""
	
	strands = reactant.strands
	structure = reactant.structure
	
	# We first take the parts apart
	strands1 = strands[0:split_start[0]]
	structure1 = structure[0:split_start[0]]
	
	strands2 = strands[split_start[0]:split_end[0]+1]
	structure2 = structure[split_start[0]:split_end[0]+1]
	
	strands3 = strands[split_end[0]+1:]
	structure3 = structure[split_end[0]+1:]
	
	out1_strands = strands1 + strands3
	out2_strands = strands2
		
		
	# These offsets are the changes to the strand numbers needed
	structure2_offset = len(strands1)
	structure3_offset = len(strands2)

	# Apply offset to structure2
	new_structure1 = []
			
	for strand in structure1:
		strand_out = []
		for struct_element in strand:
			if struct_element == None:
				strand_out.append(None)
			elif (struct_element[0] > split_end[0]):
				strand_out.append((struct_element[0] - structure3_offset, struct_element[1]))
			else:
				strand_out.append(struct_element)
		new_structure1.append(strand_out)
	
			
	# Apply offset to structure2
	new_structure2 = []
		
	for strand in structure2:
		strand_out = []
		for struct_element in strand:
			if struct_element == None:
				strand_out.append(None)
			else:
				strand_out.append((struct_element[0] - structure2_offset, struct_element[1]))
		new_structure2.append(strand_out)
		
	# Apply offset to structure3
	new_structure3 = []
	
	
	for strand in structure3:
		strand_out = []
		for struct_element in strand:
			if struct_element == None:
				strand_out.append(None)
			elif (struct_element[0] < split_start[0]):
				strand_out.append(struct_element)
			else:
				strand_out.append((struct_element[0] - structure3_offset, struct_element[1]))
		new_structure3.append(strand_out)
		
	out1_structure = new_structure1 + new_structure3
	out2_structure = new_structure2
	
	
	out1 = Complex(get_auto_name(), out1_strands, out1_structure)
	out2 = Complex(get_auto_name(), out2_strands, out2_structure)
	
	return [out1, out2]

def branch_3way(reactant):
	"""
	Returns a list of reaction pathways that can be created through one 
	iteration of a 3 way branch migration reaction (more than one molecule may 
	be produced by a reaction because branch migration can liberate strands and 
	complexes).
	"""
	
	output_sets = []
	structure = reactant.structure
	
	# We iterate through all the domains
	for (strand_index, strand) in enumerate(reactant.strands):
		for (domain_index, domain) in enumerate(strand.domains):
			# We will search in both directions
			# First direction:
			
			# The starting domain must be anchored
			if (structure[strand_index][domain_index] == None):
				continue
			
			# The starting domain must have another (displacing) domain 
			# next to it
			if ((domain_index + 1) == reactant.strands[strand_index].length):
				continue
			
			# The displacing domain must be free			
			if (structure[strand_index][domain_index + 1] != None):
				continue
			
			displacing_domain = strand.domains[domain_index + 1]
			
			# We now follow the external loop from the starting pair
			# searching for a strand to displace
			bound_loc = structure[strand_index][domain_index]
			bound_loc_orig = structure[strand_index][domain_index]
			# Follow the external loop to the end
			while True:				
				bound_loc = (bound_loc[0], bound_loc[1] - 1)
				# Check if we've reached the end of the external loop
				if (bound_loc[1] == -1):
					# Reached the end
					break
				
				
				
				# Check if this domain is unbound
				elif (structure[bound_loc[0]][bound_loc[1]] == None):
					continue
					
				# Check to see if this domain is complementary
				elif (reactant.strands[bound_loc[0]].domains[bound_loc[1]]\
							.can_pair(displacing_domain)):
					# We have found a displacement reaction
					output_sets.append(do_3way_migration(reactant,
														 (strand_index, 
														 domain_index + 1),
														 bound_loc))
					
					
				# follow the structure
				
				bound_loc = structure[bound_loc[0]][bound_loc[1]]
				if bound_loc == bound_loc_orig:
					# Caught in a loop
					break
			

	for (strand_index, strand) in enumerate(reactant.strands):
		for (domain_index, domain) in enumerate(strand.domains):			
			# Second direction:

			# The starting domain must be anchored
			if (structure[strand_index][domain_index] == None):
				continue
						
			# The starting domain must have another (displacing) domain 
			# next to it
			if ((domain_index - 1) == -1):
				continue
				
			# The displacing domain must be free			
			if (structure[strand_index][domain_index - 1] != None):
				continue
				
			displacing_domain = strand.domains[domain_index - 1]
			
			# We now follow the external loop from the starting pair
			# searching for a strand to displace
			bound_loc = structure[strand_index][domain_index]
			bound_loc_orig = structure[strand_index][domain_index]
			# Follow the external loop to the end
			while True:
				bound_loc = (bound_loc[0], bound_loc[1] + 1)
				
				# Check if we've reached the end of the external loop
				if (bound_loc[1] == reactant.strands[bound_loc[0]].length):
					# Reached the end
					break
					
				# Check if this domain is unbound
				elif (structure[bound_loc[0]][bound_loc[1]] == None):
					continue
					
				# Check to see if this domain is complementary
				elif (reactant.strands[bound_loc[0]].domains[bound_loc[1]]\
							.can_pair(displacing_domain)):
					# We have found a displacement reaction
					output_sets.append(do_3way_migration(reactant,
														 (strand_index, 
														 domain_index - 1),
														 bound_loc))
					
					
				
				bound_loc = structure[bound_loc[0]][bound_loc[1]]
				if bound_loc == bound_loc_orig:
					break
				
	output = []
	for output_set in output_sets:
		output.append(ReactionPathway('branch_3way', [reactant], output_set))		

	# Remove any duplicate reactions
	if (len(output) == 0):
		return output
		
	output.sort()
	last = output[-1]
	for i in range(len(output) - 2, -1, -1):
		if last == output[i]:
			del output[i]
		else:
			last = output[i]
	
	return output


def do_3way_migration(reactant, displacing_loc, new_bound_loc):
	"""
	Returns the product set which is the result of a 3-way branch migration
	reaction where the domain at displacing_loc displaces the domain bound to
	the domain at new_bound_loc.
	"""

	# out_reactant = copy.deepcopy(reactant)

	out_reactant_structure = copy.deepcopy(reactant.structure)

	out_reactant_structure[displacing_loc[0]][displacing_loc[1]] = new_bound_loc
	displaced_loc = out_reactant_structure[new_bound_loc[0]][new_bound_loc[1]]
	out_reactant_structure[new_bound_loc[0]][new_bound_loc[1]] = displacing_loc
	out_reactant_structure[displaced_loc[0]][displaced_loc[1]] = None
	
	out_reactant = Complex(get_auto_name(),reactant.strands[:], out_reactant_structure)
	
	global UNZIP
	
	# Check to see if an adjacent displacement is possible
	if (UNZIP):
		dstrand = displacing_loc[0]
		ddomain = displacing_loc[1]
		bstrand = new_bound_loc[0]
		bdomain = new_bound_loc[1]
		if (ddomain+1 < len(out_reactant.strands[dstrand].domains)) and (out_reactant.structure[dstrand][ddomain+1] == None) and (bdomain-1 >= 0) and (out_reactant.structure[bstrand][bdomain-1] != None) and (out_reactant.strands[bstrand].domains[bdomain-1].can_pair(out_reactant.strands[dstrand].domains[ddomain+1])):
			return do_3way_migration(out_reactant, (dstrand, ddomain+1), (bstrand, bdomain-1))
		elif (ddomain-1 >= 0) and (out_reactant.structure[dstrand][ddomain-1] == None) and (bdomain+1 < len(out_reactant.strands[bstrand].domains)) and (out_reactant.structure[bstrand][bdomain+1] != None) and (out_reactant.strands[bstrand].domains[bdomain+1].can_pair(out_reactant.strands[dstrand].domains[ddomain-1])):
			return do_3way_migration(out_reactant, (dstrand, ddomain-1), (bstrand, bdomain+1))
		else:
			return find_releases(out_reactant)
	else:
		return find_releases(out_reactant)


def branch_4way(reactant):
	"""
	Returns a list of complex sets that can be created through one iteration of
	a 4 way branch migration reaction (each set consists of the molecules that
	result from the iteration; more than one molecule may result because branch
	migration can liberate strands and complexes).	
	"""
	
	structure = reactant.structure
	output_sets = []
	
	# We loop through all domains
	for (strand_index, strand) in enumerate(reactant.strands):
		for (domain_index, domain) in enumerate(strand.domains):
			# Unbound domains can't participate in branch migration
			if (structure[strand_index][domain_index] == None):
				continue
			
			# This domain can't be at the end of a strand
			if (domain_index + 1 == reactant.strands[strand_index].length):
				continue
			
			
			# Check the 4 locations for a 4 way reaction
			
			# Displacing domain
			loc1 = (strand_index, domain_index + 1)
			dom1 = reactant.strands[loc1[0]].domains[loc1[1]]

			# Displaced domain
			loc2 = structure[strand_index][domain_index + 1]
			if (loc2 == None):
				continue
			dom2 = reactant.strands[loc2[0]].domains[loc2[1]]
			
			# Template domain (replaces displaced domain, binds loc1)
			loc3 = structure[strand_index][domain_index]
			if (loc3 == None):
				continue
			
			loc3 = (loc3[0], loc3[1] - 1)
			if (loc3[1] < 0):
				continue
			dom3 = reactant.strands[loc3[0]].domains[loc3[1]]
			
			# Displaced from template domain (replaces displacing domain, binds
			#								  loc2)
			loc4 = structure[loc3[0]][loc3[1]]
			if (loc4 == None):
				continue
			dom4 = reactant.strands[loc4[0]].domains[loc4[1]]
			
			# Confirm that the domains can in fact pair
			if not (dom1.can_pair(dom3) and dom2.can_pair(dom4)):
				continue
				
			# Confirm that this is a four way migration
			if loc2 == loc3:
				continue
				
			# If we are here, then we have found a candidate reaction
			output_sets.append(do_4way_migration(reactant, 
												 loc1, loc2, loc3, loc4))
	
	output = []
	for output_set in output_sets:
		output.append(ReactionPathway('branch_4way', [reactant], output_set))
	
	# remove any duplicate reactions	
	if (len(output) == 0):
		return output
		
	output.sort()
	last = output[-1]
	for i in range(len(output) - 2, -1, -1):
		if last == output[i]:
			del output[i]
		else:
			last = output[i]
	
	return output
	
def do_4way_migration(reactant, loc1, loc2, loc3, loc4):
	"""
	Performs a 4 way branch migration on a copy of reactant, with loc1 as the
	displacing domain, loc2 as the domain displaced from loc1, loc3 as the
	template domain, and loc4 as the domain displaced from loc3. Returns the
	set of complexes produced by this reaction (may be one or more complexes).
	"""
	
	new_struct = copy.deepcopy(reactant.structure)
	new_struct[loc1[0]][loc1[1]] = loc3
	new_struct[loc3[0]][loc3[1]] = loc1
	new_struct[loc2[0]][loc2[1]] = loc4
	new_struct[loc4[0]][loc4[1]] = loc2
	
	out = Complex(get_auto_name(), reactant.strands[:], new_struct)
	
	return find_releases(out)

