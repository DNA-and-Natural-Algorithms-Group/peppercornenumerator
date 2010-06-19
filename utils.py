#
#  utils.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/18/10.
#

import copy

SHORT_DOMAIN_LENGTH = 6
LONG_DOMAIN_LENGTH = 12

class Domain(object):
	"""
	Represents a single domain. We allow several options for specifying domain
	properties. Domains might have an explicit integer (bp) length, or may be 
	designated as short or long. If the latter method is used, the code will use
	the relevant constant as the integer domain length.
	"""
	
	def __init__(self, name, length, is_complement=False, sequence=None):
		"""
		Default constructor. Takes a domain name, a length (positive integer or 
		"short" or "long", and optionally a base sequence.
		"""
		self._name = name
		self._length = length
		self._sequence = None
		self._complement_sequence = None
		self._is_complement = is_complement
		if (sequence):
			self._sequence = sequence.upper()
			self._complement_sequence = self._sequence
			self._complement_sequence = list(self._complement_sequence[::-1])
			for i, l in enumerate(self._complement_sequence):
				if (l == 'A'):
					self._complement_sequence[i] = 'T'
				elif (l == 'T'):
					self._complement_sequence[i] = 'A'
				elif (l == 'C'):
					self._complement_sequence[i] = 'G'
				elif (l == 'G'):
					self._complement_sequence[i] = 'C'
			self._complement_sequence = ''.join(self._complement_sequence)
	
	
	def __repr__(self):
		return "Domain(%s)" % (self.name)
	
	def __str__(self):
		return self.name
	
	def __eq__(self, other):
		return (self.name == other.name)
		
	def __len__(self):
		return self.length
		
	def can_pair(self, other):
		"""
		Returns True if this domain is complementary to the argument.
		"""
		return ((self.identity == other.identity) and \
			   (self.is_complement or other.is_complement) and \
			   (not (self.is_complement and other.is_complement)))
		
	@property
	def identity(self):
		"""
		Returns the identity of this domain, which is its name without a
		complement specifier (i.e. A and A* both have identity A).
		"""
		return self._name
		
	@property
	def name(self):
		"""
		The name of this domain.
		"""
		return (self._name + ("" if not self.is_complement else "*"))
	
	@property
	def length(self):
		"""
		The length of this domain. Either uses the integer length previously 
		specified or the constant associated with "short" and "long" domains as 
		appropriate.
		"""
				
		if (type(self._length) == type(0)):
			return self._length
		elif (self._length == "short"):
			return SHORT_DOMAIN_LENGTH
		elif (self._length == "long"):
			return LONG_DOMAIN_LENGTH
	
	@property
	def sequence(self):
		"""
		This domain's sequence. May be None if no sequence was specified.
		"""
		if (not self.is_complement):
			return self._sequence
		else:
			return self._complement_sequence
			
	@property
	def is_complement(self):
		"""
		Returns true if this domain is a complement or not.
		"""		
		return self._is_complement
					
class Strand(object):
	"""
	Represents a strand, which is an ordered sequence of domains.
	"""
	
	def __init__(self, name, domains):
		"""
		Constructor for Strand objects. Takes a strand name and the ordered list
		of Domain objects from 5' to 3'.
		"""
		self._name = name
		self._domains = domains

	def __eq__(self, other):
		return (self.name == other.name) and (self.domains == other.domains)

	def __len__(self):
		return len(self.domains)

	@property
	def name(self):
		return self._name
		
	@property
	def length(self):
		return len(domains)
		
	@property
	def domains(self):
		return self._domains

class Complex(object):
	"""
	Represents a complex, which is a set of connected strands with a particular
	secondary structure.
	"""
	
	def __init__(self, name, strands, structure):
		"""
		Constructor for complex objects. Takes a name, the ordered list of
		strands, and the appropriate structure. Complexes are assumed to be
		nonpseudoknotted, and the strand list should be in the order that yields
		a nonpseudoknotted structure.
		
		A complex is in 'canonical form' if its strands are in a non-pseudoknotted
		circular permutation with the strand with the minimum name lexicographically
		as the first in the list of strands.
		
		This constructor assumes that the structure is unpseudoknotted, and will
		rotate this complex automatically until it is in canonical form.
		
		name = string holding complex name
		strands = list of Strand objects in order
		structure = list of lists of tuples indicating pairing of domains in 
					complex with None indicating unpaired --
					tuple: (strand, domain)
					Ex:
					[[(0, 2) None (0, 1)]] 
					indicates one strand with 3 domains with the first one bound 
					to the last one, and the middle one free. 
					
					[[None (1, 0) (1, 1)] [(0, 1) (0, 2) None]]
					indicates 2 strands with 3 domains each -- the first two
					domains of the second strand are bound to the last two
					of the first.
					
					Lists should be in the same order as the strands in the
					second argument, with each strand's domains from 5' to 3'
		"""
		
		# Holds the complex name
		self._name = name
		
		# Holds the list of strands
		self._strands = strands
		
		# Holds the complex structure
		self._structure = structure
		
		# Holds the number of domains total in this complex
		self._ndoms = 0
		
		strandNames = []
		
		for strand in strands:
			self._ndoms += len(strand.domains)
			strandNames.append(strand.name)
			
		# Find the minimum name to determine canonical form	
		strandNames.sort()
		
		# A list of domains on external loops, in a tuple with format
		# (domain, strand number, domain in strand number)
		self._available_domains = []
		
		# Keeps track of whether or not the _available_domains variable is valid
		# or if it needs to be updated
		self._valid_available_domains = False				
		
		# Rotate until we're in canonical form
		# TODO: this could be optimized by writing a rotation method that
		#		rotates N strands at a time
		while (strandNames[0] != self.strands[0].name):
			self._rotate_strands()
			
	def __eq__(self, other):
		"""
		Tests two complexes for equality. Complexes are equal if and only if
		they have the same strands (in the same order) and the same structure.
		This works because of the uniqueness of nonpseudoknotted representations
		of nonpseudoknotted structures (assuming everything is in canonical form
		and strands have unique names).
		"""
		return ((self._strands == other._strands) and 
			   (self._structure == other._structure))
	
	def clone(self):
		"""
		Returns a deep copy of this complex.
		"""
		return copy.deepcopy(self)
	
	@property
	def name(self):
		return self._name
	
	@property
	def strands(self):
		return self._strands
	
	@property
	def structure(self):
		return self._structure
	
	@property
	def available_domains(self):
		if not (self._valid_available_domains):
			self.update_available_domains()			
		return self._available_domains
	
	def update_available_domains(self):
		"""
		Updates the internal list of available domains. Any domains on an
		external loop are added to this list. Sets _valid_available_domains to
		true on completion.
		"""
		self._available_domains = []
		
		# We loop through every domain in the complex to check if it is free
		for (strand_num, strand) in enumerate(self.strands):
			for (dom_num, dom) in enumerate(strand.domains):
			
				# This variable will be set to true if the domain is free
				available = False
				
				# Domains cannot be on an ext. loop if bound
				if self.structure[strand_num][dom_num] == None:				
					# First we iterate 'to the left' and see if we're on an
					# external loop in that direction.
					checking_strand_num = strand_num
					checking_dom_num = dom_num - 1				
					while not available:
						# If we ever get to the other side of the original
						# domain by following the structure then the original 
						# domain may be in a hairpin, and so is not free
						if ((checking_strand_num > strand_num) or
						   ((checking_strand_num == strand_num) and 
						    (checking_dom_num >= dom_num))):
							break
						
						# If we ever get to the left end of a strand without
						# first being linked by structure to another part of the
						# complex then this strand is "dangling" off the left
						# side, and so this is a free domain
						elif checking_dom_num == -1:
							available = True
							
						# If the current domain is unbound, then we keep moving
						# to the left	
						elif self.structure[checking_strand_num]\
										   [checking_dom_num] == None:	
							checking_dom_num -= 1
							
						# If we reach some structure, then we jump to the other
						# end of the that structure (we know this is a non-
						# pseudoknotted complex, so we know we can skip any of
						# the inside structure since none of that will bind
						# with anything on the other side of the original
						# domain	
						else:
							checking_strand_num, checking_dom_num = \
								self.structure[checking_strand_num]\
								[checking_dom_num]
							checking_dom_num -= 1
							
					
					# Next we iterate 'to the right' and see if we're on an
					# external loop in that direction.
					checking_strand_num = strand_num
					checking_dom_num = dom_num + 1
					while not available:
						# If we ever get to the other side of the original
						# domain by following the structure then the original 
						# domain may be in a hairpin, and so is not free
						if ((checking_strand_num < strand_num) or
						   ((checking_strand_num == strand_num) and 
						    (checking_dom_num <= dom_num))):
							break
						
						# If we ever get to the right end of a strand without
						# first being linked by structure to another part of the
						# complex then this strand is "dangling" off the right
						# side, and so this is a free domain
						elif checking_dom_num == len(
												    self.strands
													   [checking_strand_num]
														   .domains):							
							available = True
						
						# If the current domain is unbound, then we keep moving
						# to the right
						elif self.structure[checking_strand_num]\
										   [checking_dom_num] == None:
							checking_dom_num += 1
						
						# If we reach some structure, then we jump to the other
						# end of the that structure (we know this is a non-
						# pseudoknotted complex, so we know we can skip any of
						# the inside structure since none of that will bind
						# with anything on the other side of the original
						# domain	
						else:
							checking_strand_num, checking_dom_num = \
								self.structure[checking_strand_num]\
								[checking_dom_num]
							checking_dom_num -= 1
				if available:
					self._available_domains.append((self.strands[strand_num]\
									                .domains[dom_num], \
													strand_num, dom_num))

		self._available_domains.sort(key=lambda dom: dom[0].name)
		self._valid_available_domains = True

	def _rotate_strands(self):
		"""
		Alters this complex so that all of the strands have been
		rotated to the right by one (that is, the strand formerly at index 1
		is now at index 0, and the strand formerly at index 0 is at index
		n_strands - 1). The structure is updated so that strand references are
		correct. Note that this operation will produce a non-pseudoknotted
		complex so long as the original complex was non-pseudoknotted.
		"""
		new_strands = self.strands[1:] + [self.strands[0]]
		new_struct = []
		n_strands = len(new_strands)
		
		for list in self.structure:
			new_list = []
			for el in list:
				if el == None:
					new_list.append(None)
				else:
					(strand, domain) = el
					if strand > 0:
						new_list.append((strand-1, domain))
					else:
						new_list.append((n_strands-1, domain))
			new_struct.append(new_list)
		
		new_struct = new_struct[1:] + [new_struct[0]]
			
		self._avaliable_domains = []
		self._valid_available_domains = False
		self._structure = new_struct
		self._strands = new_strands	
		
	def rotate_strands(self):
		"""
		Returns a deep copy of this complex with its strands rotated once.
		"""
		out = copy.deepcopy(self)
		out._rotate_strands()
		return out
		
	def dot_paren_string(self):
		"""
		Returns the dot paren representation of this complex.
		"""
		out = []
		for strand_num, strand in enumerate(self.structure):
			for dom_num, el in enumerate(strand):
				if el == None:
					out.append('.')
				else:
					(b_strand, b_domain) = el
					if ((b_strand > strand_num) or
						((b_strand == strand_num) and (b_dom > dom_num))):
						out.append('(')
					else:
						out.append(')')
			if strand_num != (len(self.structure) - 1):
				out.append('+')
		
		return ''.join(out)
		
	def is_unpseudoknotted(self):
		"""
		Returns True if this complex is unpseudoknotted.
		"""
		dparen = self.dot_paren_string()
		list = []
		for char in dparen:
			if (char == '.'):
				pass
			elif (char == '+'):
				pass
			elif (char == '('):
				list.append(0)
			elif (char == ')'):
				if len(list) != 0:
					list = list[1:]
				else:
					return False
		if len(list) == 0:
			return True
		return False