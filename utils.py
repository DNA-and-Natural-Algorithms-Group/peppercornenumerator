#
#  utils.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/18/10.
#

import copy
import re

SHORT_DOMAIN_LENGTH = 6
LONG_DOMAIN_LENGTH = 12

import re

def natural_sort(l): 
	"""
	Sorts a collection in the order humans would expect. Implementation from
	http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
	"""
	convert = lambda text: int(text) if text.isdigit() else text.lower() 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', str(key)) ] 
	return sorted(l, key = alphanum_key)

def find(f, seq, default=None):
	"""
	Return first item in sequence where f(item) == True.
	"""
	for item in seq:
		if f(item): 
			return item
	return default

def parse_basewise_dot_paren(structure_line, strands):
	parts = [x.strip() for x in structure_line.split("+")]
	assert len(parts) == len(strands), "Structure '%s' has %d parts, but corresponds to %d strands" % \
		(structure_line,len(parts),len(strands))

	segment_struct = ["" for s in strands]
	for (i,s) in enumerate(strands):
		strand_part = parts[i]
		for d in s.domains:
			assert len(strand_part) >= len(d), "Not enough characters for domain %s" % str(d) 
			domain_part, strand_part = strand_part[:len(d)], strand_part[len(d):]

			assert all(c == domain_part[0] for c in domain_part), "Not all parts of structure for %s are the same" % str(d)
			segment_struct[i] += domain_part[0]


	return parse_dot_paren("+".join(segment_struct))

def parse_dot_paren(structure_line):
	"""
	Parses a dot-parenthesis structure into the list of lists representeation
	used elsewhere in the enumerator.

	                         0,0  0,1  0,2   0,3   0,4     1,0   1,1
	Example: "...((+))" -> [[None,None,None,(1,0),(1,1)],[(0,4),(0,3)]] 

	"""	
	complex_structure = []
	dot_paren_stack = []			
	strand_index = 0
	domain_index = 0
	curr_strand = []
	complex_structure.append(curr_strand)
	for part in structure_line:
		# stand break
		if (part == "+"):
			strand_index += 1
			domain_index = 0
			curr_strand = []
			complex_structure.append(curr_strand)
			continue

		# unpaired
		if (part == "."):
			curr_strand.append(None)
			domain_index += 1

		# paired to later domain
		elif (part == "("):
			curr_strand.append(None)
			dot_paren_stack.append((strand_index, domain_index))
			domain_index += 1

		# paired to earlier domain
		elif (part == ")"):
			loc = dot_paren_stack.pop()
			curr_strand.append(loc)
			complex_structure[loc[0]][loc[1]] = (strand_index, domain_index)
			domain_index += 1
	return complex_structure


def index_parts(enum):
	"""
	Testing tool. Accepts an enumerator, produces a tuple 
	(domains, strands, complexes) where each element is a dict mapping names 
	of those objects to the objects in the enumerator. For instance, domains
	maps the name of each domain in the enumerator to the Domain object
	"""
	domains = {}
	strands = {}
	complexes = {}
	
	for domain in enum.domains:
		domains[domain.name] = domain
	
	for strand in enum.strands:
		strands[strand.name] = strand
	
	for complex in enum.initial_complexes:
		complexes[complex.name] = complex
	
	return (domains,strands,complexes)

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
		"short" or "long"), and optionally a base sequence.
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
			
	def __cmp__(self, other):
		out = cmp(self.name, other.name)
		if (out != 0):
			return out
			
		out = self.length - other.length
		if (out != 0):
			return out
		
		if (self.is_complement == other.is_complement):
			return 0			
		elif self.is_complement:
			return 1
		else:
			return 0
		
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
	Represents a strand---an ordered sequence of domains.
	"""
	
	def __init__(self, name, domains):
		"""
		Constructor for Strand objects. Takes a strand name and the ordered list
		of Domain objects from 5' to 3'.
		"""
		self._name = name
		self._domains = domains
		self._hash = None # assigned lazily by __hash__

	def __eq__(self, other):
		return (self.name == other.name) and (self.domains == other.domains)
	
	def __hash__(self):
		if(self._hash == None):
			self._hash = hash((self.name,tuple(self.domains)))
		
		return self._hash
	
	def __cmp__(self, other):
		return cmp(self.name, other.name)

	def __len__(self):
		return len(self.domains)

	def __repr__(self):
		return "Strand(%s)" % (self.name)
	
	def __str__(self):
		return self.name

	@property
	def name(self):
		return self._name
		
	"""
	Returns the number of domains in this strand.
	"""	
	@property
	def length(self):
		return len(self._domains)
		
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
		
		:param name: string holding complex name
		:param strands: list of Strand objects in order
		:param structure: list of lists of tuples indicating pairing of domains in the
					complex. ``None`` indicates the domain is unpaired, while a
					``(strand, domain)`` tuple indicates the domain is paired to ``domain``
					on ``strand``.

					Examples:

					*	``[[(0, 2) None (0, 0)]]`` indicates one strand with 3 domains 
						with the first one bound to the last one, and the middle one free. 
					*	``[[None (1, 0) (1, 1)], [(0, 1) (0, 2) None]]`` indicates 2 
						strands with 3 domains each -- the first two domains of the second 
						strand are bound to the last two of the first.
					
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
		perms = [tuple(strands[x:] + strands[:x]) for x in range(len(strands))] # calculate all circular permutations
		min_perm = min(range(len(strands)),key=lambda i : perms[i]) # select lexicographically minimum permutation
		self._rotate_strands_n(min_perm) # rotate until we're in that form
		
		# Holds a unique hash identifying this complex (computed lazily by self.__hash__)
		self._hash = None
	
	def __hash__(self):
		"""
		Computes a unique hash to represent this complex. Uses the tuple of 
		strands and structure
		"""
		if (self._hash == None):
			strands = tuple(self._strands)
			struct = tuple([tuple(s) for s in self._structure])
			self._hash = hash( (strands, struct) )
			
		return self._hash
	
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
	
	def __cmp__(self, other):
		out = cmp(self.strands, other.strands)
		if (out != 0):
			return out
			
		out = cmp(self.structure, other.structure)
		if (out != 0):
			return out
	
		return cmp(self.name, other.name)
	
	def clone(self):
		"""
		Returns a deep copy of this complex.
		"""
		return copy.deepcopy(self)

	def __repr__(self):
#		return self.full_string()
		return "Complex(%s)" % (self.name)
	
	def __str__(self):
		return self.name
		
	def full_string(self):
		return "Complex(%s): %s %s" % (self.name, str(self.strands), str(self.structure))
	
	@property
	def name(self):
		return self._name
	
	@property
	def strands(self):
		return self._strands[:]
	
	@property
	def structure(self):
		return self._structure[:]
	
	@property
	def available_domains(self):
		if not (self._valid_available_domains):
			self.update_available_domains()			
		return self._available_domains
	
	def get_domain(self,loc):
		"""
		Returns the domain at the given (strand,domain) index in this complex (loc is a (strand,domain) tuple)
		"""
		if(loc != None):
			return self.strands[loc[0]].domains[loc[1]]
		return None
		
	def strand_index(self, strand_name):
		"""
		Returns the index of the strand with the specified name in this
		complex, or -1 if there is no such strand.
		"""
		for (i, strand) in enumerate(self.strands):
			if strand.name == strand_name:
				return i
		
		return -1
	
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
#							checking_dom_num -= 1
							checking_dom_num += 1
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
	
	def _rotate_strands_n(self,n):
		new_strands = self.strands[:]
		new_struct = self.structure[:]
		old_struct = self.structure[:]
		
		for x in range(n):
			new_strands = new_strands[1:] + [new_strands[0]]
		
			new_struct = []
			n_strands = len(new_strands)
			
			for list in old_struct:
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
			old_struct = new_struct[:]
			
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
	
	def check_structure(self):
		"""
		Determines whether the structure includes pairs only between complementary domains. 
		Returns True if all paired domains are complementary, raises an Exception otherwise
		"""
		for (strand_index,strand_struct) in enumerate(self.structure):
			for (domain_index,target) in enumerate(strand_struct):
				
				source_domain = self.get_domain((strand_index,domain_index))
				target_domain = self.get_domain(target)
				
				if((source_domain is None) and (target_domain is None)):
					if(not source_domain.can_pair(target_domain)):
						raise Exception("In complex %s, domain %s is paired with domain %s, but the domains are not complementary." % self.name,
									source_domain.name,target_domain.name)
					
					if(self.structure[target[0]][target[1]] != (strand_index,domain_index)):
						raise Exception("In complex %s, incoherent structure at (%d, %d) and (%d, %d)" % self.name, strand_index, domain_index, target[0], target[1])
		return True
				
	def dot_paren_string(self):
		"""
		Returns the segment-wise dot paren representation of this complex.
		"""
		out = []
		for strand_num, strand in enumerate(self.structure):
			for dom_num, el in enumerate(strand):
				if el == None:
					out.append('.')
				else:
					(b_strand, b_domain) = el
					if ((b_strand > strand_num) or
						((b_strand == strand_num) and (b_domain > dom_num))):
						out.append('(')
					else:
						out.append(')')
			if strand_num != (len(self.structure) - 1):
				out.append('+')
		
		return ''.join(out)
	
	def dot_paren_string_full(self):
		"""
		Returns the base-wise dot paren representation of this complex.
		"""
		out = []
		for strand_num, strand in enumerate(self.structure):
			for dom_num, el in enumerate(strand):
				if el == None:
					out.append('.' * len(self.strands[strand_num].domains[dom_num]))
				else:
					(b_strand, b_domain) = el
					if ((b_strand > strand_num) or
						((b_strand == strand_num) and (b_domain > dom_num))):
						out.append('(' * len(self.strands[strand_num].domains[dom_num]))
					else:
						out.append(')' * len(self.strands[strand_num].domains[dom_num]))
			if strand_num != (len(self.structure) - 1):
				out.append('+')
		
		return ''.join(out)
		
		
class RestingState(object):
	"""
	Represents a resting state, which is a collection of complexes.
	"""
	
	def __init__(self, name, complexes):
		"""
		Constructor for RestingState objects. Takes a (unique) name and list of
		complexes.
		"""
		complexes.sort()
		self._complexes = complexes
		self._name = name
		self._canonical = find(lambda s: not str(s).isdigit(),sorted(complexes),str(complexes[0]))
		self._hash = None
		
	@property
	def name(self):
		return self._name
		
	@property
	def complexes(self):
		return self._complexes[:]
	
	@property
	def canonical_name(self):
		return str(self._canonical)
	
	@property
	def canonical(self):
		return self._canonical

	def __hash__(self):
		if self._hash == None:
			self._hash = hash(tuple(self.complexes))
		return self._hash
	
	def __eq__(self, other):
		return (self.complexes == other.complexes)
		
	def __cmp__(self, other):
		return cmp(self.complexes, other.complexes)

	def __str__(self):
		return self.canonical_name

	def __repr__(self):
		return "RestingState(%s: %s)" % (self.name, str(self.complexes))
