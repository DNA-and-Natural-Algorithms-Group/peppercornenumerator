#
#  input.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/26/10.
#

import logging
import json
from utils import *
from enumerator import Enumerator
from reactions import ReactionPathway
import reactions

def input_standard(filename):
	"""
	Initializes and returns an enumerator from a standard input file.
	"""
	fin = open(filename, 'r')
	
	domains = {}
	strands = {}
	complexes = {}
	
	line_counter = 1
	line = fin.readline()
	
	# We loop over all the lines in the file
	while (line != ""):
		line = line.strip()
		
		# This was an empty line
		if line == "":
			line = fin.readline()
			continue
			
		# This was a comment
		elif line.startswith("#"):
			line = fin.readline()
			continue
			
		# This is the beginning of a domain
		elif line.startswith("domain"):
			# e.g.: 
			#       "domain a : 6"
			# parts: 0      1 2 3
			
			parts = line.split()
						
			domain_name = parts[1]
			if domain_name in domains:
				logging.error("Duplicate domain name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not domain_name.isalnum():
				logging.warn("Non-alphanumeric domain name %s encountered in input line %d"
								% (domain_name, line_counter))
			
			# The domain length could be either short or long or it could be
			# an exact number
			domain_length = parts[3]
			if not ((domain_length == 'short') or (domain_length == 'long')):
				domain_length = int(domain_length)
				if domain_length <= 0:
					logging.warn("Domain of length %d found in input line %d" 
									% (domain_length, line_counter))
				
			# Check to see if a sequence is specified
			if len(parts) > 4:
				domain_sequence = parts[4]
			else:
				domain_sequence = None
			
			# Create the new domains
			new_dom = Domain(domain_name, domain_length, 
							 sequence=domain_sequence)
			new_dom_comp = Domain(domain_name, domain_length, 
								  sequence=domain_sequence,
								  is_complement=True)
			
			domains[domain_name] = new_dom
			domains["%s*" % domain_name] = new_dom_comp
		
		# This is the beginning of a strand	
		elif line.startswith("strand"):
			# e.g.: 
			#       "strand A : a x b y z* c* y* b* x*"
			# parts: 0      1 2 3 4 5 6 ...
			
			parts = line.split()
			
			strand_name = parts[1]
			if strand_name in strands:
				logging.error("Duplicate strand name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not strand_name.isalnum():
				logging.warn("Non-alphanumeric strand name %s encountered in input line %d"
								% (strand_name, line_counter))
			
			strand_doms = []
			for domain_name in parts[3:]:
				if not domain_name in domains:
					logging.error("Invalid domain name %s encountered in input line %d"
									% (domain_name, line_counter))
					raise Exception()
								
				strand_doms.append(domains[domain_name])
			
			if len(strand_doms) == 0:
				logging.warn("0-length strand encountered in input line %d")			
				
			new_strand = Strand(strand_name, strand_doms)
			strands[strand_name] = new_strand
			
		# This is the beginning of a complex	
		elif line.startswith("complex"):
			# e.g.:
			# complex A :\n
			# A\n				<- strands_line
			# .(((..)))\n		<- structure_line
			
			
			parts = line.split()
			
			complex_name = parts[1]
			if complex_name in complexes:
				logging.error("Duplicate complex name encountered in input line %d"
								% line_counter)
				raise Exception()
			
			if not complex_name.isalnum():
				logging.warn("Non-alphanumeric complex name %s encountered in input line %d"
								% (complex_name, line_counter))
			
			complex_strands = []
			
			strands_line = fin.readline()
			strands_line = strands_line.strip()
			strands_line_parts = strands_line.split()
			for strand_name in strands_line_parts:
				if not strand_name in strands:
					logging.error("Invalid strand name %s encountered in input line %d"
									% (strand_name, line_counter))
					raise Exception()
				else:
					complex_strands.append(strands[strand_name])
				
			structure_line = fin.readline()
			structure_line = structure_line.strip()
			
			complex_structure = parse_dot_paren(structure_line)
			
			
			struct_length = sum(map(len,complex_structure))
			domains_length = sum(map(len,complex_strands))
			if(struct_length != domains_length):
				logging.error("Complex %(name)s has %(doms)d domains but structure size %(struct)d."
								% {"name":complex_name,"doms":domains_length,"struct":struct_length})
				raise Exception()
			
			complex = Complex(complex_name, complex_strands, complex_structure)
			complexes[complex_name] = complex			
			
		else:
			logging.error("Unexpected characters encountered in input line %d"
							% line_counter)
			raise Exception()
		line = fin.readline()
		line_counter += 1
		
	domains = domains.values()
	strands = strands.values()
	complexes = complexes.values()
	
	enumerator = Enumerator(domains, strands, complexes)
	return enumerator

			
def load_json(filename):
	"""
	Loads a saved enumerator from a JSON output file at filename.
	"""
	
	fin = open(filename, 'r')
	saved = json.load(fin)
	saved_domains = saved['domains']
	
	domains = {}
	for saved_domain in saved_domains:
		if not 'sequence' in saved_domain:
			saved_domain['sequence'] = None
		if (saved_domain['is_complement']):
			saved_domain['name'] = saved_domain['name'][:-1]
		new_dom = Domain(saved_domain['name'], saved_domain['length'], is_complement=saved_domain['is_complement'], sequence=saved_domain['sequence'])
		domains[new_dom.name] = new_dom
	
	saved_strands = saved['strands']
	
	strands = {}
	for saved_strand in saved_strands:
		doms = []
		for domain in saved_strand['domains']:
			doms.append(domains[domain])
		
		new_strand = Strand(saved_strand['name'], doms)
		strands[saved_strand['name']] = new_strand
	
	complexes = {}
	resting_complexes = {}
	
	saved_resting_complexes = saved['resting_complexes']
	for saved_complex in saved_resting_complexes:
		c_strands = []
		for strand in saved_complex['strands']:
			c_strands.append(strands[strand])
		new_structure = []
		for strand in saved_complex['structure']:
			new_strand = []
			for tup in strand:
				if (tup == None):
					new_strand.append(None)
				else:
					new_strand.append(tuple(tup))
			new_structure.append(new_strand)
		new_complex = Complex(saved_complex['name'], c_strands, new_structure)
		resting_complexes[saved_complex['name']] = new_complex
		complexes[saved_complex['name']] = new_complex
	
	transient_complexes = {}	
	
	saved_transient_complexes = saved['transient_complexes']
	for saved_complex in saved_transient_complexes:
		c_strands = []
		for strand in saved_complex['strands']:
			c_strands.append(strands[strand])
		new_structure = []
		for strand in saved_complex['structure']:
			new_strand = []
			for tup in strand:
				if (tup == None):
					new_strand.append(None)
				else:
					new_strand.append(tuple(tup))
			new_structure.append(new_strand)
		new_complex = Complex(saved_complex['name'], c_strands, new_structure)
		transient_complexes[saved_complex['name']] = new_complex
		complexes[saved_complex['name']] = new_complex
	
	saved_reactions = saved['reactions']
	
	reactions = []
	
	for saved_reaction in saved_reactions:
		reactants = []
		for reactant in saved_reaction['reactants']:
			reactants.append(complexes[reactant])
		products = []
		for product in saved_reaction['products']:
			products.append(complexes[product])
			
		reaction = ReactionPathway(saved_reaction['name'], reactants, products)
		reactions.append(reaction)
	
	resting_states = []
	
	for resting_state in saved['resting_states']:
		comps = []
		for complex in resting_state['complexes']:
			comps.append(complexes[complex])
		resting_states.append(RestingState(resting_state['name'], comps))

	initial_complexes = {}
	for saved_complex in saved['initial_complexes']:
		c_strands = []
		for strand in saved_complex['strands']:
			c_strands.append(strands[strand])
		new_structure = []
		for strand in saved_complex['structure']:
			new_strand = []
			for tup in strand:
				if (tup == None):
					new_strand.append(None)
				else:
					new_strand.append(tuple(tup))
			new_structure.append(new_strand)
		new_complex = Complex(saved_complex['name'], c_strands, new_structure)
		initial_complexes[saved_complex['name']] = new_complex
		
			
	enumerator = Enumerator(domains.values(), strands.values(), initial_complexes.values())
	enumerator._complexes = complexes.values()
	enumerator._resting_states = resting_states
	enumerator._transient_complexes = transient_complexes.values()
	enumerator._resting_complexes = resting_complexes.values()
	enumerator._reactions = reactions

	
	return enumerator



new_input_functions = {
						'standard': input_standard
					  }
					  
load_input_functions = {
						 'json': load_json
					   }