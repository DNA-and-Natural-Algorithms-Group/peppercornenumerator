#
#  output.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/21/10.
#

import copy
import utils
import reactions
import json


text_output_functions = {
						'legacy': output_legacy,
						'json': output_json
   					    }
graph_output_functions = {

						 }

def condense_resting_states(enumerator):
	"""
	Outputs a set of 'resting state complexes' and reactions between them
	which represent a version of the reaction graph in which transient
	complexes have been removed and resting states are represented by a single
	node.
	"""
	complexes = enumerator.complexes
	reactions = enumerator.reactions
	resting_states = enumerator.resting_states
	
	new_complexes = []
	new_reactions = reactions[:]
	
	
	# First we need to mark resting state complexes
	for complex in complexes:
		complex._is_resting_state = False
	
	for resting_state in resting_states:
		for complex in resting_state.complexes:
			complex._is_resting_state = True
			complex._resting_state = resting_state
			new_complexes.append(complex)
			
	
	# Now we iterate through the complexes and remove the transient states
	for complex in complexes:
	
		# We only want to eliminate transient states
		if (complex._is_resting_state):
			continue
			
		reagent_reactions = []
		product_reactions = []
		
		# We loop through the reactions and segment them into those that have
		# the complex as either a reagent or a product (can only be at most
		# one because the reactions are normalized)
		for reaction in reactions:
			reaction.normalize()
			if complex in reaction.reagents:
				reagent_reactions.append(reaction)
			else if complex in reaction.products:
				product_reactions.append(reaction)
		
		
		# We now loop over the reactions with complex as a reagent
		for reaction in reagent_reactions:
			# Count how many times complex appears
			counter = 0
			for reagent in reaction.reagents:
				if reagent == complex:
					counter += 1
			
			# We need to make a new reaction for every reaction producing complex
			for reaction2 in product_reactions:
				new_react_reagents = []
				new_react_reagents.extend(reaction.reagents)
				new_react_products = []
				new_react_products.extend(reaction.products)
				
				for i in range(counter):
					new_react_reagents.remove(complex)
					new_react_reagents.extend(reaction2.reagents)
					new_react_products.extend(reaction2.products)
		
			new_react = ReactionPathway(str(auto_name), new_react_reagents,
										new_react_products)
			new_reactions.append(new_react)
			# Now that this reaction has been replaced, we remove it from the
			# reaction list
			new_reactions.remove(reaction)
		
		# Now we do the same for the reactions with complex as a product
		for reaction in product_reactions:
			# Count how many times complex appears
			counter = 0
			for reagent in reaction.products:
				if reagent == complex:
					counter += 1
			
		# We now loop over the reactions with complex as a reagent
		for reaction in reagent_reactions:
			# Count how many times complex appears
			counter = 0
			for reagent in reaction.reagents:
				if reagent == complex:
					counter += 1
			
			# We need to make a new reaction for every reaction producing complex
			for reaction2 in product_reactions:
				new_react_reagents = []
				new_react_reagents.extend(reaction.reagents)
				new_react_products = []
				new_react_products.extend(reaction.products)
				
				for i in range(counter):
					new_react_reagents.remove(complex)
					new_react_reagents.extend(reaction2.reagents)
					new_react_products.extend(reaction2.products)
		
			new_react = ReactionPathway(str(auto_name), new_react_reagents,
										new_react_products)
			new_reactions.append(new_react)
			# Now that this reaction has been replaced, we remove it from the
			# reaction list
			new_reactions.remove(reaction)
			
		reactions = new_reactions[:]
		new_reactions = []
	
	new_reactions = []
	for reaction in reactions:
		new_reagents = []
		for reagent in reaction.reagents:
			new_reagents.append(reagent._resting_state)
			
		new_products = []
		for product in reaction.products:
			new_products.append(product._resting_state)
			
		new_reactions.append(ReactionPathway(str(auto_name), new_reagents, new_products)
		
	output = {
				'resting_states': resting_states[:],
				'reactions': new_reactions
			 }
	
	return output


def output_legacy(enumerator, filename, output_condensed = False):
	"""
	Legacy text-based output scheme similar to that used in the original 
	enumerator. Designed to be simultaneously parsable and human-readable.
	Supports output of condensed graph in addition to the full graph. Does
	not support strands.
	"""
	complexes = enumerator.complexes
	transient_complexes = enumerator.transient_complexes
	resting_complexes = enumerator.resting_complexes
	reactions = enumerator.reactions
	resting_states = enumerator.resting_states
	
	output_file = open(filename, 'w')
	output_file.write("###### Enumerated Output ######\n")
	output_file.write("\n\n# Domains \n")
	for domain in enumerator.domains:
		output_file.write("sequence " + domain.name + " = : " + str(domain.length) + "\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Resting-state Complexes \n")
	for complex in resting_complexes:
		output_file.write(str(complex) + "\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Resting-state sets \n")
	for resting_state in resting_states:
		output_file.write("state " + resting_state.name + " = : ")
		output_file.write(resting_state.complexes[0])
		for complex in resting_state.complexes[1:]:
			output_file.write(" + complex")
		output_file.write("\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Fast (Transition) Complexes \n")
	for complex in transient_complexes:
		output_file.write(str(comp) + "\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Reactions \n")
	for reaction in reactions:
		reactants = reaction.reactants
		products = reaction.products
		reac_string_list = [reactants[0].name]
		for reactant in reactants[1:]:
			reac_string_list.append(" + " + reactant.name)
		reac_string_list.append(" -> ")
		reac_string_list.append(products[0].name)
		for product in products:
			reac_string_list.append(" + " + product.name)
		reac_string_list.append("\n")
		reac_string = ''.join(reac_string_list)
		output_file.write(reac_string)
	output_file.write("###############################\n")
	if (output_condensed):
		output_file.write("\n\n# Condensed Reactions \n")
		condensed = condense_resting_states(enumerator)
		new_reactions = condensed['reactions']
		for reaction in new_reactions:
			reactants = reaction.reactants
			products = reaction.products
			reac_string_list = [reactants[0].name]
			for reactant in reactants[1:]:
				reac_string_list.append(" + " + reactant.name)
			reac_string_list.append(" -> ")
			reac_string_list.append(products[0].name)
			for product in products:
				reac_string_list.append(" + " + product.name)
			reac_string_list.append("\n")
			reac_string = ''.join(reac_string_list)
			output_file.write(reac_string)
	output_file.write("###############################\n")
	output_file.close()

output_legacy.supports_condensed = True		

def output_json(enumerator, filename, output_condensed = False):
	"""
	JSON-based output schema intended to be easily machine parsable. Uses
	python's JSON serialization libraries.
	"""
	
	object_out = {}
	
	domains = enumerator.domains
	domains_out = []
	for domain in domains:
		temp_domain = {}
		temp_domain['name'] = domain.name
		temp_domain['length'] = domain.length
		temp_domain['is_complement'] = domain.is_complement
		if domain.sequence != None:
			temp_domain['sequence'] = domain.sequence
		domains_out.append(temp_domain)
		
	object_out['domains'] = domains_out
		
	strands = enumerator.strands
	strands_out = []
	for strand in strands:
		temp_strand = {}
		temp_strand['name'] = strand.name
		temp_domains = []
		for domain in strand.domains:
			temp_domains.append(domain.name)
		temp_strand['domains'] = temp_domains
		strands_out.append(temp_strand)
		
	object_out['strands'] = strands_out
		
	resting_complexes = enumerator.resting_complexes
	resting_complexes_out = []
	for complex in resting_complexes:
		temp_complex = {}
		temp_complex['name'] = complex.name
		temp_strands = []
		for strand in complex.strands:
			temp_strands.append(strand.name)
		temp_complex['strands'] = temp_strands
		temp_complex['structure'] = complex.structure
		resting_complexes_out.append(temp_complex)
		
	object_out['resting_complexes'] = resting_complexes
		
	transient_complexes = enumerator.transient_complexes
	transient_complexes_out = []
	for complex in transient_complexes:
		temp_complex = {}
		temp_complex['name'] = complex.name
		temp_strands = []
		for strand in complex.strands:
			temp_strands.append(strand.name)
		temp_complex['strands'] = temp_strands
		temp_complex['structure'] = complex.structure
		transient_complexes_out.append(temp_complex)
	
	object_out['transient_complexes'] = transient_complexes_out
	
	resting_states = enumerator.resting_states
	resting_states_out = []
	for resting_state in resting_states:
		temp_resting_state = {}
		temp_complexes = []
		for complex in resting_state.complexes:
			temp_complexes.append(complex.name)
		temp_resting_state['name'] = resting_state.name
		temp_resting_state['complexes'] = resting_state.complexes
		resting_states_out.append(temp_resting_state)
		
	object_out['resting_states'] = resting_states_out
		
	reactions = enumerator.reactions
	reactions_out = []
	for reaction in reactions:
		temp_reaction = {}
		temp_reaction['name'] = reaction.name
		temp_reactants = []
		for reactant in reaction.reactants:
			temp_reactants.append(reactant.name)
		temp_reaction['reactants'] = temp_reactants
		temp_products = []
		for product in temp_products:
			temp_products.append(product.name)
		temp_reaction['products'] = temp_products
		reactions_out.append(temp_reaction)
		
	object_out['reactions'] = reactions_out
	
	if output_condensed:
		condensed = condense_resting_states(enumerator)
		reactions = condensed.reactions
		reactions_out = []
		for reaction in reactions:
			temp_reaction = {}
			temp_reaction['name'] = reaction.name
			temp_reactants = []
			for reactant in reaction.reactants:
				temp_reactants.append(reactant.name)
			temp_reaction['reactants'] = temp_reactants
			temp_products = []
			for product in temp_products:
				temp_products.append(product.name)
			temp_reaction['products'] = temp_products
			reactions_out.append(temp_reaction)
		object_out['condensed_reactions'] = reactions_out
		
	fout = open(filename, 'w')
	json.dump(object_out, fout, indent=4)
	fout.close()

output_json.supports_condensed = True