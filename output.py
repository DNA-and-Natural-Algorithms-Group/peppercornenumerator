#
#  output.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/21/10.
#

import copy
import utils
from reactions import ReactionPathway, auto_name
import reactions
import json
import subprocess
import sets



def condense_resting_states(enumerator):
	"""
	Outputs a set of 'resting state complexes' and reactions between them
	which represent a version of the reaction graph in which transient
	complexes have been removed and resting states are represented by a single
	node.
	
	More specifically, for every complex C, if C is a transient complex,
	for every reaction R1 containing C, we first ensure that C only appears on 
	one side of the reaction by normalizing the reaction. Then, if C is on
	the left hand side of the reaction, for every reaction R2 producing C,
	we create a new reaction with reagent set equal to the union of the reagents
	of R1 and R2 with C removed, and a new reaction with product set equal to 
	the union of the products of R1 and R2 (Note that if more than one C is
	required, then we multiply reaction R2 by the number of Cs required).
	
	We perform a similar operation for reactions R1 with C on the right hand
	side.
	"""
	complexes = enumerator.complexes
	reactions = enumerator.reactions[:]
	resting_states = enumerator.resting_states
	
	new_complexes = []
	new_reactions = []
	
	print "Reactions: ", len(reactions)
	
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
		
		reactions_copy = reactions[:]
		# We loop through the reactions and segment them into those that have
		# the complex as either a reagent or a product (can only be at most
		# one because the reactions are normalized)
		for reaction in reactions:
			reaction.normalize()
			if complex in reaction.reactants:
				reagent_reactions.append(reaction)
				reactions_copy.remove(reaction)
			elif complex in reaction.products:
				product_reactions.append(reaction)
				reactions_copy.remove(reaction)

		reactions = reactions_copy
		

	
		
		# We now loop over the reactions with complex as a reagent
		for reaction in reagent_reactions:
			# Count how many times complex appears
			counter = 0
			for reagent in reaction.reactants:
				if reagent == complex:
					counter += 1
			
			# We need to make a new reaction for every reaction producing complex
			for reaction2 in product_reactions:
				new_react_reagents = []
				new_react_reagents.extend(reaction.reactants)
				new_react_products = []
				new_react_products.extend(reaction.products)
				
				for i in range(counter):
					new_react_reagents.extend(reaction2.reactants)
					new_react_products.extend(reaction2.products)
					new_react_products.remove(complex)
					new_react_reagents.remove(complex)
				
				new_react = ReactionPathway(str(reactions.auto_name), new_react_reagents,
										new_react_products)
				new_reactions.append(new_react)
			# Now that this reaction has been replaced, we remove it from the
			# reaction list
			#new_reactions.remove(reaction)
			

		# Now we do the same for the reactions with complex as a product
		for reaction in product_reactions:
			# Count how many times complex appears
			counter = 0
			for reagent in reaction.products:
				if reagent == complex:
					counter += 1
			
			# We need to make a new reaction for every reaction taking complex
			for reaction2 in reagent_reactions:
				new_react_reagents = []
				new_react_reagents.extend(reaction.reactants)
				new_react_products = []
				new_react_products.extend(reaction.products)				
				for i in range(counter):
					new_react_reagents.extend(reaction2.reactants)
					new_react_products.extend(reaction2.products)
					new_react_products.remove(complex)
					new_react_reagents.remove(complex)
				
				new_react = ReactionPathway(str(reactions.auto_name), new_react_reagents,
										new_react_products)
				new_reactions.append(new_react)
			# Now that this reaction has been replaced, we remove it from the
			# reaction list
			#new_reactions.remove(reaction)

		#reactions = new_reactions[:]
		reactions.extend(new_reactions)
		new_reactions = []
		
		for reaction in reactions:
			if reaction in reagent_reactions:
				print "error2"
			if reaction in product_reactions:
				print "error3"
			if complex in reaction.reactants:
				print 'error!'
			if complex in reaction.products:
				print 'error!'
					
					
	new_reactions = sets.Set()
	for reaction in reactions:
		if (len(reaction.reactants) == 0) and (len(reaction.products) == 0):
			continue
		new_reagents = []
		for reagent in reaction.reactants:
			new_reagents.append(reagent._resting_state)
			
		new_products = []
		for product in reaction.products:			
			new_products.append(product._resting_state)
			
		new_reactions.add(ReactionPathway('condensed', new_reagents, new_products))


	print "New reactions: ", len(new_reactions)

	for reaction in new_reactions:
		print reaction, reaction.reactants, reaction.products
	output = {
				'resting_states': resting_states[:],
				'reactions': list(new_reactions)
			 }
	assert False
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
	for domain in sorted(enumerator.domains):
		output_file.write("sequence " + domain.name + " = : " + str(domain.length) + "\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Resting-state Complexes \n")
	for complex in sorted(resting_complexes):
		output_file.write(str(complex) + "\n")
		names = []
		for strand in complex.strands:
			names.append(strand.name)
		strands_string = " + ".join(names)
		output_file.write(strands_string + "\n")
		output_file.write(str(complex.dot_paren_string()) + "\n")
		output_file.write("\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Resting-state sets \n")
	for resting_state in sorted(resting_states):
		output_file.write("state " + resting_state.name + " = : ")
		output_file.write(str(resting_state.complexes[0]))
		for complex in resting_state.complexes[1:]:
			output_file.write(" + complex")
		output_file.write("\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Fast (Transition) Complexes \n")
	for complex in sorted(transient_complexes):
		output_file.write(str(complex) + "\n")
		names = []
		for strand in complex.strands:
			names.append(strand.name)
		strands_string = " + ".join(names)
		output_file.write(strands_string + "\n")
		output_file.write(str(complex.dot_paren_string()) + "\n")
		output_file.write("\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Reactions \n")
	for reaction in sorted(reactions):
		reactants = reaction.reactants
		products = reaction.products
		reac_string_list = [reactants[0].name]
		for reactant in reactants[1:]:
			reac_string_list.append(" + " + reactant.name)
		reac_string_list.append(" -> ")
		reac_string_list.append(products[0].name)
		for product in products[1:]:
			reac_string_list.append(" + " + product.name)
		reac_string_list.append("\n")
		reac_string = ''.join(reac_string_list)
		output_file.write(reac_string)
	output_file.write("###############################\n")
	if (output_condensed):
		output_file.write("\n\n# Condensed Reactions \n")
		condensed = condense_resting_states(enumerator)
		new_reactions = condensed['reactions']
		for reaction in sorted(new_reactions):
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
		
	object_out['resting_complexes'] = resting_complexes_out
		
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
		temp_resting_state['complexes'] = temp_complexes
		resting_states_out.append(temp_resting_state)
		
	object_out['resting_states'] = resting_states_out

	initial_complexes = enumerator.initial_complexes
	initial_complexes_out = []
	for complex in initial_complexes:
		temp_complex = {}
		temp_complex['name'] = complex.name
		temp_strands = []
		for strand in complex.strands:
			temp_strands.append(strand.name)
		temp_complex['strands'] = temp_strands
		temp_complex['structure'] = complex.structure
		initial_complexes_out.append(temp_complex)

	object_out['initial_complexes'] = initial_complexes_out

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
		for product in reaction.products:
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

def output_graph(enumerator, filename, output_condensed=False):
	if not output_condensed:
		output_full_graph(enumerator, filename)
	else:
		output_condensed_graph(enumerator, filename)
	
output_graph.supports_condensed = True

def output_full_graph(enumerator, filename):
	"""
	Products graphical output representing the full reaction graph, with all
	reactions and complexes. Transient and resting states are differentiated
	by color.
	"""
	fout = open(filename + ".dot", "w")
	fout.write("digraph G {\n")
	fout.write('size="7,10"\n')
	fout.write('page="8.5,11"\n')
	fout.write('node[width=0.25,height=0.375,fontsize=9]\n')
	
	# We need to determine the complex clusters for the graph. Complexes with
	# the same cyclic permutation of strands are placed in the same cluster.
	
	strand_cyclic_permutations = []
	clusters = []
	
	# We loop through all complexes
	for complex in enumerator.transient_complexes:
		complex._resting = False
		flag = False
		
		# Check to see if we've already seen this strand ordering
		for (i, perm) in enumerate(strand_cyclic_permutations):
			if perm == complex.strands:
				clusters[i].append(complex)
				flag = True
				break
				
		# If not, we add it
		if not flag:
			strand_cyclic_permutations.append(complex.strands)
			clusters.append([complex])
			
	for complex in enumerator.resting_complexes:
		complex._resting = True
		flag = False
		
		# Check to see if we've already seen this strand ordering
		for (i, perm) in enumerate(strand_cyclic_permutations):
			if perm == complex.strands:
				clusters[i].append(complex)
				flag = True
				break
				
		# If not, we add it
		if not flag:
			strand_cyclic_permutations.append(complex.strands)
			clusters.append([complex])
	
	# We now draw the clusters on the graph
	for (i, cluster) in enumerate(clusters):
		fout.write("subgraph cluster%d {\n" % i)
		strands = [cluster[0].strands[0].name]
		for strand in cluster[0].strands[1:]:
			strands.append(" + ")
			strands.append(strand.name)
		strands_string = ''.join(strands)
		fout.write('label="%s"\n' % strands_string)
		fout.write('fontsize=6\n')
		for complex in cluster:
			extra_params = ""
			if complex._resting:
				extra_params = ",style=filled,color=gold1"
			fout.write('%s [label="%s:%s"%s];\n' % (complex.name, 
													complex.name, 
													complex.dot_paren_string(), 
													extra_params
													))
		fout.write("}\n")
		
	# We now draw the reactions. If there is one product and one reagent, then
	# we just draw an edge between the two. Otherwise we create a reaction node.
	for (i, reaction) in enumerate(enumerator.reactions):
		if (len(reaction.products) == 1) and (len(reaction.reactants) == 1):
			fout.write("%s -> %s\n" % (reaction.reactants[0].name,
									   reaction.products[0].name))
		else:
			reaction_label = "R_%d" % i
			# We label unimolecular reactions blue, and other reactions red
			reaction_color = "red"
			if len(reaction.reactants) == 1:
				reaction_color = "blue"
			fout.write('%s [label="",shape=circle,height=0.12,width=0.12,fontsize=1,style=filled,color=%s];\n' %
							(reaction_label, reaction_color))
			
			# We now make all the edges needed
			for reagent in reaction.reactants:
				fout.write("%s -> %s\n" % (reagent.name, reaction_label))
				
			for product in reaction.products:
				fout.write("%s -> %s\n" % (reaction_label, product.name))
				
				
	fout.write("}\n")
	fout.close()
	
	# Create the output file.
	# TODO: make 'pdf' configurable
	subprocess.call(["dot", "-O", "-Teps", "%s.dot" % filename])
	
def output_condensed_graph(enumerator, filename):
	"""
	Products graphical output representing the condensed reaction graph, with
	condensed reactions and resting states aggregated into single nodes.
	"""
	fout = open(filename + ".dot", "w")
	fout.write("digraph G {\n")
	fout.write('size="7,10"\n')
	fout.write('page="8.5,11"\n')
	fout.write('node[width=0.25,height=0.375,fontsize=9]\n')
	
	condensed_graph = condense_resting_states(enumerator)
	
	# We loop through all resting states, drawing them on the graph
	for state in condensed_graph['resting_states']:
		fout.write('%s [label="%s"]\n' % (state.name, state.name))
	
	# We now add reactions. We always create a reaction node.
	for (i, reaction) in enumerate(condensed_graph['reactions']):
		reaction_label = "R_%d" % i
		fout.write('%s [label="",shape=circle,height=0.12,width=0.12,fontsize=1,style=filled,color=red];\n' % reaction_label)
		
		for reagent in reaction.reagents:
			fout.write("%s -> %s\n" % (reagent.name, reaction_label))
			
		for product in reaction.products:
			fout.write("%s -> %s\n" % (reaction_label, product.name))
			
	fout.write("}\n")
	fout.close()
	
	# Create the output file.
	# TODO: make 'pdf' configurable
	subprocess.call(["dot", "-O", "-Tpdf", "%s.dot" % filename])
	

	
text_output_functions = {
	'legacy': output_legacy,
	'json': output_json
}

graph_output_functions = {
	'graph': output_graph
}