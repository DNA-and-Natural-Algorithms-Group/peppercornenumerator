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
from condense import condense_resting_states



def output_legacy(enumerator, filename, output_condensed = False):
	"""
	Legacy text-based output scheme similar to that used in the original 
	enumerator. Designed to be simultaneously parsable and human-readable.
	Supports output of condensed graph in addition to the full graph. Does
	not support strands.
	"""
	
	def write_complex(output_file,complex):
		output_file.write(str(complex) + "\n")
		names = []
		for strand in complex.strands:
			names.append(" ".join(map(str,strand.domains)))
		strands_string = " + ".join(names)
		output_file.write(strands_string + "\n")
		output_file.write(str(complex.dot_paren_string()) + "\n")
		output_file.write("\n")
	
	def write_reaction(output_file,reaction):
		reactants = map(str,reaction.reactants)
		products = map(str,reaction.products)
		reac_string_list = [" + ".join(reactants),"->"," + ".join(products),"\n"]
		reac_string = ' '.join(reac_string_list)
		output_file.write(reac_string)
	
#	def write_reaction(output_file,reaction):
#		reactants = reaction.reactants
#		products = reaction.products
#		reac_string_list = [reactants[0].name]
#		for reactant in reactants[1:]:
#			reac_string_list.append(" + " + reactant.name)
#		reac_string_list.append(" -> ")
#		reac_string_list.append(products[0].name)
#		for product in products[1:]:
#			reac_string_list.append(" + " + product.name)
#		reac_string_list.append("\n")
#		reac_string = ''.join(reac_string_list)
#		output_file.write(reac_string)
		
	complexes = enumerator.complexes
	transient_complexes = enumerator.transient_complexes
	resting_complexes = enumerator.resting_complexes
	reactions = enumerator.reactions
	resting_states = enumerator.resting_states
	
	output_file = open(filename, 'w')
	output_file.write("###### Enumerated Output ######\n")
	output_file.write("\n\n# Domains \n")
	for domain in sorted(enumerator.domains):
		if(not domain.is_complement):
			output_file.write("sequence " + domain.name + " = : " + str(domain.length) + "\n")
	
	output_file.write("\n\n# End-state Complexes \n")
	for complex in sorted(resting_complexes):
		write_complex(output_file,complex)
		
	# Not part of the original output; omitting
#	output_file.write("###############################\n")
#	output_file.write("\n\n# Resting-state sets \n")
#	for resting_state in sorted(resting_states):
#		output_file.write("state " + resting_state.name + " = : ")
#		output_file.write(str(resting_state.complexes[0]))
#		for complex in resting_state.complexes[1:]:
#			output_file.write(" %s" % complex)
#		output_file.write("\n")
	output_file.write("###############################\n")
	output_file.write("\n\n# Fast (Transition) Complexes \n")
	for complex in sorted(transient_complexes):
		write_complex(output_file,complex)
	output_file.write("###############################\n")
	output_file.write("\n\n# Reactions \n")
	for reaction in sorted(reactions):
		write_reaction(output_file,reaction)
		
	if (output_condensed):
		output_file.write("###############################\n")
		output_file.write("\n\n# Condensed Reactions \n")
		condensed = condense_resting_states(enumerator)
		new_reactions = condensed['reactions']
		for reaction in sorted(new_reactions):
			write_reaction(output_file,reaction)
				
	output_file.close()

output_legacy.supports_condensed = True		

def output_pil(enumerator, filename, output_condensed = False):
	"""
	Text-based output using the Pepper Intermediate Language (PIL)
	"""
	
	def write_complex(output_file,complex):
		output_file.write("structure " + str(complex) + " = ")
		names = map(lambda strand: strand.name, complex.strands)
		strands_string = " + ".join(names)
		output_file.write(strands_string + " : ")
		output_file.write(str(complex.dot_paren_string()) + "\n")
	
	def write_reaction(output_file,reaction):
		reactants = map(str,reaction.reactants)
		products = map(str,reaction.products)
		reac_string_list = ["kinetic"," + ".join(reactants),"->"," + ".join(products),"\n"]
		reac_string = ' '.join(reac_string_list)
		output_file.write(reac_string)
		
	complexes = enumerator.complexes
	transient_complexes = enumerator.transient_complexes
	resting_complexes = enumerator.resting_complexes
	reactions = enumerator.reactions
	resting_states = enumerator.resting_states
	
	output_file = open(filename, 'w')
	output_file.write("###### Enumerated Output ######\n")
	output_file.write("\n# Domains \n")
	
	def seq(dom):
		if(dom.sequence != None):
			return dom.sequence
		else:
			return "N" * len(dom)
	
	for domain in utils.natural_sort(enumerator.domains):
		if(not domain.is_complement):
			output_file.write("sequence " + domain.name + " = " + seq(domain) + " : " + str(domain.length) + "\n")
	
	output_file.write("\n# Strands \n")
	for strand in utils.natural_sort(enumerator.strands):
		output_file.write("strand " + strand.name + " = " + \
						" ".join(map(lambda dom: dom.name, strand.domains)) + "\n")
	
	output_file.write("\n# Resting-state Complexes \n")
	for complex in utils.natural_sort(resting_complexes):
		write_complex(output_file,complex)
		
	

	
	if (output_condensed):
		condensed = condense_resting_states(enumerator)

		output_file.write("\n# Resting-state sets \n")
		resting_states = condensed['resting_states']
		for resting_state in utils.natural_sort(resting_states):
			output_file.write("# state " + str(resting_state) + " = { " + " ".join(map(str,resting_state.complexes)) + " }\n")

		output_file.write("\n# Condensed Reactions \n")
		new_reactions = condensed['reactions']
		for reaction in sorted(new_reactions):
			write_reaction(output_file,reaction)
	else:	
		output_file.write("\n# Transient Complexes \n")
		for complex in utils.natural_sort(transient_complexes):
			write_complex(output_file,complex)
		output_file.write("\n# Detailed Reactions \n")
		for reaction in sorted(reactions): #utils.natural_sort(reactions):
			write_reaction(output_file,reaction)
		
			
	output_file.close()


def output_json(enumerator, filename, output_condensed = False):
	"""
	JSON-based output schema intended to be easily machine parsable. Uses
	python's JSON serialization libraries.
	"""
	
	def serializeComplex(complex):
		temp_strands = []
		for strand in complex.strands:
			temp_strands.append(strand.name)
		return { 
			'name':complex.name,
			'strands': temp_strands,
			'structure': complex.structure,
			'dot-paren': complex.dot_paren_string(),
			'dot-paren-full': complex.dot_paren_string_full(),
		}
	
	def serializeReaction(reaction):
		temp_reactants = []
		for reactant in reaction.reactants:
			temp_reactants.append(reactant.name)
		temp_products = []
		for product in reaction.products:
			temp_products.append(product.name)
		return {
			"name":reaction.name,
			"reactants":temp_reactants,
			"products":temp_products
		}
	
	def serializeDomain(domain):
		temp_domain = {}
		temp_domain['name'] = domain.name
		temp_domain['length'] = domain.length
		temp_domain['is_complement'] = domain.is_complement
		if domain.sequence != None:
			temp_domain['sequence'] = domain.sequence
		return temp_domain
		
	def serializeStrand(strand):
		temp_strand = {}
		temp_strand['name'] = strand.name
		temp_domains = []
		for domain in strand.domains:
			temp_domains.append(domain.name)
		temp_strand['domains'] = temp_domains
		return temp_strand
		
	def serializeRestingState(resting_state):
		temp_resting_state = {}
		temp_complexes = []
		for complex in resting_state.complexes:
			temp_complexes.append(complex.name)
		temp_resting_state['name'] = resting_state.name
		temp_resting_state['complexes'] = temp_complexes
		return temp_resting_state
		
	object_out = {}
	
	object_out['domains'] = map(serializeDomain,enumerator.domains)
	object_out['strands'] = map(serializeStrand,enumerator.strands)
	object_out['resting_complexes'] = map(serializeComplex,enumerator.resting_complexes)
	object_out['transient_complexes'] = map(serializeComplex,enumerator.transient_complexes)
	object_out['resting_states'] = map(serializeRestingState,enumerator.resting_states)
	object_out['initial_complexes'] = map(serializeComplex,enumerator.initial_complexes)
	object_out['reactions'] = map(serializeReaction,enumerator.reactions)
	
	if output_condensed:
		condensed = condense_resting_states(enumerator)
		object_out['condensed_reactions'] = map(serializeReaction,condensed['reactions'])
		
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
			fout.write('%s [label="%s: %s"%s];\n' % (str(complex), 
													str(complex), 
													complex.dot_paren_string(), 
													extra_params
													))
		fout.write("}\n")
		
	# We now draw the reactions. If there is one product and one reagent, then
	# we just draw an edge between the two. Otherwise we create a reaction node.
	for (i, reaction) in enumerate(enumerator.reactions):
		if (len(reaction.products) == 1) and (len(reaction.reactants) == 1):
			fout.write("%s -> %s\n" % (str(reaction.reactants[0]),
									   str(reaction.products[0])))
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
				fout.write("%s -> %s\n" % (str(reagent), reaction_label))
				
			for product in reaction.products:
				fout.write("%s -> %s\n" % (reaction_label, str(product)))
				
				
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
		fout.write('%s [label="%s"]\n' % (str(state), str(state)))
	
	# We now add reactions. We always create a reaction node.
	for (i, reaction) in enumerate(condensed_graph['reactions']):
		reaction_label = "R_%d" % i
		fout.write('%s [label="",shape=circle,height=0.12,width=0.12,fontsize=1,style=filled,color=red];\n' % reaction_label)
		
		for reagent in reaction.reactants:
			fout.write("%s -> %s\n" % (str(reagent), reaction_label))
			
		for product in reaction.products:
			fout.write("%s -> %s\n" % (reaction_label, str(product)))
			
	fout.write("}\n")
	fout.close()
	
	# Create the output file.
	# TODO: make 'pdf' configurable
	subprocess.call(["dot", "-O", "-Teps", "%s.dot" % filename])
	

def output_sbml(enumerator,filename, output_condensed = False):
	import xml.dom.minidom
	header = '<?xml version="1.0" encoding="UTF-8"?>'
	out = [header,
		'<sbml level="2" version="3" xmlns="http://www.sbml.org/sbml/level2/version3">',
		'<model name="%s">' % filename,
		'<listOfCompartments>',
		'<compartment id="reaction" />',
		'</listOfCompartments>',
		'<listOfSpecies>']
	
	if(output_condensed):
		condensed = condense_resting_states(enumerator)
		complexes = condensed['resting_states']
		reactions = condensed['reactions']
	else:
		complexes = enumerator.complexes
		reactions = enumerator.reactions
	
	for complex in complexes:
		out.append('<species compartment="reaction" id="s_%(name)s" name="%(name)s"/>' % {"name": complex.name})
	out.extend(['</listOfSpecies>','<listOfReactions>']);
	for (i, reaction) in enumerate(reactions):
		out.extend(['<reaction id="r_%d">' % i,
                '<listOfReactants>'])
		for species in reaction.reactants:
			out.append('<speciesReference species="%s"/>' % species.name)
		out.extend(['</listOfReactants>',
                '<listOfProducts>'])
		for species in reaction.products:
			out.append('<speciesReference species="%s"/>' % species.name)	
		out.extend(['</listOfProducts>','</reaction>'])

	out.extend(['</listOfReactions>','</model>','</sbml>']);

	
	doc = xml.dom.minidom.parseString("".join(out))
	
	fout = open(filename, "w")
	fout.write(header+'\n'+doc.documentElement.toprettyxml(indent="\t"))
	fout.close()
	
text_output_functions = {
	'standard': output_legacy,
	'legacy': output_legacy,
	'pil': output_pil,
	'json': output_json,
	'sbml': output_sbml
}

graph_output_functions = {
	'graph': output_graph
}