#
#  output.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/21/10.
#  modifications by Stefan Badelt.
#

import json
import subprocess

from peppercornenumerator import __version__
from peppercornenumerator.utils import natural_sort
from peppercornenumerator.condense import ReactionGraph

def write_kernel(enumerator, pil, condensed = False):
    """Write the contents of :obj:`Enumerator()` into a KERNEL file.

    Args:
      enumerator (:obj:`Enumerator()`): The enumertor object.
      pil (filehandle): The filehandle to which the output is written to.
      concensed (bool, optional): Print condensed form. Defaults to False.
    """

    pil.write("# File generated by peppercorn-{}\n".format(__version__))

    # Print Domains
    seen = set()
    pil.write("\n# Domain Specifications \n")
    for dom in natural_sort(enumerator.domains):
        if dom.is_complement : 
            dom = ~dom
        if dom not in seen :
            pil.write("length {:s} = {:d}\n".format(dom.name, len(dom)))
            seen.add(dom)

    # Print Resting-set Complexes
    pil.write("\n# Resting-set Complexes \n")
    for cplx in natural_sort(enumerator.resting_complexes):
        pil.write("{:s} = {:s}\n".format(cplx.name, cplx.kernel_string))
 
    if condensed :
        # Print Resting Sets
        enumCG = ReactionGraph(enumerator)
        enumCG.condense()
        pil.write("\n# Resting-sets \n")
        for resting in natural_sort(enumerator.resting_sets):
            pil.write("state {:s} = [{}]\n".format(resting, 
                ', '.join(map(str,resting.complexes))))

        # Print Reactions
        pil.write("\n# Condensed Reactions \n")
        for reaction in natural_sort(enumCG.condensed_reactions):
            pil.write("reaction {:s}\n".format(reaction.full_string))
    else :
        # Print Transient Complexes
        pil.write("\n# Transient Complexes \n")
        for cplx in natural_sort(enumerator.transient_complexes):
            pil.write("{:s} = {:s}\n".format(cplx.name, cplx.kernel_string))

        # Print Reactions
        pil.write("\n# Detailed Reactions \n")
        for reaction in natural_sort(enumerator.reactions):
            pil.write("reaction {:s}\n".format(reaction.full_string))

def output_graph(enumerator, filename, output_condensed=False,
                 output_rates=False):
    if not output_condensed:
        output_full_graph(enumerator, filename)
    else:
        output_condensed_graph(enumerator, filename)

def output_dotfile(enumerator, filename,
                   output_condensed=False, output_rates=False):
    if not output_condensed:
        output_full_dotfile(enumerator, filename)
    else:
        output_condensed_dotfile(enumerator, filename)

def output_full_dotfile(enumerator, filename):
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
            fout.write(
                '%s [label="%s: %s"%s];\n' %
                (str(complex),
                 str(complex),
                    complex.dot_paren_string(),
                    extra_params))
        fout.write("}\n")

    # We now draw the reactions. If there is one product and one reagent, then
    # we just draw an edge between the two. Otherwise we create a reaction
    # node.
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
            fout.write(
                '%s [label="",shape=circle,height=0.12,width=0.12,fontsize=1,style=filled,color=%s];\n' %
                (reaction_label, reaction_color))

            # We now make all the edges needed
            for reagent in reaction.reactants:
                fout.write("%s -> %s\n" % (str(reagent), reaction_label))

            for product in reaction.products:
                fout.write("%s -> %s\n" % (reaction_label, str(product)))

    fout.write("}\n")
    fout.close()

def output_full_graph(enumerator, filename):
    """
    Products graphical output representing the full reaction graph, with all
    reactions and complexes. Transient and resting states are differentiated
    by color.
    """

    output_full_dotfile(enumerator, filename)

    # Create the output file.
    # TODO: make 'pdf' configurable
    subprocess.call(["dot", "-O", "-Teps", "%s.dot" % filename])

def output_condensed_dotfile(enumerator, filename):
    fout = open(filename + ".dot", "w")
    fout.write("digraph G {\n")
    fout.write('size="7,10"\n')
    fout.write('page="8.5,11"\n')
    fout.write('node[width=0.25,height=0.375,fontsize=9]\n')

    condensed_graph = condense_resting_states(enumerator, **condense_options)

    # We loop through all resting states, drawing them on the graph
    for state in condensed_graph['resting_states']:
        fout.write('%s [label="%s"]\n' % (str(state), str(state)))

    # We now add reactions. We always create a reaction node.
    for (i, reaction) in enumerate(condensed_graph['reactions']):
        reaction_label = "R_%d" % i
        fout.write(
            '%s [label="",shape=circle,height=0.12,width=0.12,fontsize=1,style=filled,color=red];\n' %
            reaction_label)

        for reagent in reaction.reactants:
            fout.write("%s -> %s\n" % (str(reagent), reaction_label))

        for product in reaction.products:
            fout.write("%s -> %s\n" % (reaction_label, str(product)))

    fout.write("}\n")
    fout.close()

def output_condensed_graph(enumerator, filename):
    """
    Products graphical output representing the condensed reaction graph, with
    condensed reactions and resting states aggregated into single nodes.
    """

    output_condensed_dotfile(enumerator, filename)

    # Create the output file.
    # TODO: make 'pdf' configurable
    subprocess.call(["dot", "-O", "-Teps", "%s.dot" % filename])

def output_sbml(enumerator, filename, output_condensed=False,
                output_rates=True, condense_options={}):
    # # default initial concentration of all species is 100 nM
    # initial_concentration = 10e-7

    import xml.dom.minidom
    header = '<?xml version="1.0" encoding="UTF-8"?>'
    out = [
        header,
        '<sbml level="2" version="3" xmlns="http://www.sbml.org/sbml/level2/version3">',
        '<model name="%s">' %
        filename,
        '<listOfUnitDefinitions>',
        '<unitDefinition id="per_second">',
        '<listOfUnits>',
        '<unit kind="second" exponent="-1"/>',
        ' </listOfUnits>',
        '</unitDefinition>',
        '<unitDefinition id="litre_per_mole_per_second">',
        '<listOfUnits>',
        '<unit kind="mole"   exponent="-1"/>',
        '<unit kind="litre"  exponent="1"/>',
        '<unit kind="second" exponent="-1"/>',
        '</listOfUnits>',
        '</unitDefinition>',
        '</listOfUnitDefinitions>',
        '<listOfCompartments>',
        '<compartment id="reaction" size="1e-3" />',
        '</listOfCompartments>',
        '<listOfSpecies>']

    if(output_condensed):
        condensed = condense_resting_states(enumerator, **condense_options)
        complexes = condensed['resting_states']
        reactions = condensed['reactions']
    else:
        complexes = enumerator.complexes
        reactions = enumerator.reactions

    def id(species):
        return "s_" + species.name

    # build elements for each species
    if(output_condensed):
        for resting_state in complexes:
            # is_initial = any(c in enumerator.initial_complexes for c in resting_state.complexes)
            initial_concentration = sum(
                (c.concentration if c.concentration is not None else 0.0) for c in resting_state.complexes)
            out.append(
                '<species compartment="reaction" id="%(id)s" name="%(name)s" initialConcentration="%(initial)g"/>' %
                {
                    "name": resting_state.name,
                    "id": id(resting_state),
                    "initial": initial_concentration})
    else:
        for complex in complexes:
            # is_initial = (complex in enumerator.initial_complexes)
            initial_concentration = complex.concentration if complex.concentration is not None else 0.0
            out.append(
                '<species compartment="reaction" id="%(id)s" name="%(name)s" initialConcentration="%(initial)g"/>' %
                {
                    "name": complex.name,
                    "id": id(complex),
                    "initial": initial_concentration})

    out += ['</listOfSpecies>', '<listOfReactions>']

    # list reactions
    for (i, reaction) in enumerate(reactions):
        out += ['<reaction id="r_%d" reversible="false">' % i,
                '<listOfReactants>'] + \
            ['<speciesReference species="%s"/>' % id(species) for species in reaction.reactants] + \
            ['</listOfReactants>',
             '<listOfProducts>'] + \
            ['<speciesReference species="%s"/>' % id(species) for species in reaction.products] + \
            ['</listOfProducts>']

        # unimolecular rate constants have units 1/s, bimolecular rate
        # constants have units 1/M/s
        units = 'per_second' if reaction.arity[0] == 1 else 'litre_per_mole_per_second'

        out += ['<kineticLaw>',
                '<math xmlns="http://www.w3.org/1998/Math/MathML">',
                '<apply>',
                '<times />',
                '<ci>k</ci>'] + \
            ['<ci>' + id(s) + '</ci>' for s in reaction.reactants] + \
            ['</apply>',
             '</math>',
             '<listOfParameters>',
             '<parameter id="k"  value="%.10f" units="%s"/>' % (
                 reaction.rate, units),
             '</listOfParameters>',
             '</kineticLaw>',
             '</reaction>']

    out.extend(['</listOfReactions>', '</model>', '</sbml>'])

    doc = xml.dom.minidom.parseString("".join(out))

    fout = open(filename, "w")
    fout.write(header + '\n' + doc.documentElement.toprettyxml(indent="\t"))
    fout.close()


