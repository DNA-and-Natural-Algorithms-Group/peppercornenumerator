#
#  input.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/26/10.
#


import re
import logging
import json

import peppercornenumerator
import peppercornenumerator.utils as utils
import peppercornenumerator.reactions as reactions

from peppercornenumerator.objects import PepperDomain, PepperComplex, PepperRestingState
from dsdobjects.parser import parse_kernel_string, parse_kernel_file

def input_enum(filename):
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
                logging.error( "Duplicate domain name encountered in input line %d" % line_counter)
                raise Exception()

            if not re.match(r'\w+$', domain_name):
                logging.warn( "Non-alphanumeric domain name %s encountered in input line %d" % (domain_name, line_counter))

            # The domain length could be either short or long or it could be
            # an exact number
            domain_length = parts[3]
            if (domain_length == 'short') :
                domain_type = 'short'
                domain_length = None
            elif (domain_length == 'long') :
                domain_type = 'long'
                domain_length = None
            else :
                domain_length = int(domain_length)
                domain_type = None
                if domain_length <= 0: logging.warn("Domain of length %d found in input line %d" % (domain_length, line_counter))

            # # Check to see if a sequence is specified
            # if len(parts) > 4:
            #     domain_sequence = list(parts[4])
            # else:
            #     domain_sequence = list('N' * domain_length)

            # Create the new domains
            new_dom = PepperDomain(domain_name, dtype = domain_type, length = domain_length)
            new_dom_comp = ~new_dom

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
                logging.error(
                    "Duplicate strand name encountered in input line %d" %
                    line_counter)
                raise Exception()

            if not re.match(r'\w+$', strand_name):
                logging.warn(
                    "Non-alphanumeric strand name %s encountered in input line %d" %
                    (strand_name, line_counter))

            strand_doms = []
            for domain_name in parts[3:]:
                if domain_name not in domains:
                    logging.error( "Invalid domain name %s encountered in input line %d" % (domain_name, line_counter))
                    raise Exception()

                strand_doms.append(domains[domain_name])

            if len(strand_doms) == 0:
                logging.warn("0-length strand encountered in input line %d")

            strands[strand_name] = strand_doms

        # This is the beginning of a complex
        elif line.startswith("complex"):
            # e.g.:
            # complex A :\n
            # A\n				<- strands_line
            # .(((..)))\n		<- structure_line

            parts = line.split()

            complex_name = parts[1]
            if complex_name in complexes:
                logging.error( "Duplicate complex name encountered in input line %d" % line_counter)
                raise Exception()

            if not re.match(r'\w+$', complex_name):
                logging.warn( "Non-alphanumeric complex name %s encountered in input line %d" % (complex_name, line_counter))

            complex_sequence = []

            strands_line = fin.readline()
            strands_line = strands_line.strip()
            strands_line_parts = strands_line.split()
            for strand_name in strands_line_parts:
                if strand_name not in strands:
                    logging.error( "Invalid strand name %s encountered in input line %d" % (strand_name, line_counter))
                    raise Exception()
                else:
                    complex_sequence.extend(strands[strand_name] + ['+'])

            complex_sequence = complex_sequence[:-1] # remove last strandbreak

            structure_line = fin.readline()
            structure_line = structure_line.strip()

            complexes[complex_name] = PepperComplex(complex_sequence, list(structure_line.replace(' ','')), name=complex_name)
        else:
            logging.error("Unexpected characters encountered in input line %d" % line_counter)
            raise Exception()
        line = fin.readline()
        line_counter += 1

    return peppercornenumerator.Enumerator(complexes.values())

def complexes_from_kernel(kernel_string, ddlen = 15):
    parsed_kernel = parse_kernel_string(kernel_string)
    complexes = build_pepper_complexes(parsed_kernel, ddlen)
    return complexes

def build_pepper_complexes(parsed_kernel, ddlen = 15):
    def resolve_loops(loop):
        """ Return a sequence, structure pair from kernel format with parenthesis. """
        sequen = []
        struct = []
        for dom in loop :
            if isinstance(dom, str):
                sequen.append(dom)
                if dom == '+' :
                    struct.append('+')
                else :
                    struct.append('.')
            elif isinstance(dom, list):
                struct[-1] = '('
                old = sequen[-1]
                se, ss = resolve_loops(dom)
                sequen.extend(se)
                struct.extend(ss)
                sequen.append(old + '*' if old[-1] != '*' else old[:-1])
                struct.append(')')
        return sequen, struct

    # Do domains first, just in case...
    domains = {}
    for line in parsed_kernel:
        if line[0] == 'domain':
            domains[line[1]] = PepperDomain(line[1], length = int(line[2]))

    complexes = {}
    for line in parsed_kernel:
        if line[0] == 'complex':
            name = line[1]
            sequence, structure = resolve_loops(line[2])

            # Replace names with domain objects.
            # If Domain is not in domains, it has default-domain length = ddlen
            for e in range(len(sequence)):
              dname = sequence[e]
              if dname == '+':
                  continue

              if dname[-1] == '*' : # complement domain
                  if dname[:-1] in domains :
                      dom = ~domains[dname[:-1]]
                  else :
                      dom = PepperDomain(dname, length = ddlen)
                      domains[dname[:-1]] = ~dom
              else : # non-complement domain
                  if dname in domains :
                      dom = domains[dname]
                  else :
                      dom = PepperDomain(dname, length = ddlen)
                      domains[dname] = dom
              sequence[e] = dom
            complexes[name] = PepperComplex(sequence, structure, name=name)
    return complexes

def from_kernel(lines, ddlen = 15):
    """ Tranlsate a list of kernel strings. """
    print DeprecationWarning('use new function: complexes_from_kernel') 

    # split string into lines if necessary
    if isinstance(lines, basestring):
        lines = lines.split("\n")

    # remove blank lines
    lines = filter(None, lines)

    parsed_kernel = []
    for l in lines:
        parsed_kernel.extend(parse_kernel_string(l))

    complexes = build_pepper_complexes(parsed_kernel, ddlen)

    return (PepperDomain.MEMORY, None, complexes)

def input_kernel(filename, ddlen=15):
    parsed_kernel = parse_kernel_file(filename)
    complexes = build_pepper_complexes(parsed_kernel, ddlen)
    return peppercornenumerator.Enumerator(complexes.values())

def input_pil(filename):
    raise DeprecationWarning('no support for regular pil format')
    """
    Initializes and returns an enumerator from an input file in the Pepper Intermediate Language (PIL)
    """
    fin = open(filename, 'r')

    domains = {}
    strands = {}
    complexes = {}

    # maps domain-wise strand structures to auto-generated strand names
    structures_to_strands = {}

    # We loop over all the lines in the file
    for (line_counter, line) in enumerate(fin, start=1):
        line = line.strip()

        # This was an empty line
        if line == "":
            continue

        # This was a comment
        elif line.startswith("#"):
            continue

        elif line.startswith("length"):
            # e.g.:
            #       "length a = 6"
            # parts:        0   1
            parts = re.match(r"length\s*([\w-]+)\s*=\s*(\d+)\s*", line)
            if parts is None:
                logging.error("Invalid syntax on input line %d"
                              % line_counter)
                logging.error(line)
                raise Exception()

            domain_name, domain_length = parts.groups()
            if domain_name in domains:
                logging.error(
                    "Duplicate domain name encountered in input line %d" %
                    line_counter)
                raise Exception()

            if not re.match(r'[\w-]+$', domain_name):
                logging.warn(
                    "Non-alphanumeric domain name %s encountered in input line %d" %
                    (domain_name, line_counter))

            domain_length = int(domain_length)

            # Create the new domains
            new_dom = PepperDomain(domain_name, length = domain_length)
            new_dom_comp = ~new_dom

            domains[domain_name] = new_dom
            domains["%s*" % domain_name] = new_dom_comp

        # This is the beginning of a domain
        elif line.startswith("sequence"):
            # e.g.:
            #       "sequence a = 6 : 6"
            # parts: 0        1 2 3

            #                  sequence    a         =   NNNNN   :     6
            parts = re.match(
                r"sequence\s*([\w-]+)\s*=\s*(\w+)\s*:?\s*(\d?)\s*", line)
            if parts is None:
                logging.error("Invalid syntax on input line %d"
                              % line_counter)
                logging.error(line)
                raise Exception()

            domain_name, domain_sequence, length = parts.groups()
            if domain_name in domains:
                logging.error(
                    "Duplicate domain name encountered in input line %d" %
                    line_counter)
                raise Exception()

            if not re.match(r'[\w-]+$', domain_name):
                logging.warn(
                    "Non-alphanumeric domain name %s encountered in input line %d" %
                    (domain_name, line_counter))

            # The sequence specification
            # domain_sequence = parts[1]
            domain_length = len(domain_sequence)

            new_dom = PepperDomain(domain_name, length = domain_length)
            new_dom_comp = ~new_dom

            domains[domain_name] = new_dom
            domains["%s*" % domain_name] = new_dom_comp

        elif line.startswith("sup-sequence"):
            # e.g.
            #       "sequence a = b c d e : 6"
            #                 0   1         2
            parts = re.match(
                r"sup-sequence\s*([\w-]+)\s*=\s*((?:[\w-]+\s*)+):?(\d?)", line)
            if parts is None:
                logging.error("Invalid syntax on input line %d"
                              % line_counter)
                logging.error(line)
                raise Exception()

            domain_name, sequence_names, length = parts.groups()

            # domain name
            if domain_name in domains:
                logging.error(
                    "Duplicate domain name encountered in input line %d" %
                    line_counter)
                raise Exception()

            if not re.match(r'[\w-]+$', domain_name):
                logging.warn(
                    "Non-alphanumeric domain name %s encountered in input line %d" %
                    (domain_name, line_counter))

            # subsequences
            sequence = ""
            for sequence_name in sequence_names.split():
                sequence_name = sequence_name.strip()
                if sequence_name == "":
                    continue

                # make sure each subsequence is defined
                if sequence_name not in domains:
                    logging.error(
                        "Unknown domain name '%s' in super-sequence on input line %d" %
                        (sequence_name, line_counter))
                    logging.error(line)
                    raise Exception()

                # build up the full sequence
                sequence += ''.join(domains[sequence_name].sequence)

            # check for correctness
            if length:
                if int(length) != len(sequence):
                    logging.error(
                        "Sequence length for super-sequence %s is %d, not equal to expected value %d on input line %d" %
                        (domain_name, len(sequence), length, line_counter))
                    raise Exception()

            # The sequence specification
            domain_sequence = sequence
            domain_length = len(sequence)

            # Create the new domains
            new_dom = PepperDomain(domain_name, length = domain_length)
            new_dom_comp = ~new_dom

            domains[domain_name] = new_dom
            domains["%s*" % domain_name] = new_dom_comp

        elif line.startswith("equal"):
            parts = line.split()
            if len(parts) < 3:
                logging.error(
                    "'equal' statement does not specify at least 2 domains on input line %d" %
                    line_counter)
                logging.error(line)
                raise Exception()

            source_domain_name = parts[1]
            target_domain_names = parts[2:]

            if source_domain_name not in domains:
                logging.error(
                    "Unknown domain name '%s' in 'equals' statement on input line %d" %
                    (source_domain_name, line_counter))
                logging.error(line)
                raise Exception()

            source_domain = domains[source_domain_name]

            for target_domain_name in target_domain_names:
                new_dom = PepperDomain(target_domain_name, length = source_domain.length)
                new_dom_comp = ~new_dom

                domains[target_domain_name] = new_dom
                domains["%s*" % target_domain_name] = new_dom_comp

        # This is the beginning of a strand
        elif line.startswith("strand"):
            # e.g.:
            #       "strand A = a x b y z* c* y* b* x*"
            # parts: 0      1 2 3 4 5 6 ...

            parts = re.match(
                r"strand\s*([\w-]+)\s*=\s*((?:[\w*-]+\s*)+):?(\d?)", line)
            if parts is None:
                logging.error("Invalid syntax on input line %d"
                              % line_counter)
                logging.error(line)
                raise Exception()

            strand_name, strand_dom_names, length = parts.groups()

            if strand_name in strands:
                logging.error(
                    "Duplicate strand name encountered in input line %d" %
                    line_counter)
                raise Exception()

            if not re.match(r'\w+$', strand_name):
                logging.warn(
                    "Non-alphanumeric strand name %s encountered in input line %d" %
                    (strand_name, line_counter))

            strand_doms = []
            for domain_name in filter(None, strand_dom_names.split()):
                if domain_name not in domains:
                    logging.error(
                        "Invalid domain name %s encountered in input line %d" %
                        (domain_name, line_counter))
                    logging.error(line)
                    print domains
                    raise Exception()

                strand_doms.append(domains[domain_name])

            if len(strand_doms) == 0:
                logging.warn("0-length strand encountered in input line %d")

            new_strand = utils.Strand(strand_name, strand_doms)
            strands[strand_name] = new_strand

        # This is the beginning of a complex
        elif line.startswith("structure"):
            # parse `structure` line:
            # e.g.:
            # structure A = S1 : .(((..)))
            # structure    [  1nt  ]      name      =   s1 s2 s3 + s4         :
            # ....((+))...((..))....
            parts = re.match(
                r"structure\s+(\[[^\]]+\])?\s*([\w-]+)\s*=\s*((?:[\w-]+\s*\+?\s*)+):\s*([().+\s]+)",
                line)

            if parts is None:

                # parse `structure` line:
                # e.g.:
                # structure A = S1 : .(((..)))
                # structure    name      =   s1 s2 s3 + s4         :
                # ....((+))...((..))....
                parts = re.match(
                    r"structure\s+([\w-]+)\s*=\s*((?:[\w-]+\s*\+?\s*)+):\s*([().+\s]+)", line)

                if parts is None:

                    logging.error("Invalid syntax on input line %d"
                                  % line_counter)
                    logging.error(line)
                    raise Exception()
                else:
                    complex_name, strands_line, structure_line = parts.groups()
                    parameters = ""
            else:
                parameters, complex_name, strands_line, structure_line = parts.groups()

            # parse parameters
            if parameters is None:
                parameters = ""
            params = utils.parse_parameters(parameters)

            # check for duplicate complex name
            if complex_name in complexes:
                logging.error(
                    "Duplicate complex name encountered in input line %d" %
                    line_counter)
                raise Exception()

            # check for non-alphanumeric complex name
            if not re.match(r'\w+$', complex_name):
                logging.warn(
                    "Non-alphanumeric complex name %s encountered in input line %d" %
                    (complex_name, line_counter))

            # get strand names, allowing optional '+' characters
            complex_strands = []
            strands_line_parts = [
                name for name in strands_line.split() if name != "+"]
            for strand_name in strands_line_parts:
                if strand_name not in strands:
                    logging.error(
                        "Invalid strand name %s encountered in input line %d" %
                        (strand_name, line_counter))
                    raise Exception()
                else:
                    complex_strands.append(strands[strand_name])

            # parse dot-paren structure, then do some horrible magic to guess
            # if it's basewise or segmentwise...
            complex_structure = utils.parse_dot_paren(structure_line)
            struct_length = sum(map(len, complex_structure))
            # sum([ len(d) for c in complex_strands for d in c.domains ])
            domains_length = sum(map(len, complex_strands))

            if(struct_length > domains_length):
                complex_structure = utils.parse_basewise_dot_paren(
                    structure_line, complex_strands)
                struct_length = sum(map(len, complex_structure))
                if(struct_length != domains_length):
                    logging.error(
                        "Complex %(name)s has %(doms)d domains but structure size %(struct_length)d. (structure was '%(struct)s')" % {
                            "name": complex_name,
                            "doms": domains_length,
                            "struct_length": struct_length,
                            "struct": structure_line})
                    raise Exception()

            elif(struct_length != domains_length):

                logging.error(
                    "Complex %(name)s has %(doms)d domains but structure size %(struct_length)d. (structure was '%(struct)s')" % {
                        "name": complex_name,
                        "doms": domains_length,
                        "struct_length": struct_length,
                        "struct": structure_line})
                raise Exception()

            complex = utils.Complex(
                complex_name, complex_strands, complex_structure)
            complex.check_structure()

            # apply parameters
            if params['concentration'] is not None:
                complex.concentration = params['concentration']

            complexes[complex_name] = complex
        elif line.startswith("kinetic"):
            continue
        elif line.strip() == "":
            continue
        else:
            try:
                resolve_kernel([line], domains, strands,
                               structures_to_strands, complexes)
            except Exception as e:
                logging.error(
                    "Unexpected characters encountered in input line %d; tried to parse as Kernel statement but got error: %s" %
                    (line_counter, str(e)))
                raise Exception()

        # line = fin.readline()
        # line_counter += 1

    fin.close()
    domains = domains.values()
    strands = strands.values()
    complexes = complexes.values()

    enumerator = peppercornenumerator.Enumerator(complexes, strands, domains)
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
        if 'sequence' not in saved_domain:
            saved_domain['sequence'] = 'N' * int(saved_domain['length'])
        else:
            saved_domain['sequence'] = str(''.join(saved_domain['sequence']))

        if (saved_domain['is_complement']):
            #saved_domain['name'] = str(saved_domain['name'][:-1])
            new_dom = utils.PepperDomain(list(saved_domain['sequence']), 
                    name=str(saved_domain['name']), is_complement=True)
            #new_dom = tmp_dom.get_ComplementDomain(list(saved_domain['sequence']))

        else :
            new_dom = utils.PepperDomain(list(saved_domain['sequence']), 
                    name=str(saved_domain['name']))

        #new_dom = utils.Domain(
        #    saved_domain['name'],
        #    saved_domain['length'],
        #    is_complement=saved_domain['is_complement'],
        #    sequence=saved_domain['sequence'])
        domains[new_dom.name] = new_dom

    saved_strands = saved['strands']

    strands = {}
    for saved_strand in saved_strands:
        doms = []
        for domain in saved_strand['domains']:
            doms.append(domains[domain])

        new_strand = utils.Strand(saved_strand['name'], doms)
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
                if (tup is None):
                    new_strand.append(None)
                else:
                    new_strand.append(tuple(tup))
            new_structure.append(new_strand)
        new_complex = utils.Complex(
            saved_complex['name'], c_strands, new_structure)
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
                if (tup is None):
                    new_strand.append(None)
                else:
                    new_strand.append(tuple(tup))
            new_structure.append(new_strand)
        new_complex = utils.Complex(
            saved_complex['name'], c_strands, new_structure)
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

        reaction = PepperReaction(reactants, products, saved_reaction['name'])
        reactions.append(reaction)

    resting_states = []

    for resting_state in saved['resting_states']:
        comps = []
        for complex in resting_state['complexes']:
            comps.append(complexes[complex])
        resting_states.append(PepperRestingState(comps, name=resting_state['name']))

    initial_complexes = {}
    for saved_complex in saved['initial_complexes']:
        c_strands = []
        for strand in saved_complex['strands']:
            c_strands.append(strands[strand])
        new_structure = []
        for strand in saved_complex['structure']:
            new_strand = []
            for tup in strand:
                if (tup is None):
                    new_strand.append(None)
                else:
                    new_strand.append(tuple(tup))
            new_structure.append(new_strand)
        new_complex = utils.Complex(
            saved_complex['name'], c_strands, new_structure)
        initial_complexes[saved_complex['name']] = new_complex

    enumerator = peppercornenumerator.Enumerator(
        initial_complexes.values(), strands.values(), domains.values())
    enumerator._complexes = complexes.values()
    enumerator._resting_states = resting_states
    enumerator._transient_complexes = transient_complexes.values()
    enumerator._resting_complexes = resting_complexes.values()
    enumerator._reactions = reactions

    return enumerator


text_input_functions = {
    'enum': input_enum,
    'pil': input_kernel
}

load_input_functions = {
    'json': load_json
}
