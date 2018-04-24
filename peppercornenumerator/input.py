#
#  input.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/26/10.
#


import re
import logging
import json

from dsdobjects.parser import parse_kernel_string, parse_kernel_file
from dsdobjects.parser import parse_pil_string, parse_pil_file
from dsdobjects.parser import ParseException

from peppercornenumerator.objects import PepperDomain, PepperComplex, PepperReaction

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

def read_reaction(line):
    rtype = line[1][0][0] if line[1] != [] and line[1][0] != [] else None
    rate = float(line[1][1][0]) if line[1] != [] and line[1][1] != [] else None
    units = line[1][2][0] if line[1] != [] and line[1][2] != [] else None
    if rate is None or rtype is None or rtype == 'condensed' :
        r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
        logging.warning("Ignoring input reaction without a rate: {}".format(r))
        return None, None, None, None, None, None
    else :
        r = "[{} = {:12g} {}] {} -> {}".format(
                rtype, rate, units, ' + '.join(line[2]), ' + '.join(line[3]))

    return line[2], line[3], rtype, rate, units, r


def read_pil(data, is_file = False, alpha = ''):
    """ Old input standard. 

    Supports a variety of formats, including enum and pil. Sequences and sanity
    checks for domain-length are not supported.  

    Args:
        data (str): Is either the PIL file in string format or the path to a file.
        is_file (bool): True if data is a path to a file, False otherwise
        alpha (str, optiona): Is a prefix that can be applied to rename every domain
            in the system. Useful to convert non-alphanumeric domain names into 
            the alphanumeric kernel standard.
    """
    if is_file :
        parsed_file = parse_pil_file(data)
    else :
        parsed_file = parse_pil_string(data)

    domains = {'+' : '+'} # saves some code
    sequences = {}
    complexes = {}
    reactions = []
    for line in parsed_file :
        name = line[1]
        if line[0] == 'dl-domain':
            if line[2] == 'short':
                (dtype, dlen) = ('short', None)
            elif line[2] == 'long':
                (dtype, dlen) = ('long', None)
            else :
                (dtype, dlen) = (None, int(line[2]))
            if name not in domains:
                domains[name] = PepperDomain(name, dtype = dtype, length = dlen)
            logging.info('Domain {} with length {}'.format(domains[name], len(domains[name])))
            cname = name[:-1] if domains[name].is_complement else name + '*'
            if cname in domains:
                assert domains[cname] == ~domains[name]
            else :
                domains[cname] = ~domains[name]

        elif line[0] == 'sl-domain':
            logging.warning("Discarding sequence information for domain {}.".format(name))
            if len(line) == 4:
                domains[name] = PepperDomain(name, length = int(line[3]))
            else :
                domains[name] = PepperDomain(name, length = len(line[2]))
            cname = name[:-1] if domains[name].is_complement else name + '*'
            if cname in domains:
                assert domains[cname] == ~domains[name]
            else :
                domains[cname] = ~domains[name]

        elif line[0] == 'composite-domain':
            sequences[name] = map(lambda x: domains[alpha + x], line[2])

        elif line[0] == 'strand-complex':
            sequence = []
            for strand in line[2]:
                sequence += sequences[alpha + strand] + ['+']
            sequence = sequence[:-1]
            structure = line[3].replace(' ','')
            complexes[name] = PepperComplex(sequence, list(structure), name=name)

        elif line[0] == 'kernel-complex':
            sequence, structure = resolve_loops(line[2])

            # Replace names with domain objects.
            try :
                sequence = map(lambda d : domains[alpha + d], sequence)
            except KeyError:
                for e, d in enumerate(sequence):
                    if alpha + d not in domains :
                        logging.warning("Assuming {} is a long domain.".format(alpha + d))
                        domains[alpha + d] = PepperDomain(alpha + d, 'long')
                        cdom = ~domains[alpha + d]
                        domains[cdom.name] = cdom
                    sequence[e] = domains[alpha + d]

            if len(line) > 3 :
                logging.warning("Ignoring complex concentration for {}".format(name))

            complexes[name] = PepperComplex(sequence, structure, name=name)

        elif line[0] == 'reaction':
            reactants, products, rtype, rate, units, r = read_reaction(line)
            if r is None: continue

            try :
                reactants = map(lambda c : complexes[c], reactants)
                products  = map(lambda c : complexes[c], products)
            except KeyError:
                logging.warning("Ignoring input reaction with undefined complex: {}".format(r))
                continue
            
            reaction = PepperReaction(reactants, products, rtype=rtype, rate=rate)
            if reaction.rateunits != units:
                logging.warning("Rate units do not match: {} vs {}".format(
                    units, reaction.rateunits))
            reactions.append(reaction)

        elif line[0] == 'resting-macrostate':
            logging.warning("Ignoring resting-macrostate specification: {}".format(name))

        else :
            raise NotImplementedError('cannot interpret keyword:', line[0])

    return complexes, reactions

def read_kernel(data, is_file = False):
    """ New input standard, kernel notation.

    including state and reaction, ignores concetrations
    """

    if is_file :
        parsed_file = parse_kernel_file(data)
    else :
        parsed_file = parse_kernel_string(data)

    # Do domains first, just in case...
    domains = {'+' : '+'} # saves some code
    for line in parsed_file:
        name = line[1]
        if line[0] == 'dl-domain':
            if line[2] == 'short':
                (dtype, dlen) = ('short', None)
            elif line[2] == 'long':
                (dtype, dlen) = ('long', None)
            else :
                (dtype, dlen) = (None, int(line[2]))
            if name not in domains:
                domains[name] = PepperDomain(name, dtype = dtype, length = dlen)
            logging.info('Domain {} with length {}'.format(domains[name], len(domains[name])))
            cname = name[:-1] if domains[name].is_complement else name + '*'
            if cname in domains:
                assert domains[cname] == ~domains[name]
            else :
                domains[cname] = ~domains[name]

    complexes = {}
    reactions = []
    for line in parsed_file:
        name = line[1]
        if line[0] == 'dl-domain':
            pass
        elif line[0] == 'complex':
            sequence, structure = resolve_loops(line[2])

            # Replace names with domain objects.
            try :
                sequence = map(lambda d : domains[d], sequence)
            except KeyError:
                for e, d in enumerate(sequence):
                    if d not in domains :
                        logging.warning("Assuming {} is a long domain.".format(d))
                        domains[d] = PepperDomain(d, 'long')
                        cdom = ~domains[d]
                        domains[cdom.name] = cdom
                    sequence[e] = domains[d]

            if len(line) > 3 :
                logging.warning("Ignoring complex concentration for {}".format(name))

            complexes[name] = PepperComplex(sequence, structure, name=name)

        elif line[0] == 'reaction':
            reactants, products, rtype, rate, units, r = read_reaction(line)
            if r is None: continue

            try :
                reactants = map(lambda c : complexes[c], reactants)
                products  = map(lambda c : complexes[c], products)
            except KeyError:
                logging.warning("Ignoring input reaction with undefined complex: {}".format(r))
                continue

            reaction = PepperReaction(reactants, products, rtype=rtype, rate=rate)
            if reaction.rateunits != units:
                logging.warning("Rate units do not match: {} vs {}".format(
                    units, reaction.rateunits))
            reactions.append(reaction)

        elif line[0] == 'resting-macrostate':
            logging.warning("Ignoring resting-macrostate specification: {}".format(name))
        else :
            raise NotImplementedError('cannot interpret keyword:', line[0])

    return complexes, reactions

def from_kernel(lines):
    """ Tranlsate a list of kernel strings. """
    logging.warn('deprecated function: read_kernel') 

    # split string into lines if necessary
    if isinstance(lines, basestring):
        lines = lines.split("\n")

    # remove blank lines
    lines = filter(None, lines)

    # reading pil in case of non-alphanumeric names
    complexes, _ = read_pil('\n'.join(lines))

    return (PepperDomain.MEMORY, None, complexes)

