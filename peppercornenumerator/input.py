#
#  input.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/26/10.
#


from __future__ import print_function
import logging

from dsdobjects.parser import parse_pil_string, parse_pil_file
from dsdobjects.parser import ParseException, PilFormatError

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
    error = float(line[1][1][1]) if line[1] != [] and line[1][1] != [] and len(line[1][1]) == 2 else None
    units = line[1][2][0] if line[1] != [] and line[1][2] != [] else None

    if rate is None:
        r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
        logging.warning("Ignoring input reaction without a rate: {}".format(r))
        return None, None, None, None, None, None
    elif rtype is None or rtype == 'condensed' or rtype not in PepperReaction.RTYPES:
        r = "{} -> {}".format(' + '.join(line[2]), ' + '.join(line[3]))
        logging.warning("Ignoring input reaction of with rtype='{}': {}".format(rtype, r))
        return None, None, None, None, None, None
    else :
        r = "[{} = {:12g} {}] {} -> {}".format(
                rtype, rate, units, ' + '.join(line[2]), ' + '.join(line[3]))

    return line[2], line[3], rtype, rate, units, r


def read_pil(data, is_file = False, composite = False):
    """ Old input standard. 

    Supports a variety of formats, including enum and pil. Sequences and sanity
    checks for domain-length are not supported.  

    Args:
        data (str): Is either the PIL file in string format or the path to a file.
        is_file (bool): True if data is a path to a file, False otherwise
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
            logging.info("Ignoring sequence information for domain {}.".format(name))
            if len(line) == 4:
                if int(line[3]) != len(line[2]):
                    logging.error("Sequence/Length information inconsistent {} vs ().".format(
                        line[3], len(line[2])))
                domains[name] = PepperDomain(name, length = int(line[3]))
            else :
                domains[name] = PepperDomain(name, length = len(line[2]))

            domains[name].nucleotides = line[2]
            cname = name[:-1] if domains[name].is_complement else name + '*'
            if cname in domains:
                assert domains[cname] == ~domains[name]
            else :
                domains[cname] = ~domains[name]

        elif line[0] == 'composite-domain':
            # This could be a strand definition or a composite domain.
            assert name[-1] != '*'
            sequences[name] = map(lambda x: domains[x], line[2])

            def comp(name):
                return name[:-1] if name[-1] == '*' else name + '*'

            sequences[comp(name)] = map(lambda x: domains[comp(x)], reversed(line[2]))

        elif line[0] == 'strand-complex':
            sequence = []
            for strand in line[2]:
                sequence += sequences[strand] + ['+']
            sequence = sequence[:-1]
            structure = line[3].replace(' ','')
            complexes[name] = PepperComplex(sequence, list(structure), name=name)

        elif line[0] == 'kernel-complex':
            sequence, structure = resolve_loops(line[2])

            # Replace names with domain objects.
            try :
                sequence = map(lambda d : domains[d], sequence)
            except KeyError:
                for e, d in enumerate(sequence):
                    if isinstance(d, PepperDomain):
                        # Happens with composite domains, see next statement e+1
                        continue
                    if d in sequences :
                        for i, c in enumerate(sequences[d]):
                            assert c.name in domains
                            if i == 0:
                                sequence[e] = c
                            else :
                                sequence.insert(e+i, c)
                                structure.insert(e+i, structure[e])

                    elif d not in domains :
                        logging.warning("Assuming {} is a long domain.".format(d))
                        domains[d] = PepperDomain(d, 'long')
                        cdom = ~domains[d]
                        domains[cdom.name] = cdom
                        sequence[e] = domains[d]
                    else :
                        sequence[e] = domains[d]

            complexes[name] = PepperComplex(sequence, structure, name=name)

            if len(line) > 3 :
                assert len(line[3]) == 3
                complexes[name]._concentration = tuple(line[3])


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
                logging.error("Rate units must be given in {}, not: {}.".format(reaction.rateunits, units))
                raise SystemExit
            reactions.append(reaction)

        else :
            logging.warning("Ignoring {} specification: {}".format(line[0], name))
        
    if composite :
        return complexes, reactions, sequences
    else :
        return complexes, reactions

def read_kernel(data, is_file = False):
    """ New input standard, kernel notation.

    including state and reaction, ignores concetrations
    """
    logging.warning("deprecated function: use read_pil")
    return read_pil(data, is_file)

def from_kernel(lines):
    """ Tranlsate a list of kernel strings. """
    logging.warning('deprecated function: use read_pil') 

    # split string into lines if necessary
    if isinstance(lines, basestring):
        lines = lines.split("\n")

    # remove blank lines
    lines = filter(None, lines)

    # reading pil in case of non-alphanumeric names
    complexes, _ = read_pil('\n'.join(lines))

    return (PepperDomain.MEMORY, None, complexes)

def load_pil_crn(data):
    """ Input for pilsimulator
    """
    parsed_file = parse_pil_string(data)

    sysunit = None

    species = dict()
    macrostates = dict()
    reactions = []
    for line in parsed_file:
        if line[0] == 'kernel-complex':
            name = line[1]
            conc = 0
            if len(line) > 3 :
                init = line[3][0]
                conc = float(line[3][1])
                if sysunit is None:
                    sysunit = line[3][2]
                else :
                    if sysunit != line[3][2]:
                        raise PilFormatError(
                                'Conflicting units: {} vs. {}'.format(sysunit, line[3][2]))
                if init[0] != 'i':
                    raise NotImplementedError('concentrations must be specified as *initial*')
            species[name] = ('initial', conc)
        elif line[0] == 'resting-macrostate':
            name = line[1]
            conc = ['initial', 0]
            for sp in line[2]:
                assert sp in species
                if species[sp][0][0] != 'i':
                    raise NotImplementedError('concentrations must be specified as *initial*')
                conc[1] += species[sp][1]
            macrostates[name] = tuple(conc)
        elif line[0] == 'reaction':
            info = line[1]
            reactants = line[2]
            products = line[3]
            assert len(info) == 3
            rate = float(info[1][0])
            if sysunit and filter(lambda x : x != sysunit, info[2][0].split('/')[1:-1]):
                raise PilFormatError('Conflicting units: {} vs. {}'.format(sysunit, info[2][0]))
            reactions.append([reactants, products, [rate]])
        elif line[0] == 'dl-domain' :
            pass
        else :
            print('# Ignoring Keyword: {}'.format(line[0]))

    detailed = None
    for rxn in reactions:
        reac = rxn[0]
        prod = rxn[1]
        rate = rxn[2]
        if any(map(lambda r: r in macrostates, reac + prod)):
            assert all(map(lambda r: r in macrostates, reac + prod))
            d = False
        else :
            d = True
        if detailed is None:
            detailed = d
        else :
            if detailed != d:
                raise PilFormatError('Need to provide either detailed or condensed CRN for simulation. Not both!')

    return reactions, species if detailed else macrostates

