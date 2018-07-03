#
#  input.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 6/26/10.
#


from __future__ import print_function
import logging

from dsdobjects.parser import parse_pil_string, parse_pil_file
from dsdobjects.parser import parse_seesaw_string, parse_seesaw_file
from dsdobjects.parser import ParseException, PilFormatError

from peppercornenumerator.objects import PepperDomain, PepperComplex, PepperReaction

def InputFormatError(Exception):
    pass

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
    """ Peppercorn standard input.

    Supports a variety of pil-style dialects, including kernel and enum.

    Args:
        data (str): Is either the PIL file in string format or the path to a file.
        is_file (bool): True if data is a path to a file, False otherwise
        composite (bool, optional): Returns an additional dictionary that maps
            names of composite domains (or strands) to a list of domains.
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


def read_seesaw(data, 
        is_file = False, 
        conc = 100e-9, 
        explicit = True):
    """Translates a seesaw circuit into PepperComplexes.

    Args:
        is_file (bool, optional): Pareses a file if True, parses a string if False
        conc (float, optional): The 1x concentration value in M.    
            Defaults to 100 nM
        explicit (bool, optional): Translate to the exact explicit notation, or to 
            the easier shorthand-notation. Defaults to True: explicit.

    The explicit notation is more accurate, but also harder to enumerate.
    Typical domain lengths:
        gate-b: 5 + 15 + 5
        thld-b: 5 + 5 + 15
        long-domain: 10 + 5
        toehold: 5
        gate-top: 15 + 5 + 15 != 32!!
 
    """
    if is_file :
        parsed_file = parse_seesaw_file(data)
    else :
        parsed_file = parse_seesaw_string(data)

    domains = {'+' : '+'}
    def assgn_domain(name, l):
        if name not in domains:
            domains[name] = PepperDomain(name, length = l)
            cname = name[:-1] if domains[name].is_complement else name + '*'
            domains[cname] = ~domains[name]
        return domains[name]

    if explicit:
        toe = assgn_domain('toe', 3) # toehold TCT - AGA
        clp = assgn_domain('clp', 2) # clamp   CA  - TG
        dlen = 6
        tlen = 5
    else :
        toe = assgn_domain('toe', 5) # toehold CATCT - AGATG
        dlen = 10
        tlen =  5

    complexes = {}
    def make_wire(w, name=None):
        """Translates a wire into a PepperComplex.

        Every complex gets 0 M concentration as default, i.e. it is excluded
        from enumeration.

        Args:
            name (str, optional): Provida a name for the complex that is
                different from the internal wire representation.

        explicit:
            ['w', ['1', '2']] -> w1_2 = c s2 c t c s1 c
        shorthand:
            ['w', ['1', '2']] -> w1_2 = s2 t s1

        """
        if w[0] != 'w': 
            raise InputFormatError('Does not look like a wire: {}'.format(w))

        [nI, nO] = w[1]
        wire = 'w{}_{}'.format(nI, nO)
        if wire not in complexes:
            dIa = assgn_domain('s{}a'.format(nI), l=dlen)
            dIt = assgn_domain('s{}t'.format(nI), l=tlen)
            dOa = assgn_domain('s{}a'.format(nO), l=dlen)
            dOt = assgn_domain('s{}t'.format(nO), l=tlen)
            if explicit:
                sequence = [clp, dOt, dOa, clp, toe, clp, dIt, dIa, clp]
                structure = list('.........')
            else:
                sequence = [dOt, dOa, toe, dIt, dIa]
                structure = list('.....')
            complexes[wire] = PepperComplex(sequence, structure, 
                    name=name if name else wire)
            complexes[wire]._concentration = ('i', 0, 'M')
        return complexes[wire]

    def make_thld(t):
        """ ['t', [inp, name]] """
        if t[0] != 't': 
            raise InputFormatError('Does not look like a threshold: {}'.format(t))
        [nI, nO] = t[1]
        thld = 'T{}_{}'.format(nI, nO)
        if thld not in complexes:
            dIa = assgn_domain('s{}a'.format(nI), l=dlen)
            dIt = assgn_domain('s{}t'.format(nI), l=tlen)
            dOa = assgn_domain('s{}a'.format(nO), l=dlen)
            dOt = assgn_domain('s{}t'.format(nO), l=tlen)
            if explicit:
                sequence = [clp, dOt, dOa, clp, '+', ~dIt, ~clp, ~toe, ~clp, ~dOa, ~dOt, ~clp]
                structure = list('((((+...))))')
            else :
                sequence = [dOt, dOa, '+', ~dIt, ~toe, ~dOa, ~dOt]
                structure = list('((+..))')
            complexes[thld] = PepperComplex(sequence, structure, name=thld)
            complexes[thld]._concentration = ('i', 0, 'M')

            # Name the waste as well, for convenience...
            waste = '{}_c'.format(thld)
            assert waste not in complexes
            if explicit:
                sequence = [clp, dOt, dOa, clp, toe, clp, dIt, dIa, clp, '+', ~dIt, ~clp, ~toe, ~clp, ~dOa, ~dOt, ~clp]
                structure = list('(((((((..+)))))))')
            else :
                sequence = [dOt, dOa, toe, dIt, dIa, '+', ~dIt, ~toe, ~dOa, ~dOt]
                structure = list('((((.+))))')
            complexes[waste] = PepperComplex(sequence, structure, name=waste)
            complexes[waste]._concentration = ('i', 0, 'M')

        return complexes[thld]

    def make_gate(g):
        """ Two options: 
        1) ['g', [['w', ['31', '25']], '25']]
        2) ['g', ['31', ['w', ['31', '25']]]]
        """
        if g[0] != 'g': 
            raise InputFormatError('Does not look like a gate: {}'.format(g))

        if isinstance(g[1][0], list):
            # ['g', [['w', ['1', '2']], '2']]
            [[_, [nI, nO]], nG] = g[1]
            assert _ == 'w' and nG == nO
        elif isinstance(g[1][1], list):
            # ['g', ['1', ['w', ['1', '2']]]]
            [nG, [_, [nI, nO]]] = g[1]
            assert _ == 'w' and nG == nI

        if nI == nG : # this is a producing gate:
            assert nO != nG
            dGa = assgn_domain('s{}a'.format(nG), l=dlen)
            dGt = assgn_domain('s{}t'.format(nG), l=tlen)
            dOa = assgn_domain('s{}a'.format(nO), l=dlen)
            dOt = assgn_domain('s{}t'.format(nO), l=tlen)

            wire = 'w{}_{}'.format(nG, nO) # output wire
            gate = 'G_g{}_{}'.format(nG, wire)
            if explicit:
                sequence  = [clp, dOt, dOa, clp, toe, clp, dGt, dGa, clp, '+', 
                            ~clp, ~toe, ~clp, ~dGa, ~dGt, ~clp, ~toe, ~clp]
                structure = list('...((((((+..))))))')
            else :
                sequence  = [dOt, dOa, toe, dGt, dGa, '+', ~toe, ~dGa, ~dGt, ~toe]
                structure = list('..(((+.)))')

        else : # this is a consumed gate:
            assert nO == nG and nI != nG
            dIa = assgn_domain('s{}a'.format(nI), l=dlen)
            dIt = assgn_domain('s{}t'.format(nI), l=tlen)
            dGa = assgn_domain('s{}a'.format(nG), l=dlen)
            dGt = assgn_domain('s{}t'.format(nG), l=tlen)

            wire = 'w{}_{}'.format(nI, nG) # output wire
            gate = 'G_{}_g{}'.format(wire, nG)
            if explicit:
                sequence  = [clp, dGt, dGa, clp, toe, clp, dIt, dIa, clp, '+', 
                            ~clp, ~toe, ~clp, ~dGa, ~dGt, ~clp, ~toe, ~clp]
                structure = list('((((((...+))))))..')
            else :
                sequence  = [dGt, dGa, toe, dIt, dIa, '+', ~toe, ~dGa, ~dGt, ~toe]
                structure = list('(((..+))).')
 
        if gate not in complexes:
            complexes[gate] = PepperComplex(sequence, structure, name=gate)
            complexes[gate]._concentration = ('i', 0, 'M')
        return complexes[gate]

    for line in parsed_file :
        if line[0] == 'INPUT':
            # use the given name if it is not a digit...
            name = None if line[1][0].isdigit() else line[1][0]
            w = make_wire(line[2], name)
            w._concentration = None

        elif line[0] == 'OUTPUT':
            # use the given name if it is not a digit...
            name = None if line[1][0].isdigit() else line[1][0]
            if line[2][0] == 'w':
                make_wire(line[2], name)
            elif line[2][0] == 'Fluor':
                nF = line[2][1]
                fluor = "F_{}".format(line[2][1])
                name = fluor if name is None else name
                dFa = assgn_domain('s{}a'.format(nF), l=dlen)
                dFt = assgn_domain('s{}t'.format(nF), l=tlen)
                if fluor not in complexes:
                    if explicit:
                        complexes[fluor] = PepperComplex([clp, dFt, dFa, clp], ['.','.','.','.'], name=name)
                    else :
                        complexes[fluor] = PepperComplex([dFt, dFa], ['.','.'], name=name)
                    complexes[fluor]._concentration = ('i', 0, 'M')
            else:
                raise InputFormatError('Unknown output format: {}'.format(line))

        elif line[0] == 'seesaw':
            # seesaw[2,{1,4},{3}]
            [gate, inp, out] = line[1]
            for d in inp:
                make_wire(['w', [d, gate]])
                make_thld(['t', [d, gate]])
                make_gate(['g', [['w', [d, gate]], gate]])
            for d in out:
                make_wire(['w', [gate, d]])
                make_gate(['g', [gate, ['w',[gate, d]]]])

        elif line[0] in ('seesawOR', 'seesawAND'):
            # seesaw[2,5,{1,4},{3,9,10}]
            [g1, g2, inp, out] = line[1]
            n = len(inp)
            m = len(out)

            for d in inp:
                make_wire(['w', [d, g1]])
                make_thld(['t', [d, g1]])
                make_gate(['g', [['w', [d, g1]], g1]])

            # connection
            make_wire(['w',[g1, g2]])
            cx = make_thld(['t',[g1, g2]]) 
            if line[0] == 'seesawOR':
                cx._concentration = ('i', 1.1 * 0.6 * conc, 'M')
            else:
                cx._concentration = ('i', 1.1 * (n-1 + 0.2) * conc, 'M')
            cx = make_gate(['g', [g1, ['w', [g1, g2]]]])
            cx._concentration = ('i', n * conc, 'M')
            make_gate(['g', [['w',[g1, g2]], g2]])

            for d in out + ['f']:
                cxw = make_wire(['w', [g2, d]])
                cxg = make_gate(['g', [g2, ['w',[g2, d]]]])
                if d == 'f':
                    cxw._concentration = ('i', 2 * m * conc, 'M')
                else :
                    cxg._concentration = ('i', conc, 'M')

        elif line[0] == 'inputfanout':
            [g1, inp, out] = line[1]
            n = len(inp)
            m = len(out)

            make_wire(['w', [inp, g1]])
            cx = make_thld(['t', [inp, g1]])
            cx._concentration = ('i', 1.1 * 0.2 * conc, 'M')
            make_gate(['g', [['w',[inp, g1]], g1]])

            for d in out + ['f']:
                cxw = make_wire(['w', [g1, d]])
                cxg = make_gate(['g', [g1, ['w', [g1, d]]]])
                if d == 'f':
                    cxw._concentration = ('i', 2 * m * conc, 'M')
                else :
                    cxg._concentration = ('i', conc, 'M')

        elif line[0] == 'reporter':
            #['reporter', ['25', '31']]
            [nR, nI] = line[1]
            make_wire(['w', [nI, nR]])

            rep = 'R_{}'.format(nR)
            assert rep not in complexes

            dRa = assgn_domain('s{}a'.format(nR), l=dlen)
            dRt = assgn_domain('s{}t'.format(nR), l=tlen)
            if explicit:
                sequence  = [clp, dRt, dRa, clp, '+', ~clp, ~toe, ~clp, ~dRa, ~dRt, ~clp]
                structure = list('((((+..))))')
            else: 
                sequence  = [dRt, dRa, '+', ~toe, ~dRa, ~dRt]
                structure = list('((+.))')
            complexes[rep] = PepperComplex(sequence, structure, name=rep)
            complexes[rep]._concentration = ('i', 1.5 * conc, 'M')

            flr = "F_{}".format(nR)
            if flr not in complexes:
                if explicit:
                    complexes[flr] = PepperComplex([clp, dRt, dRa, clp], ['.','.','.','.'], name=flr)
                else :
                    complexes[flr] = PepperComplex([dRt, dRa], ['.','.'], name=flr)
                complexes[flr]._concentration = ('i', 0, 'M')

        elif line[0] == 'conc':
            if line[1][0] == 'w':
                cx = make_wire(line[1])
            elif line[1][0] == 'g':
                cx = make_gate(line[1])
            elif line[1][0] == 'th':
                cx = make_thld(line[1])
            else:
                print(line)
                raise NotImplementedError
            cx._concentration = ('i', conc*float(line[2][0]), 'M')

        else :
            print('WARNING: keyword not supported: {}'.format(line[0]))

    return complexes, []
