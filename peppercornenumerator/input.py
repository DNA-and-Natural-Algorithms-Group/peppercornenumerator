#
#  peppercornenumerator/input.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

from dsdobjects.objectio import set_io_objects, clear_io_objects
from dsdobjects.objectio import read_pil as dsd_read_pil
from dsdobjects.objectio import read_pil_line as dsd_read_pil_line
from dsdobjects import (parse_pil_string,
                        parse_seesaw_string,
                        parse_seesaw_file)
from .objects import (PepperDomain, 
                      PepperStrand,
                      PepperComplex, 
                      PepperReaction, 
                      PepperMacrostate)


class InputFormatError(Exception):
    pass

class PilFormatError(Exception):
    pass

def read_pil_line(*args, **kwargs):
    set_io_objects(D = PepperDomain, S = PepperStrand, C = PepperComplex, R = PepperReaction, M = PepperMacrostate)
    out = dsd_read_pil_line(*args, **kwargs)
    clear_io_objects()
    return out

def read_pil(data, is_file = False, ignore = None):
    set_io_objects(D = PepperDomain, S = PepperStrand, C = PepperComplex, R = PepperReaction, M = PepperMacrostate)
    out = dsd_read_pil(data, is_file, ignore)
    clear_io_objects()

    if len(out['con_reactions']):
        log.warning('Ignoring condensed reaction input.')
    if len(out['macrostates']):
        log.warning('Ignoring complex macrostate input.')

    # Make sure that strands are not within initial complexes!
    for strand in out['strands'].values():
        strand.concentration = ('constant', 0, 'M')
        out['complexes'][('s', strand.name)] = strand
    out['domains'].clear()
    return out['complexes'], out['det_reactions']

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
            mode = 'initial'
            if len(line) > 3 :
                init = line[3][0]
                conc = float(line[3][1])
                if sysunit is None:
                    sysunit = line[3][2]
                else :
                    if sysunit != line[3][2]:
                        raise PilFormatError('Conflicting units: {} vs. {}'.format(sysunit, line[3][2]))
                if init[0] == 'c':
                    mode = 'constant'
            species[name] = (mode, conc)
        elif line[0] == 'resting-macrostate':
            name = line[1]
            conc = 0
            mode = None
            for sp in line[2]:
                assert sp in species
                conc += species[sp][1]
                if mode: 
                    if mode[0] != species[sp][0][0]:
                        raise PilFormatError('Constant and initial species in the same resting macrostate.')
                else:
                    mode = 'initial' if species[sp][0][0] == 'i' else 'constant'
            macrostates[name] = (mode, conc)
        elif line[0] == 'reaction':
            info = line[1]
            reactants = line[2]
            products = line[3]
            assert len(info) == 3
            rate = float(info[1][0])
            if sysunit and [x for x in info[2][0].split('/')[1:-1] if x != sysunit]:
                raise PilFormatError('Conflicting units: {} vs. {}'.format(sysunit, info[2][0]))
            reactions.append([reactants, products, [rate]])
        elif line[0] == 'dl-domain' :
            pass
        else :
            logging.warning('Ignoring Keyword: {}'.format(line[0]))

    crnspecies = dict()
    condensed = bool(macrostates)
    for [reac,prod,rate] in reactions:
        for sp in reac + prod:
            if condensed:
                # If macrostates are specified, then they must contain all
                # species of the CRN.. otherwise, we'd have to make sure that
                # we don't return a mix of detailed and condensed CRN.
                if sp not in macrostates:
                    raise InputFormatError('Undefined macrostate in condensed input CRN: {}'.format(sp))
                crnspecies[sp] = macrostates[sp]
            elif sp not in crnspecies:
                crnspecies[sp] = species.pop(sp, ('initial', 0))

    if not condensed and bool(species):
        log.warning("Some species do not appear in the CRN: {}".format(list(species.keys())))

    return reactions, crnspecies

def read_seesaw(data, is_file = False, conc = 100e-9, explicit = True, reactions = ''):
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
        dlen = 10 # length of the domain
        tlen =  5 # lenth of the threshold

    complexes = {}
    def make_wire(w, name = None):
        """ Translates a wire into a PepperComplex.

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
            dIa = assgn_domain(f's{nI}a', l=dlen)
            dIt = assgn_domain(f's{nI}t', l=tlen)
            dOa = assgn_domain(f's{nO}a', l=dlen)
            dOt = assgn_domain(f's{nO}t', l=tlen)
            if explicit:
                sequence = [clp, dOt, dOa, clp, toe, clp, dIt, dIa, clp]
                structure = list('.........')
            else:
                sequence = [dOt, dOa, toe, dIt, dIa]
                structure = list('.....')
            complexes[wire] = PepperComplex(sequence, structure, 
                    name = name if name else wire)
            complexes[wire].concentration = ('i', 0, 'M')
        return complexes[wire]

    def make_thld(t):
        """ ['t', [inp, name]] """
        if t[0] != 't': 
            raise InputFormatError('Does not look like a threshold: {}'.format(t))
        [nI, nO] = t[1]
        thld = 'T{}_{}'.format(nI, nO)
        waste1 = 'Waste{}_{}'.format(nI, nO)
        waste2 = "Waste_{}".format(nO)
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
            complexes[thld].concentration = ('i', 0, 'M')

            # Name the waste as well, for convenience...
            assert waste1 not in complexes
            if explicit:
                sequence = [clp, dOt, dOa, clp, toe, clp, dIt, dIa, clp, '+', ~dIt, ~clp, ~toe, ~clp, ~dOa, ~dOt, ~clp]
                structure = list('(((((((..+)))))))')
            else :
                sequence = [dOt, dOa, toe, dIt, dIa, '+', ~dIt, ~toe, ~dOa, ~dOt]
                structure = list('((((.+))))')
            complexes[waste1] = PepperComplex(sequence, structure, name=waste1)
            complexes[waste1].concentration = ('i', 0, 'M')

            if waste2 not in complexes:
                if explicit:
                    raise NotImplementedError
                else :
                    complexes[waste2] = PepperComplex([dOt, dOa], ['.','.'], name=waste2)
                    complexes[waste2].concentration = ('i', 0, 'M')

        return complexes[thld], complexes[waste1], complexes[waste2]

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
            complexes[gate].concentration = ('i', 0, 'M')
        return complexes[gate]

    for line in parsed_file :
        if line[0] == 'INPUT':
            # use the given name if it is not a digit ...
            name = None if line[1][0].isdigit() else line[1][0]
            w = make_wire(line[2], name)
            assert w.concentration[1] == 0
            w.concentration = None

        elif line[0] == 'OUTPUT':
            # use the given name if it is not a digit ...
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
                    complexes[fluor].concentration = ('i', 0, 'M')
            else:
                raise InputFormatError('Unknown output format: {}'.format(line))

        elif line[0] == 'seesaw':
            # seesaw[5,{2},{6,7}]
            # rxn w(2-5) + g(5,5-6) <=> w(5-6) + g(2-5,5)
            #
            # produces all species for a custom seesaw gate, but with ZERO
            # concentration. If no concentrations are specified later via conc
            # statements), then this species is initially not present.
            #
            [gate, inp, out] = line[1]
            for d in inp:
                make_wire(['w', [d, gate]])                 # input wire
                make_thld(['t', [d, gate]])                 # threshold cplx (consuming input)
                make_gate(['g', [['w', [d, gate]], gate]])  # consumed gate (producing input)

            for d in out:
                make_wire(['w', [gate, d]])                 # output wire
                make_gate(['g', [gate, ['w',[gate, d]]]])   # active gate (producing output)

        elif line[0] in ('seesawOR', 'seesawAND'):
            # seesaw[2,5,{1,4},{3,9,10}]
            [g1, g2, inp, out] = line[1]
            n = len(inp)
            m = len(out)

            # make gates for all different inputs
            for d in inp:
                Iw = make_wire(['w', [d, g1]])
                Ig = make_gate(['g', [g1, ['w', [g1, g2]]]])

                Ow = make_wire(['w', [g1, g2]])
                Og = make_gate(['g', [['w', [d, g1]], g1]]) # consumed gate

            # Thresholding reaction
            Iw = make_wire(['w',[g1, g2]])
            Th,w1,w2 = make_thld(['t',[g1, g2]]) 
            if line[0] == 'seesawOR':
                assert Th.concentration[1] == 0
                Th.concentration = ('i', 1.1 * 0.6 * conc, 'M')
            else:
                assert Th.concentration[1] == 0
                Th.concentration = ('i', 1.1 * (n-1 + 0.2) * conc, 'M')

            # Connection reaction (assign intial concentration)
            Ig = make_gate(['g', [g1, ['w', [g1, g2]]]])
            assert Ig.concentration[1] == 0
            Ig.concentration = ('i', n * conc, 'M')

            # fanout for all outputs
            for d in out + ['f']:
                Iw = make_wire(['w',[g1, g2]])
                Ig = make_gate(['g', [g2, ['w',[g2, d]]]])

                Ow = make_wire(['w', [g2, d]])
                Og = make_gate(['g', [['w', [g1, g2]], g2]]) #consumed gate

                if d == 'f':
                    assert Ow.concentration[1] == 0
                    fmod = Ow.concentration[0]
                    fcon = Ow.concentration[1]
                    funi = Ow.concentration[2]
                    fcon += 2 * m * conc
                    Ow.concentration = (fmod, fcon, funi)
                else :
                    Ig.concentration = ('i', conc, 'M')

        elif line[0] == 'inputfanout':
            [g1, inp, out] = line[1]
            n = len(inp)
            m = len(out)

            make_wire(['w', [inp, g1]])
            cx,_,_ = make_thld(['t', [inp, g1]])
            assert cx.concentration[1] == 0
            cx.concentration = ('i', 1.1 * 0.2 * conc, 'M')
            make_gate(['g', [['w',[inp, g1]], g1]])

            for d in out + ['f']:
                cxw = make_wire(['w', [g1, d]])
                cxg = make_gate(['g', [g1, ['w', [g1, d]]]])
                if d == 'f':
                    assert cxw.concentration[1] == 0
                    cxwmod = Ow.concentration[0]
                    cxwcon = Ow.concentration[1]
                    cxwuni = Ow.concentration[2]
                    cxwcon += 2 * m * conc
                    cxw.concentration = (cxwmod, cxwcon, cxwuni)
                else :
                    assert cxg.concentration[1] == 0
                    cxg.concentration = ('i', conc, 'M')

        elif line[0] == 'reporter':
            #['reporter', ['25', '31']]
            [nR, nI] = line[1]
            sw = make_wire(['w', [nI, nR]])

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
            complexes[rep].concentration = ('i', 1.5 * conc, 'M')

            # TODO: swap F and Q!!!
            flr = "F_{}".format(nR)
            if flr not in complexes:
                if explicit:
                    complexes[flr] = PepperComplex([clp, dRt, dRa, clp], 
                                                   ['.','.','.','.'], name=flr)
                else :
                    complexes[flr] = PepperComplex([dRt, dRa], ['.','.'], name=flr)
                complexes[flr].concentration = ('i', 0, 'M')

            que = "Q_{}".format(nR)
            if que not in complexes:
                dIa = assgn_domain('s{}a'.format(nI), l=dlen)
                dIt = assgn_domain('s{}t'.format(nI), l=tlen)
                if explicit:
                    raise NotImplementedError
                else :
                    sequence  = [dRt, dRa, toe, dIt, dIa, '+', ~toe, ~dRa, ~dRt]
                    structure = list('(((..+)))')
                    complexes[que] = PepperComplex(sequence, structure, name=que)
                    complexes[que].concentration = ('i', 0, 'M')

        elif line[0] == 'conc':
            if line[1][0] == 'w':
                cx = make_wire(line[1])
            elif line[1][0] == 'g':
                cx = make_gate(line[1])
            elif line[1][0] == 'th':
                tmp  = line[1][1][0]
                tmp[0] = 't'
                cx,_,_ = make_thld(tmp)
            else:
                raise NotImplementedError('Cannot interpret line {}'.format(line))
            # input wires have None at default
            assert cx.concentration is None or cx.concentration[1] == 0
            cx.concentration = ('i', conc*float(line[2]), 'M')

        else :
            log.warning('Keyword not supported: {}'.format(line[0]))

    if reactions: # e.g. 'seesaw-T25-utbr-leak-reduced'
        for term in reactions.split('-'):
            assert term in ('seesaw', 'T20', 'T25', 'utbr', 'leak', 'reduced')
        T = 25 if 'T25' in reactions else 20
        utbr = True if 'utbr' in reactions else False
        leak = True if 'leak' in reactions else False
        reduced = True if 'reduced' in reactions else False
        reactions = get_seesaw_compiler_reactions(complexes, T, utbr, leak, reduced)
    else :
        reactions = []

    return complexes, reactions

def get_seesaw_compiler_reactions(complexes, T = 20, utbr = False, leak = False, global_species = False):
    """Write the standard seesaw reactions. 

    Takes a dictionary of seesaw complexes and finds all designed and side
    reactions described in Qian & Winfree 2011, Supplemental Section S13

    Args:
        complexes (dict, optional): Maps seesaw names to complexes. Note, if
            complexes got renamed due to INPUT/OUTPUT statements, then
            complex.name is not the same as the dictionary key.
        T (int, optional): Temperature, can be either 20 or 25 *C.
        utbr (bool, optional): Return universal toehold binding reactions.
        leak (bool, optional): Return leak reactions.
        global_species (bool, optional): Use global Threshold, Wire and Gate
            species to reduce the number of explicit toehold binding reactions.
    
    """
    PepperReaction.RTYPES.add('seesaw')      # desinged seesawing
    PepperReaction.RTYPES.add('seesaw-thld') # designed thresholding
    PepperReaction.RTYPES.add('seesaw-rep')  # designed reporting
    PepperReaction.RTYPES.add('side-utbr')   # side reaction: universal toehold binding 
    PepperReaction.RTYPES.add('side-leak')   # side reaction: leak

    ks  = 5e4 # /M/s
    kf  = 2e6 # /M/s
    if T == 20:
        krs = 0.5 # /s
        krf = 10  # /s
        kl  = 1   # /M/s
    elif T == 25:
        krs = 1.3 # /s
        krf = 26  # /s
        kl  = 10  # /M/s
    else :
        raise NotImplementedError('Please use temperature 20 or 25.')

    # complexes[k] = PepperComplex ... k != complex.name 
    wires = []
    gates = []
    thlds = []
    freps = []
    for k,x in list(complexes.items()):
        if k[0] == 'w':
            x.seesawname = k
            wires.append(x)
        elif k[0]=='G': 
            gates.append(x)
        elif k[0]=='R':
            freps.append(x)
        elif k[0]=='T':
            thlds.append(x)

    count = 0
    reactions = []
    if utbr and global_species:
        # Make the universal wire W
        # This is the fuel value... 
        conc = sum([x.concentration[1] if x.concentration else 0 for x in wires])
        if conc:
            uW = PepperDomain('dummy_{}'.format('W'), dtype='long')
            uW = PepperComplex(sequence=[uW], structure=['.'], name='W')
            uW.concentration = ('initial', conc, 'M')
            complexes[uW.name] = uW

        # Make the universal gate G (Gates and Reporters)
        uG = PepperDomain('dummy_{}'.format('G'), dtype='long')
        uG = PepperComplex(sequence=[uG], structure=['.'], name='G')
        conc = sum([x.concentration[1] if x.concentration else 0 for x in gates + freps])
        uG.concentration = ('initial', conc, 'M')
        complexes[uG.name] = uG

        # Make the universal threshold TH 
        uT = PepperDomain('dummy_{}'.format('TH'), dtype='long')
        uT = PepperComplex(sequence=[uT], structure=['.'], name='TH')
        conc = sum([x.concentration[1] if x.concentration else 0 for x in thlds])
        uT.concentration = ('initial', conc, 'M')
        complexes[uT.name] = uT

    for w in wires:
        wi, wo = w.seesawname[1:].split('_')

        # Find all reactions with gates:
        tot_G_conc = 0
        for g in gates:

            G, gn, gi, go = str(g).split('_')
            if gn[0] == 'g': # a forward gate G5:5,6
                assert gi[0] == 'w'
                gw = complexes[gi+'_'+go]
                G = '+'
            else : # a reverse gate G5,6:6
                assert gn[0] == 'w'
                assert go[0] == 'g'
                [gi, go, gn] = [gn, gi, go]
                gw = complexes[gi+'_'+go]
                G = '-'

            gwi, gwo = gw.seesawname[1:].split('_')
            gn = gn[1:]

            if G == '+' and wo == gn:
                # a regular interaction
                og = complexes['G_' + w.seesawname + '_' + 'g' + gn]
                rxn = PepperReaction([w, g], [gw, og], rtype = 'seesaw')
                rxn.rate_constant = (ks, '/M/s')
                reactions.append(rxn)
                continue

            if G == '-' and wi == gn:
                # a regular interaction
                og = complexes['G_' + 'g' + gn + '_' + w.seesawname]
                rxn = PepperReaction([w, g], [gw, og], rtype = 'seesaw')
                rxn.rate_constant = (ks, '/M/s')
                reactions.append(rxn)
                continue

            if leak and G == '+' and wi == gwi:# and wo != gwo:
                og = complexes['G_' + 'g' + gn + '_' + w.seesawname]
                rxn = PepperReaction([w, g], [gw, og], rtype = 'side-leak')
                rxn.rate_constant = (kl, '/M/s')
                reactions.append(rxn)

            if leak and G == '-' and wo == gwo:# and wi != gwi:
                og = complexes['G_' + w.seesawname + '_' + 'g' + gn]
                rxn = PepperReaction([w, g], [gw, og], rtype = 'side-leak')
                rxn.rate_constant = (kl, '/M/s')
                reactions.append(rxn)

            if utbr and global_species:
                # TODO: check if that makes sense...?
                if w.concentration and w.concentration[1] != 0:
                    #try : # all Wires together
                    gWn = g.name + '_W'
                    gW = PepperDomain('dummy_{}'.format(gWn), dtype='long')
                    gW = PepperComplex(sequence=[gW], structure=['.'], name=gWn)
                    gW.concentration = ('initial', 0, 'M')
                    complexes[gWn] = gW

                    rxn = PepperReaction([g, uW], [gW], rtype = 'side-utbr')
                    rxn.rate_constant = (kf, '/M/s')
                    reactions.append(rxn)
                    rxn = PepperReaction([gW], [g, uW], rtype = 'side-utbr')
                    rxn.rate_constant = (krf, '/s')
                    reactions.append(rxn)
                    #except DSDDuplicationError as err:
                    #    pass
                else :
                    #try: # all Gates together
                    wGn = w.name + '_G'
                    wG = PepperDomain('dummy_{}'.format(wGn), dtype='long')
                    wG = PepperComplex(sequence=[wG], structure=['.'], name=wGn)
                    wG.concentration = ('initial', 0, 'M')
                    complexes[wGn] = wG

                    rxn = PepperReaction([w, uG], [wG], rtype = 'side-utbr')
                    rxn.rate_constant = (kf, '/M/s')
                    reactions.append(rxn)
                    rxn = PepperReaction([wG], [w, uG], rtype = 'side-utbr')
                    rxn.rate_constant = (krf, '/s')
                    reactions.append(rxn)
                    #except DSDDuplicationError as err:
                    #    pass

            elif utbr :
                count += 1
                a = PepperDomain('dummy_{}'.format(count), dtype='long')
                a = PepperComplex(sequence=[a], structure=['.'], prefix='s')
                complexes[a.name] = a
                rxn = PepperReaction([w, g], [a], rtype = 'side-utbr')
                rxn.rate_constant = (kf, '/M/s')
                reactions.append(rxn)
                rxn = PepperReaction([a], [w, g], rtype = 'side-utbr')
                rxn.rate_constant = (krf, '/s')
                reactions.append(rxn)
        
        for g in freps:
            # R_6
            if '_{}'.format(wo) in str(g):
                # a regular interaction
                of = complexes['F_'+wo]
                oq = complexes['Q_'+wo]
                rxn = PepperReaction([w, g], [of, oq], rtype = 'seesaw-rep')
                rxn.rate_constant = (ks, '/M/s')
                reactions.append(rxn)
                continue

            if utbr and global_species:
                # TODO: check if that makes sense...?
                if w.concentration and w.concentration[1] != 0:
                    #try: # all-wires product
                    gWn = g.name + '_W'
                    gW = PepperDomain('dummy_{}'.format(gWn), dtype='long')
                    gW = PepperComplex(sequence=[gW], structure=['.'], name=gWn)
                    gW.concentration = ('initial', 0, 'M')
                    complexes[gWn] = gW

                    rxn = PepperReaction([g, uW], [gW], rtype = 'side-utbr')
                    rxn.rate_constant = (kf, '/M/s')
                    reactions.append(rxn)
                    rxn = PepperReaction([gW], [g, uW], rtype = 'side-utbr')
                    rxn.rate_constant = (krf, '/s')
                    reactions.append(rxn)
                    #except DSDDuplicationError as err:
                    #    pass
                else:
                    #try: # all Gates together
                    wGn = w.name + '_G'
                    wG = PepperDomain('dummy_{}'.format(wGn), dtype='long')
                    wG = PepperComplex(sequence=[wG], structure=['.'], name=wGn)
                    wG.concentration = ('initial', 0, 'M')
                    complexes[wGn] = wG

                    rxn = PepperReaction([w, uG], [wG], rtype = 'side-utbr')
                    rxn.rate_constant = (kf, '/M/s')
                    reactions.append(rxn)
                    rxn = PepperReaction([wG], [w, uG], rtype = 'side-utbr')
                    rxn.rate_constant = (krf, '/s')
                    reactions.append(rxn)
                    #except DSDDuplicationError as err:
                    #    pass

            elif utbr:
                count+=1
                a = PepperDomain('dummy_{}'.format(count), dtype='long')
                a = PepperComplex(sequence=[a], structure=['.'], prefix='s')
                complexes[a.name] = a
                rxn = PepperReaction([w, g], [a], rtype = 'side-utbr')
                rxn.rate_constant = (kf, '/M/s')
                reactions.append(rxn)
                rxn = PepperReaction([a], [w, g], rtype = 'side-utbr')
                rxn.rate_constant = (krf, '/s')
                reactions.append(rxn)

        for g in thlds:
            # T2_5
            if 'T{}_{}'.format(wi, wo) == str(g):
                tc = complexes['Waste{}_{}'.format(wi, wo)]
                tw = complexes['Waste_'+wo]
                rxn = PepperReaction([w, g], [tc, tw], rtype = 'seesaw-thld')
                rxn.rate_constant = (kf, '/M/s')
                reactions.append(rxn)
                continue

            if utbr and global_species:
                # TODO: check if that makes sense...?
                if w.concentration and w.concentration[1] != 0:
                    #try: # all Wires together
                    gWn = g.name + '_W'
                    gW = PepperDomain('dummy_{}'.format(gWn), dtype='long')
                    gW = PepperComplex(sequence=[gW], structure=['.'], name=gWn)
                    gW.concentration = ('initial', 0, 'M')
                    complexes[gWn] = gW

                    rxn = PepperReaction([g, uW], [gW], rtype = 'side-utbr')
                    rxn.rate_constant = (kf, '/M/s')
                    reactions.append(rxn)
                    rxn = PepperReaction([gW], [g, uW], rtype = 'side-utbr')
                    rxn.rate_constant = (krs, '/s')
                    reactions.append(rxn)
                    #except DSDDuplicationError as err:
                    #    pass
                else:
                    #try: # all Gates together
                    wGn = w.name + '_TH'
                    wG = PepperDomain('dummy_{}'.format(wGn), dtype='long')
                    wG = PepperComplex(sequence=[wG], structure=['.'], name=wGn)
                    wG.concentration = ('initial', 0, 'M')
                    complexes[wGn] = wG

                    rxn = PepperReaction([w, uT], [wG], rtype = 'side-utbr')
                    rxn.rate_constant = (kf, '/M/s')
                    reactions.append(rxn)
                    rxn = PepperReaction([wG], [w, uT], rtype = 'side-utbr')
                    rxn.rate_constant = (krs, '/s')
                    reactions.append(rxn)
                    #except DSDDuplicationError as err:
                    #    pass

            elif utbr:
                count+=1
                a = PepperDomain('dummy_{}'.format(count), dtype='long')
                a = PepperComplex(sequence=[a], structure=['.'], prefix='s')
                complexes[a.name] = a
                rxn = PepperReaction([w, g], [a], rtype = 'side-utbr')
                rxn.rate_constant = (kf, '/M/s')
                reactions.append(rxn)
                rxn = PepperReaction([a], [w, g], rtype = 'side-utbr')
                rxn.rate_constant = (krs, '/s')
                reactions.append(rxn)

    return reactions
