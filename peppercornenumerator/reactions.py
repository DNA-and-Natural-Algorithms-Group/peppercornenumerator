#
#  reactions.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/18/2010.
#  Modifications by Casey Grun and Erik Winfree 8/15/2014.
#  Modifications by Stefan Badelt 09/2019

import copy

from peppercornenumerator.utils import Loop, wrap

from peppercornenumerator.objects import PepperComplex, PepperReaction, DSDDuplicationError
from peppercornenumerator.objects import make_pair_table, pair_table_to_dot_bracket


auto_name = 0

def get_auto_name(prefix=''):
    """
    Returns a new unique name
    """
    global auto_name
    auto_name += 1
    return prefix + str(auto_name)

# Rate constant formulae
# ----------------------------------------------------------------------------


def opening_rate(length):
    """
    Rate constant formula for opening a duplex of a given `length`.
    """
    # use k_open = k_hybrid exp( (length * dG_bp + dG_assoc) / RT )
    # where k_hybrid = 3x10^6 /M/s   from Zhang&Winfree 2009 and Srinivas et al 2013
    #       dG_bp = -1.7 kcal/mol
    #       dG_assoc = +1.9 kcal/mol
    #       R = 0.001987 kcal/mol/K
    #       T = (273.15 + 25) K
    # return 7.41e7 * (0.0567 ** length)
    #
    # instead, use k_hybrid = L * 3 * 10^5, which matches the above for L=10.
    # this is to be consistent with the bimolecular binding rate.
    return length * 7.41e6 * (0.0567 ** length)


def polymer_link_length(before, after):
    """
    Effective length estimate for (ss+ds) linkers between two domains, 
    one or both of which may be open.
    """
    L_stem = 2.0 / 0.43  # rough equivalent number of single-stranded nucleotides to span a stem
    if not before.is_open:
        L_before = 1 + before.bases + L_stem * before.stems
    if not after.is_open:
        L_after = 1 + after.bases + L_stem * after.stems
    # for both closed & open cases, assume shorter length matters most
    if not after.is_open and not before.is_open:
        return min(L_before, L_after)
    if not before.is_open:
        return L_before
    if not after.is_open:
        return L_after
    assert False, "should not have reached this case -- how can both sides be open?"
    raw_input("bad bad bad -- computing polymer lengths in disconnected complex!")


def polymer_link_rate(linker_length):
    """
    Unimolecular hairpin closing rate, as a function of effective linker length. 
    """
    # a = 2.5e7    # per second; fit from data in Bonnet et al 1998 only
    a = 1e6        # per second; Kuznetsov et al 2008, Nayak et al 2012, Tsukanov et al 2013ab, all say at least 10x slower
    # fit from data in Bonnet et al 1998, consistent with Kuznetsov et al 2008
    b = -2.5
    # hairpin closing adapted (simplified) from data in Bonnet et al 1998, modified as above
    k = a * linker_length**b
    # hairpins shorter than about 4 nt can't close any faster.
    return min(k, 33000)


def binding_rate(length, before, after):
    """
    Rate constant formula for unimolecular binding of a domain of the given length.
    Could be zippering, hairpin closing, bubble closing, bulge formation, multiloop formation,
    depending on the flanking loops, which may be open or closed.
    """
    if not before.is_open and before.stems == 1 and before.bases == 0:
        if not after.is_open and after.stems == 1 and after.bases == 0:
            return 1e4  # bubble closing rate from Altan-Bonnet 2003
        # zippering from Wetmur&Davidson 1968, Gueron&Leroy 1995, Srinivas et al 2013, low end
        return (1e6) / length
    if not after.is_open and after.stems == 1 and after.bases == 0:
        # zippering from Wetmur&Davidson 1968, Gueron&Leroy 1995, Srinivas et al 2013, low end
        return (1e6) / length

    # bulge closing assumed to be similar to faster of two hairpin closings
    L = polymer_link_length(before, after)
    # hairpin closing adapted (simplified) from data in Bonnet et al 1998
    return polymer_link_rate(L)

# Diagram for 3-way branch migration, general case.
# Loops could be listed either 5'->3' or 3'->5, but they always go from
# invading domain to bound stem (non-inclusive).
#
#                before
#                _______  x (bound domain)
#               /       \____
# (invading) x |         ____
#               \_______/ x*
#
#                 after


def show_loops(before, after, message):
    """
    Debugging help: show (partial) loops returned by find_on_loop().  

    ! indicates a stem, | indicates open loop break.  
    """
    print "before: [ ",
    for step in before.parts:
        print " | " if step is None else step[0].name + ("!" if step[1] is not None else "") + " ",
    print " ] is_open = %r, stems = %d, bases = %d" % (before.is_open, before.stems, before.bases)
    print "after: [ ",
    for step in after.parts:
        print " | " if step is None else step[0].name + ("!" if step[1] is not None else "") + " ",
    print " ] is_open = %r, stems = %d, bases = %d" % (after.is_open, after.stems, after.bases)
    raw_input(message)


def branch_3way_remote_rate(length, before, after, debug = False):
    """
    Rate constant formula for 3-way branch migration, possibly with a remote toehold
    """
    # step = 0.1e-3  # 100us branch migration step time from Srinivas et al 2013 (not relevant)
    # k_init = k_uni exp(-dGsp / RT) with k_uni = 7.5e7, dGsp = 7.3 kcal/mol,
    # T = 273.15 + 25 K, R = 0.001987 kcal/mol/K
    init = 3.0e-3  # sec, = 1/k_init from Srinivas et al 2013

    if debug :
        show_loops(before, after, "...before & after loops for 3-way branch migration...")

    # "standard" 3-way bm initiation (plus "before" being closed)
    if not after.is_open and after.stems == 1 and after.bases == 0:
        # each initiation has probability 1/length of succeeding.  
        # how long it takes doesn't matter.
        return 1.0 / init / length
    
    if debug:
        show_loops(before, after, 
                "run_tests.py should not have any remote toeholds for 3-way branch migration")

    # consider a slowdown analogous to Genot et al 2011 (remote) supp info derivation
    # bulge closing assumed to be similar to faster of two hairpin closings
    L = polymer_link_length(before, after)
    # how much slower than our (admittedly slow) zippering rate is this?
    ratio = polymer_link_rate(L) / (1e6)
    # we slow down initiation and thus success probability (note: ratio < 1/30)
    return ratio / init / length


def branch_4way_remote_rate(length, before, after, debug=False):
    """
    Rate constant formula for 4-way branch migration, possibly with a remote toehold
    """
    # rates recalculated from Nadine Dabby, PhD Thesis, 2013, based on
    # assuming last 6 bp dissociate faster than 4way bm step
    open_step = 107   # sec, = 1/k_first  (this is for open 4-way case only)
    # sec, = 1/k_bm     (this is used for initiating closed 4-way case;
    # consistent with Panyutin&Hsieh 1993)
    closed_step = 3.0

    # open_step = 200 # fudge !
    
    if debug:
        show_loops(before, after, "before & after loops for 4-way branch migration")

    if not before.is_open and not after.is_open:
        init = closed_step
        if before.bases == 0 and before.stems == 1 and after.bases == 0 and after.stems == 1:
            return 1.0 / init / length   # closed & ready-to-rock-and-roll 4 way initiation
    if before.is_open:
        init = open_step
        if after.bases == 0 and after.stems == 1:
            # we care about probability of taking this path, not actual time
            return 1.0 / init / length
    if after.is_open:
        init = open_step
        if before.bases == 0 and before.stems == 1:
            # we care about probability of taking this path, not actual time
            return 1.0 / init / length

    if debug: 
        show_loops(before, after, 
                "run_tests.py should not have any remote toeholds for 4-way branch migration")

    # consider a slowdown analogous to Genot et al 2011 (remote) supp info derivation
    # bulge closing assumed to be similar to faster of two hairpin closings
    L = polymer_link_length(before, after)
    # how much slower than our (admittedly slow) zippering rate is this?
    ratio = polymer_link_rate(L) / (1e6)
    # we slow down initiation and thus success probability (note: ratio < 1/30)
    return ratio / init / length


def bimolecular_binding_rate(length):
    """
    Rate constant formula for bimolecular association (binding).
    """
    # use k_hybrid = 3x10^6 /M/s   from Zhang&Winfree 2009
    # return 3.0e6
    #
    # instead, use k_hybrid = L * 3 * 10^5, which matches the above for L=10.
    # see Wetmur 1976 review, and Srinivas et al 2013 AEL model.
    # another motivation is to have binding rate approx = if a domain is
    # divided into two domains.
    return length * 3e5


# Reaction functions
# ----------------------------------------------------------------------------

def bind11(reactant, max_helix=True, remote=None):
    """
    Returns a list of reaction pathways which can be produced by 1-1 binding
    reactions of the argument complex. The 1-1 binding reaction is the
    hybridization of two complementary unpaired domains within a single complex
    to produce a single unpseudoknotted product complex.

    Note: Remote is ineffective, but may be set for convencience
    """
    
    # Filter function for find_on_loop()
    def filter_bind11(dom1, struct1, loc1, dom2, struct2, loc2, freed):
        return struct1 is None and struct2 is None and dom2.can_pair(dom1)

    reactions = []
    structure = reactant.pair_table

    # We iterate through all the domains
    for (strand_index, strand) in enumerate(structure):
        for (domain_index, domain) in enumerate(strand):

            # The displacing domain must be free
            if structure[strand_index][domain_index] is not None :
                continue

            loc1 = (strand_index, domain_index)

            # search both directions around the loop for a bound domain that
            # has the same sequence (and therefore can be displaced)
            locs = find_on_loop(reactant, loc1, -1, filter_bind11, max_helix=max_helix) + \
                   find_on_loop(reactant, loc1, +1, filter_bind11, max_helix=max_helix)  

            # build products
            for (loc1s, loc2s, before, after) in locs:
                product = do_bind11(reactant, loc1s.locs, loc2s.locs)

                try :
                    reaction = PepperReaction([reactant], [product], 'bind11')
                except DSDDuplicationError, e :
                    reaction = e.existing

                # length of invading domain
                length = len(loc1s)

                # calculate reaction constant
                reaction.rate = binding_rate(length, before, after)

                reactions.append(reaction)

    # remove any duplicate reactions
    return sorted(list(set(reactions)))

def do_bind11(reactant, loc1s, loc2s):
    """ Returns PepperComplex after the bind11 reaction. """
    struct = reactant.pair_table
    for loc1, loc2 in zip(loc1s, loc2s):
        struct[loc1[0]][loc1[1]] = loc2
        struct[loc2[0]][loc2[1]] = loc1
    newstr = pair_table_to_dot_bracket(struct)
    try:
        product = PepperComplex(reactant.sequence, newstr)
        product.pair_table = struct
    except DSDDuplicationError, e:
        product = e.existing
    return product

def bind21(reactant1, reactant2, max_helix = True, remote=None):
    """
    Returns a list of reaction pathways which can be produced by 2-1 binding
    reactions of the argument complexes. The 2-1 binding reaction is the
    hybridization of two complementary unpaired domains, each in a different
    complex, to produce a single, unpseudoknotted product complex containing
    all of the strands contained in either of the original complexes.

    Note: remote is ineffective, but may be set for convencience
    """
    r1_doms = reactant1.available_domains
    r2_doms = reactant2.available_domains

    def filter_bind21(dom1, struct1, loc1, dom2, struct2, loc2, freed):
        return struct1 is None and struct2 is None and dom1.can_pair(dom2)

    reactions = []

    # Iterate through all the free domains in reactant1
    for (dom1, strand_num1, dom_num1) in r1_doms:
        # For each, find any domains in reactant2 that could bind
        for (dom2, strand_num2, dom_num2) in r2_doms:
            # If it can pair, this is one possible reaction (this kind of
            # reaction cannot possibly produce a pseudoknotted structure)
            if (dom1.can_pair(dom2)):
                # combine the two complexes into one, but do not perform the association
                reactions.append(join_complexes_21(
                    reactant1, (strand_num1, dom_num1),
                    reactant2, (strand_num2, dom_num2)))
    
    output = []
    for complex, location1, location2 in reactions:
        # build "before" and "after" loop structures via find_on_loop ...
        out = find_on_loop(complex, location1, 1,
            lambda dom1, struct1, loc1, dom2, struct2, loc2, freed: 
                loc1 == location1 and loc2 == location2, max_helix=False)

        [(loc1s, loc2s, before, after)] = out

        # zipper for max-helix semantics
        if max_helix :
            (loc1s, loc2s, before, after) = zipper(complex, location1, location2, 
                    before.parts, after.parts, 1, None, None, filter_bind21)

        product = do_bind11(complex, loc1s.locs, loc2s.locs)

        try :
            reaction = PepperReaction(sorted([reactant1, reactant2]), [product], 'bind21')
        except DSDDuplicationError, e :
            #assert opening_rate(length) == PepperReaction.dictionary[e.solution].rate
            reaction = e.existing

        length = len(loc1s)
        reaction.rate = bimolecular_binding_rate(length)

        output.append(reaction)

    return sorted(list(set(output)))

def join_complexes_21(complex1, location1, complex2, location2):
    """
    Combines two complexes to form one complex, binding the domain in
    complex1 at location1 to the domain in complex2 at location2.

    Returns the new complex.
    """

    # make sure that this value is larger than any possible loc in the pair_table
    maxlen = complex1.size + complex2.size + 1

    if complex1.get_loop_index(location1) == 0 :
        seq1 = complex1.sequence
        ptb1 = complex1.pair_table
        loc1 = location1
        ptb1[loc1[0]][loc1[1]] = (maxlen,maxlen) # add an additional '(' 
    else :
        seen = False # don't break, otherwise you might lose the original strand ordering...
        for e, rot in enumerate(complex1.rotate(),1):
            tmp = complex1.rotate_location(location1, -e)
            if not seen and complex1.get_loop_index(tmp) == 0 :
                seq1 = complex1.sequence
                ptb1 = complex1.pair_table
                loc1 = tmp
                ptb1[loc1[0]][loc1[1]] = (maxlen,maxlen) # add an additional '(' 
                seen = True

    if complex2.get_loop_index(location2) == 0 :
        seq2 = complex2.sequence
        ptb2 = complex2.pair_table
        loc2 = location2
        ptb2[loc2[0]][loc2[1]] = (-1,-1) # add an additional ')' 
    else :
        seen = False # don't break, otherwise you might lose the original strand ordering...
        for e, rot in enumerate(complex2.rotate(),1):
            tmp = complex2.rotate_location(location2, -e)
            if not seen and complex2.get_loop_index(tmp) == 0 :
                seq2 = complex2.sequence
                ptb2 = complex2.pair_table
                loc2 = tmp
                ptb2[loc2[0]][loc2[1]] = (-1,-1) # add an additional ')' 
                seen = True

    # build the new sequence and structure *including* the new pair
    newseq = seq1 + ['+'] + seq2
    newstr = pair_table_to_dot_bracket(ptb1) + ['+'] + pair_table_to_dot_bracket(ptb2)

    # update loc2 from the new structure
    combined = make_pair_table(newstr)
    loc2 = combined[loc1[0]][loc1[1]]

    # remove the new pair again ...
    combined[loc1[0]][loc1[1]] = None
    combined[loc2[0]][loc2[1]] = None

    # update the structure to the unpaired (disconnected) version.
    newstr = pair_table_to_dot_bracket(combined)

    try :
        new_complex = PepperComplex(newseq, newstr)
    except DSDDuplicationError, e:
        new_complex = e.existing
        # strands may be re-ordered in new complex, so we need to back
        # out where the new strands ended up
        loc1 = new_complex.rotate_location(loc1, e.rotations)
        loc2 = new_complex.rotate_location(loc2, e.rotations)

    if loc1 > loc2:
        (loc1, loc2) = (loc2,loc1)

    return new_complex, loc1, loc2

def open(reactant, max_helix = True, release_11=6, release_1N=6):
    """
    Returns a list of reaction product sets that can be produced by the
    'open' reaction, in which a short helix dissociates. Each product
    set are the results of one particular dissociation; each strand in the
    reactant occurs exactly once in one of the complexes in the product set.

    A dissociation can happen to any helix under the threshold length

    """

    # remember the larger release cutoff; don't enumerate any reactions
    # for helices longer than this
    max_release = max(release_11, release_1N)

    reactions = []
    structure = reactant.pair_table

    # for no-max-helix mode:
    if not max_helix:
        # We iterate through all the domains
        for (strand_index, strand) in enumerate(structure):
            for (domain_index, domain) in enumerate(strand):

                # The bound domain must be... bound
                if reactant.get_structure((strand_index, domain_index)) is None :
                    continue

                bound_loc = (strand_index, domain_index)
                bound_domain = reactant.get_domain(bound_loc)

                release_reactant = do_single_open(reactant, bound_loc)
                product_set = release_reactant.split()
                reactions.append((product_set, len(bound_domain)))

    # for max-helix mode:
    else:
        # We loop through all stands, domains
        for (strand_index, strand) in enumerate(structure):
            for (domain_index, domain) in enumerate(strand):
                # If the domain is unpaired, skip it
                if (structure[strand_index][domain_index] is None):
                    continue

                # A: Strand/domain position on "top" strand - CG 5/21
                helix_startA = [strand_index, domain_index]

                # B: Strand/domain position on "bottom" strand - CG 5/21
                helix_startB = list(structure[strand_index][domain_index])

                # If the domain is bound to an earlier domain, then we have
                # already considered it, so skip it
                if (helix_startB < helix_startA):
                    continue

                helix_endA = helix_startA[:]
                helix_endB = helix_startB[:]

                bound_loc = (strand_index, domain_index)
                bound_domain = reactant.get_domain(bound_loc)
                helix_length = len(bound_domain)

                # Now iterate through the whole helix to find the other end of
                # this one (The helix ends at the first strand break from
                # either direction)
                while True:
                    # Strands run in opposite directions, so A must be incremented
                    # and B decremented in order that both pointers move "right"
                    # along the helix- CG 5/21
                    helix_endA[1] += 1
                    helix_endB[1] -= 1

                    # If one of the strands has broken, the helix has ended
                    if helix_endA[1] >= reactant.strand_length(helix_endA[0]) :
                        break
                    elif helix_endB[1] < 0 :
                        break

                    # If these domains aren't bound to each other, the helix has ended
                    if (tuple(helix_endA) != structure[helix_endB[0]][helix_endB[1]]):
                        break

                    # Add the current domain to the current helix
                    temp_bl = (helix_endA[0], helix_endA[1])
                    helix_length += len(reactant.get_domain(temp_bl))

                # We must also iterate in the other direction
                while True:
                    helix_startA[1] -= 1
                    helix_startB[1] += 1

                    # If one of the strands has broken, the helix has ended
                    if helix_startA[1] < 0 :
                        break
                    elif helix_startB[1] >= reactant.strand_length(helix_startB[0]):
                        break

                    # If these domains aren't bound to each other, the helix has ended
                    if tuple(helix_startA) != structure[helix_startB[0]][helix_startB[1]] :
                        break

                    # Add the current domain to the current helix
                    temp_bl = (helix_startA[0], helix_startA[1])
                    helix_length += len(reactant.get_domain(temp_bl))

                # Move start location to the first domain in the helix
                helix_startA[1] += 1
                helix_startB[1] -= 1

                # If the helix is short enough, we have a reaction
                if (helix_length <= max_release):

                    new_struct = reactant.pair_table

                    # Delete all the pairs in the released helix
                    for dom in range(helix_startA[1], helix_endA[1]):
                        bound_loc = reactant.get_structure((helix_startA[0], dom))
                        new_struct[helix_startA[0]][dom] = None
                        new_struct[bound_loc[0]][bound_loc[1]] = None

                    newstr = pair_table_to_dot_bracket(new_struct)

                    try:
                        release_reactant = PepperComplex(reactant.sequence, newstr)
                    except DSDDuplicationError, e:
                        release_reactant = e.existing

                    product_set = release_reactant.split()
                    reactions.append((product_set, helix_length))

    output = []
    for product_set, length in reactions:
        try :
            reaction = PepperReaction([reactant], sorted(product_set), 'open')
        except DSDDuplicationError, e :
            reaction = e.existing

        # discard reactions where the release cutoff is greater than the threshold
        if len(reaction.products) == 1 and length > release_11:
            continue
        elif len(reaction.products) > 1 and length > release_1N:
            continue

        reaction.rate = opening_rate(length)
        output.append(reaction)

    return sorted(list(set(output)))

def do_single_open(reactant, loc):
    struct = reactant.pair_table
    loc1 = loc
    loc2 = struct[loc1[0]][loc1[1]]
    assert struct[loc2[0]][loc2[1]] == loc1
    struct[loc1[0]][loc1[1]] = None
    struct[loc2[0]][loc2[1]] = None
    newstr = pair_table_to_dot_bracket(struct)
    try:
        product = PepperComplex(reactant.sequence, newstr)
        product.pair_table = struct
    except DSDDuplicationError, e:
        product = e.existing
    return product

def branch_3way(reactant, max_helix = True, remote=True):
    """
    Returns a list of reaction pathways that can be created through one
    iteration of a 3 way branch migration reaction (more than one molecule may
    be produced by a reaction because branch migration can liberate strands and
    complexes).
    """

    def filter_3way(dom1, struct1, loc1, dom2, struct2, loc2, freed):
        return (struct1 is None or struct1 in freed) \
                and struct2 is not None and dom1.can_pair(dom2)

    reactions = []
    structure = reactant.pair_table

    # We iterate through all the domains
    for (strand_index, strand) in enumerate(structure):
        for (domain_index, domain) in enumerate(strand):

            # The displacing domain must be free
            if (structure[strand_index][domain_index] is not None):
                continue

            displacing_loc = (strand_index, domain_index)

            # search 5'->3' and 3'->5' directions around the loop for a bound
            # domain that is complementary (and therefore can be displaced)

            bound_doms = (find_on_loop(reactant, displacing_loc, -1, filter_3way, 
                                max_helix=max_helix) +
                          find_on_loop(reactant, displacing_loc, +1, filter_3way, 
                                max_helix=max_helix))

            for (displacing, bound, before, after) in bound_doms:
                # Sometimes, the direction of 3-way branch migration matters, 
                # and we don't want to allow invalid intermediate states:
                # x( y( x( y x + ) ) ) -> x y x( y( x( + ) ) )
                try :
                    products = do_3way_migration(reactant, displacing.locs, bound.locs)
                except AssertionError:
                    products = do_3way_migration(reactant, displacing.revlocs, bound.revlocs)


                try :
                    reaction = PepperReaction([reactant], products, 'branch-3way')
                except DSDDuplicationError, e :
                    #assert opening_rate(length) == PepperReaction.dictionary[e.solution].rate
                    reaction = e.existing

                # length of invading domain
                length = len(displacing)

                # calculate reaction constant
                (after, before) = (before, after)
                reaction.rate = branch_3way_remote_rate(length, before, after)

                # skip remote toehold reactions if directed
                if not remote :
                    if not (not after.is_open and after.stems == 1 and after.bases == 0):
                        # print "Rejecting... " + reaction.kernel_string
                        # import pdb; pdb.set_trace()
                        continue

                reactions.append(reaction)

    # Remove any duplicate reactions
    return sorted(list(set(reactions)))

def do_3way_migration(reactant, displacing_locs, bound_locs):
    """
    Each location in displacing_locs will end up bound to the corresponding
    location in bound_locs. The stuff bound to bound_locs will end up un-bound
    """

    def update_structure(struct, displacing_loc, new_bound_loc):
        """
        displacing_loc will be bound to new_bound_loc; whatever new_bound_loc
        was bound to will be unbound.
        """
        displaced_loc = struct[new_bound_loc[0]][new_bound_loc[1]]

        assert struct[displacing_loc[0]][displacing_loc[1]] is None
        assert struct[new_bound_loc[0]][new_bound_loc[1]] is not None
        assert struct[displaced_loc[0]][displaced_loc[1]] is not None
        assert struct[displaced_loc[0]][displaced_loc[1]] == new_bound_loc
    
        struct[displacing_loc[0]][displacing_loc[1]] = new_bound_loc
        struct[new_bound_loc[0]][new_bound_loc[1]] = displacing_loc
        struct[displaced_loc[0]][displaced_loc[1]] = None

    struct = reactant.pair_table
    for displacing_loc, new_bound_loc in zip(displacing_locs, bound_locs):
        update_structure(struct, displacing_loc, new_bound_loc)
    newstr = pair_table_to_dot_bracket(struct)
    try:
        product = PepperComplex(reactant.sequence, newstr)
        product.pair_table = struct
    except DSDDuplicationError, e:
        product = e.existing
  
    return product.split()

def branch_4way(reactant, max_helix = False, remote=True):
    """
    Returns a list of complex sets that can be created through one iteration of
    a 4 way branch migration reaction (each set consists of the molecules that
    result from the iteration; more than one molecule may result because branch
    migration can liberate strands and complexes).
    """

    def filter_4way(dom1, struct1, loc1, dom2, struct2, loc2, freed):
        """ A filter function for *find_on_loop()* """
        # struct1 is necessary within the zipper function...
        return struct1 is not None and struct2 is not None and dom1 == dom2

    reactions = []
    structure = reactant.pair_table

    # We loop through all domains
    for (strand_index, strand) in enumerate(structure):
        for (domain_index, domain) in enumerate(strand):

            # Unbound domains can't participate in branch migration
            if (structure[strand_index][domain_index] is None):
                continue

            displacing_loc = (strand_index, domain_index)

            # searches only 5'->3' direction around the loop for a bound domain that
            # has the same sequence (and therefore can be displaced)
            #
            #   z  _~_  z* (displacing)
            #  ___/   \___>
            #
            # <___     ___
            #     \_ _/
            #   z*  ~   z
            #

            bound_doms = find_on_loop(reactant, displacing_loc, +1, 
                    filter_4way, max_helix=max_helix, b4way = True)

            # build products
            for (displacing, displaced, before, after) in bound_doms:

                products = do_4way_migration(reactant, displacing.locs,
                            (structure[displacing_loc[0]][displacing_loc[1]] 
                                for displacing_loc in displacing.locs),
                            (structure[bound_loc[0]][bound_loc[1]]
                                for bound_loc in displaced.locs), 
                            displaced.locs)

                try :
                    reaction = PepperReaction([reactant], products, 'branch-4way')
                except DSDDuplicationError, e :
                    #assert opening_rate(length) == PepperReaction.dictionary[e.solution].rate
                    reaction = e.existing

                # length of invading domain
                length = len(displacing)

                # calculate reaction constant
                reaction.rate = branch_4way_remote_rate(length, before, after)

                # skip remote toehold reactions
                if not remote:
                    if not (not after.is_open and after.stems == 1 and after.bases == 0 and
                            not before.is_open and before.stems == 1 and before.bases == 0):
                        continue

                reactions.append(reaction)

    # remove any duplicate reactions
    return sorted(list(set(reactions)))

def do_4way_migration(reactant, loc1s, loc2s, loc3s, loc4s):
    """
    Perform a sequence of max_helix 4-way branch migration reactions.
    """
    def update_structure(struct, loc1, loc2, loc3, loc4):
        """
        Performs a 4 way branch migration on a copy of reactant, with loc1 as the
        displacing domain, loc2 as the domain displaced from loc1, loc3 as the
        template domain, and loc4 as the domain displaced from loc3. Returns the
        set of complexes produced by this reaction (may be one or more complexes).
    
        loc1:loc2, loc3:loc4 -> loc1:loc3, loc2:loc4
        """
        struct[loc1[0]][loc1[1]] = loc3
        struct[loc3[0]][loc3[1]] = loc1
        struct[loc2[0]][loc2[1]] = loc4
        struct[loc4[0]][loc4[1]] = loc2

    struct = reactant.pair_table
    for loc1, loc2, loc3, loc4 in zip(loc1s, loc2s, loc3s, loc4s):
        update_structure(struct, loc1, loc2, loc3, loc4)
    newstr = pair_table_to_dot_bracket(struct)
    try:
        product = PepperComplex(reactant.sequence, newstr)
        product.pair_table = struct
    except DSDDuplicationError, e:
        product = e.existing
    return product.split()

def find_on_loop(reactant, start_loc, direction, filter, max_helix=True, b4way = False):
    """
    Finds the next domain within `reactant` that's on the same inner loop as
    `start_loc` and matches the passed `filter` function. Looks in either the
    5'->3' (+1) or 3'->5' (-1) `direction`.

    Filter should accept the following arguments and return True or False:
      -	dom (utils.Domain) : the domain at `loc`
      -	struct (tuple or None): a (strand index, domain index) pair
              indicating what `dom` is bound to, or None if `dom` is unpaired.
      -	loc (tuple) : a (strand index, domain index) pair
      -       Note that while every single-stranded domain is tested,
              only the "first" domain of a stem helix (in the direction of
              search) will be passed to the filter.


    Returns an array of tuples: `(loc, before, after)`, where:
      -	`loc` is a (strand index, domain index) pair indicating the
              position of the matched domain
      -	`before` is a list of (domain, struct, loc) triples giving the
              domains after `start_loc` but before the matched domain on the loop
              (or None instead of triple where there is a break in the loop)
      -	`after` is a list of (domain, struct, loc) triples giving the
              domains after the matched domain but before `start_loc` on the loop
              (or None instead of triple where there is a break in the loop)

    Where a loop involves stems, only one of the complementary domains will be
    listed in the array of tuples, specifically, the "first" one in the search
    direction. Thus, a multiloop with n unpaired domains and m stems will
    result, for closed loops, in `len(before+after) == n+m-2`, as the match
    location and `start_loc` are omitted.

    `before` and `after` are converted to Loop objects (see utils.py) prior
    to being returned, so that the number of bases and number of stems and
    open/closed status is readily accessible.

    Note 1: `before` and `after` refer to the partial loops between `start_loc`
    and each of the results, _in the `direction`_ of the search. For example:

       A
          ____
         /    \ 
     x  |     |  x*
        |
         \____> 3'

            B

    If `start_loc` pointed to `x` and `direction` is +1, then `before` would
    be `A` and `after` would be `B`. If instead `direction` is -1, then
    `before` is `B` and `after` is `A`.

    Note 2: If the domain passed to `start_loc` is a duplex, the results may
    be unexpected:

           ___  x  ___
    5' ___/   \___/   \___
    3' ___  A  ___  B  ___)
          \___/   \___/
                x*

    Notice that the duplex x() participates in two internal loops (A and B).
    By convention, the internal loop considered is the _internal loop which
    encloses this domain_. That means if you pass domain x and +1, you'll get
    loop A, whereas if you pass x and -1, you'll get loop B. This is in an
    attempt to be consistent with the case where you pass an unpaired domain
    (and therefore the internal loop searched is the one which encloses the
    unpaired domain).

    Returns a list of [(4xLoop() Objects), ...] for every reaction matching
    the filter function. 4xLoop = The 


    """
    results = []
    loop = []
    freed = set()

    def triple(loc):
        return (reactant.get_domain(loc), reactant.get_structure(loc), loc)

    # We now follow the external loop from the starting pair
    # searching for a bound domain to displace
    bound_loc = start_loc

    # Avoid getting stuck inside an internal loop enclosed by this domain,
    # if the starting domain is a duplex.
    #
    #   1      2
    #  ___________
    #  ____  _____
    #   1*  /  2*
    #
    #  If we start at domain 1, going in the - direction, then
    #  immediately continue to the next domain, we'll go to 1*
    if reactant.get_structure(bound_loc) is not None:
        bound_loc = reactant.get_structure(bound_loc)

    # Follow the external loop to the end
    while True:
        # move to the next domain in the indicated direction
        # (+1 = 5' -> 3', -1 = 3' -> 5')
        bound_loc = (bound_loc[0], bound_loc[1] + direction)

        # if we've reached the end of the strand (5')
        if (bound_loc[1] == -1):

            # Continue to next strand
            bound_loc = (wrap(bound_loc[0] - 1, reactant.size),)
            bound_loc = (bound_loc[0], reactant.strand_length(bound_loc[0]))
            loop.append(None)  # EW
            continue

        # if we've reached the end of the strand (3')
        elif (bound_loc[1] == reactant.strand_length(bound_loc[0])):

            # Continue to next strand
            bound_loc = (wrap(bound_loc[0] + 1, reactant.size), -1)
            loop.append(None)  # EW
            continue

        if bound_loc == start_loc:
            # We've returned to the original location of the
            # displacing domain
            break

        # try to match the filter function
        elif (filter(
                reactant.get_domain(start_loc),
                reactant.get_structure(start_loc),
                start_loc,
                reactant.get_domain(bound_loc),
                reactant.get_structure(bound_loc),
                bound_loc, freed)):

            # store the id(s) of the displaced strand
            disp_strands = [None, None]
            if reactant.get_structure(start_loc) is not None:
                disp_strands[0] = reactant.get_structure(start_loc)[0]
            if reactant.get_structure(bound_loc) is not None:
                disp_strands[1] = reactant.get_structure(bound_loc)[0]
                freed.add(bound_loc)

            # append the location
            results.append((bound_loc, len(loop), disp_strands))

        # store unpaired domains and "first" domain of each stem
        loop.append(triple(bound_loc))  # EW
        #if filtered and not b4way and reactant.structure[bound_loc[0]][bound_loc[1]] is not None:
        #    loop.append(triple(reactant.structure[bound_loc[0]][bound_loc[1]]))
        
        # if the domain at bound_loc is unbound
        if (reactant.get_structure(bound_loc) is None):
            # look to the next domain
            continue

        # so it's bound to something: follow the structure to stay on the same loop
        bound_loc = reactant.get_structure(bound_loc)

    if max_helix: 

        zipped_results = []
        for (bound_loc, i, disp_strands) in results:

            if b4way :
                # for 4-way we don't zip along before/after, but along top/bottom...
                # both top and bottom have to mach -- i.d. eht filter function has
                # to compare top and bottom.

                # 4-way zipper
                top_loop = []
                top_loc = (start_loc[0], start_loc[1] + direction)
                while True:
                    if (top_loc[1] == -1):
                        top_loop.append(None)
                        break
                    elif (top_loc[1] == reactant.strand_length(top_loc[0])):
                        top_loop.append(None)
                        break
                    else :
                        top_loop.append(triple(top_loc))
                        top_loc = (top_loc[0], top_loc[1] + direction)

                bottom_loop = []
                bottom_loc = (bound_loc[0], bound_loc[1] + direction)
                while True:
                    if bottom_loc[1] == -1 :
                        break
                    elif bottom_loc[1] == reactant.strand_length(bottom_loc[0]) :
                        break
                    else :
                        bottom_loop.append(triple(bottom_loc))
                        bottom_loc = (bottom_loc[0], bottom_loc[1] + direction)

                b4way_loop = top_loop + bottom_loop[::-1] # reversed
            else :
                b4way_loop = None
 
            zipped_results.append(zipper(
                reactant, start_loc, bound_loc, loop[:i], loop[i + 1:], 
                direction, disp_strands, freed, filter, b4way=b4way_loop))
        return zipped_results
    else:
        return [(Loop([triple(start_loc)]), Loop([triple(bound_loc)]), 
            Loop(loop[:i]), Loop(loop[i + 1:])) for (bound_loc, i, _) in results]

def zipper(reactant, start_loc, bound_loc, before, after, direction, disp_strands, freed, 
        filter, b4way=None):
    """
    Takes a result from `find_on_loop` and "zips" it inwards (in the given
    `direction`); that is, given some `start_loc` and some `bound_loc`, tries
    to find as many adjacent domains as possible such that the `filter`
    function still returns True.

    For example, if `start_loc` was b1 and `bound_loc` was b1*, and the filter
    function specified that the domain at `start_loc` must be complementary to
    the domain at `bound_loc`, then the function would return [b1,b2] as
    start_locs and [b1*, b2*] as bound_locs

    Note: if you pass disp_strands = None, you will turn of the new max-helix semantics
    and enable Caseys max helix semantics


            b1* b2*
            ______
         __/      \__>
        <__        __
           \______/
            b1  b2


    return start_locs, bound_locs, before, after

    """
    if freed is None:
        freed = set()

    def triple(loc):
        return (reactant.get_domain(loc), reactant.get_structure(loc), loc)

    def move_towards_middle(middle, direction):

        dstrand, ddomain = start_loc
        bstrand, bdomain = bound_loc

        while True:
            # move domain pointers "inwards" towards each other
            ddomain += direction
            if b4way :
                bdomain += direction
            else :
                bdomain -= direction

            # if ddomain is still on dstrand
            if ((ddomain < reactant.strand_length(dstrand) and ddomain >= 0) and

                    # and bdomain is still on bstrand
                    (bdomain < reactant.strand_length(bstrand) and bdomain >= 0) and

                    # and ddomain hasn't passed bound_loc
                    (cmp((dstrand, ddomain), bound_loc) == cmp(start_loc, bound_loc)) and

                    # and bdomain hasn't passed start_loc
                    (cmp((bstrand, bdomain), start_loc) == cmp(bound_loc, start_loc)) and

                    # If we are displacing, we are still displacing the same strand.
                    # Delete these conditions if you want to enable the previous max-helix semantics
                    ((reactant.get_structure((dstrand, ddomain)) is None) or 
                        (disp_strands is None) or 
                        # disp_strands is none! but it was just freed during this move!
                        ((disp_strands[0] is None) and 
                            reactant.get_structure((dstrand, ddomain)) in freed) or
                        (reactant.get_structure((dstrand, ddomain))[0] == disp_strands[0])) and

                    ((reactant.get_structure((bstrand, bdomain)) is None) or 
                        (disp_strands is None) or 
                        # disp_strands is none! but it was just freed during this move!
                        ((disp_strands[1] is None) and 
                            reactant.get_structure((bstrand, bdomain)) in freed) or
                        (reactant.get_structure((bstrand, bdomain))[0] == disp_strands[1])) and
                    # ~~~

                    # and filter condition still applies
                    filter(reactant.get_domain((dstrand, ddomain)),
                        reactant.get_structure((dstrand, ddomain)),
                        (dstrand, ddomain),
                        reactant.get_domain((bstrand, bdomain)),
                        reactant.get_structure((bstrand, bdomain)),
                        (bstrand, bdomain), freed)):

                # add new positions to list
                if direction == 1:
                    start_locs.append(triple((dstrand, ddomain)))
                    bound_locs.append(triple((bstrand, bdomain)))
                elif direction == -1:
                    start_locs[:0] = [triple((dstrand, ddomain))]
                    bound_locs[:0] = [triple((bstrand, bdomain))]

                # remove zipped positions from `middle` loop
                displacing_index = (direction - 1) / 2
                bound_index = (-direction - 1) / 2

                freed.add((bstrand, bdomain))

                if middle and middle[displacing_index] is not None and \
                        middle[displacing_index][2] == (dstrand, ddomain):
                    del middle[displacing_index]

                if middle and middle[bound_index] is not None and \
                        middle[bound_index][2] == (bstrand, bdomain):
                    del middle[bound_index]

            else:
                break

    start_locs = [triple(start_loc)]
    bound_locs = [triple(bound_loc)]

    if b4way :
        for (d, middle) in [(direction, b4way)]:
            move_towards_middle(middle, d)
    else :
        for (d, middle) in [(direction, before), (-direction, after)]:
            move_towards_middle(middle, d)

    start_locs = Loop(start_locs)
    bound_locs = Loop(bound_locs)
    before = Loop(before)
    after = Loop(after)

    return start_locs, bound_locs, before, after

def breathing(reactant, ends = 1):
    # every complex reacts by opening one dummy nucleotide at the end
    # add an universal dummy-domain to each helix-end.
    pass

