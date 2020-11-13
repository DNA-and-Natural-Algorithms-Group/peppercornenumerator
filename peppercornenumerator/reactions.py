#
#  peppercornenumerator/reactions.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

from dsdobjects.complex_utils import (make_pair_table,
                                      pair_table_to_dot_bracket, 
                                      make_strand_table,
                                      strand_table_to_sequence, 
                                      make_loop_index)
from .utils import wrap
from .objects import (SingletonError,
                      PepperComplex, 
                      PepperReaction,
                      Loop)
from .ratemodel import (unimolecular_binding_rate,
                        bimolecular_binding_rate,
                        opening_rate,
                        branch_3way_remote_rate,
                        branch_4way_remote_rate)

def bind11(reactant, max_helix = True):
    """
    Returns a list of reaction pathways which can be produced by 1-1 binding
    reactions of the argument complex. The 1-1 binding reaction is the
    hybridization of two complementary unpaired domains within a single complex
    to produce a single unpseudoknotted product complex.
    """
    reactions = set()
    structure = list(reactant.pair_table)
    for (strand_index, strand) in enumerate(structure):
        for (domain_index, domain) in enumerate(strand):
            # The displacing domain must be free
            if structure[strand_index][domain_index] is not None :
                continue
            start_loc = (strand_index, domain_index)
            # search (one direction) around the loop for an open domain that can be bound.
            results = find_on_loop(reactant, start_loc, filter_bind11)
            assert len(results) == len(find_on_loop(reactant, start_loc, filter_bind11, direction = -1))
            for e, (invader, before, target, after) in enumerate(results):
                if max_helix:
                    invader, before, target, after = zipper(
                            reactant, invader[0], before, target[0], after, filter_bind11)
                results[e] = list(map(Loop, [invader, before, target, after]))
            # build products
            for (loc1s, before, loc2s, after) in results:
                # Should be reversed loc2s right?
                assert [x == ~y for x,y in zip(loc1s.domains, loc2s.domains)]
                product = do_bind11(reactant, loc1s.domain_locs, loc2s.domain_locs)
                reaction = PepperReaction([reactant], [product], 'bind11')
                if reaction.rate_constant[0] is None:
                    reaction.rate_constant = (unimolecular_binding_rate(loc1s.dlength, before, after), '/s')
                reactions.add(reaction)
    return sorted(reactions)

def do_bind11(reactant, loc1s, loc2s):
    """ Returns PepperComplex after the bind11 reaction. """
    news = list(reactant.pair_table)
    for loc1, loc2 in zip(loc1s, loc2s):
        assert news[loc1[0]][loc1[1]] is None
        assert news[loc2[0]][loc2[1]] is None
        news[loc1[0]][loc1[1]] = loc2
        news[loc2[0]][loc2[1]] = loc1
    newstr = pair_table_to_dot_bracket(news)
    try:
        new = PepperComplex(list(reactant.sequence), newstr)
    except SingletonError as err:
        new = err.existing
    return new

def bind21(reactant1, reactant2, max_helix = True, pkwarning = False):
    """
    Returns a list of reaction pathways which can be produced by 2-1 binding
    reactions of the argument complexes. The 2-1 binding reaction is the
    hybridization of two complementary unpaired domains, each in a different
    complex, to produce a single, unpseudoknotted product complex containing
    all of the strands contained in either of the original complexes.

    """
    r1_doms = reactant1.available_domains
    r2_doms = reactant2.available_domains
    reactions = []
    # Iterate through all the free domains in reactant1
    for (dom1, s1, d1) in r1_doms:
        # For each, find any domains in reactant2 that could bind
        for (dom2, s2, d2) in r2_doms:
            # If it can pair, this is one possible reaction (this kind of
            # reaction cannot possibly produce a pseudoknotted structure)
            assert (dom1 is ~dom2) is (dom2 is ~dom1)
            if dom1 is ~dom2:
                # combine the two complexes into one, but do not perform the association
                reactions.append(join_complexes_21(
                    reactant1, (s1, d1),
                    reactant2, (s2, d2)))

    if pkwarning: # check if potential pseudoknots are there:
        pk1_doms = reactant1.pk_domains
        pk2_doms = reactant2.pk_domains
        # Iterate through all the free domains in reactant1
        for (dom1, strand_num1, dom_num1) in r1_doms + pk1_doms:
            # For each, find any domains in reactant2 that could bind
            for (dom2, strand_num2, dom_num2) in r2_doms + pk2_doms:
                if (dom1, strand_num1, dom_num1) in r1_doms and \
                    (dom2, strand_num2, dom_num2) in r2_doms:
                    # Exclude the non-pseudoknotted interactions
                    continue
                if dom1 is ~dom2:
                    log.warning("potential pk-interaction: {} and {}".format(reactant1, reactant2))
   
    output = set()
    for cplx, loc1, loc2 in reactions:
        def findloc(trip1, trip2):
            (dom1, _, dloc1) = trip1 
            (dom2, _, dloc2) = trip2
            return dloc1 == loc1 and dloc2 == loc2
        # build "before" and "after" loop structures via find_on_loop ...
        [(loc1s, before, loc2s, after)] = find_on_loop(cplx, loc1, findloc)
        if max_helix:
            loc1s, before, loc2s, after = zipper(cplx, loc1s[0], before, 
                                                 loc2s[0], after, filter_bind11)
        [loc1s, before, loc2s, after] = list(map(Loop, [loc1s, before, loc2s, after]))
        product = do_bind11(cplx, loc1s.domain_locs, loc2s.domain_locs)
        reaction = PepperReaction([reactant1, reactant2], [product], 'bind21')
        if reaction.rate_constant[0] is None:
            assert [x == ~y for x,y in zip(loc1s.domains, loc2s.domains)]
            reaction.rate_constant = (bimolecular_binding_rate(loc1s.dlength), '/M/s')
        output.add(reaction)
    return sorted(output)

def join_complexes_21(complex1, loc1, complex2, loc2):
    """ Joins two complexes into one disconnected, base-pair-compatible complex.

    Returns the disconnected complex and the loci for pairing.
    """
    # a value larger than any possible location on the pair_table.
    maxlen = complex1.size + complex2.size + 1

    for e, (st, pt) in enumerate(complex1.rotate_pt()):
        l1 = complex1.rotate_pairtable_loc(loc1, e)
        li, ext = make_loop_index(pt)
        if li[l1[0]][l1[1]] == 0:
            seq1 = strand_table_to_sequence(st)
            pt[l1[0]][l1[1]] = (maxlen, maxlen) # add an additional '(' 
            ptb1 = pair_table_to_dot_bracket(pt)
            break
    for e, (st, pt) in enumerate(complex2.rotate_pt()):
        l2 = complex2.rotate_pairtable_loc(loc2, e)
        li, ext = make_loop_index(pt)
        if li[l2[0]][l2[1]] == 0:
            seq2 = strand_table_to_sequence(st)
            pt[l2[0]][l2[1]] = (-1, -1) # add an additional ')' 
            ptb2 = pair_table_to_dot_bracket(pt)
            break

    # build the new sequence and structure *including* the new pair
    newseq = seq1 + ['+'] + seq2
    newstr = ptb1 + ['+'] + ptb2
    # update l2 from the new structure
    combined = make_pair_table(newstr)
    l2 = combined[l1[0]][l1[1]]
    # remove the new pair again ...
    combined[l1[0]][l1[1]] = None
    combined[l2[0]][l2[1]] = None
    # update the structure to the unpaired (disconnected) version.
    newstr = pair_table_to_dot_bracket(combined)
    try:
        new_complex = PepperComplex(newseq, newstr)
    except SingletonError as err:
        new_complex = err.existing
    # strands may be rotated in the new complex ...
    for e, (st, pt) in enumerate(new_complex.rotate()):
        if st == newseq and pt == newstr:
            rotate = e
            break
    else:
        raise ValueError(f'Joining of complexes {complex1} and {complex2} failed.')
    loc1 = new_complex.rotate_pairtable_loc(l1, -rotate)
    loc2 = new_complex.rotate_pairtable_loc(l2, -rotate)
    if loc1 > loc2:
        (loc1, loc2) = (loc2, loc1)
    return new_complex, loc1, loc2

def open1N(reactant, max_helix = True, release_11 = 6, release_1N = 6, dG_bp = -1.7):
    """ Returns a list of open reactions.

    Args:
        reactant (PepperComplex): The reactant complex
        max_helix (bool, optional): Use max-helix notion. Defaults to True.
        release_11 (int, optional): Threshold length for a open11 reaction. Defaults to 6.
        release_1N (int, optional): Threshold length for a open12 reaction. Defaults to 6.

    Returns:
        [PepperReactions]
    """

    def get_max_helix(loc, structure):
        # A: Strand/domain position on "top" strand - CG 5/21
        helix_startA = list(loc)
        helix_length = len(reactant.get_domain(loc))
        # B: Strand/domain position on "bottom" strand - CG 5/21
        helix_startB = list(structure[loc[0]][loc[1]])
        # If the domain is bound to an earlier domain, then we have
        # already considered it, so skip it
        if (helix_startB < helix_startA):
            return None, None, None
        helix_endA = helix_startA[:]
        helix_endB = helix_startB[:]
        # Now iterate through the whole helix to find the other end of
        # this one (The helix ends at the first strand break from
        # either direction)
        ext_fw, ext_bw = False, False
        while True:
            # Strands run in opposite directions, so A must be incremented
            # and B decremented in order that both pointers move "right"
            # along the helix - CG 5/21
            helix_endA[1] += 1
            helix_endB[1] -= 1
            # If one of the strands has broken, the helix has ended
            if helix_endA[1] >= reactant.strand_length(helix_endA[0]):
                break
            elif helix_endB[1] < 0:
                break
            # If these domains aren't bound to each other, the helix has ended
            if tuple(helix_endA) != structure[helix_endB[0]][helix_endB[1]]:
                break
            # Add the current domain to the current helix
            temp_bl = (helix_endA[0], helix_endA[1])
            helix_length += len(reactant.get_domain(temp_bl))
            ext_fw = True
        helix_endA[1] -= 1
        helix_endB[1] += 1
        # We must also iterate in the other direction
        while True:
            helix_startA[1] -= 1
            helix_startB[1] += 1
            # If one of the strands has broken, the helix has ended
            if helix_startA[1] < 0:
                break
            elif helix_startB[1] >= reactant.strand_length(helix_startB[0]):
                break
            # If these domains aren't bound to each other, the helix has ended
            if tuple(helix_startA) != structure[helix_startB[0]][helix_startB[1]]:
                break
            # Add the current domain to the current helix
            temp_bl = (helix_startA[0], helix_startA[1])
            helix_length += len(reactant.get_domain(temp_bl))
            ext_bw = True
        # Move start location back to the first domain in the helix
        helix_startA[1] += 1
        helix_startB[1] -= 1
        if ext_fw and ext_bw :
            # we only want to allow moves that start at helix ends!
            return None, None, None
        return helix_startA, helix_endA, helix_length

    # remember the larger release cutoff to avoid reactions opening longer helices.
    max_release = max(release_11, release_1N) if release_11 and release_1N else 0

    reactions = []
    structure = list(reactant.pair_table)
    if not max_helix:
        # Iterate through all the domains
        for (strand_index, strand) in enumerate(structure):
            for (domain_index, domain) in enumerate(strand):
                # The bound domain must be... bound
                if reactant.get_paired_loc((strand_index, domain_index)) is None:
                    continue
                loc = (strand_index, domain_index)
                domain = reactant.get_domain(loc)
                products = do_open(reactant, loc, loc)
                reactions.append((products, len(domain)))
    else: # max-helix mode:
        for (strand_index, strand) in enumerate(structure):
            for (domain_index, domain) in enumerate(strand):
                # If the domain is unpaired, skip it
                if reactant.get_paired_loc((strand_index, domain_index)) is None:
                    continue
                loc = (strand_index, domain_index)
                domain = reactant.get_domain(loc)
                length = len(domain)
                if max_release and length > max_release:
                    continue
                loc1, loc2, length = get_max_helix(loc, structure)
                if loc1 is None:
                    continue
                if max_release and length > max_release:
                    continue
                # If the helix is short enough, we have a reaction
                products = do_open(reactant, loc1, loc2)
                reactions.append((products, length))
    output = set()
    for products, length in reactions:
        if release_11 and len(products) == 1 and length > release_11:
            continue
        elif release_1N and len(products) > 1 and length > release_1N:
            continue
        reaction = PepperReaction([reactant], products, 'open')
        if reaction.rate_constant[0] is None:
            reaction.rate_constant = (opening_rate(length, dG_bp = dG_bp), '/s')
        output.add(reaction)
    return sorted(output)

def do_open(reactant, sloc, eloc):
    assert sloc[0] == eloc[0]
    news = list(reactant.pair_table)
    for dom in range(sloc[1], eloc[1]+1):
        loc = news[sloc[0]][dom]
        news[sloc[0]][dom] = None
        news[loc[0]][loc[1]] = None
    newstr = pair_table_to_dot_bracket(news)
    try:
        new = PepperComplex(list(reactant.sequence), newstr)
    except SingletonError as err:
        new = err.existing
    return list(new.split())

def branch_3way(reactant, max_helix = True, remote = True):
    """
    Returns a list of reaction pathways that can be created through one
    iteration of a 3 way branch migration reaction (more than one molecule may
    be produced by a reaction because branch migration can liberate strands and
    complexes).
    """

    reactions = set()
    structure = list(reactant.pair_table)
    for (strand_index, strand) in enumerate(structure):
        for (domain_index, domain) in enumerate(strand):
            # The displaced domain must not be free
            if (structure[strand_index][domain_index] is not None):
                continue
            # search 5'->3' and 3'->5' directions around the loop for a bound
            # domain that is complementary (and therefore can be displaced)
            start_loc = (strand_index, domain_index)

            # build products
            fwresults = find_on_loop(reactant, start_loc, filter_3way, direction = 1)
            bwresults = find_on_loop(reactant, start_loc, filter_3way, direction = -1)

            results = []
            for (invader, before, target, after) in fwresults:
                if max_helix:
                    invader, before, target, after = zipper(
                            reactant, invader[0], before, target[0], after, filter_3way)
                after += invader
                results.append(list(map(Loop, [invader[::-1], before, target[::-1], after])))

            for (invader, before, target, after) in bwresults:
                if max_helix:
                    invader, before, target, after = zipper(
                            reactant, invader[0], before, target[0], after, filter_3way)
                after += invader[::-1]
                results.append(list(map(Loop, [invader, before, target, after])))

            for (displacing, before, bound, after) in results:
                if remote is False:
                    if not (not before.is_open and before.stems == 1 and before.bases == 0):
                        continue
                products = do_3way_migration(reactant, list(displacing.domain_locs), list(bound.domain_locs))
                reaction = PepperReaction([reactant], products, 'branch-3way')
                if reaction.rate_constant[0] is None:
                    reaction.rate_constant = (branch_3way_remote_rate(displacing.dlength, before, after), '/s')
                reactions.add(reaction)
    return sorted(reactions)

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

    struct = list(reactant.pair_table)
    for displacing_loc, new_bound_loc in zip(displacing_locs, bound_locs):
        update_structure(struct, displacing_loc, new_bound_loc)
    newstr = pair_table_to_dot_bracket(struct)
    try:
        product = PepperComplex(list(reactant.sequence), newstr)
    except SingletonError as err:
        product = err.existing
    return list(product.split())

def branch_4way(reactant, max_helix = False, remote=True):
    """
    Returns a list of complex sets that can be created through one iteration of
    a 4 way branch migration reaction (each set consists of the molecules that
    result from the iteration; more than one molecule may result because branch
    migration can liberate strands and complexes).
    """

    reactions = set()
    structure = list(reactant.pair_table)
    for (strand_index, strand) in enumerate(structure):
        for (domain_index, domain) in enumerate(strand):
            # Unbound domains can't participate in 4-way branch migration
            if (structure[strand_index][domain_index] is None):
                continue
            start_loc = (strand_index, domain_index)
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
            results = find_on_loop(reactant, start_loc, filter_4way)
            for e, (invader, before, target, after) in enumerate(results):
                if max_helix:
                    invader, _, target, _ = zipper(
                            reactant, invader[0], None, target[0], None, filter_4way)
                results[e] = list(map(Loop, [invader, before, target, after]))
            for (displacing, before, displaced, after) in results:
                if remote is False:
                    if not ((not after.is_open and after.stems == 1 and after.bases == 0) or
                            (not before.is_open and before.stems == 1 and before.bases == 0)):
                        continue
                products = do_4way_migration(reactant, displacing.domain_locs,
                            (structure[dl[0]][dl[1]] for dl in displacing.domain_locs),
                            (structure[bl[0]][bl[1]] for bl in displaced.domain_locs), 
                            displaced.domain_locs)

                reaction = PepperReaction([reactant], products, 'branch-4way')
                if reaction.rate_constant[0] is None:
                    reaction.rate_constant = (branch_4way_remote_rate(displacing.dlength, before, after), '/s')
                reactions.add(reaction)
    return sorted(reactions)

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
        assert None not in (loc1, loc2, loc3, loc4)
        assert struct[loc1[0]][loc1[1]] == loc2
        assert struct[loc3[0]][loc3[1]] == loc4
        assert struct[loc2[0]][loc2[1]] == loc1
        assert struct[loc4[0]][loc4[1]] == loc3

        struct[loc1[0]][loc1[1]] = loc3
        struct[loc3[0]][loc3[1]] = loc1
        struct[loc2[0]][loc2[1]] = loc4
        struct[loc4[0]][loc4[1]] = loc2

    struct = list(reactant.pair_table)
    for loc1, loc2, loc3, loc4 in zip(loc1s, loc2s, loc3s, loc4s):
        update_structure(struct, loc1, loc2, loc3, loc4)
    newstr = pair_table_to_dot_bracket(struct)
    try:
        product = PepperComplex(list(reactant.sequence), newstr)
    except SingletonError as err:
        product = err.existing
    return list(product.split())

# Filter functions for find_on_loop()
def filter_bind11(trip1, trip2):
    (dom1, struct1, loc1) = trip1
    (dom2, struct2, loc2) = trip2
    return (struct1 is None) and (struct2 is None) and (dom2 is ~dom1)

def filter_3way(trip1, trip2):
    (dom1, struct1, loc1) = trip1 
    (dom2, struct2, loc2) = trip2 
    return (struct1 is None) and (struct2 is not None) and (dom2 is ~dom1)

def filter_4way(trip1, trip2):
    (dom1, struct1, loc1) = trip1
    (dom2, struct2, loc2) = trip2
    return (struct1 is not None) and (struct2 is not None) and (dom1 is dom2)

def find_on_loop(reactant, start_loc, pattern, direction = 1):
    r""" Find a reaction pattern within a loop.
    
    Starts at a particular locus and searches for every possible pattern that
    preseves secondary structure, i.e. matches within the "loop".  For a loop
    that involves stems, only one of the complementary domains will be listed
    in the array of tuples, specifically, the "first" one in the search
    direction.

    Note 1: `before` and `after` refer to the partial loops between `start_loc`
    and each of the results, _in the `direction`_ of the search. For example:
        A ____
         /    \ 
     x  |     |  x*
        |
         \____> 3'
        B
    If `start_loc` pointed to `x` and `direction` is +1, then `before` would be
    `A` and `after` would be `B`. If instead `direction` is -1, then `before`
    is `B` and `after` is `A`.  

    Note 2: If the domain passed to `start_loc` is a duplex, the results may be
    unexpected:
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
  
    Args:
        reactant (:obj:) = The reactant complex object.
        start_loc ((int, int)) = A tuple pointing to the respective strand and
            domain indices for reactant.pair_table
        pattern (function) = A pattern matching function which determines if
            there exists a valid result when given the start locus and a target locus. 
        directions (-1, +1, optional): Searching for binding partners in 
            5'->3' (+1) or 3'->5' (-1) direction.

    Returns: 
       Four loops: start_loc, after, target_loc, before
    """
    results = []
    loop = []
    def triple(loc):
        return (reactant.get_domain(loc), reactant.get_paired_loc(loc), loc)
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
    if reactant.get_paired_loc(bound_loc) is not None:
        bound_loc = reactant.get_paired_loc(bound_loc)
    
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
            loop.append(None)
            continue

        # if we've reached the end of the strand (3')
        if (bound_loc[1] == reactant.strand_length(bound_loc[0])):
            # Continue to next strand
            bound_loc = (wrap(bound_loc[0] + 1, reactant.size), -1)
            loop.append(None)
            continue

        if bound_loc == start_loc:
            # We've returned to the original location of the displacing domain
            break

        # try to match the pattern function
        if pattern(triple(start_loc), triple(bound_loc)):
            # append the location
            results.append((bound_loc, len(loop)))

        # store unpaired domains and "first" domain of each stem
        loop.append(triple(bound_loc))
 
        # if the domain at bound_loc is unbound
        if reactant.get_paired_loc(bound_loc) is None:
            # look to the next domain
            continue
        else :
            # follow the structure to stay on the same loop
            bound_loc = reactant.get_paired_loc(bound_loc)

    return [([triple(start_loc)],   # invading domain
             loop[:i],              # first linker 
             [triple(bound_loc)],   # target domain
             loop[i+1:]) for (bound_loc, i) in results]

def zipper(reactant, start_trp, before, bound_trp, after, pattern):
    r""" Max-helix mode zipping to extend a given move type.

    Takes a result from `find_on_loop` and finds as many adjacent domains as
    possible such that the `pattern` function still returns True.

    In the following example, if `start_loc` was b1 and `bound_loc` was b1*,
    and the pattern function specified that the domain at `start_loc` must be
    complementary to the domain at `bound_loc`, then the function would return
    [b1, b2] as start_locs and [b1*, b2*] as bound_locs:

            b1* b2*
            ______
         __/      \__>
        <__        __
           \______/
            b1  b2

    Note that 4-way branch migration max-helix never extends into (changes) the 
    before / after loops, so it is safe to pass None instead.

    Returns:
        start_locs (extended list of triples)
        before (reduced list of triples)
        bound_locs (extended list of triples)
        after (reduced list of triples)
    """
    assert pattern(start_trp, bound_trp)

    def triple(loc):
        return (reactant.get_domain(loc), reactant.get_paired_loc(loc), loc)

    def extend_match(start_loc, bound_loc, sdir, bdir):
        # Move upstream or downstream from start locus and bound locus and
        # check if the (global) pattern condition still applies.
        (sstrand, sdomain) = start_loc
        (bstrand, bdomain) = bound_loc
        start_pair = reactant.get_paired_loc(start_loc)
        bound_pair = reactant.get_paired_loc(bound_loc)
        sloc = (sstrand, sdomain + sdir)
        bloc = (bstrand, bdomain + bdir)
        try:
            spair = reactant.get_paired_loc(sloc)
            bpair = reactant.get_paired_loc(bloc)
        except IndexError:
            # sdomain is no longer on sstrand or 
            # bdomain is no longer on bstrand
            return
        if (start_pair is None) != (spair is None): 
            return
        if (bound_pair is None) != (bpair is None): 
            return
        if sloc == start_pair or bloc == bound_pair: 
            return
            
        if (# sdomain hasn't passed bound_loc
            #(cmp(sloc, bound_loc) == cmp(start_loc, bound_loc))
            ((sloc > bound_loc) - (sloc < bound_loc) == \
            (start_loc > bound_loc) - (start_loc < bound_loc)) and
            # and bdomain hasn't passed start_loc
            #cmp(bloc, start_loc) == cmp(bound_loc, start_loc)
            ((bloc > start_loc) - (bloc < start_loc) == \
            (bound_loc > start_loc) - (bound_loc < start_loc)) and
            # if paired, then also the paired domain must be adjacent
            ((start_pair is None) and (spair is None) or
                start_pair[0] == spair[0] and         # same strand!
                start_pair[1] == spair[1] + sdir) and # adjacent!
            ((bound_pair is None) and (bpair is None) or
                bound_pair[0] == bpair[0] and         # same strand!
                bound_pair[1] == bpair[1] + bdir) and # adjacent!
            # and, finally, the pattern condition still applies.
            pattern(triple(sloc), triple(bloc))):

            # add new positions to list
            if sdir == 1:
                start_locs.append(triple(sloc))
                if before and before[0] == triple(sloc):
                    before.pop(0)
            else:
                start_locs[:0] = [triple(sloc)]
                if after and after[-1] == triple(sloc):
                    after.pop(-1)
            if bdir == 1:
                bound_locs[:0] = [triple(bloc)]
                if after and after[0] == triple(bloc):
                    after.pop(0)
            else:
                bound_locs.append(triple(bloc))
                if before and before[-1] == triple(bloc):
                    before.pop(-1)
            extend_match(sloc, bloc, sdir, bdir)
        return

    start_locs = [start_trp]
    bound_locs = [bound_trp]
    if filter_4way(start_trp, bound_trp):
        extend_match(start_trp[2], bound_trp[2], sdir = 1, bdir = 1)
        bound_locs = bound_locs[::-1]
    else :
        extend_match(start_trp[2], bound_trp[2], sdir = 1, bdir = -1)
        extend_match(start_trp[2], bound_trp[2], sdir = -1, bdir = 1)
    return start_locs, before, bound_locs, after

