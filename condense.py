import itertools
import operator
import logging
from reactions import get_auto_name,ReactionPathway
from utils import Complex, Domain, Strand, RestingState

class ReachableRestingStates(object):
    """
    Corresponds to the set of multisets that can be reached from a given 
    complex. Thin wrapper around the frozenset, mostly for prettier printing.
    """
    
    def __init__(self,states):
        # self._states = frozenset(states) 
        self._states = frozenset([ tuple(sorted(s)) for s in states ])
        self._hash = None
        
    @property
    def states(self):
        return self._states
    
    def is_singleton(self):
        x = list(self._states)[0]
        return (len(self._states) == 1) and (len(x) == 1) 
    
    def get_singleton(self):
        x = list(self._states)[0]
        return x[0]
    
    def __iter__(self):
        return iter(self.states)
    
    def __hash__(self):
        if(self._hash is None):
            self._hash = hash(self._states)
        return self._hash
    
    def __eq__(self,other):
        return sorted(list(self.states)) == sorted(list(other.states))    
    
    def __str__(self):
        return "{ " + ", ".join(["{" + ", ".join(map(str,inner)) + "}" for inner in self.states]) + " }"
    
    def __len__(self):
        return len(self._states)
    
    def __repr__(self):
        return "ReachableRestingStates([ " + ", ".join(["(" + ", ".join(map(repr,inner)) + ")" for inner in self.states]) + " ])"
        #return 'ReachableRestingStates('+repr(self.states)+")"

def is_fast(reaction):
    """
    Current heuristic to determine if reaction is fast: unimolecular in reactants
    """
    return reaction.arity == (1,1) or reaction.arity == (1,2)

def is_outgoing(reaction,SCC_set):
    """
    Determines if the passed ReactionPathway leads out of the passed SCC
    """
    if SCC_set is not set:
        SCC_set = set(SCC_set)

    return not (set(reaction.products) <= SCC_set)

def tuple_sum(iterable):
    """
    Given an iterable of iterables, concatenates each element of the 
    internal iterables into a tuple. 

    e.g.: [ [a,b,c], [d], [e,f] ] -> 
    """
    return reduce(operator.concat, iterable, tuple())

def cartesian_multisets(reactions):
    """
    Accepts a list of lists of iterables 
    """
    cartesian = []
    for prods in reactions:
        x = list(itertools.product(*prods))
        y = map(tuple_sum,x)
        z = sorted(y)
        cartesian.append(z)
    return cartesian

def get_reactions_consuming(complexes,reactions):
    """
    Returns a dict mapping each passed complex C to a list of reactions for which 
    C is a reactant.  
    """
    return dict((c,[r for r in reactions if (c in r.reactants)]) for c in complexes)

def get_reactions_producing(complexes,reactions):
    """
    Returns a dict mapping each passed complex C to a list of reactions for which 
    C is a product.  
    """
    return dict((c,[r for r in reactions if (c in r.products)]) for c in complexes)



def tarjans(complexes,reactions,reactions_consuming):
    """
    Implementation of tarjans algorithm which segments SCCs of fast reactions
    """
    
    def strongconnect(complex,index):
        """
        Traverse a complex and all outgoing reactions
        """
        complex._index = index
        complex._lowlink = index
        # print "<"+str(complex)+" i='"+str(complex._index)+"'>"
        index += 1
        S.append(complex)
        
        for reaction in reactions_consuming[complex]:
            if reaction.arity == (1,1):
            # if is_fast(reaction):
                # print '<rxn>' + repr(reaction) + '</rxn>'
                for product in reaction.products:
                    
                    # Product hasn't been traversed; recurse
                    if(product._index is None):
                        index = strongconnect(product,index)
                        complex._lowlink = min(product._lowlink,complex._lowlink)
                    
                    # Product is in the current neighborhood
                    elif product in S:
                        complex._lowlink = min(product._index,complex._lowlink)
                # print "<ll complex='"+str(complex)+"'>"+str(complex._lowlink)+"</ll>"
        
        if(complex._lowlink == complex._index):
            scc = []
            while(True):
                next = S.pop()
                scc.append(next)
                if(next == complex): break
        
            SCCs.append(scc)
            # print "<scc>" + str(scc) + "</scc>"

        # print "<index i='"+str(index)+"' />"
        # print "</"+str(complex)+">"

        return index
    
    index = 0
    S = []
    SCCs = []
    
    # Just to make sure we don't try to access an undeclared property, initialize 
    # all indicies to None
    for complex in complexes:
        complex._index = None
    
    # Outer loop ensures that no complexes get left behind by the DFS
    for complex in complexes:
        if(complex._index is None):
            index = strongconnect(complex,index)
    
    return SCCs


def condense_graph(enumerator):

    # Approach: compute SCCs using Tarjan's algorithm, including only fast 
    # reactions as edges. Compute R(Xn)---the set of destinies for each 
    # complex Xn---for each complex in each SCC. Once R(Xn) is computed for 
    # each Xn, compute the condensed reactions from R(Xn).
    # 
    # 
    # let R(Xn) = the set of multisets of resting states reachable from 
    #   detailed species Xn via fast reactions. That is, R(Xn) is the set of 
    #   possible resting-state destinies for Xn (via fast reactions); each of 
    #   these "destinies" is a multiset, because Xn might e.g. break apart 
    #   into multiple complexes, each of which is a resting state.  
    #
    # let S(Xn) = the resting state to which Xn belongs, or the value 
    #   `undefined` if Xn is transient. 
    # 
    # 
    # To compute R(Xn):
    # 
    #   for each SCC:
    #       let outgoing_reactions = [r for r in reactions if (SCC < r.products and r is fast)]  
    #         (SCC < r.products means SCC is a proper subset of _and not equal to_ r.products)
    #       
    #          If 0 (fast) outgoing_reactions: for each complex Xn in the SCC, 
    #             R(Xn) = { { SCC } }
    #          
    #          If 1+ (fast) outgoing_reactions: for each complex Xn in the SCC,
    #             let outgoing_1_1_products = [r.products for r in outgoing_reactions if r.arity == (1,1)]
    #             let outgoing_1_n_products = [r.products for r in outgoing_reactions if r.arity[1] > 1]
    #             R(Xn) = union( [R(A) for A in outgoing_1_1_products], 
    #                           [cartesian product of [R(A), R(B), R(C) ... R(N)] for [A, B, C ... N] in outgoing_1_n_products]
    #   
    # Compute R(Xn) for Xn in complexes by iterating through each SCC, 
    # classifying it as above, and recursing when necessary to compute R(Xn) 
    # for inner values.
    #
    #
    # Once R(Xn) is computed for all Xn in complexes, compute condensed 
    # reactions as follows:
    #
    # For each reaction r, generate all possible combinations of reachable
    # resting state multisets (elements of R(Xn)) for each of the reactants Xn,
    # then do the same for  the products. For the reactants, R(Xn) is a
    # singleton (because slow reactions only happen between resting state
    # complexes). For the products, this means choose each possible combination
    # of destinies for each of the products.   For each combination of reactants
    # and each combination of products, generate a new reaction rc.   Prune
    # duplicate and trivial reactions from condensed reactions.
    
    # This stores mappings between X1 and R(X1) for each X1 in complexes
    resting_state_targets = {}
    
    # This stores mappings between an SCC and the associated resting state
    resting_states = {}
    
    # This stores mappings between each complex and the SCC which contains it
    SCC_containing = {}
    
    # This remembers which SCCs we've processed
    processed_SCCs = set([])
    
    def get_reachable_resting_states(complex):
        """
        Returns the resting state associated with the passed `complex`, or 
        computes the resting state if not already calculated
        """
        if(complex not in resting_state_targets):
            process_scc(SCC_containing[complex])
        return resting_state_targets[complex]
    
    def get_resting_state(complex):
        scc = SCC_containing[complex]
        scc_set = frozenset(scc)
        if(scc_set not in resting_states):
            process_scc(scc)
        return resting_states[scc_set]
    
    def process_scc(scc):
        """
        Processes a single SCC neighborhood, generating resting state multisets 
        for each complex, and storing the mappings in the outer-scope 
        `resting_state_targets` dict 
        """
        # print "Begin process SCC:" + str(scc)
        
        # Convert to a set for fast lookup
        scc_set = frozenset(scc)
        
        # Check that we haven't been here already
        if(scc_set in processed_SCCs): 
            return
        
        outgoing_reactions = [r for c in scc for r in reactions_consuming[c] if (is_fast(r) and is_outgoing(r,scc_set)) ]
        # print "\tOutgoing reactions: " + str(outgoing_reactions)
                
        # Remember that we've processed this neighborhood
        processed_SCCs.add(scc_set)

        
        # If this SCC is a resting state:  
        if(len(outgoing_reactions) == 0):
            
            resting_state = RestingState(get_auto_name(),scc)
            resting_states[scc_set] = resting_state
            
            for c in scc:
                if c in resting_state_targets: 
                    raise Exception()
                
                resting_state_targets[c] = frozenset([ (resting_state,) ]) #frozenset([ (get_resting_state(c),) ])
                
        # Otherwise, if there are outgoing fast reactions:
        else:
            # Compute the various combinations of resting states (the so-called 
            # 'reachable resting state multisets') which can be reached from this SCC:
            
            # NOTE: Sets are represented as frozensets, and multisets are represented as tuples, 
            # so { } notation is used to represent frozensets and tuple notation is used to 
            # represent multisets. Therefore this: { {A, B, B}, {C, D}, {E, E, F}, {G} }
            # will be represented as { (A, B, B), (C, D), (E, E, F), (G,) }
            
            
            # First, compute the reachable resting state multisets due to outgoing 1-1 reactions: 
            # all these reactions should have 1 product
            outgoing_1_1 = frozenset().union(*[get_reachable_resting_states(A.products[0]) for A in outgoing_reactions if A.arity == (1,1)])
            
            
            # outgoing_1_n is a list (reactions) of lists (products) of frozensets (resting states) of tuples (multisets)  
            # That is, for each 1-n reaction where n>1, outgoing_1_n contains a list of the 
            # reachable resting state multisets for each of that reaction's products.
            # e.g.: outgoing 1-n reactions: [A1 -> B1+C1, D1 -> E1+F1+G1] =>
            #         outgoing_1_n = [ [R(B1), R(C1)], [R(E1), R(F1), R(G1)] ]
            #         where R(B1), R(C1), etc. are each something like { (B,), (H, J), ... }
            outgoing_1_n = [map(get_reachable_resting_states,B.products) for B in outgoing_reactions if B.arity[1] != 1]
            

            
            # Now we need to find the reachable resting state multisets due to outgoing 1-n reactions.
            # for each (outgoing 1-n) reaction r, let prods = the list of reachable 
            # resting state multisets for each of r's products i.e. prods = [R(Xn) for Xn in r.products]. 
            # Each element in prods is therefore a frozenset (reachable resting states) of tuples (multisets).  
            # Take all cartesian products of prods. This will result in all possible sets containing one element from
            # each frozenset in prods. 
            # For example: if prods = [ { (A,B), (C,) }, { (D,), (E,F) } ], cartesian(prods) = [ ((A,B),(D,)), ((A,B), (E,F)), ((C,),(D,)), ((C,),(E,F)) ]. 
            # We then flatten the outer tuples, so we get [ (A,B,D), (A,B,E,F), (C,D), (C,E,F) ]
            # by summing each inner tuple. This is now a list of reachable resting state multisets for
            # each outgoing 1-n reaction from this SCC.
            
            # outgoing_1_n = map(lambda prods: sorted(tuple_sum( itertools.product(*prods) )),outgoing_1_n)
            # outgoing_1_n = map(lambda prods: sorted( itertools.chain(*itertools.product(*prods)) ),outgoing_1_n)
            
            # outgoing_1_n becomes a list (reactions) of lists (products) of tuples (multisets)
            outgoing_1_n = cartesian_multisets(outgoing_1_n)
            
            # We want now to generate the set of reachable resting state multsets 
            # outgoing_1_n is a list of lists of tuples (multisets); we want to make a frozenset
            # of tuples, by concatenating the inner lists and converting the result to a frozenset.
            outgoing_1_n = frozenset().union(*outgoing_1_n)
 
            
            
            # FINALLY, the set of reachable resting state multisets for this SCC is the union or the 
            # set of reachable resting state multisets from outgoing 1-1 reactions and outgoing 1-n 
            # reactions
            resting_state_target = outgoing_1_1 | outgoing_1_n
            
            for c in scc:
                if c in resting_state_targets: raise Exception()
                resting_state_targets[c] = resting_state_target
        
        
        
        # print "End process SCC: " + str(scc)
    
    
    # Cache some things for speed
    #complex_set = set(enumerator.complexes)
    #reactions_producing = get_reactions_producing(enumerator.complexes,enumerator.reactions) 
    reactions_consuming = get_reactions_consuming(enumerator.complexes,enumerator.reactions) 
    
    # Compute list of SCCs
    SCCs = tarjans(enumerator.complexes,enumerator.reactions,reactions_consuming)
    
    # Map each complex to the SCC which contains the complex (each complex 
    # should be in 1 SCC)
    SCC_containing.update({ c: scc for scc in SCCs for c in scc })
    
    # Generate resting state multisets by processing each SCC. Outer loop
    # guarantees that all SCCs will be processed even if DFS recursion
    # misses some
    for scc in SCCs: 
        process_scc(scc)
        
    # Map each complex to the set of multisets of resting states which it can reach    
    resting_state_targets = { complex: ReachableRestingStates(resting_state_targets[complex]) for complex in resting_state_targets }
        
    # Generate condensed reactions
    condensed_reactions = set()
    
    
    for reaction in enumerator.reactions:
        
        # Filter reactions which are not fast
        if( is_fast(reaction) ):
            continue
        
        # Filter reactions with no products/reactants
        if (len(reaction.reactants) == 0) and (len(reaction.products) == 0):
            continue
        
        
        reactant_resting_states = map(get_reachable_resting_states,reaction.reactants)
        product_resting_states = map(get_reachable_resting_states,reaction.products)
        
        # Take cartesian products of reaction product reachable resting state sets multisets to get all 
        # combinations of resting states
        # e.g. for X1 -> A1 + B1, R(A1) = { {A}, {A', A''} }, R(B1) = { {B, B'}, {B''} },
        #    we get [ ({A},{B,B'}), ({A',A''},{B,B'}), ({A},{B''}), ({A',A''},{B''}) ]
        
        # Then, sum each individual tuple in the above list, so we get:
        #    [ {A,B,B'}, {A',A'',B,B'}, {A,B''}, {A',A'',B''} ]
        # This is the list of possible product sets for the new reaction
        
        # Note: { multisets } are actually represented as tuples in this implementation, so 
        # the above lists are more like [((A,), (B, B')), ((A',A''),(B,B')), ...] and
        # [ (A,B,B'), (A',A'',B,B'), ... ], respectively.
        
        # Note 2: product_resting_states is a *list* of frozensets of tuples, so we use the
        # extended call notation itertools.product(*product_resting_states), which is like
        # writing itertools.product(product_resting_states[0], product_resting_states[1], 
        # ..., product_resting_states[n]).
        new_product_combinations = map(tuple_sum,
                                       itertools.product(*product_resting_states))    
        
        # We could do the same thing for reactants; in principle, all resting states of reactants 
        # in bimolecular reactions should be of the form R(X1) = { {X} }, so the results of
        # this call should be trivial:
        #     new_reactant_combinations = map(tuple_sum,
        #                                    itertools.product(*reactant_resting_states))
        # Therefore we'll skip the itertools call and just map the reactants to their 
        # resting states. 
        # new_reactant_combinations = [tuple(r_xn.get_singleton() for r_xn in reactant_resting_states)]
        
        
        new_reactant_combinations = map(tuple_sum,
                                        itertools.product(*reactant_resting_states))
        
        # TODO: throw exception if R(X1) has multiple elements
        for (i, r_xn) in enumerate(reactant_resting_states):
            if(not r_xn.is_singleton()):
                logging.error("Cannot condense reactions: R(%s) has multiple elements: %s"
                                % ( str(reaction.reactants[i]), str(r_xn) ))
                raise Exception()
            
        
        # Now, we'll take all of the combinations of reactants and products
        new_reactant_product_combinations = itertools.product(new_reactant_combinations,
                                                              new_product_combinations)
        
        # And generate new reactions with each combination which is not trivial (reactants == products)
        for (reactants,products) in new_reactant_product_combinations:
            reactants = sorted(reactants)
            products = sorted(products)
            
            # Prune trivial reactions
            if(reactants != products):
                condensed_reactions.add(ReactionPathway('condensed', list(reactants), list(products)))
        
#        print "End reaction "+str(reaction)
    
    return {
                'resting_states': resting_states.values(),
                'resting_state_map': resting_states,
                'resting_state_targets':resting_state_targets,
                'condensed_reactions':list(condensed_reactions),
                'reactions': list(condensed_reactions)
             }
    
    #return (resting_states,resting_state_targets,condensed_reactions)


condense_resting_states = condense_graph