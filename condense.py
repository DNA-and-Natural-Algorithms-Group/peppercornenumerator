import itertools as it
import operator
import logging
import numpy as np
from reactions import get_auto_name,ReactionPathway
from utils import Complex, Domain, Strand, RestingState

class SetOfFates(object):
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
        return "SetOfFates([ " + ", ".join(["(" + ", ".join(map(repr,inner)) + ")" for inner in self.states]) + " ])"
        #return 'SetOfFates('+repr(self.states)+")"

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

def times(iterable):
    return reduce(operator.mul,iterable)

def tuple_sum(iterable):
    """
    Given an iterable of tuples, concatenates each of the tuples
    and returns one big tuple. Precisely, for some list of tuples 
    P = [p1, p2, p3 ...]:

        tuple_sum(P) = p1 + p2 + p3 + ... 

    where + represents the concatenation of two tuples, e.g. 
    (a1, a2) + (a3, a4) = (a1, a2, a3, a4).

    Example: 
        tuple_sum([(1,2,3),(4,),(5,6)]) == (1, 2, 3, 4, 5, 6)
    """
    return reduce(operator.concat, tuple(iterable), tuple())

def cartesian_sum(iterable):
    """
    Takes the cartesian sum over a set of multisets (counters). We define the 
    cartesian sum X [+] Y between two sets X and Y as:

        X [+] Y = {x + y : x is in X and y is in Y}

    where + represents the sum of the multisets.

    Equivalently, in terms of the cartesian product X x Y, we can write X [+] Y as:

        X [+] Y = { (tuple_sum_{a in p} a) : p in X x Y }
    """
    return map(tuple_sum,it.product(*iterable))

def cartesian_multisets(reactions):
    """
    Accepts an iterable of iterables of tuples. The outer iterable 
    represents a sequence of reactions, the inner tuple represents 
    products of those reactions.
    """
    # cartesian = []
    # for prods in reactions:
    #     # x = list(it.product(*prods))
    #     # y = map(tuple_sum,x)
        
    #     y = cartesian_sum(prods)
    #     z = sorted(y)
    #     cartesian.append(z)
    # return cartesian
    
    return list(it.chain( sorted(cartesian_sum(prods)) for prods in reactions ))

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


def condense_graph(enumerator, compute_rates=True):

    # Approach: compute SCCs using Tarjan's algorithm, including only fast 
    # 1-1 reactions as edges. For each complex in the SCC, compute the 
    # set of _fates_ of that complex (a fate is a multiset of resting states
    # that are reachable from that complex by fast reactions). For each 
    # detailed reaction, generate one condensed reaction for each combination
    # of the fates of the products.
    # 
    # Detailed description:
    # A fate of a complex is the multiset of resting states that can be reached
    # from the complex by a series of fast reactions. For instance, in the 
    # network A -> B -> C, C -> D, D -> E + F where all reactions are fast,
    # there are three resting states: ^D = { D }, ^E = { E }, ^F = { F }, and
    # the complex A has two fates: { ^D } and { ^E, ^F }. We define F(x) for 
    # some complex x to be the _set_ of fates of x; therefore 
    # F(A) = { { ^D }, { ^E, ^F } }. 
    # 
    # Relatedly, for some reaction r = (A, (b1, b2, ...)) where A is the set of 
    # reactants and B is the set of products, let the fate of the reaction
    # R(r) = F(b1) [+] F(b2) [+] ..., where [+] is the cartesian sum. Finally,
    # let S(x) be the SCC containing complex x. Note that S(x) may or may not
    # be a resting state. 
    # 
    # With these definitions, we can calculate F(x) by a recursion;
    # 
    # F(x) = { { S(x) }            : if S(x) is a resting state
    #        { U_{r in R_o} R(r)   : if S(x) is not a resting state and R_o 
    #                                represents the set of outgoing reactions
    #                                from S(x)
    # 
    # Once we have calculated F(x) (which can be done by a DFS on the DAG formed
    # by the SCCs of the fast, detailed, 1-1 reactions---see the paper for
    # details on why this forms a DAG and therefore F(x) can be calculated in
    # finite time), the set of condensed reactions can be easily computed. We
    # generate one condensed reaction for each combination of fates of the
    # products.
    # 
    # ------------------------------------------------------------------------
    
    # Stores mapping between a complex x and the set of its possible
    # fates, F(x)
    complex_fates = {}
    
    # Stores mappings between an SCC and the associated resting state
    resting_states = {}
    
    # Stores mappings between each complex x and the SCC S(x) which 
    # contains x
    SCC_containing = {}
    
    # Remembers which SCCs we've processed (for which we've computed
    # the fates)
    processed_SCCs = set()

    # The following dicts map SCCs to the various matrices used to 
    # calculate the condensed reaction rates:
    
    # For resting state SCCs:
    # stationary_distributions[SCC] =  s^ : maps complexes to stationary probabilities
    stationary_distributions = {}

    # For transient SCCs:
    # exit_probabilities[SCC] = B : maps (incoming complex, outgoing reaction) tuples to exit 
    # probabilities
    exit_probabilities = {}
    
    # decay_probabilities[(x,F)] = P( x decays to F ), where F is a fate of x
    decay_probabilities = {}

    # Define helper functions for calculating condensed reaction rates
    def get_stationary_distribution(scc):
        scc_set = frozenset(scc)
        scc_list = sorted(scc)
        complex_indices = { c:i for (i,c) in enumerate(scc_list) }

        reactions = [r for c in scc for r in reactions_consuming[c] if (r.arity == (1,1) and not is_outgoing(r,scc_set))]
        L = len(scc)
        T = np.zeros((L,L))

        # add transition rates for each reaction
        for r in reactions:
            a = r.reactants[0]
            b = r.products[0]

            T[complex_indices[a]][complex_indices[b]] = r.rate()

        # calculate eigenvalues
        (w,v) = np.linalg.eig(T)

        # find eigenvector corresponding to eigenvalue zero
        # http://stackoverflow.com/questions/9542738/python-find-in-list
        s = next((v[i,:] for (i, x) in enumerate(w) if x == 0)) # np.ones(L)/L

        # return dict mapping complexes to stationary probabilities
        return { c: s[i] for (c,i) in complex_indices.iteritems() }

    def get_exit_probabilities(scc):
        # build set and list of elements in SCC; assign a numerical index to each complex
        scc_set = frozenset(scc)
        scc_list = sorted(scc)
        complex_indices = { c:i for (i,c) in enumerate(scc_list) }

        # find all fast reactions involving complexes in this SCC
        reactions = [r for c in scc for r in reactions_consuming[c] if is_fast(r)]

        # sort reactions into internal and outgoing; assign a numerical index 
        # to each reaction
        r_outgoing = sorted([r for r in reactions if is_outgoing(r, scc_set)])
        r_internal = sorted([r for r in reactions if not is_outgoing(r, scc_set)])
        exit_indices = { r:i for (i,r) in enumerate(r_outgoing) }

        # L = # of complexes in SCC
        L = len(scc)

        # e = # of exit pathways
        e = len(r_outgoing)

        # add transition rates for each internal reaction
        T = np.zeros((L,L))
        for r in r_internal:
            a = r.reactants[0]
            b = r.products[0]
            T[complex_indices[a]][complex_indices[b]] = r.rate()

        # add transition rates for each outgoing reaction
        Te = np.zeros((L,e))
        for r in r_outgoing:
            a = r.reactants[0]
            Te[complex_indices[a]][exit_indices[r]] = r.rate()            
        

        #  P = (T_{LxL} Te_{Lxe})
        P = np.hstack((T, Te))

        # normalize P along each row, to get the overall transition probabilities
        P = P / np.sum(P,1)[:,np.newaxis]

        # extract the interior transition probabilities
        Q = P[:,0:L]

        # extract the exit probabilities
        R = P[:,L:L+e]

        # calculate the fundamental matrix
        N = np.linalg.inv(np.eye(L) - Q)

        # calculate the exit probabilities
        B = np.dot(N,R)

        # return dict mapping tuples of (incoming complex, outgoing reaction) to exit probabilities
        return { (c, r):B[i,j] for (c,i) in complex_indices.iteritems() for (r,j) in exit_indices.iteritems()  }

    # Define helper functions
    def get_fates(complex):
        """
        Returns the set of possible fates for this complex, calling 
        compute_fates if necessary
        """
        if(complex not in complex_fates):
            compute_fates(SCC_containing[complex])
        return complex_fates[complex]
    
    def get_resting_state(complex):
        """
        Returns the resting state associated with this complex, if 
        one exists.
        """
        scc = SCC_containing[complex]
        scc_set = frozenset(scc)
        if(scc_set not in resting_states):
            compute_fates(scc)
        return resting_states[scc_set]
    
    def compute_fates(scc):
        """
        Processes a single SCC neighborhood, generating resting state multisets 
        for each complex, and storing the mappings in the outer-scope 
        `complex_fates` dict 
        """
        
        # Convert to a set for fast lookup
        scc_set = frozenset(scc)
        
        # Check that we haven't been here already
        if(scc_set in processed_SCCs): 
            return
        
        outgoing_reactions = [r for c in scc for r in reactions_consuming[c] if (is_fast(r) and is_outgoing(r,scc_set)) ]
                
        # Remember that we've processed this neighborhood
        processed_SCCs.add(scc_set)

        
        # If this SCC is a resting state:  
        if(len(outgoing_reactions) == 0):
            
            # build new resting state
            resting_state = RestingState(get_auto_name(),scc)
            resting_states[scc_set] = resting_state
            
            # calculate stationary distribution
            stationary_distributions[resting_state] = get_stationary_distribution(scc)

            # assign fates to each complex in the SCC
            for c in scc:
                if c in complex_fates: raise Exception()
                complex_fates[c] = frozenset([ (resting_state,) ]) #frozenset([ (get_resting_state(c),) ])
                
        # Otherwise, if there are outgoing fast reactions:
        else:

            # Compute the fates of each of the outgoing reactions
            # This is a list, with each element corresponding to the frozenset of fates for one reaction.
            reaction_fates = [ cartesian_sum(map(get_fates,r.products)) for r in outgoing_reactions ]
            

            # calculate the exit probabilities
            exit_probabilities[scc_set] = get_exit_probabilities(scc)

            # calculate the probability of each fate
            # for each outgoing reaction, `r`
            for (i,r) in enumerate(outgoing_reactions):
                # for each possible `fate` of that reaction `r`
                for fate in reaction_fates[i]:
                    # the decay probability will be different for different complexes in the SCC
                    for c in scc:
                        # P(x decays to F) = P(SCC exits via r | SCC was entered via c)
                        decay_probabilities[(c, fate)] = exit_probabilities[scc_set][(c,r)]


            # The set of fates for the complexes in this SCC is the union of 
            # the fates for all outgoing reactions. 
            
            # Note that frozenset().union(*X) === X[0] U X[1] U X[2] ... 
            # where X[i] is a frozenset and a U b represents the union of frozensets a and b 
            set_of_fates = frozenset().union( *reaction_fates )


            for c in scc:
                if c in complex_fates: raise Exception()
                complex_fates[c] = set_of_fates
            
    
    # Cache some things for speed
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
        compute_fates(scc)
        
    # Map each complex to the set of multisets of resting states which it can reach    
    complex_fates = { complex: SetOfFates(complex_fates[complex]) for complex in complex_fates }
        
    # Generate condensed reactions
    condensed_reactions = set()
    for reaction in enumerator.reactions:
        
        # Filter reactions which are fast (these have been handled by the fates)
        if( is_fast(reaction) ):
            continue
        
        # Filter reactions with no products/reactants
        if (len(reaction.reactants) == 0) and (len(reaction.products) == 0):
            continue
        
        
        reactant_fates = map(get_fates,reaction.reactants)
        product_fates = map(get_fates,reaction.products)


        # Take cartesian sum of reactant fates
        new_reactant_combinations = cartesian_sum(reactant_fates)
        
        # Take cartesian sum of reaction product fates to get all 
        # combinations of resting states
        new_product_combinations = cartesian_sum(product_fates)

        
        # make sure each reactant is a resting state (F(xn) is a singleton for each reactant xn)
        for (i, f_xn) in enumerate(reactant_fates):
            if(not f_xn.is_singleton()):
                logging.error("Cannot condense reactions: F(%s) has multiple elements: %s"
                                % ( str(reaction.reactants[i]), str(f_xn) ))
                raise Exception()
            
        
        # Now, we'll take all of the combinations of reactants and products
        new_reactant_product_combinations = it.product(new_reactant_combinations,
                                                              new_product_combinations)
        
        # And generate new reactions with each combination which is not trivial (reactants == products)
        for (reactants,products) in new_reactant_product_combinations:
            reactants = sorted(reactants)
            products = sorted(products)
            
            # Prune trivial reactions
            if(reactants != products):
                reaction = ReactionPathway('condensed', list(reactants), list(products))

                # calculate reaction rate by summing over all representative detailed reactions
                detailed_reactions_consuming = set([ r for reactant in reactants for c in reactant.complexes for r in reactions_consuming[c] ])
                reaction_rate = 0
                for r in detailed_reactions_consuming:
                    
                    # this obnoxious routine is necessary to compute the probability that any possible 
                    detailed_products = sorted(tuple(r.products))
                    condensed_products = products

                    # represents all possible ways the products of the detailed 
                    # reaction could decay into various fates.
                    # list of tuples, each containing several (product,fate) tuples:
                    # ex: [  ((p1,f11), (p2,f21)), ((p1,f12), (p2,f21)), ((p1,f11), (p2,f22)), ... ]
                    fate_combinations = it.product(*[it.izip(it.repeat(product), get_fates(product)) for product in detailed_products])
                    
                    # fate_combinations, but filtered to only include combinations which yield the products of the condensed reaction
                    decay_possibilities = [ combination for combination in fate_combinations if sorted(tuple_sum( f for (p,f) in combination )) == condensed_products]

                    # if there's any way the products of this reaction can decay to the condensed reaction products, this reaction will contribute to the rate of the condensed reaction
                    if len(decay_possibilities) > 0:

                        # probability that the products of this detailed reaction decay into fates that yield the condensed reaction
                        # is the sum of the joint probabilities of each of the decay_possibilities
                        product_probability = sum([  times((decay_probabilities[pair] if pair in decay_probabilities else 0) for combination in decay_possibilities for pair in combination) ])

                        reactant_probabilities = times(stationary_distributions[resting_states[frozenset(SCC_containing[a])]][a] for a in r.reactants)
                        k = r.rate()

                        # overall contribution of detailed reaction r to rate of the condensed reaction ^r = 
                        # P(reactants of ^r are present as reactants of r) * k_r * P(products of r decay to products of ^r)
                        reaction_rate += reactant_probabilities * k * product_probability

                reaction._const = reaction_rate
                condensed_reactions.add(reaction)
            
    return {
                'resting_states': resting_states.values(),
                'resting_state_map': resting_states,
                'resting_state_targets':complex_fates,
                'condensed_reactions':list(condensed_reactions),
                'reactions': list(condensed_reactions)
             }
    
    #return (resting_states,complex_fates,condensed_reactions)


condense_resting_states = condense_graph
