#
#  peppercornenumerator/utils.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

import sys
sys.setrecursionlimit(10**6) 

class PeppercornUsageError(Exception):
    pass

def wrap(x, m):
    """
    Mathematical modulo; wraps x so that 0 <= wrap(x,m) < m. x can be negative.
    """
    return (x % m + m) % m

def tarjans(complexes, products):
    """ Tarjans algorithm to find strongly connected components. """
    stack, SCCs = [], []
    def strongconnect(at, index):
        stack.append(at)
        at.index = index
        at.llink = index
        index += 1
        for to in products[at]:
            if to.index is None:
                # Product hasn't been traversed; recurse
                strongconnect(to, index)
            if to in stack:
                # Product is in the current neighborhood
                at.llink = min(at.llink, to.llink)

        if at.index == at.llink:
            # Back to the start, get the SCC.
            scc = [] 
            while True:
                nc = stack.pop()
                nc.llink == at.index
                scc.append(nc)
                if nc == at:
                    break
            SCCs.append(scc)

    for cplx in complexes: 
        cplx.index = None

    for i, cplx in enumerate(complexes):
        if cplx.index is None:
            strongconnect(cplx, i)
    return SCCs

