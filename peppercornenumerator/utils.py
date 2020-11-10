#
#  peppercornenumerator/utils.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

import re

class PeppercornUsageError(Exception):
    pass

def wrap(x, m):
    """
    Mathematical modulo; wraps x so that 0 <= wrap(x,m) < m. x can be negative.
    """
    return (x % m + m) % m

def natural_sort(l):
    """
    Sorts a collection in the order humans would expect. Implementation from
    http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    """
    def convert(text): 
        return int(text) if text.isdigit() else text.lower()
    def alphanum_key(key): 
        return [convert(c) for c in re.split('([0-9]+)', str(key))]
    return sorted(l, key=alphanum_key)

def tarjans(complexes, products):
    """ Tarjans algorithm to find strongly connected components. """
    def strongconnect(cplx, index):
        cplx.tarjans_index = index
        cplx.tarjans_lowlink = index
        index += 1
        S.append(cplx)
        for prod in products[cplx]:
            # Product hasn't been traversed; recurse
            if prod.tarjans_index is None :
                index = strongconnect(prod, index)
                cplx.tarjans_lowlink = min(prod.tarjans_lowlink, cplx.tarjans_lowlink)
            # Product is in the current neighborhood
            elif prod in S:
                cplx.tarjans_lowlink = min(prod.tarjans_index, cplx.tarjans_lowlink)
        if cplx.tarjans_lowlink == cplx.tarjans_index:
            scc = []
            while True:
                n = S.pop()
                scc.append(n)
                if n == cplx:
                    break
            SCCs.append(scc)
        return index

    index = 0
    S, SCCs = [], []
    for cplx in complexes:
        cplx.tarjans_index = None
    for cplx in complexes:
        if cplx.tarjans_index is None:
            index = strongconnect(cplx, index)
    return SCCs

