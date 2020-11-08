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

