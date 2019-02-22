#
#  test_utils.py
#  EnumeratorProject
#
#  Created by Karthik Sarma on 4/21/10.
#

import unittest

import logging
logging.disable(logging.CRITICAL)

from peppercornenumerator.utils import *
from peppercornenumerator.input import read_pil
from peppercornenumerator.objects import clear_memory

class MiscTests(unittest.TestCase):
    def testWrap(self):
        assert wrap(0, 3) == 0
        assert wrap(2, 3) == 2
        assert wrap(3, 3) == 0
        assert wrap(4, 3) == 1
        assert wrap(-1, 3) == 2
        assert wrap(-2, 3) == 1

    def testNaturalSort(self):
        assert natural_sort(['c10', 'b', 'c2', 'c1', 'a']) == [
            'a', 'b', 'c1', 'c2', 'c10']

class LoopTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def testLoop(self):
        complexes, _ = read_pil("C = 1 2 3() + 4")

        parts = []
        domain_length = 0
        complex = complexes['C']

        # skip the closing helical domain
        for loc in [(0, 0), (0, 1), (0, 2), None, (1, 0)]:
            if loc is None:
                parts.append(None)
                continue

            domain = complex.get_domain(loc)
            parts.append((domain, complex.get_structure(loc), loc))

            # we'll naively assume all domains have the same length for the
            # sake of this example
            domain_length = len(domain)

        # construct a loop from parts
        loop = Loop(parts)

        assert loop.bases == domain_length * 3
        assert loop.stems == 1
        assert loop.is_open

    def testLoop_fail(self):
        complexes, _ = read_pil("C = 1 2 3() + 4")

        parts = []
        domain_length = 0
        complex = complexes['C']

        # do NOT skip the closing helical domain
        for loc in [(0, 0), (0, 1), (0, 2), (0,3), None, (1, 0)]:
            if loc is None:
                parts.append(None)
                continue

            domain = complex.get_domain(loc)
            parts.append((domain, complex.get_structure(loc), loc))

        # construct a loop from parts
        with self.assertRaises(PeppercornUsageError):
            loop = Loop(parts)

if __name__ == '__main__':
    unittest.main()
