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
from peppercornenumerator.input import from_kernel
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

    def testParseDotParen(self):
        #						012345678
        assert parse_dot_paren('(((...)))') == [
            [(0, 8), (0, 7), (0, 6), None, None, None, (0, 2), (0, 1), (0, 0)]]

        #						0	1
        #						012 012345
        assert parse_dot_paren(
            '(((+...)))') == [[(1, 5), (1, 4), (1, 3)], [None, None, None, (0, 2), (0, 1), (0, 0)]]

        #						0    1
        #						0123 01234
        assert parse_dot_paren(
            '((((+))).)') == [[(1, 4), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), None, (0, 0)]]

    def testParseConcentration(self):
        assert parse_concentration('0.05 uM') == 5e-8
        assert parse_concentration('17e-4 M') == 17e-4
        assert parse_concentration('2e-5pM') == 2e-17
        assert parse_concentration('1pM') == 1e-12

    def testParseParameters(self):
        assert parse_parameters('[1nt]') == {'concentration': None}
        assert parse_parameters('[@1uM]') == {'concentration': 1e-6}
        assert parse_parameters('1nt@1pM') == {'concentration': 1e-12}

    def testFormatSI(self):
        assert format_si(1e-6) == (1, 'u')
        assert format_si(5e-9) == (5, 'n')
        assert format_si(50e-12) == (50, 'p')

class LoopTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_memory()

    def testLoop(self):
        (domains, strands, complexes) = from_kernel(["C = 1 2 3() + 4"])

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
        (domains, strands, complexes) = from_kernel(["C = 1 2 3() + 4"])

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
