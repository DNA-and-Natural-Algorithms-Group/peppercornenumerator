#!/usr/bin/env python
#
#  test_utils.py
#  EnumeratorProject
#
import unittest

from peppercornenumerator.utils import wrap, natural_sort

class MiscTests(unittest.TestCase):
    def test_wrap(self):
        assert wrap(0, 3) == 0
        assert wrap(2, 3) == 2
        assert wrap(3, 3) == 0
        assert wrap(4, 3) == 1
        assert wrap(-1, 3) == 2
        assert wrap(-2, 3) == 1

    def test_natural_sort(self):
        assert natural_sort(['c10', 'b', 'c2', 'c1', 'a']) == [
            'a', 'b', 'c1', 'c2', 'c10']

if __name__ == '__main__':
    unittest.main()
