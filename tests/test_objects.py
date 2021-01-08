#!/usr/bin/env python
#
#  test_objects.py
#  EnumeratorProject
#
import unittest

from peppercornenumerator.objects import (SingletonError, 
                                          clear_memory,
                                          ObjectInitError,
                                          PepperDomain,
                                          PepperComplex,
                                          Loop)

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class TestPepperDomain(unittest.TestCase):
    def tearDown(self): 
        clear_memory()

    def test_domain_length(self):
        foo = PepperDomain('foo', dtype = 'short')
        bar = PepperDomain('bar', dtype = 'long')
        assert foo.length < bar.length

    def test_new_constructor(self):
        with self.assertRaises(SingletonError):
            bar = PepperDomain('foo')

        foo = PepperDomain('foo', dtype = 'short')
        bar = PepperDomain('foo', dtype = 'short')
        self.assertTrue(foo is bar)

        bar = PepperDomain('foo')
        self.assertTrue(foo is bar)

        with self.assertRaises(SingletonError):
            bar = PepperDomain('foo', dtype = 'long')

class TestLoop(unittest.TestCase):
    def setUp(self):
        self.d1 = PepperDomain('d1', 15)
        self.d2 = PepperDomain('d2', 15)
        self.d3 = PepperDomain('d3', 15)
        self.d4 = PepperDomain('d4', 15)

    def tearDown(self):
        clear_memory()

    def _triple(self, cplx, loc):
        return None if loc is None else (cplx.get_domain(loc), cplx.get_paired_loc(loc), loc)

    def testLoop_01(self):
        d1, d2, d3, d4 = self.d1, self.d2, self.d3, self.d4
        cplx = PepperComplex([d1, d2, d3, ~d3, '+', d4], list('..()+.'))

        llocs = [(0, 0), (0, 1), (0, 2), None, (1, 0)]
        parts = [self._triple(cplx, l) for l in llocs]

        loop = Loop(parts)
        assert loop.is_open
        assert loop.bases == 15 * 3
        assert loop.stems == 1
        with self.assertRaises(AssertionError):
            loop.dlength == 15 * 4

        llocs = [(0, 0), (0, 1), (0, 2), (0, 3), None, (1, 0)]
        parts = [self._triple(cplx, l) for l in llocs]
        with self.assertRaises(ObjectInitError):
            loop = Loop(parts)

if __name__ == '__main__':
  unittest.main()
