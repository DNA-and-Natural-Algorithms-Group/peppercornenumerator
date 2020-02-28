#!/usr/bin/env python
#
#  test_objects.py
#  EnumeratorProject
#
import unittest

import peppercornenumerator.objects as objects

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class PepperDomainTest(unittest.TestCase):

    def setUp(self): 
        pass

    def tearDown(self): 
        objects.clear_memory()

    def test_domain_length(self):
        foo = objects.PepperDomain('foo', dtype = 'short')
        bar = objects.PepperDomain('bar', dtype = 'long')

        self.assertEqual(len(foo), objects.DL_Domain.SHORT_DOM_LEN)
        self.assertEqual(len(bar), objects.DL_Domain.LONG_DOM_LEN)

    def test_new_constructor(self):
        with self.assertRaises(objects.DSDObjectsError):
            bar = objects.PepperDomain('foo')

        foo = objects.PepperDomain('foo', dtype = 'short')
        bar = objects.PepperDomain('foo', dtype = 'short')
        self.assertTrue(foo is bar)

        bar = objects.PepperDomain('foo')
        self.assertTrue(foo is bar)

        with self.assertRaises(objects.PepperObjectsError):
            bar = objects.PepperDomain('foo', dtype = 'long')


if __name__ == '__main__':
  unittest.main()

