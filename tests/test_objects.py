
import unittest

import peppercornenumerator.objects as objects

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class PepperDomain_Test(unittest.TestCase):
    def setUp(self): 
        pass

    def tearDown(self): 
        objects.clear_memory()

    def test_domain_length(self):
        foo = objects.PepperDomain('foo', dtype = 'short')
        with self.assertRaises(objects.DSDDuplicationError):
            bar = objects.PepperDomain('foo', dtype = 'long')
        bar = objects.PepperDomain('bar', dtype = 'long')

        self.assertEqual(len(foo), objects.DL_Domain.SHORT_DOM_LEN)
        self.assertEqual(len(bar), objects.DL_Domain.LONG_DOM_LEN)

if __name__ == '__main__':
  unittest.main()

