#!/usr/bin/env python
#
#  test_utils.py
#  EnumeratorProject
#
import unittest

from peppercornenumerator.utils import wrap, tarjans

class MiscTests(unittest.TestCase):
    def test_wrap(self):
        assert wrap(0, 3) == 0
        assert wrap(2, 3) == 2
        assert wrap(3, 3) == 0
        assert wrap(4, 3) == 1
        assert wrap(-1, 3) == 2
        assert wrap(-2, 3) == 1

class mock:
    def __init__(self, n):
        self.name = n
    def __repr__(self):
        return f'{self.name}'

class TarjansTest(unittest.TestCase):
    def test_tarjans(self):
        # https://www.youtube.com/watch?v=wUgWX0nc4NY
        a = mock('a')
        b = mock('b')
        c = mock('c')
        d = mock('d')
        e = mock('e')
        f = mock('f')
        g = mock('g')
        h = mock('h')
        complexes = [a, b, c, d, e, f, g, h]
        products = {a: [b, e],
                    b: [f],
                    c: [b, g, d],
                    d: [g],
                    e: [a, f],
                    f: [c, g],
                    g: [h],
                    h: [d]}

        sccs = tarjans(complexes, products)
        outp = [[d, h, g], [c, f, b], [e, a]]
        ss = sorted([sorted(s, key = repr) for s in sccs], key = repr)
        so = sorted([sorted(s, key = repr) for s in outp], key = repr)
        assert ss == so

if __name__ == '__main__':
    unittest.main()
