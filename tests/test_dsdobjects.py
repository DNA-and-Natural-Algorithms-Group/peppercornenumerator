
import os
import sys
import unittest
import peppercornenumerator.dsdobjects as objects

SKIP = True

#@unittest.skipIf(SKIP, "skipping tests")
class UtilityTests(unittest.TestCase):
    def setUp(self):
        pass

    def test_make_pair_table(self):
        res = objects.make_pair_table('(((...)))')
        exp = [[(0, 8), (0, 7), (0, 6), None, None, None, (0, 2), (0, 1), (0, 0)]]
        self.assertEqual(res, exp)

        res = objects.make_pair_table('(((+...)))')
        exp = [[(1, 5), (1, 4), (1, 3)], [None, None, None, (0, 2), (0, 1), (0, 0)]]
        self.assertEqual(res, exp)

        res = objects.make_pair_table('((((+))).)')
        exp = [[(1, 4), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), None, (0, 0)]]
        self.assertEqual(res, exp)

        res = objects.make_pair_table('((((&))).)', strand_break='&')
        exp = [[(1, 4), (1, 2), (1, 1), (1, 0)], [(0, 3), (0, 2), (0, 1), None, (0, 0)]]
        self.assertEqual(res, exp)

        with self.assertRaises(objects.DSDObjectsError):
            res = objects.make_pair_table('((((+)).)')

        with self.assertRaises(objects.DSDObjectsError):
            res = objects.make_pair_table('((((&))).)')


#@unittest.skipIf(SKIP, "skipping tests")
class Test_SequenceConstraint(unittest.TestCase):
    def setUp(self):
        pass

    def test_add_constraints(self):
        self.seq1 = 'ACGTAGTTC'
        self.seq2 = 'NNNNNNNNN'
        self.seq3 = 'RRRRRRRRR'
        self.seq4 = 'HHHHHHHHH'
        self.seq5 = 'AAAAAAAAA'
        self.seq6 = 'TGYRTYRRG'

        con1 = objects.SequenceConstraint(self.seq1)
        con2 = objects.SequenceConstraint(self.seq2)
        con3 = objects.SequenceConstraint(self.seq3)
        con4 = objects.SequenceConstraint(self.seq4)
        con5 = objects.SequenceConstraint(self.seq5)
        con6 = objects.SequenceConstraint(self.seq6)

        self.assertEqual(con1, con1 + con2)

        self.assertEqual('AAAAAAAAA', str(con3 + con4))
        self.assertEqual(con5, con3 + con4)

        with self.assertRaises(objects.DSDObjectsError):
            con1 + con3

        self.assertEqual(str(con2), con2.complement)
        self.assertEqual('TGYRTYRRG', con1.complement)
        self.assertEqual(con6, ~con1)
        self.assertEqual(len(con1), 9)


#@unittest.skipIf(SKIP, "skipping tests")
class DSD_DomainObjectTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        objects.DSD_Domain.id_counter = 0
        objects.DSD_Domain.dictionary = dict()

    def test_DSD_DomainInit(self):
        doodle = objects.DSD_Domain(list('Y' * 5), prefix='doodle')

        self.assertIsInstance(doodle, objects.DSD_Domain, "doodle is a DSD_Domain")
        self.assertEqual(str(doodle), '{}'.format(doodle.name), "print DSD_Domain")
        self.assertEqual(len(doodle), 5, "DSD_Domain length")
        self.assertEqual(doodle.sequence, list('Y' * 5), "DSD_Domain sequence")

        moodle = objects.DSD_Domain(list('Y' * 5))
        self.assertEqual(str(moodle), '{}'.format(moodle.name),
                         "Automatic DSD_Domain Name")

    def test_DSD_ComplementDomainInit(self):
        foo = objects.DSD_Domain(list('Y' * 5))

        # Conflicting Constraints
        with self.assertRaises(objects.DSDObjectsError):
            bar = foo.get_ComplementDomain(list('R' * 3))
        with self.assertRaises(objects.DSDObjectsError):
            bar = foo.get_ComplementDomain(list('Y' * 5))
        with self.assertRaises(objects.DSDObjectsError):
            foo.update_constraints(list('R' * 6))

        bar = foo.get_ComplementDomain(list('R' * 5))

        moo = ~foo
        self.assertTrue(bar == moo, "bar is moo")
        self.assertTrue(bar is moo, "bar is moo")
        self.assertTrue(foo == ~bar, "foo is complement of bar")
        self.assertFalse(moo == ~bar, "moo is not complement of bar")
        self.assertFalse(foo == bar)
        self.assertTrue(bar.is_complement, "bar is complement")
        self.assertFalse(foo.is_complement, "foo is complement")

    def test_domains_of_domains(self):
        d1aa = objects.DSD_Domain(list('N'))
        d1ab = objects.DSD_Domain(list('Y'))
        d1a = objects.DSD_Domain([d1aa, d1ab])

        d1b = objects.DSD_Domain(list('RR'))
        d1c = objects.DSD_Domain(list('NN'))
        d1 = objects.DSD_Domain([d1a, d1b, d1c])

        self.assertIsInstance(d1, objects.DSD_Domain, "d1 is a DSD_Domain")

        for d in d1.sequence:
            self.assertIsInstance(d, objects.DSD_Domain, "Sequence of DSD_Domains")

        self.assertEqual(len(d1), 3, "Length of DSD_Domain Sequence")
        self.assertListEqual(d1.base_sequence, list('NYRRNN'))
        self.assertEqual(d1.base_length, 6)

    def test_complement_domains_of_domains_of_domains(self):
        d1aa = objects.DSD_Domain(list('N'))
        d1ab = objects.DSD_Domain(list('Y'))
        d1a = objects.DSD_Domain([d1aa, d1ab])

        d1b = objects.DSD_Domain(list('RR'))
        d1c = objects.DSD_Domain(list('NN'))
        d1 = objects.DSD_Domain([d1a, d1b, d1c])
        with self.assertRaises(NotImplementedError):
            d2 = d1.get_ComplementDomain(list('R' * 6))

        with self.assertRaises(NotImplementedError):
            d1.update_constraints(list('R' * 6))

    def dont_test_complementarity(self):
        # TODO: This is not properly implemented yet.
        # e.g. foo is in ~bar
        pass


@unittest.skipIf(SKIP, "skipping tests")
class DSD_ComplexObjectTest(unittest.TestCase):

    def setUp(self):
        self.d1 = objects.DSD_Domain(list('Y' * 5))
        self.d2 = objects.DSD_Domain(list('N' * 5))
        self.d3 = objects.DSD_Domain(list('R' * 5))
        self.d1c = self.d1.get_ComplementDomain(list('N' * 5))
        self.d2c = self.d2.get_ComplementDomain(list('D' * 5))
        self.d3c = self.d3.get_ComplementDomain(list('H' * 5))

    def tearDown(self):
        objects.DSD_Domain.id_counter = 0
        objects.DSD_Complex.id_counter = 0

    def test_DSD_ComplexInit(self):
        # NOTE: There is no particular reason for this Error, so it might
        # change!
        with self.assertRaises(objects.DSDObjectsError):
            foo = objects.DSD_Complex()

        foo = objects.DSD_Complex(sequence=list('RNNNY'), structure=list('(...)'))
        self.assertIsInstance(foo, objects.DSD_Complex)

        self.assertEqual(foo.sequence, list('RNNNY'))
        self.assertEqual(foo.structure, list('(...)'))
        self.assertEqual(foo.lol_sequence, [list('RNNNY')])
        self.assertEqual(foo.nucleotide_sequence, list('RNNNY'))
        self.assertEqual(foo.rotate_once, foo)
        for r in foo.rotate:
            self.assertEqual(r.sequence, list('RNNNY'))
            self.assertEqual(r.structure, list('(...)'))

    def test_DSD_ComplexDomains(self):
        foo = objects.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+', self.d1,
                                        '+', self.d1c, self.d3c, self.d1c, self.d2], structure=list('..(+(+))..'))

        self.assertEqual(foo.sequence,
                         [self.d1,
                          self.d2,
                          self.d3,
                          '+',
                          self.d1,
                          '+',
                          self.d1c,
                          self.d3c,
                          self.d1c,
                          self.d2])
        self.assertEqual(foo.lol_sequence, [[self.d1, self.d2, self.d3], [
                         self.d1], [self.d1c, self.d3c, self.d1c, self.d2]])
        self.assertEqual(foo.nucleotide_sequence, list(
            'YYYYYNNNNNRRRRR+YYYYY+RRRRRYYYYYRRRRRNNNNN'))

        bar = objects.DSD_Complex(sequence=[self.d1c, self.d3c, self.d1c, self.d2, '+',
                                        self.d1, self.d2, self.d3, '+', self.d1], structure=list('((..+..)+)'))

        self.assertEqual(foo, bar)
        self.assertTrue(foo == bar)

    def test_names(self):
        objects.reset_names()
        foo = objects.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+', self.d1,
                                        '+', self.d1c, self.d3c, self.d1c, self.d2], structure=list('..(+(+))..'))
        self.assertEqual(foo.name, 'cplx0')
        with self.assertRaises(objects.DSDObjectsError):
            foo.name = 'bar'

        objects.reset_names()

        foo.name = 'foo'
        self.assertEqual(foo.name, 'foo')
        self.assertSetEqual(objects.DSD_Complex.names, set(['foo']))

        foo.name = 'bar'
        self.assertEqual(foo.name, 'bar')
        self.assertSetEqual(objects.DSD_Complex.names, set(['bar']))

        objects.reset_names()

    def test_rotations(self):
        foo = objects.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+', self.d1,
                                        '+', self.d1c, self.d3c, self.d1c, self.d2], structure=list('..(+(+))..'))
        self.assertEqual(foo.rotate_once, foo)
        self.assertEqual(foo.sequence,
                         [self.d1,
                          '+',
                          self.d1c,
                          self.d3c,
                          self.d1c,
                          self.d2,
                          '+',
                          self.d1,
                          self.d2,
                          self.d3])
        self.assertEqual(foo.structure, list('(+)(..+..)'))
        self.assertEqual(foo.nucleotide_sequence, list(
            'YYYYY+RRRRRYYYYYRRRRRNNNNN+YYYYYNNNNNRRRRR'))
        for r in foo.rotate:
            self.assertEqual(r, foo)
        self.assertEqual(foo.sequence,
                         [self.d1,
                          '+',
                          self.d1c,
                          self.d3c,
                          self.d1c,
                          self.d2,
                          '+',
                          self.d1,
                          self.d2,
                          self.d3])
        self.assertEqual(foo.structure, list('(+)(..+..)'))
        self.assertEqual(foo.nucleotide_sequence, list(
            'YYYYY+RRRRRYYYYYRRRRRNNNNN+YYYYYNNNNNRRRRR'))


@unittest.skipIf(SKIP, "skipping tests")
class TestTubeTests(unittest.TestCase):
    def setUp(self):
        self.d1 = objects.DSD_Domain(list('Y' * 5))
        self.d2 = objects.DSD_Domain(list('N' * 5))
        self.d3 = objects.DSD_Domain(list('R' * 5))
        self.d1c = self.d1.get_ComplementDomain(list('N' * 5))
        self.d2c = self.d2.get_ComplementDomain(list('D' * 5))
        self.d3c = self.d3.get_ComplementDomain(list('H' * 5))

        self.cplx1 = objects.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+', self.d1,
                                               '+', self.d1c, self.d3c, self.d1c, self.d2], structure=list('..(+(+))..'))

        self.cplx2 = objects.DSD_Complex(sequence=[self.d1, self.d2, self.d3, '+',
                                               self.d3c, self.d2c, self.d2], structure=list('.((+)).'))

    def tearDown(self):
        objects.reset_names()

    def test_TestTubeInit(self):
        foo = objects.TestTube()
        foo.add_complex(self.cplx1, (None, None))
        self.assertTrue(foo.has_complex(self.cplx1))

    def test_TestTubeSum(self):
        foo = objects.TestTube()
        bar = objects.TestTube()
        foo.add_complex(self.cplx1, (None, True))
        bar.add_complex(self.cplx2, (float('inf'), None))

        # Both versions should work
        foobar = foo + bar
        foobar = sum([foo, bar])

        # Check if original TestTubes remain unchanged
        self.assertTrue(foo.has_complex(self.cplx1))
        self.assertFalse(foo.has_complex(self.cplx2))

        self.assertTrue(bar.has_complex(self.cplx2))
        self.assertFalse(bar.has_complex(self.cplx1))

        # Check if new TestTube has everything
        self.assertTrue(foobar.has_complex(self.cplx1))
        self.assertTrue(foobar.has_complex(self.cplx2))

        # Check if attributes were copied correctly
        self.assertEqual(
            foo.get_complex_concentration(
                self.cplx1), foobar.get_complex_concentration(
                self.cplx1))
        self.assertEqual(
            bar.get_complex_concentration(
                self.cplx2), foobar.get_complex_concentration(
                self.cplx2))


@unittest.skipIf(SKIP, "skipping tests")
class TestTubeIOTest(unittest.TestCase):
    def setUp(self):
        self.t0 = objects.DSD_Domain(list('N' * 5), prefix='t')
        self.t1 = objects.DSD_Domain(list('N' * 5), prefix='t')
        self.t2 = objects.DSD_Domain(list('N' * 5), prefix='t')
        self.d3 = objects.DSD_Domain(list('N' * 15), prefix='d')
        self.d4 = objects.DSD_Domain(list('N' * 15), prefix='d')
        self.d5 = objects.DSD_Domain(list('N' * 15), prefix='d')
        self.d6 = objects.DSD_Domain(list('N' * 15), prefix='d')
        self.d7 = objects.DSD_Domain(list('N' * 15), prefix='d')
        self.t0c = self.t0.get_ComplementDomain(list('N' * 5))
        self.t1c = self.t1.get_ComplementDomain(list('N' * 5))
        self.t2c = self.t2.get_ComplementDomain(list('N' * 5))
        self.d3c = self.d3.get_ComplementDomain(list('N' * 15))
        self.d4c = self.d4.get_ComplementDomain(list('N' * 15))
        self.d5c = self.d5.get_ComplementDomain(list('N' * 15))
        self.d6c = self.d6.get_ComplementDomain(list('N' * 15))
        self.d7c = self.d7.get_ComplementDomain(list('N' * 15))

    def test_IO_dna(self):
        # NOTE: This function writes to a file, so we need to compare it to
        # pre-written output-file. Doesn't feel necessary at this point ...
        t0 = self.t0
        t1 = self.t1
        t2 = self.t2
        d3 = self.d3
        d4 = self.d4
        d5 = self.d5
        d6 = self.d6
        d7 = self.d7
        t0c = self.t0c
        t1c = self.t1c
        t2c = self.t2c
        d3c = self.d3c
        d4c = self.d4c
        d5c = self.d5c
        d6c = self.d6c
        d7c = self.d7c

        sequence = [
            d4,
            t0,
            '+',
            d6,
            t2,
            '+',
            d3,
            '+',
            d3,
            t2,
            '+',
            t0,
            d5,
            t1,
            '+',
            t0,
            d7,
            t1,
            '+',
            t1c,
            d7c,
            t1c,
            d5c,
            t2c,
            '+',
            d3c,
            d3c,
            t2c,
            d6c,
            t0c,
            d4c,
            t0c]
        structure = [
            '(',
            '(',
            '+',
            '(',
            '(',
            '+',
            '(',
            '+',
            '(',
            '(',
            '+',
            '.',
            '(',
            '(',
            '+',
            '.',
            '(',
            '(',
            '+',
            ')',
            ')',
            ')',
            ')',
            ')',
            '+',
            ')',
            ')',
            ')',
            ')',
            ')',
            ')',
            '.']

        foo = objects.DSD_Complex(sequence=sequence, structure=structure)
        fooIO = objects.TestTubeIO(
            objects.TestTube(
                complexes={
                    foo.name: (
                        foo,
                        float("inf"),
                        None)}))

        # fooIO.write_dnafile(sys.stdout)
        f = open(os.devnull, "w")
        fooIO.write_dnafile(f)

    def test_IO_kernel(self):
        t0 = self.t0
        t1 = self.t1
        t2 = self.t2
        d3 = self.d3
        d4 = self.d4
        d5 = self.d5
        d6 = self.d6
        d7 = self.d7
        t0c = self.t0c
        t1c = self.t1c
        t2c = self.t2c
        d3c = self.d3c
        d4c = self.d4c
        d5c = self.d5c
        d6c = self.d6c
        d7c = self.d7c

        sequence = [
            d4,
            t0,
            '+',
            d6,
            t2,
            '+',
            d3,
            '+',
            d3,
            t2,
            '+',
            t0,
            d5,
            t1,
            '+',
            t0,
            d7,
            t1,
            '+',
            t1c,
            d7c,
            t1c,
            d5c,
            t2c,
            '+',
            d3c,
            d3c,
            t2c,
            d6c,
            t0c,
            d4c,
            t0c]
        structure = [
            '(',
            '(',
            '+',
            '(',
            '(',
            '+',
            '(',
            '+',
            '(',
            '(',
            '+',
            '.',
            '(',
            '(',
            '+',
            '.',
            '(',
            '(',
            '+',
            ')',
            ')',
            ')',
            ')',
            ')',
            '+',
            ')',
            ')',
            ')',
            ')',
            ')',
            ')',
            '.']

        foo = objects.DSD_Complex(sequence=sequence, structure=structure)
        fooIO = objects.TestTubeIO(
            objects.TestTube(
                complexes={
                    foo.name: (
                        foo,
                        float("inf"),
                        None)}))

        f = open(os.devnull, "w")
        fooIO.write_pil_kernel(f)


if __name__ == '__main__':
    unittest.main()
