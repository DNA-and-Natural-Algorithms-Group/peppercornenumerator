#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

version = __import__('peppercornenumerator').__version__

setup(
    name='peppercornenumerator',
    version=version,
    description='Domain-level nucleic acid reaction enumerator)',
    long_description=readme,
    author='Karthik Sarma, Casey Grun, and Erik Winfree',
    author_email='winfree@caltech.edu',
    # url='http://www.dna.caltech.edu/peppercorn/',
    license=license,
    test_suite='tests',
    install_requires=['argparse>=1.2.1', 'nose', 'pyparsing', 'numpy', 'dsdobjects'],
    packages=['peppercornenumerator'],
    scripts=['scripts/peppercorn']
)
