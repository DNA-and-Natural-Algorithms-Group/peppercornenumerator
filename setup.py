#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'peppercornenumerator',
    version = '1.1',
    description = 'Domain-level nucleic acid reaction enumerator',
    long_description = LONG_DESCRIPTION,
    long_description_content_type = 'text/markdown',
    author = 'Stefan Badelt, Casey Grun, Karthik Sarma, Brian Wolfe, Seung Woo Shin and Erik Winfree',
    author_email = 'winfree@caltech.edu',
    maintainer = 'Stefan Badelt',
    maintainer_email = 'bad-ants-fleet@posteo.eu',
    url = 'http://www.github.com/DNA-and-Natural-Algorithms-Group/peppercorn/',
    license = 'MIT',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        ],
    python_requires = '>=3.7',
    install_requires = [
        'numpy',
        'pandas',
        'natsort',
        'dsdobjects>=0.8',
        'crnsimulator>=0.9'],
    packages = ['peppercornenumerator'],
    test_suite = 'tests',
    entry_points = {
        'console_scripts': [
            'peppercorn=peppercornenumerator.peppercorn:main',
            'pilsimulator=peppercornenumerator.pilsimulator:main'
            ],
        }
)

