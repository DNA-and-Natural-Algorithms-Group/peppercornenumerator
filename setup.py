#!/usr/bin/env python

from setuptools import setup, find_packages

LONG_DESCRIPTION="""
This package enumerates domain-level strand displacement (DSD) reaction
networks assuming low species concentrations, such that unimolecular reaction
pathways always equilibrate before bimolecular reactions initiate. The
enumerator can handle arbitrary non-pseudoknotted secondary structures and
supports a complete set of unimolecular and bimolecular domain-level reactions:
bind/unbind reactions, 3-way branch-migration and 4-way branch-migration
reactions and remote-toehold migration. For more background on reaction
semantics we refer to the publication [Grun et al. (2014)].
"""

setup(
    name='peppercornenumerator',
    version='0.7',
    description='Domain-level nucleic acid reaction enumerator',
    long_description=LONG_DESCRIPTION,
    author='Casey Grun, Stefan Badelt, Karthik Sarma, Brian Wolfe, Seung Woo Shin and Erik Winfree',
    author_email='winfree@caltech.edu',
    license='MIT',
    test_suite='tests',
    install_requires=[
        'numpy',
        'pyparsing>=1.5.5', 
        'crnsimulator>=0.5',
        'dsdobjects>=0.6.3'],
    packages=['peppercornenumerator'],
    scripts=['scripts/peppercorn',
             'scripts/pilsimulator']
)

