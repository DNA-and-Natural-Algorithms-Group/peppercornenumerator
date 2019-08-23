#!/usr/bin/env python

from setuptools import setup, find_packages

LONG_DESCRIPTION="""
This package enumerates domain-level strand displacement (DSD) reaction
networks assuming a timescale separation, such that fast reaction pathways
always equilibrate before slow reactions initiate. The enumerator can handle
arbitrary non-pseudoknotted secondary structures and supports all typical
unimolecular and bimolecular domain-level reactions: bind, unbind, 3-way
branch-migration and 4-way branch-migration including remote-toehold branch
migration. For more background on reaction semantics we refer to the README.
"""

setup(
    name='peppercornenumerator',
    version='0.7',
    description='Domain-level nucleic acid reaction enumerator',
    long_description=LONG_DESCRIPTION,
    author='Stefan Badelt, Casey Grun, Karthik Sarma, Brian Wolfe, Seung Woo Shin and Erik Winfree',
    author_email='winfree@caltech.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        ],
    test_suite='tests',
    install_requires=[
        'future',
        'pandas', # for case-studies
        'numpy',
        'crnsimulator>=0.6',
        'dsdobjects>=0.7'],
    packages=['peppercornenumerator'],
    scripts=['scripts/peppercorn',
             'scripts/pilsimulator']
)

