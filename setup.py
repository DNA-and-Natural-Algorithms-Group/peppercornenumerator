#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name='peppercornenumerator',
    version='0.9',
    description='Domain-level nucleic acid reaction enumerator',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author='Stefan Badelt, Casey Grun, Karthik Sarma, Brian Wolfe, Seung Woo Shin and Erik Winfree',
    author_email='winfree@caltech.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 2.7',
        ],
    python_requires='>=2.7',
    test_suite='tests',
    install_requires=[
        'future',
        'pandas', # for case-studies
        'numpy',
        'crnsimulator>=0.6',
        'dsdobjects>=0.7.1'],
    packages=['peppercornenumerator'],
    entry_points = {
        'console_scripts': [
            'peppercorn=peppercornenumerator.peppercorn:main',
            'pilsimulator=peppercornenumerator.pilsimulator:main'
            ],
        }
)

