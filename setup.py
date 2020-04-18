#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name='peppercornenumerator',
    version='1.0',
    description='Domain-level nucleic acid reaction enumerator',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author='Stefan Badelt, Casey Grun, Karthik Sarma, Brian Wolfe, Seung Woo Shin and Erik Winfree',
    author_email='winfree@caltech.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        ],
    python_requires='>=3.7',
    test_suite='tests',
    install_requires=[
        'pandas', # for case-studies
        'numpy',
        'crnsimulator>=0.8',
        'dsdobjects>=0.7.1'],
    packages=['peppercornenumerator'],
    entry_points = {
        'console_scripts': [
            'peppercorn=peppercornenumerator.peppercorn:main',
            'pilsimulator=peppercornenumerator.pilsimulator:main'
            ],
        }
)

