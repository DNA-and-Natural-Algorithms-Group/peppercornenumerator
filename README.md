# peppercornenumerator 

[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/dna-and-natural-algorithms-group/peppercornenumerator)](https://github.com/dna-and-natural-algorithms-group/peppercornenumerator/tags)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/dna-and-natural-algorithms-group/peppercornenumerator?include_prereleases)](https://github.com/dna-and-natural-algorithms-group/peppercornenumerator/releases)
[![GitHub](https://img.shields.io/github/license/dna-and-natural-algorithms-group/peppercornenumerator)](https://opensource.org/licenses/MIT)
[![Travis (.org)](https://api.travis-ci.com/dna-and-natural-algorithms-group/peppercornenumerator.svg)](https://travis-ci.com/github/dna-and-natural-algorithms-group/peppercornenumerator)
[![Codecov](https://img.shields.io/codecov/c/github/dna-and-natural-algorithms-group/peppercornenumerator)](https://codecov.io/gh/dna-and-natural-algorithms-group/peppercornenumerator)


This is a package for domain-level strand displacement (DSD) system analysis.
The reaction enumerator **Peppercorn** reads a file with initially present
domain-level complexes, and returns all possible reactions and products.

Peppercorn supports arbitrary non-pseudoknotted structures and the following
domain-level reactions: bind, open, proximal 3-way and 4-way branch migration,
remote 3-way and 4-way branch migration.  For more background on reaction
semantics we refer to [Badelt et al. (2020)].

Given a specification of initial species concentrations, the simulation
software **Pilsimulator** can read Peppercorn's output and simulate expected
species concentrations over time. Alternatively, you can export the results
into **SBML** format for analysis in your favorite simulator.
Note that the reaction rates assume DNA
molecules!

## Installation
Please use Python 3 if possible. Python 3 is supported starting with version
0.7.  Python 2.7 support will be dropped at version 1.0 (released upon
manuscript acceptance).


```bash
$ python setup.py install
```
Please consider testing the installation first, e.g. using any of the following
commands:
``` 
$ python setup.py test
$ pytest tests
```

## Quickstart using the executable "peppercorn"

### General workflow
After installation of the **peppercornenumerator** package, you must have
two excutable scripts, try if they work and look at the options:
```sh
$ peppercorn --help
$ pilsimulator --help
```

Use the executable **peppercorn** to load the file [example.pil] and write results to example_enum.pil:

```sh
# either using commandline flags
$ peppercorn -o example_enum.pil example.pil
# or read from STDIN and write to STDOUT:
$ cat example.pil | peppercorn > example_enum.pil
```
Your can simulate the enumerated system using the **pilsimulator** executable.
```sh
$ cat example_enum.pil | pilsimulator --t8 1800 --p0 S1=1e-7 S2=1e-7 C1=1e-9 --atol 1e-10 --rtol 1e-10
```
Note that default reaction rate constants assume 'M' concentration units, hence
we use the same units for specification of initial concentrations (--p0). Due
to those small numbers (molar concentrations), we have to specify more
sensitive realtive and absolute tolerances for the solver.
Check commandline options of peppercorn to change units, e.g. to 'nM', as well 
as to provide initial concentrations directly in the input file.

### Input/Output format

The following input format is recommended. The lengths of all domains are
defined first, then all initial complexes are defined in terms of their
sequence and secondary structure. For more details on the **kernel notation**
of complexes, see Figure 1 or Section 2 of [Badelt et al. (2020)]. 

```
# <- this is a comment
#
# Shohei Kotani and William L. Hughes (2017)
# Multi-Arm Junctions for Dynamic DNA Nanotechnology
# 
# Figure 2A: Single-layer catalytic system with three-arm junction substrates.
#

#
# Initialize domains (and their complements) and specify their lengths:
#
length a   = 22  # Domains a and a* with length = 22
length b   = 22
length c   = 22
length t1  = 6   # name = 1 in Figure
length t2  = 6   # name = 2 in Figure
length t3  = 10  # name = 3 in Figure
length T2  = 2

length d1s = 16
length d2  = 6

#
# Initialize all initial complexes using kernel notation, 
# which combines name, sequence and structure into a single line!
# Always use 5' to 3' direction of sequences.
#

# The following complex is called C1, has a single strand with 3 unpaired domains:
C1 = t1 c a

# The complex S1 has multiple strands and is connected via paired domains!
# Sequence of S1:  d1s T2 b a t2 + t2* a* c* t1* + c b* 
# Structure of S1:  .  .  ( ( (  +  )  )  (   .  + ) )
S1 = d1s T2 b( a( t2( + ) ) c*( t1* + ) ) 

S2 = t1( c( a( + t2* ) b*( d2 t3 + ) ) )

P1 = t2* a*( c*( t1*( + ) ) )
I1 = d1s T2 b( a t2 + c )
I2 = d1s T2 b( a( t2( + ) ) b*( d2 t3 + ) c*( t1* + ) )

P2 = d1s T2 b( a( t2( + ) ) ) d2 t3
P3 = b( c*( t1* + ) )

R = d1s( d2( + t3* ) )

D = d1s d2
RW = d1s( T2 b( a( t2( + ) ) ) d2( t3( + ) ) )
```

Let's use reaction condensation for a more compact representation of the
reaction network and increase the default maximum complex size to avoid
(warnings about) incomplete enumeration.
```
$ peppercorn -o system-enum.pil --max-complex-size 10 --condensed < system.pil
```

And then the output file should look something like this. The layout may vary
between different Peppercorn versions.
```
# File generated by peppercorn-v0.6

# Domain specifications 
length a = 22
length b = 22
length c = 22
length d1s = 16
length d2 = 6
length t1 = 6
length T2 = 2
length t2 = 6
length t3 = 10

# Resting complexes 
C1 = t1 c a
D = d1s d2
e48 = t3*( d2*( d1s*( + ) ) + b( c*( t1*( + ) ) a( + t2* ) ) d2 )
e51 = t3*( d2*( d1s*( + ) d2 + b( c*( t1*( + ) ) a( + t2* ) ) ) )
I1 = d1s T2 b( a t2 + c )
P1 = t2* a*( c*( t1*( + ) ) )
P2 = d1s T2 b( a( t2( + ) ) ) d2 t3
P3 = b( c*( t1* + ) )
R = d1s( d2( + t3* ) )
RW = d1s( T2 b( a( t2( + ) ) ) d2( t3( + ) ) )
S1 = d1s T2 b( a( t2( + ) ) c*( t1* + ) )
S2 = t1( c( a( + t2* ) b*( d2 t3 + ) ) )

# Resting macrostates 
macrostate C1 = [C1]
macrostate D = [D]
macrostate e51 = [e51, e48]
macrostate I1 = [I1]
macrostate P1 = [P1]
macrostate P2 = [P2]
macrostate P3 = [P3]
macrostate R = [R]
macrostate RW = [RW]
macrostate S1 = [S1]
macrostate S2 = [S2]

# Condensed reactions 
reaction [condensed      =       588645 /M/s ] e51 + I1 -> P3 + RW + D + C1
reaction [condensed      =      3083.77 /M/s ] I1 + P1 -> S1 + C1
reaction [condensed      =        3e+06 /M/s ] P2 + R -> RW + D
reaction [condensed      =    1.637e+06 /M/s ] S1 + C1 -> I1 + P1
reaction [condensed      =       588645 /M/s ] S2 + I1 -> P3 + P2 + C1
reaction [condensed      =        3e+06 /M/s ] S2 + R -> e51
```

### Tips & Tricks
  * You can specify concentrations of complexes in the input file, e.g.:

        C1 = t1 c a @initial 10 nM

    The benefit of hard-coding initial concentrations in the input file, is
    that they do not have to be specified again when simulating the output
    using the Pilsimulator software.

  * You can name complexes, even though you do not want them initially present. To do that, give them explicitly 0 initial concentration:

        I1 = d1s T2 b( a t2 + c ) @initial 0 nM

  * Every file produced by Peppercorn, can be used as an input to Peppercorn,
    although some types of input (e.g. macrostate assignments, condensed
    reactions) may be ignored with a warning.  You can use that to specify
    a reaction with a specific type and rate constant and that reaction will
    then be part of the enumerated and condensed output of peppercorn. If
    Peppercorn enumerates an already specified reaction (same type, same
    complexes), the rate will not be updated.

        reaction [branch-4way    =        1e-10 /s   ] I2 -> P3 + P2

    Of course, all reactants and products must be previously defined complexes.
    You can use the supported reaction types: 
    open, bind21, bind11, branch-3way, branch-4way ... but not condensed!

  * You can use a short-hand notation for a *super sequence* that is composed
    of multiple domains. For example, if a domain *k* is sometimes split into
    *k1*, *k2*, and *k3*, e.g. when they are not all bound at the same time,
    you can define

        sup-sequence k = k1 k2 k3

    and use the domain *k* as a shorthand for *k1 k2 k3* when applicable.
    Note, that the output will always use the extended format, explicitly
    listing the most detailed domain composition.
        
  * You can specify the nucleotide sequence of a domain, and it will remain
    present in the ouput file. Hower, note that rates (at the moment) only
    depend on domain length, and not on the nucleotide sequence. Use the keyword
    `sequence` to initialize a domain with the nucleotide sequence:

        sequence a = ACGTUGCA : 8

    The digit at the end denotes the length of the domain, it can be omitted,
    but must not be wrong.

  * Peppercorn supports a range of different input file formats, e.g. the
    specification language of the [Seesaw Compiler]. You can use the option
    `--dry-run` to translate every valid input into kernel format without
    enumeration.

        $ cat seesaw_ciruit.ssw | peppercorn --dry-run

## Version
0.9 -- Python 2.7 and Python 3 compatible, including SBML output standard.

## Authors
Stefan Badelt, Casey Grun, Karthik V. Sarma, Brian Wolfe, Seung Woo Shin and Erik Winfree.

[Badelt et al. (2020)]: <https://doi.org/10.1098/rsif.2019.0866>
[example.pil]: <https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator/blob/development/tests/examples/literature/kotani2017_F2.pil>
[Seesaw Compiler]: <http://www.qianlab.caltech.edu/SeesawCompiler/AOtoSEESAW.php>
