# Case study data 
This directory contains the data used for the publication [Badelt et al. (2020)]:

Badelt, Grun, Sarma, Wolfe, Shin, Winfree: 
*A domain-level DNA strand displacement reaction enumerator allowing arbitrary non-pseudoknotted secondary structures*

# Commandline verison
Execute the python script:

    $ ./figure_analysis.py

... you may have to make a directory called **tmp/** where temporary files are
stored. Be aware that the temporary files remain in this directory for further
analysis.

The script reads the experimental data from library files, where each of those
library files corresponds to a particular publictation, currently:

  * zhang2007.py
  * yin2008.py
  * zhang2009.py
  * zhang2010.py
  * zhang2011.py
  * genot2011.py
  * qian2011.py
  * qian2011_sqrt.py
  * kotani2017.py
  * zhang2009_rates.py
  * dabby2013_rates.py

In order to (re)produce results, you might want to adapt settings in the
library files.

# Jupyter-Lab notebooks
If you are interested in particular figures, the following notebook files have
been used to produce Figures 6 - 10. We recommend using a Python 3 kernel for
the analysis.

  * Reactionrates.ipynb (Badelt et al. submitted, Figure 6)
  * Kotani2017.ipynb (Badelt et al. submitted, Figure 7)
  * Yin2008.ipynb (Badelt et al. submitted, Figure 8)
  * Qian2011.ipynb (Badelt et al. submitted, Figure 9)
  * Many_systems.ipynb (Badelt et al. submitted, Figure 10)

[Badelt et al. (2020)]: <https://doi.org/10.1098/rsif.2019.0866>
