# Case study data 
This directory contains the data used for the publication [Badelt et al. (2020)]:

Badelt, Grun, Sarma, Wolfe, Shin, Winfree: 
*A domain-level DNA strand displacement reaction enumerator allowing arbitrary non-pseudoknotted secondary structures*

Note: The notebooks have been produced with peppercorn version 1.0. 
If you are experiencing problems at later versions of this software,
please install peppercornenumerator==1.0.

You must make a directory called **tmp/** where temporary files are stored. Be
aware that the temporary files remain in this directory for further analysis.

# Jupyter notebooks
If you are interested in particular figures, the following notebook files have
been used to produce Figures 6 - 11.

  * Condensation.ipynb ([Badelt et al. (2020)], [Figure 6])
  * Reactionrates.ipynb ([Badelt et al. (2020)], [Figure 7])
  * Kotani2017.ipynb ([Badelt et al. (2020)], [Figure 8])
  * Yin2008.ipynb ([Badelt et al. (2020)], [Figure 9])
  * Qian2011.ipynb ([Badelt et al. (2020)], [Figure 10])
  * Many_systems.ipynb ([Badelt et al. (2020)], [Figure 11])

# Commandline verison
Then execute the python script:

    $ ./figure_analysis.py

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

You might want to adapt settings in those files.

[Badelt et al. (2020)]: <https://doi.org/10.1098/rsif.2019.0866>
[Figure 6]: <https://royalsocietypublishing.org/doi/10.1098/rsif.2019.0866#RSIF20190866F6>
[Figure 7]: <https://royalsocietypublishing.org/doi/10.1098/rsif.2019.0866#RSIF20190866F7>
[Figure 8]: <https://royalsocietypublishing.org/doi/10.1098/rsif.2019.0866#RSIF20190866F8>
[Figure 9]: <https://royalsocietypublishing.org/doi/10.1098/rsif.2019.0866#RSIF20190866F9>
[Figure 10]: <https://royalsocietypublishing.org/doi/10.1098/rsif.2019.0866#RSIF20190866F10>
[Figure 11]: <https://royalsocietypublishing.org/doi/10.1098/rsif.2019.0866#RSIF20190866F11>
