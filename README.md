![Image](logo.jpg)

viQC - Visual&Intuitive quality control for bottom-up proteomics experiments
===========================================================================================


Requirements
------------
- Python 2.7 or 3.x
- matplotlib
- pyteomics
- pandas
- seaborn
- statsmodels

Installation via pip

``sudo apt-get install python-setuptools python-dev build-essential``

``sudo easy_install pip``

``sudo pip install -U lxml numpy matplotlib pyteomics pandas seaborn statsmodels ``

Before use
----------

Convert your raw files to mzML format. [MSconvert](<http://proteowizard.sourceforge.net/projects.html>) can be used for this purpose.
For the angle score calculation a list with PSMs is required (by default viQC works with [MPscore](<https://bitbucket.org/markmipt/mp-score>)/[Scavager](<https://bitbucket.org/markmipt/scavager>) results, but it can be changed by optional parameters).
For additional information about angle score see the [article](<https://www.sciencedirect.com/science/article/pii/S138738061730146X>)

How to use
----------

To install viQC, clone the repository and install with `pip` or distutils:

```
$ git clone https://bitbucket.org/lisavetasol/viqc
$ cd viqc
$ pip install .
```

or install directly from BitBucket:

```
$ pip install git+https://bitbucket.org/lisavetasol/viqc.git
```

To run viQC:

``viQC input_mzML``

For help on parameters:

``viQC --help``


Optional parameters
-------------------

``-o`` - path to save result, by default save in the same folder as an input file

``-start`` - delay time before sample elution; using for precursor intensity and injection time (MS/MS) calculation; 0 by default

``-stop`` - time of wash starting; using for precursor intensity and injection time (MS/MS) calculation. By default maximum analysis time

``-charge`` - max charge of precursor ions; by default, all charges are considered

ANGLE SCORE

``-refFile`` - mzML or mgf file for angle score calculation

``-refPSM`` - csv file with psm identifications for angle score calculation (should be based on mzML or mgf file specified in refFile parameter)

``-d`` - delimiter in csv file with psm identifications for angle score calculation; tab by default

``-cn``- column name with spectra names in csv file with psm identifications for angle score calculation; "spectrum" by default


Output
------
The output file has the name "InputFileName_viQC.png" and contains 8 quality metrics based on MS and MS/MS spectra characteristics.
The program doesn't require identifications, allowing fast and unbiased control of instrument performance.


Questions
---------
- [Create an issue](<https://bitbucket.org/lisavetasol/viqc/issues>) with BitBucket.
- Email us at pyteomics@googlegroups.com.
