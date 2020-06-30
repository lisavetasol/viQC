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
- scipy
- sklearn

Installation via pip

``sudo apt-get install python-setuptools python-dev build-essential``

``sudo easy_install pip``

``sudo pip install -U lxml numpy matplotlib pyteomics pandas seaborn statsmodels ``

Before use
----------

Convert your raw files to mzML format. [MSconvert](<http://proteowizard.sourceforge.net/projects.html>) can be used for this purpose.
For the angle score calculation a list with PSMs is required (by default viQC works with [MPscore](<https://bitbucket.org/markmipt/mp-score>)/[Scavager](<https://bitbucket.org/markmipt/scavager>) results, but it can be changed by optional parameters). 

For the additional information about viQC program and user cases check the [article](<https://link.springer.com/article/10.1134/S1061934819140119>); Russian speakers can enjoy the mother-tongue [version](<http://mass-spektrometria.ru/viqc-%D0%B1%D1%8B%D1%81%D1%82%D1%80%D1%8B%D0%B9-%D0%B8-%D0%BD%D0%B0%D0%B3%D0%BB%D1%8F%D0%B4%D0%BD%D1%8B%D0%B9-%D0%BA%D0%BE%D0%BD%D1%82%D1%80%D0%BE%D0%BB%D1%8C-%D0%BA%D0%B0%D1%87%D0%B5%D1%81%D1%82/>).

For the additional information about angle score see the [article](<https://www.sciencedirect.com/science/article/pii/S138738061730146X>)

How to use
----------

To install viQC, clone the repository and install with `pip` or distutils:

```
$ git clone https://github.com/lisavetasol/viQC
$ cd viQC
$ pip install .
```


To run viQC:

``viQC input_mzML(s)``

In case of multiple files, the program builds summarizing graphs as well as separate pictures for each analysis. For easier interpreting results start/stop parameters should be specified, otherwise, they will be calculated for each file separately and corresponding graphs (precursor intensity and MS/MS injection time) may not be comparable between files. 

For help on parameters:

``viQC --help``


Optional parameters
-------------------

``-o`` - path to save result, by default save in the same folder as an input file

``-start`` - delay time before sample elution; using for precursor intensity and injection time (MS/MS) calculation; 0 by default

``-stop`` - time of wash starting; using for precursor intensity and injection time (MS/MS) calculation. By default maximum analysis time

``-charge`` - max charge of precursor ions; by default, all charges are considered

``-pic`` - the output figure format (png or svg for vector graphic)

ANGLE SCORE

``-refFile`` - mzML or mgf file for angle score calculation

``-refPSM`` - csv file with psm identifications for angle score calculation (should be based on mzML or mgf file specified in refFile parameter)

``-d`` - delimiter in csv file with psm identifications for angle score calculation; tab by default

``-cn``- column name with spectra names in csv file with psm identifications for angle score calculation; "spectrum" by default


Output
------
The output file has the name "InputFileName_viQC.png" and contains 8 quality metrics based on MS and MS/MS spectra characteristics.
The program doesn't require identifications, allowing fast and unbiased control of instrument performance.

Citing viQC
------
Solovyeva E.M. et al. viQC: Visual and Intuitive Quality Control for Mass Spectrometry-Based Proteome Analysis. J Anal Chem 74, 1363–1370 (2019). https://doi.org/10.1134/S1061934819140119

Solovyeva E.M. et al. Semi-supervised quality control method for proteome analyses based on tandem mass spectrometry. International Journal Of Mass Spectrometry 427, 59-64 (2018).
https://doi.org/10.1016/j.ijms.2017.09.008

Questions
---------
- [Create an issue](<https://github.com/lisavetasol/viQC/issues>).
- Email us at pyteomics@googlegroups.com.
