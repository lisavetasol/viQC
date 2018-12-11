viQC - Visual&Intuitive quality control for bottom-up proteomics experiments
===========================================================================================

![Image](toc_fig.png)

Requirements
------------
- Python 2.7
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
For angle score calculation the list with PSMs is required (by default viQC works with [MPscore](<https://bitbucket.org/markmipt/mp-score>)/[Scavager](<https://bitbucket.org/markmipt/scavager>) results, but it can be changed by optional parameters). 
For additional information about angle score see the [article](<https://www.sciencedirect.com/science/article/pii/S138738061730146X>)

How to use
----------
    
``python2 viQC.py input_mzML_file optional_parameters ``


Optional parameters 
----------

``-o`` - path to save result, by default save in the same folder as an input file

``-start`` - delay time before sample actually comes to mass spec; using for precursor intensity and injection time (MS/MS) calculation; 0 by default

``-stop`` - time of wash starting; using for precursor intensity and injection time (MS/MS) calculation. By default maximum analysis time

``-charge`` - max charge of precursor ions; 4 by default

ANGLE SCORE 

``-refmzML`` - mzML file for angle score calculating

``-refPSM`` - csv file with psm identifications for angle score calculating. It should be noted that by default identification proccess is assumed to be strarted from mgf file, if mzML file was used parameter "-f" should be changed to 0 

``-d`` - delimiter in csv file with psm identifications for angle score calculating; tab by default

``-cn``- column name with spectra names in csv file with psm identifications for angle score calculating; "spectrum" by default


Output 
------
The output file has the name "viQC_results_InputFileName.png" and contains 8 quality metrics based on MS and MS/MS spectra characteristics. 
The program doesn't require identifications that allow fast and unbiasied quality control of instrument performance.  




Questions
---------
- [Create an issue](<https://bitbucket.org/levitsky/fractionoptimizer/issues>) with BitBucket.
- Email us at pyteomics@googlegroups.com or biolccc@googlegroups.com.
