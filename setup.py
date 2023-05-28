#!/usr/bin/env python

'''
setup.py file for viQC
'''

from setuptools import setup, find_packages

version = open('VERSION').readline().strip()

setup(
    name                 = 'viQC',
    version              = version,
    description          = '''Visual & intuitive quality control for bottom-up proteomics experiments.''',
    long_description     = (''.join(open('README.md').readlines())),
    author               = 'Elizaveta Solovyeva',
    author_email         = 'pyteomics@googlegroups.com',
    url                  = 'http://hg.theorchromo.ru/viQC',
    install_requires     = ['pyteomics[DF]', 'seaborn', 'statsmodels', 'lxml', 'scikit-learn', 'scipy'],
    classifiers          = ['Intended Audience :: Science/Research',
                            'Programming Language :: Python :: 2.7',
                            'Programming Language :: Python :: 3',
                            'Topic :: Education',
                            'Topic :: Scientific/Engineering :: Bio-Informatics',
                            'Topic :: Scientific/Engineering :: Chemistry',
                            'Topic :: Scientific/Engineering :: Physics',
                            'Topic :: Software Development :: Libraries'],
    license              = 'License :: OSI Approved :: Apache Software License',
    packages             = find_packages(),
    entry_points         = {'console_scripts': ['viQC = viqc:main']}
)
