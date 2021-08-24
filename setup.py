#!/usr/bin/env python

from distutils.core import setup

setup(name='wasserstein_flows',
      version='0.0.1',
      description='Algorithms for Wasserstein metric',
      author='Michał Startek',
      author_email='michal.startek@mimuw.edu.pl',
      url='https://github.com/michalsta/wasserstein_flows',
      packages=['wasserstein_flows'],
      depends=['cppyy']
     )

