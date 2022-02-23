#!/usr/bin/env python
from setuptools import setup

setup(name='ecalphisym',
      version='0.1',
      description='NanoAOD + coffea based implementation of the ECAL PhiSym calibration algorithm',
      #long_description = long_description,
      url='https://github.com/simonepigazzini/ecalphisym',
      author='Simone Pigazzini',
      author_email='simone.pigazzini@cern.ch',
      license='GPLv3',
      packages=[
          'ecalphisym'
      ],
      classifiers=[
          "Programming Language :: Python :: 3",
          "Operating System :: Linux",
      ],
      python_requires=">=3.10",
      install_requires=[
          'coffea>=0.7.12'
      ]
      )
