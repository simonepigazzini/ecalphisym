# ECAL PhiSym calibration package

This package provides tools and examples to perform phisym/eflow calibration of the CMS
ECAL detector using ZeroBias events. This version works with reconstructed events 
processed with CMSSW and saved in NanoAOD format. Please read the compiled [documentation](https://spigazzi.web.cern.ch/docs/ecalphisym/)

## Installation
This package can be installed from github using pip:

```
pip install git+https://github.com/simonepigazzini/ecalphisym.git
```

An easy way to install and test is to create a dedicated conda environment

```
conda create -n phisym python==3.10
conda activate phisym
pip install git+https://github.com/simonepigazzini/ecalphisym.git
```

