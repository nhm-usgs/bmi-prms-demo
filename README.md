# PRMS6-BMI Test Projects and python evaluation code

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nhm-usgs/bmi-test-binder/master?urlpath=https%3A%2F%2Fgithub.com%2Fnhm-usgs%2Fbmi-test-binder%2Fblob%2Fmaster%2Fnotebooks%2F03_Surf_Soil_GW_SF_Validation.ipynb)

## Contents
2. /prms - example PRMS projects for testing PRMS6 BMIs
3. /notebooks - notebooks for testing and demonstrating PRMS6 BMIs

### Create and activate conda environment
1. Create conda environment from environment.yml file
    - `conda env create -f environment.yml`
2. Activate environment
    - `conda activate prms6-bmi `


The pipestem project used in testing is found in /prms/pipestem. 

The Fortran git repositories for each BMI can be found at https://github.com/nhm-usgs

1. bmi-prms6-surface
2. bmi-prms6-soil
3. bmi-prms6-groundwater
4. bmi-prms6-streamflow

The python pymt gir repositories can be found here:
https://github.com/pymt-lab

1. pymt_prms_surface
2. pymt_prms_soil
3. pymt_prms_groundwater
3. pymt_prms_streamflow

Each BMI depends on the PRMS6 model and library found at the prms repository currently using the 6.0.0_dev_bmi branch here: https://github.com/nhm-usgs/prms

### Run on the CSDMS JupyterHub

Instead of installing locally,
the PRMS component Notebooks can be run on the CSDMS JupyterHub.
Follow these steps:

1. [Create an account](https://csdms.rc.colorado.edu/hub/signup) on the CSDMS JupyterHub, providing a username and password--they can be whatever you like
1. [Request authorization](https://github.com/csdms/help-desk/issues/new?assignees=mdpiper&labels=jupyterhub&template=new-csdms-jupyterhub-account.md&title=CSDMS+JupyterHub+account) for your new account through the CSDMS Help Desk--if you don't already have a GitHub account, you'll be asked to make one
1. Once approved, [run Jupyter Notebooks](https://csdms.rc.colorado.edu/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Fnhm-usgs%2Fbmi-prms-demo&urlpath=tree%2Fbmi-prms-demo%2Fnotebooks&branch=master)

Disclaimer
----------

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.

