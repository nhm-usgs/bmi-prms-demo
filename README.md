# PRMS6-BMI Test Projects and python evaluation code

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nhm-usgs/bmi-test-binder/master?urlpath=https%3A%2F%2Fgithub.com%2Fnhm-usgs%2Fbmi-test-binder%2Fblob%2Fmaster%2Fnotebooks%2F03_Surf_Soil_GW_SF_Validation.ipynb)

## Contents
1. /prms6bmi - Package prms6bmi contains some plotting tools for evaluating ouput of PRMS6 BMIs
2. /prms - example PRMS projects for testing PRMS6 BMIs
3. /notebooks - notebooks for testing and demonstrating PRMS6 BMIs

### Create python environment and install prms6bmi
1. Create conda environment from environment.yml file
    - `conda env create -f environment.yml`
2. Activate environment
3. Install prms6bmi into environment
    - ``` pip install prms6bmi ```

To build and install the prms6-bmi python package used in python evaluation:
1. Create a conda environment from the environment.yml file.
2. Activate the environment
3. cd to pkg/prms6-bmi
4. run "pip install ."

Or if you are developing the code inside of the prms6-bmi package -

- run "pip install -e ."

The pipestem project used in testing is found in /prms/pipestem. To test the bmi's copy the /pipestem folder and and save as:
1. pipestem_surface - control file: control_surface.simple1 
2. pipestem_soil - control file: control_soil.simple1 
3. pipestem_groundwater - control file: control_groundwater.simple1 
4. pipestem_streamflow - control file: control_streamflow.simple1 

The git repositories for each BMI can be found at https://github.com/nhm-usgs
1. bmi-prms6-surface
2. bmi-prms6-soil
3. bmi-prms6-groundwater
4. bmi-prms6-streamflow

Each BMI depends on the PRMS6 model and library found at the prms repository currently using the 6.0.0_dev_bmi branch


### Run on the CSDMS JupyterHub

Instead of installing locally,
the PRMS component Notebooks can be run on the CSDMS JupyterHub.
Follow these steps:

1. [Create an account](https://csdms.rc.colorado.edu/hub/signup) on the CSDMS JupyterHub, providing a username and password--they can be whatever you like
1. [Request authorization](https://github.com/csdms/help-desk/issues/new?assignees=mdpiper&labels=jupyterhub&template=new-csdms-jupyterhub-account.md&title=CSDMS+JupyterHub+account) for your new account through the CSDMS Help Desk--if you don't already have a GitHub account, you'll be asked to make one
1. Once approved, [run Jupyter Notebooks](https://csdms.rc.colorado.edu/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Fnhm-usgs%2Fbmi-prms-demo&urlpath=tree%2Fbmi-prms-demo%2Fnotebooks&branch=mdpiper/to-csdms-hub)

