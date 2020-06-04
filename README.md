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



