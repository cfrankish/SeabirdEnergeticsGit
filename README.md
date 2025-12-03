# ⚡Workflow for estimating energetics of North Atlantic seabirds ⚡

This repository contains the workflow needed to calculate and map estimates of non-breeding season energy expenditure for six species of pelagic seabirds in the North Atlantic from immersion and positional (GLS) data collected via the SEATRACK program (https://seatrack.net/). It uses a probabilistic approach where uncertainty is included where possible and propagated through the workflow to the final results. 

This workflow stems from another project (SeabirdEnergetics) which was developed locally in R and is now being properly seperated from this initial work. 

## 📄Simple Schematic of workflow📄

<img width="2401" height="893" alt="image" src="https://github.com/user-attachments/assets/ffca257f-b7f3-4df6-b46b-6f82fcbf25d0" />

Cirles 1 & 2 form the basis for the work conducted in a paper that is in review (details to follow).

Circle 3 is under construction. 

## Pre-requisites 
Access to a high-performance computing cluster is required. This workflow has been built to run on SAGA (https://documentation.sigma2.no/hpc_machines/saga.html) but can be modified to run on other hpc systems, or eventually on a local machine.

## 🔧How to use🔧

### Install conda environment 
conda environment (i.e. all software needed to conduct this analysis) is contained within the `cbird.yml` file. To install it, enter the following code in console: 

<pre> sbatch Create_environment.sh </pre>

Make sure `Create_environment.sh` has been personalized first! (with account name etc.)

Environment can then be activated using the following code:

<pre>module load Miniconda3/22.11.1-1
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn******/PROJECTNAME/cbird </pre>

where nn****** is the project number and PROJECTNAME is your project name. 

### Make all directories

<pre>mkdir -p data/popdata_raw data/sst data/positionsIRMA code data/wetdry_raw</pre>

### Upload data that is needed to run workflow (update once SAGA is up and running again)

- **raw immersion data**: Unless you are part of the SEATRACK project group, you will need to upload ***
- **IRMA maps ([Fauchald et al. 2019](https://www.researchgate.net/profile/Arnaud-Tarroux-2/publication/334458632_Arctic-breeding_seabirds'_hotspots_in_space_and_time_-_A_methodological_framework_for_year-round_modelling_of_environmental_niche_and_abundance_using_light-logger_data/links/5d2c2ed292851cf44085033c/Arctic-breeding-seabirds-hotspots-in-space-and-time-A-methodological-framework-for-year-round-modelling-of-environmental-niche-and-abundance-using-light-logger-data.pdf))**: 
- **Population maps ([Fauchald et al. 2021](https://www.int-res.com/articles/meps_oa/m676p255.pdf))**:
- **SST & Sea-ice rasters ([Hersbach et al., 2023](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview))**: 

### Run through steps in worflow
The workflow is created via gwf (https://gwf.app/) and is scripted within `workflow.py`. This file contains all steps running from raw data extraction to the creation of figures. I recommend running through the steps one by one. Each step specifies expected input and output files, the R script which is called (contained in `code` folder), and requested walltime and memory usage. 

To query the status of a step, enter the following code:

<pre> gwf status s1_1_wetdry_* </pre>

To run one of the steps, enter the following code:

<pre> gwf run s1_1_wetdry_* </pre>

Note this specific example would submit six jobs as this step is parallelized. 

To query the status of your job on SAGA, enter the following code:

<pre> squeue -l -u ACOUNTNAME </pre>

## Key Contributors
- [Caitlin Frankish](https://github.com/cfrankish): Lead developer and maintainer
- [Mads Reinholdt Jensen](https://github.com/MadsRJ): Co-developer (input on using GWF) 


