# ⚡Workflow for estimating activity budgets and energy expenditure for North Atlantic seabirds ⚡

This repository contains a workflow which can be used to calculate non-breeding season activity budgets and energy expenditure for six species of pelagic seabirds in the North Atlantic from immersion and positional (GLS) data collected via the SEATRACK program (https://seatrack.net/). It uses a probabilistic approach where uncertainty and variation is included where possible and propagated through the workflow to the final results. It also performs all statistical analysis contained in xxxx (paper under review). 

## 📄Simple Schematic of workflow📄

<img width="2353" height="1340" alt="image" src="https://github.com/user-attachments/assets/267f6d6c-9d27-44cf-bfbe-45fe08e16538" />

## Pre-requisites 

- **Access to a high-performance computing cluster**:  
  - This workflow has been built to run on SAGA (https://documentation.sigma2.no/hpc_machines/saga.html), but can be modified to run on other HPC systems or eventually on a local machine.

- **Tracking data** (Can be requested via the
    [SEATRACK data request form](https://seatrack.net/data/data-request-form/))
  - **Immersion/activity data from GLS:** 

    Has the following columns:

    | session_id | individ_id | date_time | conductivity | std_conductivity |
    |------------|------------|-----------|--------------|------------------|

  - **IRMA data** ([Fauchald et al. 2019](https://www.researchgate.net/profile/Arnaud-Tarroux-2/publication/334458632_Arctic-breeding_seabirds'_hotspots_in_space_and_time_-_A_methodological_framework_for_year-round_modelling_of_environmental_niche_and_abundance_using_light-logger_data/links/5d2c2ed292851cf44085033c/Arctic-breeding-seabirds-hotspots-in-space-and-time-A-methodological-framework-for-year-round-modelling-of-environmental-niche-and-abundance-using-light-logger-data.pdf)):

    Has the following columns:

    | individ_id | session_id | age_class | timestamp | loc_type | lon | lat |
    |------------|------------|-----------|-----------|----------|-----|-----|
    

## 🔧How to use : local steps🔧

### Create local project with the following folders ###

```text
project/
├── code/
├── data/
│   ├── wetdry_database/
│   ├── positionsIRMA/
│   ├── ice/
│   ├── sst/
│   └── temp/
```

where:

- `code/` – scripts (in this case `s1_1_merge_tracking_data_local.R` and `s1_2_download_env_data_local.R`)
- `data/` – input datasets used by the project.
  - `wetdry_database/` – activity data by species with following format `wetdry_database_Black-legged kittiwake_2026-08-06.rds`
  - `wetdry_raw/` 
  - `positionsIRMA/` – IRMA position data.
  - `ice/` 
  - `sst/` 
  - `temp/` 

### Run local code (steps s1_1 and s1_2)  
Step s1_1: processes the immersion data and merges it with the IRMA data.

Step s1_2: Extracts all unique dates and uses these to download environmental (SST, air temperature and sea-ice concentration rasters) from (https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview))

## 🔧How to use : cluster steps🔧

### Make all directories on cluster
<pre>mkdir -p data/wetdry_raw data/sst data/ice data/sst data/temp data/positionsIRMA data/birddata_raw/ data/birddata_ind/ scripts results results/tempPlots</pre>

### Upload local data to cluster
1. Upload all environmental data to `data/sst`, `data/ice` and `data/temp`

2. Upload all activity data produced using `s1_1` to `data/wetdry_raw`

3. Upload all IRMA positions to `data/positionsIRMA`

4. Upload all scripts to scripts folder

### Install conda environment 
conda environment (i.e. all software needed to conduct this analysis) is contained within the `cbird.yml` file. To install it, enter the following code in console: 

<pre> sbatch Create_environment_cluster.sh </pre>

Make sure `Create_environment_cluster.sh` has been personalized first! (with account name etc.)

Environment can then be activated using the following code:

<pre>module load Miniconda3/22.11.1-1
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn******/PROJECTNAME/cbird </pre>

where nn****** is the project number and PROJECTNAME is your project name. 

 ### Run through steps in worflow
The workflow is created via gwf (https://gwf.app/) and is scripted within `workflow.py`. This file contains all steps running from raw data extraction to the creation of figures. I recommend running through the steps one by one. Each step specifies expected input and output files, the R script which is called (contained in `code` folder), and requested walltime and memory usage. 

To query the status of a step, enter the following code:

<pre> gwf status s1_3_merge_wetdry_lox* </pre>

To run one of the steps, enter the following code:

<pre> gwf run s1_3_merge_wetdry_lox* </pre>

Note this specific example would submit six jobs as this step is parallelized. 

To query the status of your job on SAGA, enter the following code:

<pre> squeue -l -u ACOUNTNAME </pre>

## Key Contributors
- [Caitlin Frankish](https://github.com/cfrankish): Lead developer and maintainer
- [Mads Reinholdt Jensen](https://github.com/MadsRJ): Co-developer (input on using GWF) 


