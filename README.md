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
- **IRMA maps ([Fauchald et al. 2019](https://nva-resource-storage-755923822223.s3.eu-west-1.amazonaws.com/12b6a384-2006-4a27-98d6-f217a8542f00?response-content-type=application%2Fpdf&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEEsaCWV1LXdlc3QtMSJHMEUCIQC%2BlAhvyAE5PPWK9UxG90bHdJq1iLn9UBzW%2B8alw3QdHQIgeRntmRrnI3ttaToNdKUJsmtZ8yuLSJ3oWyqjdz3y9oMq0wMIExADGgw3NTU5MjM4MjIyMjMiDGBWOr%2FLBggMVTrMJSqwA%2FunjdJnsWdBASWkN2GN2RsAHtkVTO6OIUDtgSmcKfuxsHispGYwFsfuYBOQy3B9IkSJs3OmVN%2BwTdOXYZvj9NWxCbc%2BoYdqNRd1NiaQFiat5nxC0kD09AUDNp4V89Kydb5rrcF6XYuP5LeQiJPa%2Bf1W5AnDlreeoiBCGFJ%2B50j9lYU5dsm53vpfpSSG2ZaXGM7Tl9lW2ywavh4bBt25sFEf0EiorSvlmsE92XKm8PSx83nCIpoxpBoIqc9g%2Fct5CBa9mCD6sR%2B%2FSUgkus3m9gaizfIQ24lbjlLuTaUquitvkzbBnh7tTZFYSFNbgzHTmBv3q1wETmP5HeZL7cZwdsbMJq6cPdAB8mTwpEYgj9jWeGilHlWD55wXH5RJqC87uFVMZ%2Bq7Iyt2%2BgYi0XhATQ%2FdoPVY0G%2B%2B0%2BycFKOKHCmyoXN5Hw597VRoch2uVhc4D3UQURzof%2BkZ%2FtP4grT%2F5bbsPWXe4%2BUIBCGMuCNEyrs7wZr%2FPDE%2F0Xp9j%2B9GratBOI85apstsEm6MX0EUY2NaHmmwC0GLPZpN12YNOq4JEJsClJG4S2vQcV4CjSso2ntizDU%2BrrJBjqeAQCE28s%2BlEW%2FNvpNjhLa%2Fnm3oh7IZKun%2FK8HrJg5cBlx1xFdXgbt437vGsL0%2FwXmitW1xlmWc7pJ33dhkUemItqkPn3eAg39F9S%2F28mZzK1N4pVxpgPZlnwoBg0wn6khR8nKZ2Pl6UowvulrdtP5ByoXofMMBrkm5aS9xATtqq7WbxNlUzOEFotFpFtdzVYZLBC9qqJayEvuVVmXcMIm&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20251202T103300Z&X-Amz-SignedHeaders=host&X-Amz-Credential=ASIA3AAESE2HQZRMXTST%2F20251202%2Feu-west-1%2Fs3%2Faws4_request&X-Amz-Expires=600&X-Amz-Signature=fe8573b6e6cc221e0cb96420cbb693803def6ddd7f35086d4124961b2781cdfb))**: 
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


