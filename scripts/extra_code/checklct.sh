#!/bin/bash

#SBATCH --account=nn11080k
#SBATCH --job-name=example
#SBATCH --partition=normal
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

# it is good to have the following lines in any bash script
set -o errexit  # make bash exit on any error
set -o nounset  # treat unset variables as errors

#module restore
#module load R/4.4.2-gfbf-2024a

Rscript timeUnderLCT.R > timeUnderLCT.Rout