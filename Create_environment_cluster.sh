#!/bin/bash
#SBATCH --job-name Create_environment   ## Name of the job
#SBATCH --output slurm-%j.out   ## Name of the output-script (%j will be replaced with job number)
#SBATCH --account xxxxxxxx   ## The billed account
#SBATCH --time=24:00:00   ## Walltime of the job
#SBATCH --mem-per-cpu=15G   ## Memory allocated to each task
#SBATCH --cpus-per-task=1   ## Number of CPUs allocated for each task

set -o errexit   ## Exit the script on any error
set -o nounset   ## Treat any unset variables as an error

module --quiet purge
module load Miniconda3/22.11.1-1
source ${EBROOTMINICONDA3}/bin/activate
my_conda_storage=/cluster/projects/xxxxxxxx/conda
export CONDA_PKGS_DIRS=${my_conda_storage}/package-cache
conda env create --prefix /cluster/projects/xxxxxxxx/conda/cbird --file cbird.yml
