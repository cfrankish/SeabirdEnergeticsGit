"""
A gwf-workflow file for running Caitlins seabird data on Saga.
Before running the script, run the following command:
mkdir -p data/birddata_raw data/popdata_raw data/popdata_ind data/sst data/positionsIRMA data/metadata scripts data/wetdry_raw data/birddata_ind
1) Place all bird ".csv" files into the "data/birddata_raw/" folder
2) Place all raster files in the "data/popdata_raw/" folder
3) Place all sst data files in the "data/sst/" folder - has to be data/sst/current/HadGEM3-GC31-MM
4) Place all IRMA files in the "data/positionsIRMA" folder
5) Place all rds metadata files in the "data/metadata" folder
6) Place all scripts in the "scripts/" folder
"""

from gwf import Workflow
import os
import csv
from glob import glob

gwf = Workflow(defaults={"account": "nn11080k"}) # Remember to change account number as needed!

############################################################################
### Step 1_1: Extract wet-dry data from SEATRACK database ##################
############################################################################

import csv
import time
import re

# Get today's date in YYYY-MM-DD format
today = time.strftime("%Y-%m-%d", time.localtime())

# Create list of species to iterate through
species_list = ["Littleauk", "Northernfulmar", "Atlanticpuffin", "Blackleggedkittiwake", "Commonguillemot", "Brunnichsguillemot"]

# Generate a GWF target per species
for sp in species_list:
    output_file = f"data/{sp}_wetdry_{today}.csv"
    
    # Create GWF target for each species
    gwf.target(
        name=f"s1_1_wetdry_{sp}",
        inputs=[],
        outputs=[output_file],
        queue="bigmem",
        cores=1,
        memory="50G",
        walltime="12:00:00"
    ) << f"""
    Rscript scripts/s1_1_extract_wetdry.R {sp} {today} {output_file}
    """
    
############################################################################
### Step 1_2: Clean data & split into individual bird files ################
############################################################################

# Directory containing the CSV files
csv_dir_wetdry = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/data/wetdry_raw"

# Find all CSV files in the directory
csv_files_wetdry = glob(os.path.join(csv_dir_wetdry, "*.csv"))

for csv_file in csv_files_wetdry:
    
    # Make a target name specific to each file
    filename = csv_file.split("/")[-1]  # get the file name only
    base_name = filename.split("_")[0]  # take text before first underscore
    clean_base_name = base_name.replace("-", "") # remove hyphen
    
    # Make an output file
    output_files1 = f"data/birddata_raw/{clean_base_name}_IDS_test.txt"
    
    # Set-up the target
    gwf.target(
        name=f"s1_3_process_wetdry_{clean_base_name}",
        inputs=[],
        outputs=[output_files1],
        queue="bigmem",
        cores=1,
        memory="50G",
        walltime="30:00:00"
    ) << f"""
Rscript scripts/s1_3_merge_wetdry_lox_env.R {clean_base_name} {output_files1} 
"""    

#############################################################################################
### Step 2_1: Calculating activity & energy budgets for individual CSV files with Rscript ###
#############################################################################################

# Directory containing the CSV files
csv_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/data/birddata_ind/*/"

# Find all CSV files in the directory
csv_files = glob(os.path.join(csv_dir, "*.csv"))

# Process each CSV file
for csv_file in csv_files:
    bird_species = os.path.basename(os.path.dirname(csv_file))
    bird_name = os.path.basename(csv_file).replace('.csv', '')
    
    input_files = [csv_file]
    output_files1 = f"tmp/{bird_name}_energyDay.csv"
    
    gwf.target(
        name=f"s2_1_calculate_budgets_{bird_species}_{bird_name}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="04:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s2_1_calculatebudgets.R {csv_file} {output_files1} 
    """.format(csv_dir=csv_dir, csv_files=csv_files)

##############################################################
### Step 2_2: Sensitivity analysis - part 1 (running it)  ####
##############################################################

# Directory containing the CSV files
csv_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/data/birddata_ind/*/"

# Find all CSV files in the directory
csv_files = glob(os.path.join(csv_dir, "*.csv"))

# Process each CSV file
for csv_file in csv_files:
    bird_species = os.path.basename(os.path.dirname(csv_file))
    bird_name = os.path.basename(csv_file).replace('.csv', '')
    
    input_files = [csv_file]
    output_files1 = f"tmp2/{bird_name}_energyDay.csv"
    
    gwf.target(
        name=f"s2_2_run_sensitivity_{bird_species}_{bird_name}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="04:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s2_2_sensitivity_analysis_part1.R {csv_file} {output_files1} 
    """.format(csv_dir=csv_dir, csv_files=csv_files)    
    
##############################################################
### Step 2_3: Sensitivity analysis - part 2 (running it)  ####
##############################################################

input_files = []  # Add actual input files if necessary
output_files1 = f"./results/tables/supplementary/table0_sensitivityResults.csv"

gwf.target(
    name="s2_3_run_sensitivity_2",
    inputs=input_files,
    outputs=[output_files1],
    cores=1,
    memory="7G",
    walltime="04:00:00"
) << f"""
Rscript scripts/s2_3_sensitivity_analysis_part2.R {output_files1} 
"""

##############################################################################
### Step 3_1: Conduct analyses for supplementary part 1 - min iteration no ###
##############################################################################
 
# Determine start & end of study period # 
startDate = "2021-09-01"
endDate = "2022-04-30"
 
# Determine names of input & output files 
input_files = []  # Add actual input files if necessary
output_files1 = f"./results/tables/main/table1_idcatalogue.csv"
output_files2 = f"./results/tables/supplementary/table1_miniteration.csv"

gwf.target(
    name="s3_1_analysis_min_iteration_number",
    inputs=input_files,
    outputs=[output_files1, output_files2],
    queue="bigmem",
    cores=1,
    memory="40G",
    walltime="06:00:00"
) << f"""
Rscript scripts/s3_1_supanalysis_min_iteration_number.R {startDate} {endDate} {output_files1} {output_files2}
"""

###########################################################
### Step 3_2: Conduct analyses for supplementary part 2 ###
###########################################################
 
# Determine names of input & output files # 
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"  # Add actual input files if necessary
output_files1 = f"./results/tables/supplementary/table2_timeRandom.csv"
output_files2 = f"./results/tables/supplementary/table3_timeDark.csv"

gwf.target(
    name="s3_2_analysis_check_bout_lengths",
    inputs=input_files1,
    outputs=[output_files1, output_files2],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="20:00:00"
) << f"""
Rscript scripts/s3_2_supanalysis_check_bout_lengths.R {input_files1} {output_files1} {output_files2}
"""

##########################################################################################
### Step 4_1: Analysis MS prt1 - Maps + temporal patterns in activity & energy budgets ###
##########################################################################################
 
# Determine study period 
startDate = "2021-09-15"
endDate = "2022-04-15"

# Determine input & output files
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"  # Add actual input files if necessary 
output_files1 = f"./results/tables/main/table2_budgets_species.csv" # Table with species-specific activity and energy budgets
output_files2 = f"./results/tables/main/table3_budgets_population.csv" # Table with population-specific activity and energy budgets
output_files3 = f"./results/tables/main/table4_budgets_individual.csv" # Table with population-specific activity and energy budgets

gwf.target(
    name="s4_1_analysis_plot_maps_and_budgets",
    inputs=input_files1,
    outputs=[output_files1, output_files2, output_files3],
    queue="bigmem",
    cores=1,
    memory="60G",
    walltime="72:00:00"
) << f"""
Rscript scripts/s4_1_mainanalysis_plot_maps_and_budgets.R {input_files1} {startDate} {endDate} {output_files1} {output_files2} {output_files3}
"""

#########################################################################################
### Step 4_2: Analysis MS prt2 - Calculating metrics needed to conduct statistics #######
#########################################################################################
 
# Determine study period 
startDate = "2021-09-15"
endDate = "2022-04-15"

# Determine input & output files 
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"  # Add actual input files if necessary 
output_files1 = f"./results/tables/main/table5_migratory_distance.csv"    # Table with individual migratory distance

gwf.target(
    name="s4_2_analysis_calc_metrics",
    inputs=input_files1,
    outputs=[output_files1],
    queue="bigmem",
    cores=1,
    memory="80G",
    walltime="72:00:00"
) << f"""
Rscript scripts/s4_2_mainanalysis_calc_metrics.R {input_files1} {startDate} {endDate} {output_files1} 
"""

#######################################################################
### Step 4_3: Analysis MS prt3 - All statistics presented in MS #######
#######################################################################
 
# Determine study period 
startDate = "2021-09-15"
endDate = "2022-04-15"

# Determine intput and output files
output_files1 = f"./results/tables/main/table6_totalNBCosts.csv"
output_files2 = f"./results/tables/main/table7_stats_WEE_vs_TEE.csv"
output_files3 = f"./results/tables/main/table8_stats_WEE_vs_pred.csv"
output_files4 = f"./results/tables/main/table9_stats_TEE_vs_pred.csv"

gwf.target(
    name="s4_3_anaysis_weekly_vs_tot_energy",
    inputs=[], 
    outputs=[output_files1, output_files2, output_files3, output_files4],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="72:00:00"
) << f"""
Rscript scripts/s4_3_mainanalysis_statistics.R {startDate} {endDate} {output_files1} {output_files2} {output_files3} {output_files4} 
"""

#######################################################################
### Step 4_4: Supplementary analysis (extra) - effect sizes ###########
#######################################################################
 
# Determine study period 
startDate = "2021-09-15"
endDate = "2022-04-15"

# Determine intput and output files
input_files1 = f"./results/tables/main/table6_totalNBCosts.csv"
input_files2 = f"./results/tables/main/table7_stats_WEE_vs_TEE.csv"
input_files3 = f"./results/tables/main/table8_stats_WEE_vs_pred.csv"
input_files4 = f"./results/tables/main/table9_stats_TEE_vs_pred.csv"

gwf.target(
    name="s4_4_effectSize_vs_reps",
    inputs=[input_files1, input_files2, input_files3, input_files4], 
    outputs=[],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="72:00:00"
) << f"""
Rscript scripts/s4_3_mainanalysis_statistics.R {startDate} {endDate} {input_files1} {input_files2} {input_files3} {input_files4} 
"""