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
        walltime="80:00:00"
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
#output_files2 = f"./results/tables/main/table6_migratory_distance_sst_anomaly.csv" # Table with individual migratory distance
output_files3 = f"./results/tables/main/table7_species_mean_deviance.csv" # Table with species deviance in WEE, weighted for varying sample size
output_files4 = f"./results/tables/main/table8_population_mean_deviance.csv" # Table with population deviance in WEE

gwf.target(
    name="s4_2_analysis_calc_metrics",
    inputs=input_files1,
    outputs=[output_files1, output_files3, output_files4],
    queue="bigmem",
    cores=1,
    memory="80G",
    walltime="72:00:00"
) << f"""
Rscript scripts/s4_2_mainanalysis_calc_metrics.R {input_files1} {startDate} {endDate} {output_files1} {output_files3} {output_files4} 
"""

#######################################################################
### Step 4_3: Analysis MS prt3 - All statistics presented in MS #######
#######################################################################
 
# Determine study period 
startDate = "2021-09-15"
endDate = "2022-04-15"

# Determine intput and output files
output_files1 = f"./results/tables/main/table9_totalNBCosts.csv"
output_files2 = f"./results/tables/main/table10_stats3_dredge_nbCosts.csv"
output_files3 = f"./results/tables/main/table11_cluster_pop.csv"
output_files4 = f"./results/tables/main/table12_stats4_dredge_nbCosts_species.csv"
output_files5 = f"./results/tables/main/table13_stats5_dredge_nbCosts_deviance.csv"
output_files6 = f"./results/tables/main/table14_stats6_dredge_nbCosts_cov.csv"

gwf.target(
    name="s4_3_anaysis_weekly_vs_tot_energy",
    inputs=[], 
    outputs=[output_files1, output_files2, output_files3, output_files4, output_files5, output_files6],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="72:00:00"
) << f"""
Rscript scripts/s4_3_mainanalysis_statistics.R {startDate} {endDate} {output_files1} {output_files2} {output_files3} {output_files4} {output_files5} {output_files6}
"""

#####################################################################################
### Step 5_1: Determine distribution of activity budgets per pop (non spatial) ######
#####################################################################################

# Directory containing the CSV files
nc_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/data/popdata_raw"

# Find all CSV files in the directory
nc_files = glob(os.path.join(nc_dir, "*.nc"))

# Define secondary input file
input_file2 = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/main/table1_idcatalogue.csv"

# Process each CSV file
for nc_file in nc_files:
    
    # Get just the file name
    filename = os.path.basename(nc_file)  # "file1.nc"
    
    # Remove the extension
    name_no_ext = filename.replace(".nc", "")

    # Split by underscore
    parts = name_no_ext.split("_")

    # Extract the species (4th and 5th parts)
    bird_species = f"{parts[3]}_{parts[4]}"
    input_files = [nc_file, input_file2]
    
    # Determine intput and output files
    output_files1 = f"./results/tables/supplementary/{bird_species}_actBudget.csv"
    
    gwf.target(
        name=f"s5_1_1_prep_activity_temporal_{bird_species}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="04:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s5_1_1_prep_activity_distribution_temporal.R {nc_file} {input_file2} {output_files1}
    """   
    
#########################################################################################
### Step 5_2: Map energy expenditure (assume same behaviour regardless of location ######
#########################################################################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*monthly.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp4", exist_ok=True)

# Process each CSV file
for budget_file in budget_files:
    
    # Get just the file name
    filename = os.path.basename(budget_file)  
    
    # Remove the extension
    name_no_ext = filename.replace(".csv", "")
    
    # Remove special characters
    name_no_ext_clean = (name_no_ext
               .replace("(", "")
               .replace(")", "")
               .replace(",", "")
               .replace("'", "")
               .replace("–", "")
               .replace("ý", "")
               .replace("ð", "")
               .replace("ó", "")
               .replace("á", "")
               .replace(":", "")
               .replace("í", "")
               .replace("`", "")
               .replace("ú", "")
               .replace("é", ""))
    
    
    # Remove special characters
    budget_clean = (budget_file
               .replace("(", "")
               .replace(")", "")
               .replace("`", "")
               .replace("'", ""))
    
    # Define the input file
    input_files = [budget_file]
    
    # Determine intput and output files
    output_files1 = f"./tmp4/{name_no_ext_clean}_energyMap.csv"
    
    gwf.target(
        name=f"s5_2_1_map_energy_{name_no_ext_clean}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="04:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s5_2_1_map_energy_v1.R {budget_file} {output_files1}
    """       
    
##################################
### Step 5_3: Map activity  ######
##################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*daily.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp4", exist_ok=True)

# Process each CSV file
for budget_file in budget_files:
    
    # Get just the file name
    filename = os.path.basename(budget_file)  
    
    # Remove the extension
    name_no_ext = filename.replace(".csv", "")
    
    # Remove special characters
    name_no_ext_clean = (name_no_ext
               .replace("(", "")
               .replace(")", "")
               .replace(",", "")
               .replace("'", "")
               .replace("–", "")
               .replace("ý", "")
               .replace("ð", "")
               .replace("ó", "")
               .replace("á", "")
               .replace(":", "")
               .replace("í", "")
               .replace("`", "")
               .replace("ú", "")
               .replace("é", ""))
    
    
    # Remove special characters
    budget_clean = (budget_file
               .replace("(", "")
               .replace(")", "")
               .replace("`", "")
               .replace("'", ""))
    
    # Define the input file
    input_files = [budget_file]
    
    # Determine intput and output files
    output_files1 = f"./tmp4/{name_no_ext_clean}_activityMap.csv"
    
    gwf.target(
        name=f"s5_3_map_activity_monthly_{name_no_ext_clean}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="04:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s5_3_map_activity.R {budget_file} {output_files1}
    """    
    
##################################
### Step 5_3: Map activity - yearly  ######
##################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*daily.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp4", exist_ok=True)

# Process each CSV file
for budget_file in budget_files:
    
    # Get just the file name
    filename = os.path.basename(budget_file)  
    
    # Remove the extension
    name_no_ext = filename.replace(".csv", "")
    
    # Remove special characters
    name_no_ext_clean = (name_no_ext
               .replace("(", "")
               .replace(")", "")
               .replace(",", "")
               .replace("'", "")
               .replace("–", "")
               .replace("ý", "")
               .replace("ð", "")
               .replace("ó", "")
               .replace("á", "")
               .replace(":", "")
               .replace("í", "")
               .replace("`", "")
               .replace("ú", "")
               .replace("é", ""))
    
    
    # Remove special characters
    budget_clean = (budget_file
               .replace("(", "")
               .replace(")", "")
               .replace("`", "")
               .replace("'", ""))
    
    # Define the input file
    input_files = [budget_file]
    
    # Determine intput and output files
    output_files1 = f"./tmp4/{name_no_ext_clean}_activityMap_yearly.csv"
    
    gwf.target(
        name=f"s5_3_map_activity_yearly_{name_no_ext_clean}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="04:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s5_3_map_activity_winter.R {budget_file} {output_files1}
    """    

##############################################################
### Step 5_4: Aggregate activity - spatially mapped  #########
##############################################################

import re

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp4"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*.csv"))

# Extract species names
results = [re.findall(r"tmp4/([^_]+)_", str(p))[0] for p in budget_files]
unique_names = set(results)

# Process each CSV file
for name in unique_names:
    # Define output files
    output_file1 = f"./results/tables/{name}_activityMap_WeightedMean.csv"
    output_file2 = f"./results/tables/{name}_activityMap_WeightedSD.csv"
    
    gwf.target(
        name=f"s5_4_aggregate_maps_{name}",
        inputs=[],
        outputs=[output_file1, output_file2],
        cores=1,
        memory="8G",
        walltime="04:00:00",
    ) << f"""
    Rscript scripts/s5_4_aggregate_maps_activity_spatial.R {name} {output_file1} {output_file2}
    """     

##############################################################
### Step 5_4: Aggregate energy - non-spatial activity  #######
##############################################################

import re

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp4"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*.csv"))

# Extract species names
results = [re.findall(r"tmp4/([^_]+)_", str(p))[0] for p in budget_files]
unique_names = set(results)

# Process each CSV file
for name in unique_names:
    # Define output files
    output_file1 = f"./results/tables/{name}_energyMap_monthly.csv"
    
    gwf.target(
        name=f"s5_4_aggregate_maps_energy_{name}",
        inputs=[],
        outputs=[output_file1],
        queue="bigmem",
        cores=1,
        memory="20G",
        walltime="10:00:00",
    ) << f"""
    Rscript scripts/s5_3_1_aggregate_maps_energy_species_v1.R {name} {output_file1}
    """         
    
######################################################################
### Step 5_5: Aggregate energy - non-spatial activity - total  #######
######################################################################

# Define output files
output_file1 = f"./results/tables/allBirds_energyMap_monthly.csv"

gwf.target(
    name="s5_5_aggregate_maps_energy_all",
    inputs=[],
    outputs=[output_file1],
    queue="bigmem",
    cores=1,
    memory="20G",
    walltime="10:00:00"
) << f"""
Rscript scripts/s5_5_aggregate_maps_energy_all.R {output_file1}
"""
   

########################################################################
### Step XX: For Lila - Calculate activity budgets every day of year ###
########################################################################

import csv
import time
import re

# Get today's date in YYYY-MM-DD format
today = time.strftime("%Y-%m-%d", time.localtime())

# Function to clean species names
def clean_species_name(species_name):
    # Replace 'ü' with 'u'
    cleaned_name = species_name.replace("ü", "u")
    
    # Remove spaces, dashes, apostrophes
    cleaned_name = cleaned_name.replace(" ", "")  # Remove spaces
    cleaned_name = cleaned_name.replace("-", "")  # Remove dashes
    cleaned_name = cleaned_name.replace("'", "")  # Remove apostrophes
    
    # Remove non-alphanumeric characters (except for the allowed ones like 'u')
    cleaned_name = re.sub(r'[^a-zA-Z0-9]', '', cleaned_name)
    
    return cleaned_name

# Read and deduplicate species names
species_set = set()

# Open the CSV file and read the species column
with open("./results/tables/main/table1_idcatalogue.csv", "r") as f:
    reader = csv.DictReader(f)  # Use DictReader to read the file as a dictionary
    for row in reader:
        species = row['species'].strip()  # Extract species from the 'species' column
        if species:
            cleaned_species = clean_species_name(species)  # Clean the species name
            species_set.add(cleaned_species)

# Convert the set to a sorted list
species_list = sorted(species_set)

# Debugging: Check cleaned species names
#print("Cleaned species list:", species_list)

# Generate a GWF target per species
for sp in species_list:
    output_file = f"results/dailyactivity_{sp}_{today}.csv"
    
    # Create GWF target for each species
    gwf.target(
        name=f"forLila_{sp}",
        inputs=[],
        outputs=[output_file],
        queue="bigmem",
        cores=1,
        memory="50G",
        walltime="12:00:00"
    ) << f"""
    mkdir -p results
    Rscript scripts/forLila.R {sp} {today} {output_file}
    """
    
#######################################################################
### Step XX: Extra analysis for Seb - make Svalbard bird energetic ###
#######################################################################
 
# Determine intput and output files
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"
output_files1 = f"./results/tables/main/seb_svalbard_energy.csv"

gwf.target(
    name="forSeb",
    inputs=[input_files1], 
    outputs=[output_files1],
    cores=1,
    memory="10G",
    walltime="72:00:00"
) << f"""
Rscript scripts/s2_4_energy_seb.R {input_files1} {output_files1} 
"""