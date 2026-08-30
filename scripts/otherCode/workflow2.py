#########################################################
### Step 4: Conduct analyses for supplementary part 2 ###
#########################################################
  
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"  # Add actual input files if necessary
output_files1 = f"./results/tables/supplementary/table2_timeRandom.csv"
output_files2 = f"./results/tables/supplementary/table3_timeDark.csv"

gwf.target(
    name="s2_2_analysis_check_bout_lengths",
    inputs=input_files1,
    outputs=[output_files1, output_files2],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="20:00:00"
) << f"""
Rscript scripts/s2_2_analysis_check_bout_lengths.R {input_files1} {output_files1} {output_files2}
"""

#######################################################################################
### Step 5: Analysis MS prt1 - Maps + temporal patterns in activity & energy budgets ###
#######################################################################################
  
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"  # Add actual input files if necessary 
output_files1 = f"./results/tables/main/table2_budgets_species.csv" # Table with species-specific activity and energy budgets
output_files2 = f"./results/tables/main/table3_budgets_population.csv" # Table with population-specific activity and energy budgets
output_files3 = f"./results/tables/main/table4_budgets_individual.csv" # Table with population-specific activity and energy budgets

gwf.target(
    name="s2_3_analysis_plot_maps_and_budgets",
    inputs=input_files1,
    outputs=[output_files1, output_files2, output_files3],
    cores=1,
    memory="8G",
    walltime="07:00:00"
) << f"""
Rscript scripts/s2_3_analysis_plot_maps_and_budgets.R {input_files1} {output_files1} {output_files2} {output_files3}
"""

#######################################################################################
### Step 6: Analysis MS prt2 - Calculating migratory distance & weekly deviance #######
#######################################################################################
  
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"  # Add actual input files if necessary 
output_files1 = f"./results/tables/main/table5_migratory_distance.csv"    # Table with individual migratory distance
output_files2 = f"./results/tables/main/table6_migratory_distance_sst_anomaly.csv" # Table with individual migratory distance
output_files3 = f"./results/tables/main/table7_species_mean_deviance.csv" # Table with species deviance in WEE, weighted for varying sample size
output_files4 = f"./results/tables/main/table8_population_mean_deviance.csv" # Table with population deviance in WEE

gwf.target(
    name="s2_4_analysis_calc_migratory_dist",
    inputs=input_files1,
    outputs=[output_files1, output_files2, output_files3, output_files4],
    queue="bigmem",
    cores=1,
    memory="80G",
    walltime="13:00:00"
) << f"""
Rscript scripts/s2_4_analysis_calc_migratory_dist.R {input_files1} {output_files1} {output_files2} {output_files3} {output_files4} 
"""

#######################################################################
### Step 7: Analysis MS prt3 - Effect deviance on total NB costs ######
#######################################################################
 
# Determine study period 
startDate = "2021-09-01"
endDate = "2022-04-30"

# Determine intput and output files
 output_files1 = f"./results/tables/main/table9_totalNBCosts.csv"
output_files2 = f"./results/tables/main/table10_stats3_dredge_nbCosts.csv"
output_files3 = f"./results/tables/main/table11_cluster_pop.csv"
output_files4 = f"./results/tables/main/table12_stats4_dredge_nbCosts_species.csv"
output_files5 = f"./results/tables/main/table13_stats5_dredge_nbCosts_deviance.csv"

gwf.target(
    name="s2_5_anaysis_weekly_vs_tot_energy",
    inputs=[], 
    outputs=[output_files1, output_files2, output_files3, output_files4, output_files5],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="13:00:00"
) << f"""
Rscript scripts/s2_5_anaysis_weekly_vs_tot_energy.R {startDate} {endDate} {output_files1} {output_files2} {output_files3} {output_files4} {output_files5}
"""

#######################################################################
### Step 8: For Lila - Calculate activity budgets every day of year ###
#######################################################################

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