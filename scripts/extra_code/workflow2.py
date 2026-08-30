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
        walltime="72:00:00"
    ) << f"""
    mkdir -p results
    Rscript scripts/forLila.R {sp} {today} {output_file}
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
    
###############################################
### Step 5_1_2: Is activity spatial?  #########
###############################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*daily.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp4", exist_ok=True)

# Make a list of species to loop through:
species=["Blackleggedkittiwake", "Littleauk", "Commonguillemot", "Northernfulmar", "Atlanticpuffin", "Brunnichsguillemot"]

# Make a list of months to loop through:
months = [9, 10, 11, 12, 1, 2, 3, 4]

# Make a list of behaviours to loop through:
beh = ["Forage", "RestWater", "Active", "Land", "Flight"]

# Define secondary input file
input_file2 = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/main/table1_idcatalogue.csv"

# Process files for each species
for speciesSub in species:

    # Loop over months and behaviours
    for m in months:
        for b in beh:

            # Output name includes month + behaviour
            output_file1 = f"./tmp4/{speciesSub}_{b}_m{m}_stats.csv"
            output_file2 = f"./tmp4/{speciesSub}_{b}_m{m}_rasters_gam1.rds" # No spatial component
            output_file3 = f"./tmp4/{speciesSub}_{b}_m{m}_rasters_gam2.rds" # spatial (species-level)
            output_file4 = f"./tmp4/{speciesSub}_{b}_m{m}_rasters_gam3.rds" # spatial with population deviations

            # Define inputs
            input_files = [input_file2]

            gwf.target(
                name=f"s5_1_2_test_spatial_{speciesSub}_{b}_m{m}",
                inputs=input_files,
                outputs=[output_file1, output_file2, output_file3, output_file4],
                queue="bigmem",
                cores=1,
                memory="20G",
                walltime="72:00:00"
            ) << f"""
            mkdir -p tmp/
            Rscript scripts/s5_1_2_is_activity_spatial.R {speciesSub} {input_file2} {b} {m} {output_file1} {output_file2} {output_file3} {output_file4}
            """

###############################################
### Step 5_1_3: Is activity spatial?  #########
###############################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*daily.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp4", exist_ok=True)

# Make a list of species to loop through:
species=["Blackleggedkittiwake", "Littleauk", "Commonguillemot", "Northernfulmar", "Atlanticpuffin", "Brunnichsguillemot"]

# Make a list of months to loop through:
months = [9, 10, 11, 12, 1, 2, 3, 4]

# Make a list of behaviours to loop through:
beh = ["Forage", "RestWater", "Active", "Land", "Flight"]

# Define secondary input file
input_file2 = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/main/table1_idcatalogue.csv"

# Process files for each species
for speciesSub in species:

    # Loop over months and behaviours
    for m in months:
        for b in beh:

            # Output name includes month + behaviour
            output_file1 = f"./tmp4/{speciesSub}_{b}_m{m}_stats_v2.csv"
            output_file2 = f"./tmp4/{speciesSub}_{b}_m{m}_rasters_gam1_v2.rds" # No spatial component
            output_file3 = f"./tmp4/{speciesSub}_{b}_m{m}_rasters_gam2_v2.rds" # spatial (species-level)
            output_file4 = f"./tmp4/{speciesSub}_{b}_m{m}_rasters_gam3_v2.rds" # spatial with population deviations

            # Define inputs
            input_files = [input_file2]

            gwf.target(
                name=f"s5_1_3_test_spatial_{speciesSub}_{b}_m{m}",
                inputs=input_files,
                outputs=[output_file1, output_file2, output_file3, output_file4],
                queue="bigmem",
                cores=1,
                memory="20G",
                walltime="72:00:00"
            ) << f"""
            mkdir -p tmp/
            Rscript scripts/s5_1_3_is_activity_spatial_part2.R {speciesSub} {input_file2} {b} {m} {output_file1} {output_file2} {output_file3} {output_file4}
            """

###########################################################################################
### Step 5_2_1: Map energy expenditure (assume same behaviour regardless of location ######
###########################################################################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*monthly.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp5", exist_ok=True)

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
    output_files1 = f"./results/tables/supplementary/{name_no_ext_clean}_energyMap.csv"
    
    gwf.target(
        name=f"s5_2_1_map_energy_{name_no_ext_clean}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="24:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s5_2_1_map_energy_v1.R {budget_file} {output_files1}
    """       

###########################################################################################
### Step 5_2_2: Map energy expenditure (spatialized #######################################
###########################################################################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*monthly.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp5", exist_ok=True)

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
    output_files1 = f"./results/tables/supplementary/{name_no_ext_clean}_energyMap_v2.csv"
    
    gwf.target(
        name=f"s5_2_2_map_energy_spatial_{name_no_ext_clean}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="10:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s5_2_2_map_energy_v2.R {budget_file} {output_files1}
    """       
    
################################################################
### Step 5_3_1: Aggregate energy - non-spatial activity  #######
################################################################

import re

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp5"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*.csv"))

# Extract species names
results = [re.findall(r"tmp5/([^_]+)_", str(p))[0] for p in budget_files]
unique_names = set(results)

# Process each CSV file
for name in unique_names:
    # Define output files
    output_file1 = f"./results/tables/{name}_energyMap_monthly_map_v1.csv"
    output_file2 = f"./results/tables/{name}_energyMap_monthly_sum_v1.csv"
    
    gwf.target(
        name=f"s5_3_1_aggregate_maps_energy_{name}",
        inputs=[],
        outputs=[output_file1, output_file2],
        queue="bigmem",
        cores=1,
        memory="20G",
        walltime="10:00:00",
    ) << f"""
    Rscript scripts/s5_3_1_aggregate_maps_energy_species_v1.R {name} {output_file1} {output_file2} 
    """    

################################################################
### Step 5_3_2: Aggregate energy - spatial activity  ###########
################################################################

import re

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp5"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*.csv"))

# Extract species names
results = [re.findall(r"tmp5/([^_]+)_", str(p))[0] for p in budget_files]
unique_names = set(results)

# Process each CSV file
for name in unique_names:
    # Define output files
    output_file1 = f"./results/tables/{name}_energyMap_monthly_map_v2.csv"
    output_file2 = f"./results/tables/{name}_energyMap_monthly_sum_v2.csv"
    
    gwf.target(
        name=f"s5_3_2_aggregate_maps_energy_{name}",
        inputs=[],
        outputs=[output_file1, output_file2],
        queue="bigmem",
        cores=1,
        memory="20G",
        walltime="10:00:00",
    ) << f"""
    Rscript scripts/s5_3_2_aggregate_maps_energy_species_v2.R {name} {output_file1} {output_file2} 
    """        

################################################################
### Step 5_3_3: Compare both methods, species by species #######
################################################################

import re

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp5"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*.csv"))

# Extract species names
results = [re.findall(r"tmp5/([^_]+)_", str(p))[0] for p in budget_files]
unique_names = set(results)

# Process each CSV file
for name in unique_names:
    
    # Define output files
    input_file1 = f"./results/tables/{name}_energyMap_monthly_map_v1.csv"
    input_file2 = f"./results/tables/{name}_energyMap_monthly_map_v2.csv"
    input_file3 = f"./results/tables/{name}_energyMap_monthly_sum_v1.csv"
    input_file4 = f"./results/tables/{name}_energyMap_monthly_sum_v2.csv"
    
    output_file1 = f"./results/tables/{name}_energyMap_monthly_sums.csv"
    
    gwf.target(
        name=f"s5_3_3_compare_energy_{name}",
        inputs=[input_file1, input_file2, input_file3, input_file4],
        outputs=[output_file1],
        queue="bigmem",
        cores=1,
        memory="20G",
        walltime="5:00:00",
    ) << f"""
    Rscript scripts/s5_3_3_compare.R {input_file1} {input_file2} {input_file3} {input_file4} {output_file1} 
    """     

################################################################
### Step 5_4_1: Aggregate everything ###########################
################################################################

# Determine names of input & output files # 
output_files1 = f"./results/tables/tablex_allenergy_v1.csv"

gwf.target(
    name="s5_4_1_aggregate_all_energy_v1",
    inputs=[],
    outputs=[output_files1],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="20:00:00"
) << f"""
Rscript scripts/s5_4_1_aggregate_maps_energy_all_v1.R {output_files1} 
"""

################################################################
### Step 5_5_1: Identify energy demand hotspots ################
################################################################

# Define output files
output_file1 = f"./results/tables/energyhotspots.csv"

gwf.target(
    name="s5_5_1_identify_hotspots",
    inputs=[],
    outputs=[output_file1],
    queue="bigmem",
    cores=1,
    memory="10G",
    walltime="20:00:00"
) << f"""
Rscript scripts/s5_5_1_identify_hotspots.R {output_file1} 
"""

######################################################################
### Step 2_4_1: Extra analysis for Seb - make Svalbard bird energetic ###
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
    walltime="4:00:00"
) << f"""
Rscript scripts/s2_4_energy_seb_svalbard.R {input_files1} {output_files1} 
"""

######################################################################
### Step 2_4_2: Extra analysis for Seb - make Svalbard bird energetic ###
#######################################################################
 
# Determine intput and output files
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"
output_files1 = f"./results/tables/main/seb_svalbard_bjornoya.csv"

gwf.target(
    name="forSeb2",
    inputs=[input_files1], 
    outputs=[output_files1],
    cores=1,
    memory="10G",
    walltime="4:00:00"
) << f"""
Rscript scripts/s2_4_energy_seb_bjornoya.R {input_files1} {output_files1} 
"""

######################################################################
### Step 2_4_1: Extra analysis for Seb - make Svalbard bird energetic ###
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
    walltime="4:00:00"
) << f"""
Rscript scripts/s2_4_energy_seb_svalbard.R {input_files1} {output_files1} 
"""

######################################################################
### Step 2_4_2: Extra analysis for Seb - make Svalbard bird energetic ###
#######################################################################
 
# Determine intput and output files
input_files1 = f"./results/tables/main/table1_idcatalogue.csv"
output_files1 = f"./results/tables/main/seb_svalbard_bjornoya.csv"

gwf.target(
    name="forSeb2",
    inputs=[input_files1], 
    outputs=[output_files1],
    cores=1,
    memory="10G",
    walltime="4:00:00"
) << f"""
Rscript scripts/s2_4_energy_seb_bjornoya.R {input_files1} {output_files1} 
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
    
###############################################
### Step 5_1_2: Prep activity (spatial)  ######
###############################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*daily.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp4", exist_ok=True)

# Make a list of species to loop through:
species=["Blackleggedkittiwake", "Littleauk", "Commonguillemot", "Northernfulmar", "Atlanticpuffin", "Brunnichsguillemot"]

# Make a list of months to loop through:
months = [9, 10, 11, 12, 1, 2, 3, 4]

# Make a list of behaviours to loop through:
beh = ["Forage", "RestWater", "Active", "Land", "Flight"]

# Define secondary input file
input_file2 = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/main/table1_idcatalogue.csv"

# Process files for each species
for speciesSub in species:

    # Loop over months and behaviours
    for m in months:
        for b in beh:

            # Output name includes month + behaviour
            output_file1 = f"./tmp4/{speciesSub}_m{m}_{b}_rasters.csv"
            output_file2 = f"./tmp4/{speciesSub}_m{m}_{b}_stats.csv"

            # Define inputs
            input_files = [input_file2]

            gwf.target(
                name=f"s5_1_2_map_activity_{speciesSub}_m{m}_{b}",
                inputs=input_files,
                outputs=[output_file1, output_file2],
                queue="bigmem",
                cores=1,
                memory="20G",
                walltime="72:00:00"
            ) << f"""
            mkdir -p tmp/
            Rscript scripts/s5_1_2_prep_activity_distribution_spatial.R {speciesSub} {input_file2} {b} {m} {output_file1} {output_file2}
            """

###########################################################################################
### Step 5_2_1: Map energy expenditure (assume same behaviour regardless of location ######
###########################################################################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*monthly.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp5", exist_ok=True)

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
    output_files1 = f"./results/tables/supplementary/{name_no_ext_clean}_energyMap.csv"
    
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
    
###########################################################################################
### Step 5_2_2: Map energy expenditure (spatialized #######################################
###########################################################################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*monthly.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp5", exist_ok=True)

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
    output_files1 = f"./results/tables/supplementary/{name_no_ext_clean}_energyMap_v2.csv"
    
    gwf.target(
        name=f"s5_2_2_map_energy_spatial_{name_no_ext_clean}",
        inputs=input_files,
        outputs=[output_files1],
        cores=1,
        memory="8G",
        walltime="04:00:00"
    ) << f"""
    mkdir -p tmp/
    Rscript scripts/s5_2_2_map_energy_v2.R {budget_file} {output_files1}
    """       

################################################################
### Step 5_3_1: Aggregate energy - non-spatial activity  #######
################################################################

import re

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp5"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*.csv"))

# Extract species names
results = [re.findall(r"tmp5/([^_]+)_", str(p))[0] for p in budget_files]
unique_names = set(results)

# Process each CSV file
for name in unique_names:
    # Define output files
    output_file1 = f"./results/tables/{name}_energyMap_monthly_map_v1.csv"
    output_file2 = f"./results/tables/{name}_energyMap_monthly_sum_v1.csv"
    
    gwf.target(
        name=f"s5_3_1_aggregate_maps_energy_{name}",
        inputs=[],
        outputs=[output_file1, output_file2],
        queue="bigmem",
        cores=1,
        memory="20G",
        walltime="10:00:00",
    ) << f"""
    Rscript scripts/s5_3_1_aggregate_maps_energy_species_v1.R {name} {output_file1} {output_file2} 
    """    

################################################################
### Step 5_3_2: Aggregate energy - spatial activity  ###########
################################################################

import re

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp5"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*.csv"))

# Extract species names
results = [re.findall(r"tmp5/([^_]+)_", str(p))[0] for p in budget_files]
unique_names = set(results)

# Process each CSV file
for name in unique_names:
    # Define output files
    output_file1 = f"./results/tables/{name}_energyMap_monthly_map_v2.csv"
    output_file2 = f"./results/tables/{name}_energyMap_monthly_sum_v2.csv"
    
    gwf.target(
        name=f"s5_3_2_aggregate_maps_energy_{name}",
        inputs=[],
        outputs=[output_file1, output_file2],
        queue="bigmem",
        cores=1,
        memory="20G",
        walltime="10:00:00",
    ) << f"""
    Rscript scripts/s5_3_2_aggregate_maps_energy_species_v2.R {name} {output_file1} {output_file2} 
    """        

################################################################
### Step 5_3_3: Compare both methods, species by species #######
################################################################

import re

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp5"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*.csv"))

# Extract species names
results = [re.findall(r"tmp5/([^_]+)_", str(p))[0] for p in budget_files]
unique_names = set(results)

# Process each CSV file
for name in unique_names:
    
    # Define output files
    input_file1 = f"./results/tables/{name}_energyMap_monthly_map_v1.csv"
    input_file2 = f"./results/tables/{name}_energyMap_monthly_map_v2.csv"
    input_file3 = f"./results/tables/{name}_energyMap_monthly_sum_v1.csv"
    input_file4 = f"./results/tables/{name}_energyMap_monthly_sum_v2.csv"
    
    output_file1 = f"./results/tables/{name}_energyMap_monthly_sums.csv"
    
    gwf.target(
        name=f"s5_3_3_compare_energy_{name}",
        inputs=[input_file1, input_file2, input_file3, input_file4],
        outputs=[output_file1],
        queue="bigmem",
        cores=1,
        memory="20G",
        walltime="5:00:00",
    ) << f"""
    Rscript scripts/s5_3_3_compare.R {input_file1} {input_file2} {input_file3} {input_file4} {output_file1} 
    """     
    
################################################################
### Step 5_4_1: Aggregate everything ###########################
################################################################

# Determine names of input & output files # 
output_files1 = f"./results/tables/tablex_allenergy_v1.csv"

gwf.target(
    name="s5_4_1_aggregate_all_energy_v1",
    inputs=[],
    outputs=[output_files1],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="20:00:00"
) << f"""
Rscript scripts/s5_4_1_aggregate_maps_energy_all_v1.R {output_files1} 
"""

################################################################
### Step 5_4_2: Aggregate everything ###########################
################################################################

# Determine names of input & output files # 
output_files1 = f"./results/tables/tablex_allenergy_v2.csv"

gwf.target(
    name="s5_4_2_aggregate_all_energy_v2",
    inputs=[],
    outputs=[output_files1],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="20:00:00"
) << f"""
Rscript scripts/s5_4_2_aggregate_maps_energy_all_v2.R {output_files1} 
"""

################################################################
### Step 5_4_3: Aggregate everything ###########################
################################################################

# Define output files
input_file1 = f"./results/tables/tablex_allenergy_v1.csv"
input_file2 = f"./results/tables/tablex_allenergy_v2.csv"

gwf.target(
    name="s5_4_3_compare_all_energy",
    inputs=[input_file1, input_file2],
    outputs=[],
    queue="bigmem",
    cores=1,
    memory="50G",
    walltime="20:00:00"
) << f"""
Rscript scripts/s5_4_3_compare.R {input_file1} {input_file2}
"""

################################################################
### Step 5_5_1: Identify energy demand hotspots ################
################################################################

# Define output files
output_file1 = f"./results/tables/energyhotspots.csv"

gwf.target(
    name="s5_5_1_identify_hotspots",
    inputs=[],
    outputs=[output_file1],
    queue="bigmem",
    cores=1,
    memory="10G",
    walltime="20:00:00"
) << f"""
Rscript scripts/s5_5_1_identify_hotspots.R {output_file1}
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
        walltime="72:00:00"
    ) << f"""
    mkdir -p results
    Rscript scripts/forLila.R {sp} {today} {output_file}
    """