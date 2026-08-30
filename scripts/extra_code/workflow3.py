###############################################
### Step 5_1_3: Prep activity (spatial)  ######
###############################################

# Directory containing the CSV files
budget_dir = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3"

# Find all CSV files in the directory
budget_files = glob(os.path.join(budget_dir, "*daily.csv"))

# Ensure tmp4 exists
os.makedirs("./tmp5_spatial_activity", exist_ok=True)

# Make a list of months to loop through:
months = [9, 10, 11, 12, 1, 2, 3, 4]

# Make a list of behaviours to loop through:
beh = ["Forage", "RestWater", "Active", "Land", "Flight"]

# Define secondary input file
input_file2 = "/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/main/table1_idcatalogue.csv"

# Process each CSV file
for budget_file in budget_files:

    # Get just the file name
    filename = os.path.basename(budget_file)

    # Remove the extension
    name_no_ext = filename.replace(".csv", "")

    # Remove special characters
    name_no_ext_clean = (
        name_no_ext.replace("(", "")
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
                   .replace("é", "")
    )

    # Remove special characters in path
    budget_clean = (
        budget_file.replace("(", "")
                   .replace(")", "")
                   .replace("`", "")
                   .replace("'", "")
    )

    # Loop over months and behaviours
    for m in months:
        for b in beh:

            # Output name includes month + behaviour
            output_file = f"./tmp5_spatial_activity/{name_no_ext_clean}_m{m}_{b}.csv"

            # Define inputs
            input_files = [budget_file, input_file2]

            gwf.target(
                name=f"s5_1_3_map_activity_{name_no_ext_clean}_m{m}_{b}",
                inputs=input_files,
                outputs=[output_file],
                cores=1,
                memory="8G",
                walltime="04:00:00"
            ) << f"""
            mkdir -p tmp/
            Rscript scripts/s5_1_3_prep_activity_distribution_spatial.R {budget_file} {input_file2} {b} {m} {output_file}
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
    output_file1 = f"./results/tables/{name}_energyMap_monthly_v1.csv"
    
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