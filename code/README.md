# Code

## Scripts needed to run workflow for paper 1 (in review)

'North Atlantic seabirds exhibit strong variability in seasonal energy expenditure with implications for energetic bottlenecks and total non-breeding costs' (Frankish et al.)

### Step 1 extracts & processes immersion data:
- `S1_1_extract_wetdry.R`: this step does NOT work in parallel yet. I have to run in the front-end. By species, the script extracts immersion data from the SEATRACK database. It then retains tracks with IDs that are in the IRMA data, removes tracks with weird sampling intervals & juveniles. This step outputs species-specific files ('species_wetdry_date.csv') which are backed-up locally.
- `S1_2_downloadSST_underConstruction.R`: this step is meant to donwload SST data automatically from COPERNICUS but I can't get this working on SAGA. Instead I have uploaded monthly SST files myself from my pc. These were downloaded using the script called `s4_environmental_data.R`
- `S1_3_merge_wetdry_lox_env.R`: merges immersion, location and environmental data by individual bird. Each 10-mins of data is also assigned day/night/twilight. Files are split into individual files stored in 'birddata_raw/speciesname/id.csv'
  
### Step 2 calculates daily activity budgets and energy expenditure:
- `S2_1_calculatebudgets.R`: calculates the above for every individual bird x 100. Input file is 'birddata_raw/speciesname/id.csv', and output file is 'tmp/id_energyDay.csv'.
- `S2_2_sensitivity_analysis_part1.R`: conducts sensitivity analysis whereby activity & energy is calculated for every individual where each parameter used in our energetic approach is modified in turn according to a few values. Input file is 'birddata_raw/speciesname/id.csv', and output file is 'tmp2/id_energyDay.csv'.
- `S2_3_sensitivity_analysis_part2.R`: gathers individual files, and summarizes the effect of changing each parameter on weekly variation and total energy expenditure for the six species.

### Step 3 runs some supplementary analyses:
- `S3_1_supanalysis_min_iteration_number.R`: Conducts an analysis where it calculates total non-breeding season energy expenditure for an increasing number of iterarions so I can decide what a good number is for further analysis. It outputs Figure S11. 
- `S3_2_supanalysis_check_bout_lengths.R`: Does some checks to validate my approach for dry bout allocation. Specificially, it looks into how many dry bouts are allocated to flight vs. land at random, calculates the duration of flight bouts during darkness, and creates a distirbution plot of fulmar flight bout lengths. It outputs Figures S5-7.

### Step 4 mainly conducts the main analysis:
- `S4_1_mainanalysis_plot_maps_budgets.R`: This script maps the distribution of the study populations and their non-breeding season distributions (Figures 1, S1 & S2). It also estimates yearly activity budgets, SST & energy expenditure at a species and population-level for visual exploration.
- `S4_2_mainanalysis_calc_metrics.R`: This script calculates migratory distance and weekly variation in energy expenditure used to make Figure 3B. It also outputs all supplementary figures showing weekly deviation in different behaviours & SST (Figures S15-S20).
- `S4_3_mainanalysis_calc_metrics.R`: This script conducts all statistics in the manuscript. It outputs figures 2, 3A & all remaining supplementary figures (Figures S8, 12-14, 21-34). 
