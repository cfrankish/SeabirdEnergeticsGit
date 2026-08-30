library(ggplot2)
library(dplyr)
library(fields)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(fasterize)
library(gdistance)
library(sf)
library(sp)
library(tidyr)
library(suncalc)
library(lubridate)
library(data.table)
library(terra)
library(ncdf4)
library(gridExtra)
library(lme4)
library(MuMIn)
library(pracma)
library(lmerTest)

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract species name and date from command line arguments
species_name <- args[1]  # First argument is the species name
today_date <- args[2]  # Second argument is today's date

# Print the species name for debugging
cat("Processing species:", species_name, "\n")

# Set up minimum sample size & number of iterations
minSampleSize<-5
print(paste0("min sample size per colony is: ", minSampleSize))
reps<-50
print(paste0("min iteration number is: ", reps))

### Step 1: Subset main id file to species of interest ###

# Read the input data
input_file <- "./results/tables/main/table1_idcatalogue.csv"
data <- read.csv(input_file)
data$species2 <- data$species %>%
  gsub("-", "", .) %>%       # Remove dashes
  gsub(" ", "", .) %>%       # Remove spaces
  gsub("'", "", .) %>%       # Remove apostrophes
  gsub("ü", "u", .)          # Replace 'ü' with 'u'

# Subset the data for the specific species
energyAll <- data %>% filter(species2 == species_name)

### Step 2: Creat species-specific files with activity budgets throughout the year ### 

allResults<-list.files("./tmp", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# Make a species list to cycle through
speciesList<-unique(energyAll$species)

activityAll<-list() 

for (i in 1:length(speciesList)) {

# subset to species i
speciesSub<-speciesList[i]

# Make a list of ids to cycle through
energySub<-subset(energyAll, species==speciesSub)
ids<-unique(energySub$individ_id)

# Make a list to save temporary results in
activitySpecies<-list()

# Cycle through iterations
for (j in 1:reps) {

# make an empty list to save results in
activityRep<-list()

# cycle through ids
for (k in 1:length(ids)) {

print(paste0("Species ", i, ", Rep ", j, ", Bird ", k, "/", length(ids)))

# subset to bird k
birdSub<-fread(energyRes_day[grepl(ids[k], energyRes_day)])

# Subset to rep j
birdSubSub<-subset(birdSub, rep==j)

# Add column tActive if it doesn't exist
if ("tActive" %in% colnames(birdSubSub) == FALSE) {
		  
birdSubSub$tActive<-0

}

# Determine day of the year (doy) as we are calculting activity budgets for days 1-365. We write some code below to remove 29th of Feb in leap years. 
bird_activity<-birdSubSub %>%
ungroup() %>%
dplyr::mutate(day=as.numeric(substr(date, 9 ,10)), month=as.numeric(substr(date, 6, 7)), day_month=paste0(day, "_", month), year=substr(date, 1, 4)) %>%
dplyr::group_by(year) %>%
dplyr::mutate(leapyear=ifelse(day_month=="29_2", 1, 0)) %>%
dplyr::mutate(leapyear=cumsum(leapyear)) %>%
dplyr::filter(!day_month %in% c("29_2")) %>%
dplyr::select(rep, species, colony, individ_id, date, tActive, tForage, tFlight, tLand, tRestWater, leapyear, year) %>%
dplyr::mutate(doy=as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01")), unit=c("days"))) + 1) %>%
dplyr::mutate(doy=ifelse(max(leapyear)>0 & doy > 59, doy-1, doy)) %>%
dplyr::mutate(doy=ifelse(max(doy)==366 & doy >59, doy-1, doy)) %>%
ungroup()%>%
dplyr::select(-c(leapyear, year)) 

maxdoy<-max(bird_activity$doy)

if (maxdoy==366) {
stop(print("Error"))}

# Save result
activityRep<-rbind(activityRep, bird_activity)

}

# Summarise species-specific information for rep j
activityRep_sum1<-activityRep %>%
ungroup() %>%
dplyr::group_by(rep, species, colony, doy) %>%
dplyr::summarise(birds=n_distinct(individ_id), meanActive=mean(tActive), sdActive=sd(tActive), seActive=sdActive/sqrt(birds), lowerCI_active=meanActive-1.96*seActive, upperCI_active=meanActive + 1.96*seActive,
meanFlight=mean(tFlight), sdFlight=sd(tFlight), seFlight=sdFlight/sqrt(birds), lowerCI_flight=meanFlight-1.96*seFlight, upperCI_flight=meanFlight + 1.96*seFlight,
meanForage=mean(tForage), sdForage=sd(tForage), seForage=sdForage/sqrt(birds), lowerCI_forage=meanForage-1.96*seForage, upperCI_forage=meanForage + 1.96*seForage,
meanLand=mean(tLand), sdLand=sd(tLand), seLand=sdLand/sqrt(birds), lowerCI_land=meanLand-1.96*seLand, upperCI_land=meanLand + 1.96*seLand,
meanRestWater=mean(tRestWater), sdRestWater=sd(tRestWater), seRestWater=sdRestWater/sqrt(birds), lowerCI_RestWater=meanRestWater-1.96*seRestWater, upperCI_RestWater=meanRestWater + 1.96*seRestWater) 

#naFlight<-subset(activityRep_sum1, is.na(sdFlight))

# Save result
activitySpecies<-rbind(activitySpecies, activityRep_sum1)

	}

# Summarize this result and Save
activitySpecies_sum<-activitySpecies %>%
ungroup() %>%
dplyr::group_by(species, colony, doy) %>%
dplyr::summarise(repsTot=reps-sum(is.na(seFlight)), birdsMean=mean(birds), minbirds=min(birds), maxbirds=max(birds), meanActive=mean(meanActive, na.rm=TRUE), sdActive=mean(sdActive, na.rm=TRUE), seActive=mean(seActive, na.rm=TRUE), lowerCI_active=mean(lowerCI_active, na.rm=TRUE), upperCI_active=mean(upperCI_active, na.rm=TRUE),
meanFlight=mean(meanFlight, na.rm=TRUE), sdFlight=mean(sdFlight, na.rm=TRUE), seFlight=mean(seFlight, na.rm=TRUE), lowerCI_flight=mean(lowerCI_flight, na.rm=TRUE), upperCI_flight=mean(upperCI_flight, na.rm=TRUE),
meanForage=mean(meanForage, na.rm=TRUE), sdForage=mean(sdForage, na.rm=TRUE), seForage=mean(seForage, na.rm=TRUE), lowerCI_forage=mean(lowerCI_forage, na.rm=TRUE), upperCI_forage=mean(upperCI_forage, na.rm=TRUE),
meanLand=mean(meanLand, na.rm=TRUE), sdLand=mean(sdLand, na.rm=TRUE), seLand=mean(seLand, na.rm=TRUE), lowerCI_land=mean(lowerCI_land, na.rm=TRUE), upperCI_land=mean(upperCI_land, na.rm=TRUE),
meanRestWater=mean(meanRestWater, na.rm=TRUE), sdRestWater=mean(sdRestWater, na.rm=TRUE), seRestWater=mean(seRestWater, na.rm=TRUE), lowerCI_RestWater=mean(lowerCI_RestWater, na.rm=TRUE), upperCI_RestWater=mean(upperCI_RestWater, na.rm=TRUE))

}

# Define the output file name based on species and date
output_file <- paste0("./results/activity_Lila/dailyactivity_", species_name, "_", today_date, ".csv")

# Save output files
print("Saving output files")

# Number 2
output_file1 <- args[3]
print("Saving output file 1")
write.csv(activitySpecies_sum, file = output_file, row.names = FALSE) # species-specific budgets



