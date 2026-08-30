## This script maps activity & energy expenditure for our species ##
## To do this, by pop, the script will open the corresponding distribution of activity budgets ##
## For energy, it will extract the correct SST map & calculate energy 100 times (in line with the pop maps too?) ##
## Then it will map activity budgets somewhow too ##
## The output will be five rasters per pop ##

library(dplyr)
library(ncdf4)
library(raster)
library(sp)
library(sf)
library(lubridate)
library(data.table)
library(tidyr)
library(ggplot2)

#### Step 0: setting up basic conditions ####

# Set-up some projections
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # Flat projection
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"

# Set-up number of iterations...
overall.iterations<-3 # how many times this is calculated per individual
print(paste0("Determining activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# input_file1<-"tmp3/Northernfulmar_breidafjordur_actbudgets_daily.csv"
input_file1 <- args[1] # This will read in a population-specific activity budget  

print(input_file1)

# open input file 1
actBudgets<-fread(input_file1) # actBudgets<-read.csv("tmp3/Northernfulmar_breidafjordur_actbudgets_daily.csv")

### Step 1: Source all data needed for calculations ###

print("Sourcing matching data...")

# Find matching species raster #
rasterList<-data.frame(files=list.files("data/popdata_raw/", full.names=TRUE), species=c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot"))
rasterSub<-subset(rasterList, species==actBudgets$species[1]) # Subset to correct species
nc<-nc_open(rasterSub$files[1]) # open the file

# Create metadata of this file so I can easily pull the correct raster from it
modelcolonies <- (ncvar_get(nc,"SmcolName")) # Find list of populations
colonies <- (ncvar_get(nc,"colonyName")) # Find list of ALL colonies
colcode <- (ncvar_get(nc,"SmcolCode")) # Colony code
meta <- data.frame(modelcolonies, colonies, colcode) # Make some metadata so easier to understand raster structure
meta$species<-actBudgets$species[1] # Add species names
colonybase<-basename(input_file1) # Now we try and extract the colony from this name (first we remove the path)
parts <- strsplit(colonybase, "_")[[1]] # now we split it into sections divided by "_"
colonyPart <- parts[2] # select colony
print(colonyPart)
colonyFinal<-paste0(toupper(substr(colonyPart, 1, 1)), substr(colonyPart, 2, nchar(colonyPart))) # Change first letter to uppercase

if(colonyFinal == c("Karmoy") & actBudgets$species[1] == "Northern fulmar") {
colonyFinal="Jarsteinen"
}

if(colonyFinal == c("Kaphoegh")) {
colonyFinal="Kap_Hoegh"
}

if(colonyFinal == c("Faroeislands")) {
colonyFinal="Faroe_Islands"
}

if(colonyFinal == c("Isleofcanna")) {
colonyFinal="Isle_of_Canna"
}

if(colonyFinal == c("Isleofmay")) {
colonyFinal="Isle_of_May"
}

if(colonyFinal == c("Capekrutik")) {
colonyFinal="Cape_Krutik"
}

if(colonyFinal == c("Russkayagavan")) {
colonyFinal="Russkaya_Gavan"
}

if(colonyFinal == c("Franzjosefland")) {
colonyFinal="Franz_Josef_Land"
}

if(colonyFinal == c("Karagate")) {
colonyFinal="Kara_Gate"
}

if(colonyFinal == c("Rundeandaalesund")) {
colonyFinal="Runde_and_Aalesund"
}

if(colonyFinal == c("Skelligmichael")) {
colonyFinal="Skellig_Michael"
}

if(colonyFinal == c("Capegorodetskiy")) {
colonyFinal="Cape_Gorodetskiy"
}

if(colonyFinal == c("Littlesaltee")) {
colonyFinal="Little_Saltee"
}

if(colonyFinal == c("Coatsisland")) {
colonyFinal="Coats_Island"
}

if(colonyFinal == c("Janmayen")) {
colonyFinal="Jan_Mayen"
}

if(colonyFinal == c("Rost")) {
colonyFinal="Roest"
}

if(colonyFinal == c("Oranskieislands")) {
colonyFinal="Oranskie_Islands"
}

if(colonyFinal == c("Witlessbay")) {
colonyFinal="Witless_Bay"
}

if(colonyFinal == c("Karmoy")) {
colonyFinal="Karmøy"
}

if(colonyFinal == c("Bjornoya")) {
colonyFinal="Bjoernoeya"
}

if(colonyFinal == c("Hornoya")) {
colonyFinal="Hornoeya"
}

if(colonyFinal == c("Hjelmsoya")) {
colonyFinal="Hjelmsoeya"
}

metaSub<-subset(meta, modelcolonies==colonyFinal) # Find out model colony
print(nrow(metaSub))
metaAll<-distinct(subset(meta, modelcolonies==metaSub$modelcolonies[1])) # Subset to these...
print(nrow(metaAll))

# Determine location of SST data #
sstMonthly<-list.files("./data/sstPopMaps", full.names=TRUE)
sstMaps<-sstMonthly[grepl("tif", sstMonthly)] # Choose the tif file
sstMap<-raster::stack(sstMaps[1]) # Stack all layers (lol)
sstKey<-data.frame(layerNo=1:nlayers(sstMap)) # Make a key with dates so sstMap is easy to access
# Set start and end dates
start_date <- as.Date("2012-01-01")
end_date   <- as.Date("2023-01-01")
date_seq <- seq(from = start_date, to = end_date, by = "month") # Generate sequence for the 1st of each month
# Annotate key
sstKey$dates<-date_seq
sstKey$Year<-substr(sstKey$dates, 1, 4)
sstKey$Month<-as.numeric(substr(sstKey$dates, 6, 7))

# Determine location of air temp data #
airMonthly<-list.files("./data/temp", full.names=TRUE)
airMaps<-airMonthly[grepl("nc", airMonthly)] # Choose the tif file
airMaps_df<-data.frame(fileNames=airMaps) # Add file names 
airMaps_df$month<-sub(".*Month_(.*)\\.nc$", "\\1", airMaps_df$fileNames) # add month
airMaps_df$year<-sub(".*sst_([0-9]+)_Month.*", "\\1", airMaps_df$fileNames) # add year

# Determine model parameters to choose from #
speciesNo<-6
paramNo<-19
modelParams<-tibble(species=rep(c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot"), paramNo))
modelParams$parameter<-rep(c("L1", "Th1", "Th2", "L1_colony", "dist_colony", "pLand", "c",
"RMR", "c1", "c2", "c3", "c4", "c5", "TC_water", "Beta_active", "Beta_rest", "LCT_water", "TC_air", "LCT_air"), each=speciesNo)
modelParams$values<-list(
240, 810, 90, 134, 88, 88, # Species-specific flight bout duration (minutes)
0.95, 0.95, c(0.85, 0.9), c(0.85, 0.9), c(0.85, 0.9), c(0.85, 0.9), # Th1 % wet threshold for differenciating between behaviors
0, 0, 0, 0, 0, 0, # Th2 #% wet for determing dry vs. intermediate
c(0, 240), c(0, 810), c(0, 90), c(0, 134), c(0, 88), c(0, 88), # L1_colony: duration of dry bouts at start of tLand that can re-allocated to flight
500, 500, 500, 500, 500, 500,  # Distance to colony (km) below which it is considered possible to be on land
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, # pland: probability of dry being land or something else
c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), # coefficient for adjusting for leg-tucking
c(0), c(seq(0.96, 1.04, 0.01)), c(0), c(0), c(0), c(0), # RMR is just for a few of the species...BlKi / NoFu (Resting metabolic rate)
c(seq(2, 5.7, 0.1)), c(2.2), c(seq(106, 176, 1)), c(seq(5.4, 11.5, 0.1)), c(seq(106, 176, 1)), c(seq(106, 176, 1)), # c1 is the cost of Flight
c(seq(0.3, 3.8, 0.1)), c(seq(0.3, 3.8, 0.1)), c(0), c(0), c(0), c(0), # c2 is the cost of foraging...
c(seq(0.05, 1.1, 0.1)), c(0.8), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), # c3 is cost of being on land...
c(seq(0.1, 2.8, 0.1)), c(2), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), # c4 is the cost of resting on the water...
c(0), c(0), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), # c5 is the cost of being active on water when thermoneutral...(currently made in the code)
c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), 
c(0.0997), c(0.07), # TC is thermal conductivity (code is TC - 1.96*TCError)
c(2.75), c(2.75),
c(2.75), c(2.75),
c(0), c(0), c(118),  # Intercepts of resting metabolic rate at 0°C during different activities (active)
c(118), c(118), c(118),
c(1.87), c(1.34),# # Intercepts of resting metabolic rate at 0°C during different activities (rest)
c(72.2), c(72.2),
c(72.2), c(72.2),
12.5, 9, 14.18, 14.18, 14.18, 14.18, # LCT in water
0.0466, 0.0336, 0.0282, 0.05, 0.0282, 0.0282, # TC in air
4.5, 9, 5.72, 4.5, 2, 2) # LCt in air

### Step 2: Calculate energy & map it many times ####

# Determine number of model colonies #
modelColony_no<-nrow(metaAll)

# Make a list to save output file in
allResults<-list()

for (p in 1:nrow(metaAll)) {

# Determine name of model colony for later
colonyName<-metaAll[p,]

# Determine number of months to loop through
Months<-c(9, 10, 11, 12, 1, 2, 3, 4)

# Make an an empty list to save the results in
PopEnergyFinal_All<-list()

for (i in 1:length(Months)) {

# Subset to month i
monthSub<-Months[i]

# Prepare corresponding SST surface
sstKeySub<-subset(sstKey, Month == monthSub)
sstSub <- sstMap[[c(unique(sstKeySub$layerNo))]]
sstMean<-raster::calc(sstSub, fun = mean, na.rm = TRUE)

# Prepare correct air Temp surface 
airSub<-subset(airMaps_df, month==monthSub & year >= 2013 & year <=2023) # Subset to correct month
airStack<-raster::stack(airSub$fileNames) # Stack all these months
airMean<-raster::calc(airStack, fun = mean, na.rm = TRUE) # Average
values(airMean)<-values(airMean) - 273 # to change from kelvin to degrees C
airMean_crop<-projectRaster(airMean, sstMean) # crop to extent of SST raster

# Cut actBudgets
actMonth<-subset(actBudgets, month==monthSub)

# Extract correct population map
iMth<-monthSub
icols<-which(ncvar_get(nc, "colonyName") %in% colonyName$colonies[1]) # Determine where colony info is located
if (length(icols)>1) {icols<-icols[1]}

# Add number of birds & error
rast.mean <- raster(nc$filename, varname="PredictedAbundanceMean", band = icols, level = iMth) # Mean density
rast.ci.high<- raster(nc$filename, varname="PredictedAbundance95%CIHigh", band = icols, level = iMth) # 95% CI
rast.ci.low<- raster(nc$filename, varname="PredictedAbundance95%CILow", band = icols, level = iMth) # 95% CI

# Turn into data frames and add to energy Map file
meanDensity<-as.data.frame(rast.mean, xy=TRUE)
colnames(meanDensity)<-c("x", "y", "NoBirds_mean")
meanDensity_upper<-as.data.frame(rast.ci.high, xy=TRUE)
colnames(meanDensity_upper)<-c("x", "y", "NoBirds_ci_upper")
meanDensity_lower<-as.data.frame(rast.ci.low, xy=TRUE)
colnames(meanDensity_lower)<-c("x", "y", "NoBirds_ci_lower")

for (j in 1:overall.iterations) {

print(paste0("Mapping energy for month ", Months[i], " : iteration ", j))

### We generate activity budget from spatial maps ###

print("Randomising individual activity budget...")

# To do this we first locate where these maps are #
mapFiles<-list.files("tmp4/", full.names=TRUE)
speciesSelect<-sub("^.*/([^_]+)_.*$", "\\1", input_file1) # Determine species
mapFilesSpecies<-mapFiles[grep(speciesSelect, mapFiles)] # subset to species of interest
mapFilesRaster<-mapFilesSpecies[grep("rasters_gam3_v2", mapFilesSpecies)] # select raster files # v2 is the Gams based on gridded data and gam3 is population-specific predictions
mapFilesMonth<-mapFilesRaster[grep(paste0("_", "m", Months[i], "_"), mapFilesRaster)] # subset to specific month

# To do this, we select a behaviour value at random from the CIs and adjust the last one accordingly
speciesSub<-actBudgets$species[1]
if (speciesSub %in% c("Black-legged kittiwake", "Northern fulmar")) {
behaviors<-c("Flight", "Forage", "RestWater", "Land") 
} else {
behaviors<-c("Flight", "Active", "RestWater", "Land") 
}
behaviors_random<-sample(behaviors, 4, replace=FALSE) # Choose at random the order in which these will be determined

##### I got to there: now I need to figure out how to assign behaviors from these rasters in a good way :) #####

round<-0

for (b in behaviors_random) {

round<-round + 1

# Subset to raster of behavior b
beh_b<-mapFilesMonth[grep(b, mapFilesMonth)]

# open raster list
beh_rast<-readRDS(beh_b)

# Stack together all the iterations for a given colony #
restacked <- lapply(seq_len(nlayers(beh_rast[[1]])), function(i) {
  stack(lapply(beh_rast, raster::subset, i))
})

# Find the correct stack in this list # (i.e. the one corresponding to our population of interest )
allNames<-unlist(lapply(restacked, names)) # Determine all population names
allNames2<-gsub("\\.\\d+$", "", allNames) # Remove the number at the end of each name (corresponds to layer number in raster stack)
allNames3<-gsub("\\.", " ", allNames2) # replace dot between two words with a space
populations<-data.frame(colony=unique(allNames3)) # Create a data frame with unique population names
populations$listItem<-1:nrow(populations) # Add number in the list so I can easily extract
populationsSub<-subset(populations, colony==actBudgets$colony[1]) # Subset to colony of interest
populationStack<-restacked[[populationsSub$listItem[1]]] # Extract this stack from the list
populationRast_mean<-overlay(populationStack, fun="mean") # Calculate mean accross iterations
populationRast_sd<-overlay(populationStack, fun="sd") # Calculate sd accross iterations
populationRast_se<-populationRast_sd/sqrt(nlayers(populationStack)) # Calculate se accross iterations

# Rasterize, project, and re-sample to match Per's maps

# Mean first
beh_proj<-raster::projectRaster(populationRast_mean, sstMean)
beh_proj_df<-as.data.frame(beh_proj, xy=TRUE)

# Sd second
beh_proj2<-raster::projectRaster(populationRast_sd, sstMean)
beh_proj_df2<-as.data.frame(beh_proj2, xy=TRUE)

# Determine lower & upper limits for sampling behavior
beh_rast<-beh_proj_df 
beh_rast$sd<-beh_proj_df2$layer
beh_rast$col_lower<-beh_rast$layer-beh_rast$sd
beh_rast$col_upper<-beh_rast$layer+beh_rast$sd
 
# Sample value at random from this range for every cell 
  beh_rast[[paste0("prop", b)]] <- runif(
  nrow(beh_rast),
  min = beh_rast$col_lower,
  max = beh_rast$col_upper
)

# Make sure values are not below 0 or above 1
beh_rast[[paste0("prop", b)]]<-ifelse(beh_rast[[paste0("prop", b)]]<0, 0, beh_rast[[paste0("prop", b)]])
beh_rast[[paste0("prop", b)]]<-ifelse(beh_rast[[paste0("prop", b)]]>1, 1, beh_rast[[paste0("prop", b)]])

if (round ==1) {

beh_final<-beh_rast

} else {

beh_final[[paste0("prop", b)]]<-beh_rast[[paste0("prop", b)]]

}

}

# Select columns of interest (just xy & behaviours)
meanBeh<-beh_final[,c(1:2, 7:10)]

if(speciesSub %in% c("Black-legged kittiwake", "Northern fulmar")) {
meanBeh$propActive<-0
} else {
meanBeh$propForage<-0
}

# make sure nothing is negative
	meanBeh$propForage<-ifelse(meanBeh$propForage<0, 0, meanBeh$propForage)
	meanBeh$propFlight<-ifelse(meanBeh$propFlight<0, 0, meanBeh$propFlight)
	meanBeh$propActive<-ifelse(meanBeh$propActive<0, 0, meanBeh$propActive)
	meanBeh$propLand<-ifelse(meanBeh$propLand<0, 0, meanBeh$propLand)
	meanBeh$propRestWater<-ifelse(meanBeh$propRestWater<0, 0, meanBeh$propRestWater)

# Change value of last column so it adds up to 1
meanBeh[,6]<-1-rowSums(meanBeh[,3:5])

# Figure out if there are any negative rows
negativeRows<-subset(meanBeh, meanBeh[[paste0("prop", behaviors_random[4])]]<0)

if (nrow(negativeRows)>0) {

rowChange<-3
rowOld<-4

repeat {

meanBeh[[paste0("prop", behaviors_random[rowChange])]]<-ifelse(meanBeh[[paste0("prop", behaviors_random[rowOld])]]<0, meanBeh[[paste0("prop", behaviors_random[rowChange])]] + meanBeh[[paste0("prop", behaviors_random[rowOld])]], meanBeh[[paste0("prop", behaviors_random[rowChange])]])
meanBeh[[paste0("prop", behaviors_random[rowOld])]]<-ifelse(meanBeh[[paste0("prop", behaviors_random[rowOld])]]<0, 0, meanBeh[[paste0("prop", behaviors_random[rowOld])]])

negativeRows<-subset(meanBeh, meanBeh[[paste0("prop", behaviors_random[rowChange])]]<0)

if (nrow(negativeRows) < 1) {
break}

# Change which behaviours we are altering
rowChange<-rowChange - 1
rowOld<-rowOld - 1

if (rowChange==0) {
break}

}

}

# Calculate total hours per row & make sure none are bigger than they should be (at least they should all be equal)
meanBeh$totProp<-rowSums(meanBeh[,3:7])

# Change proportions into hours
hoursMonth<-c(30*24, 31*24, 30*24, 31*24, 31*24, 28*24, 31*24, 30*24) # determine hours in months 9-2
names(meanBeh) <- gsub("^prop", "t", names(meanBeh)) # Change prop to t to match energy functions
names(meanBeh) <- paste0(names(meanBeh), "_month")
hoursSub<-hoursMonth[i] # select length of month i 
meanBeh[,3:7]<-hoursSub*meanBeh[,3:7]
 
### Now we randomize the energy expenditure parameters ###
print("Randomizing energy parameters...")

actRes<-meanBeh # Allocate activity Sum to ActRes because that fits with previous code
actRes$species<-speciesSub
actRes$colony<-actMonth$colony[1]
actRes$session_id<-1
actRes$rep<-j
actRes$individ_id<-"Pop"
actRes$month<-Months[i]
modelParamsSub<-subset(modelParams, species==actRes$species[1]) # Subset parameter data frame to species of interest

if (actRes$species[1] %in% c("Northern fulmar", "Black-legged kittiwake")) {

actRes$c2<-mean(sample(subset(modelParamsSub, parameter=="c2")$values[[1]], 1), sample(c(seq(2, 5.7, 0.1)) ,1))
actRes$c5<-0

} else {

actRes$c2<-sample(subset(modelParamsSub, parameter=="c2")$values[[1]], 1) # choose at random from uniform distribution
actRes$c5<-mean(sample(subset(modelParamsSub, parameter=="c5")$values[[1]][[1]], 1), sample(subset(modelParamsSub, parameter=="c5")$values[[1]][[2]], 1)) # mean of two numbers

}

actRes$RMR<-sample(subset(modelParamsSub, parameter=="RMR")$values[[1]], 1) # choose at random from uniform distribution
actRes$c1<-sample(subset(modelParamsSub, parameter=="c1")$values[[1]], 1) # choose at random from uniform distribution
actRes$c3<-sample(subset(modelParamsSub, parameter=="c3")$values[[1]], 1) # choose at random from uniform distribution
actRes$c4<-sample(subset(modelParamsSub, parameter=="c4")$values[[1]], 1) # Choose at random from uniform distribution
actRes$TC_air<-subset(modelParamsSub, parameter=="TC_air")$values[[1]] # Choose at random from uniform distribution
actRes$TC_water<-subset(modelParamsSub, parameter=="TC_water")$values[[1]] # Choose at random from uniform distribution
actRes$Beta_active<-subset(modelParamsSub, parameter=="Beta_active")$values[[1]] # Choose at random from uniform distribution
actRes$Beta_rest<-subset(modelParamsSub, parameter=="Beta_rest")$values[[1]] # Choose at random from uniform distribution
actRes$LCT_air<-subset(modelParamsSub, parameter=="LCT_air")$values[[1]] # Choose the number
actRes$LCT_water<-subset(modelParamsSub, parameter=="LCT_water")$values[[1]] # Choose the number

### Calculate energetics (monthly) ####
print("Calculating energetics")

# We will stack the two environmental datas
envStack<-stack(sstMean, airMean_crop)

energyMap<-calculateEnergetics(species=actRes$species[1], data=actRes, colonySub=actRes$colony[1], sstVals=envStack, type="daily", type2="map")
energyMap$month<-actRes$month[1]

### Draw a random number of birds value ###
print("Drawing random birds number...")

NoBirds <- meanDensity %>%
  dplyr::left_join(meanDensity_lower, by = c("x", "y")) %>%
  dplyr::left_join(meanDensity_upper, by = c("x", "y")) %>%
  replace_na(list(NoBirds_ci_lower=0, NoBirds_ci_upper=0)) %>%
  dplyr::mutate(
    
    # sample random integer per row
    randomBirds = ifelse(
      NoBirds_ci_upper >= NoBirds_ci_lower,
      floor(runif(n(), NoBirds_ci_lower, NoBirds_ci_upper + 1)),
      0  # or NA_real_ if you prefer
    )
  ) %>%
  dplyr::select(x, y, NoBirds_mean, randomBirds)

#### Join the two ####

if (j==1) {

PopEnergy<-NoBirds %>%
dplyr::left_join(energyMap, by=c("x", "y")) %>%
dplyr::select(x, y, species, month, colony, weight, sst, temp, NoBirds_mean, randomBirds, DEEkJ) %>%
#dplyr::mutate(energyPopkJ=randomBirds*DEEkJ)
dplyr::mutate(energyPopkJ=NoBirds_mean*DEEkJ)

IndEnergy<-NoBirds %>%
dplyr::left_join(energyMap, by=c("x", "y")) %>%
dplyr::select(x, y, species, month, colony, weight, sst, temp, DEEkJ)

} else {

# This is a base layer
PopEnergy2<-NoBirds %>%
dplyr::left_join(energyMap, by=c("x", "y")) %>%
dplyr::select(x, y, species, colony, weight, sst, temp, NoBirds_mean, DEEkJ) %>%
#dplyr::mutate(energyPopkJ=randomBirds*DEEkJ)
dplyr::mutate(energyPopkJ=NoBirds_mean*DEEkJ)
colname <- paste0("energyPopkJ_", j)   # creates "col_2025"
PopEnergy[[colname]] <- PopEnergy2$energyPopkJ

IndEnergy2<-NoBirds %>%
dplyr::left_join(energyMap, by=c("x", "y")) %>%
dplyr::select(x, y, species, colony, weight, sst, temp, DEEkJ) 
colname <- paste0("DEEkJ_", j)   # creates "col_2025"
IndEnergy[[colname]] <- IndEnergy2$DEEkJ

}

}

# Summarize result and make one raster layer with mean energy, lower, upper # (multplied by number of birds) per month 
PopEnergyFinal<-PopEnergy %>%
dplyr::select(x, y, species, month, colony, weight, sst, temp, NoBirds_mean)

PopEnergyFinal$energyPopkJ_mean<-rowMeans(PopEnergy[,12:ncol(PopEnergy)])
PopEnergyFinal$energyPopkJ_sd <- apply(PopEnergy[, 12:ncol(PopEnergy)], 1, sd)
PopEnergyFinal$energyPopkJ_se<-PopEnergyFinal$energyPopkJ_sd/sqrt(overall.iterations)

# Summarize the same for energy without multiplying by numbers of birds...
IndEnergyFinal<-IndEnergy %>%
dplyr::select(x, y, species, month, colony, weight, sst)

IndEnergyFinal$DEEkJ_mean<-rowMeans(IndEnergy[,9:ncol(IndEnergy)])
IndEnergyFinal$DEEkJ_sd <- apply(IndEnergy[, 9:ncol(IndEnergy)], 1, sd)

# Save result in a monthly list
PopEnergyFinal$DEEkJ_mean<-IndEnergyFinal$DEEkJ_mean
PopEnergyFinal$DEEkJ_sd<-IndEnergyFinal$DEEkJ_sd
PopEnergyFinal_All<-rbind(PopEnergyFinal_All, PopEnergyFinal)

}

PopEnergyFinal_All$colony<-colonyName$colonies[1]
PopEnergyFinal_All$modelcolony<-colonyName$modelcolonies[1]

# Save results #
speciesSave<-gsub(" ", "", speciesSub)
speciesSave<-gsub("-", "", speciesSave)
speciesSave<-gsub("ü", "u", speciesSave)
speciesSave<-gsub("'", "", speciesSave)
speciesSave<-gsub("`", "", speciesSave)

# Make a special case for a population name which is cray cray
startName<-substr(PopEnergyFinal_All$colony[1], 1, 7)

if (startName =="Rodolph") {
popSave<-"rodolphisland"
} else {
popSave<-gsub("_", "", PopEnergyFinal_All$colony[1])
}
popSave<-gsub("ö", "o", popSave)
popSave<-gsub("ø", "o", popSave)
popSave<-gsub("Å", "Aa", popSave)
popSave<-gsub("æ", "ae", popSave)
popSave<-gsub(" ", "", popSave)
popSave<-gsub("-", "", popSave)
popSave<-gsub(".,", "", popSave)
popSave<-gsub(",", "", popSave)
popSave<-gsub("/", "", popSave)
popSave<-gsub("\\)", "", popSave)
popSave<-gsub("\\(", "", popSave)
popSave<-gsub("–", "", popSave)
popSave<-gsub("ý", "", popSave)
popSave<-gsub("ð", "", popSave)
popSave<-gsub("á", "", popSave)
popSave<-gsub("ó", "", popSave)
popSave<-gsub(":", "", popSave)
popSave<-gsub("í", "", popSave)
popSave<-gsub("`", "", popSave)
popSave<-gsub("ú", "", popSave)
popSave<-gsub("é", "", popSave)
popSave<-gsub("'", "", popSave)
popSave<-tolower(popSave)
print(popSave)

write.csv(PopEnergyFinal_All, file=paste0("tmp5/", speciesSave, "_", popSave, "_energyMap_v2.csv"), row.names=FALSE)

# Save metadata #
results<-data.frame(species=speciesSub, colony=PopEnergyFinal_All$colony[1], modelcolony=PopEnergyFinal_All$modelcolony[1])
allResults<-rbind(allResults, results)

}

print("Saving results...")

# Save result
output_file1 <- args[2]

# Save metadata as main output file
write.csv(allResults, file=output_file1, row.names=FALSE)
