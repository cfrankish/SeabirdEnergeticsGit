## This script maps activity & energy expenditure for our species ##
## To do this, by pop, the script will open the corresponding distribution of activity budgets ##
## For energy, it will extract the correct SST map & calculate energy 100 times (in line with the pop maps too?) ##
## Then it will map activity budgets somewhow too ##
## The output will be five rasters per pop ##

### Here we use SST projections (year 2050) from COPERNICUS
# CMIP6 climate projections
# HadGEM3-GC31-MM (UK) 
# Downloaded 27/01/26

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

# Set-up number of iterations...
overall.iterations<-5 # how many times this is calculated per individual
print(paste0("Determining activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1] # This will read in a population-specific activity budget  

print(input_file1)

# open input file 1
actBudgets<-read.csv(input_file1) # actBudgets<-read.csv("tmp3/Blackleggedkittiwake_faroeislands_actbudgets_monthly.csv")

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
print(colonyFinal)

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

# Determine location of SST data (used in Per's maps) so we can align the resolution #
sstMonthly<-list.files("./data/sstPopMaps", full.names=TRUE)
sstMaps<-sstMonthly[grepl("tif", sstMonthly)] # Choose the tif file
sstMap<-raster::stack(sstMaps[1]) # Stack all layers (lol)
sstSub<-subset(sstMap, 1)

# Find predicted layers #
sstMonthly_predict<-list.files("./data/sst_future", full.names=TRUE)
sstMaps_predict<-nc_open(sstMonthly_predict[1]) # Stack all layers (lol)

# Decode (so I now which layers correspond to which months #
t_num   <- ncvar_get(sstMaps_predict, "time")
t_units <- ncatt_get(sstMaps_predict, "time", "units")$value
t_cal   <- ncatt_get(sstMaps_predict, "time", "calendar")$value  # may be NULL/NA
origin <- as.Date(sub("days since ", "", t_units))
d <- t_num
year  <- as.integer(format(origin, "%Y")) + (d %/% 360)
d2    <- d %% 360
month <- 1 + (d2 %/% 30)
day   <- 1 + (d2 %% 30)
ymd <- sprintf("%04d-%02d-%02d", year, month, day)
ym  <- sprintf("%04d-%02d", year, month)
sstPredictKey<-data.frame(month=as.numeric(substr(ymd, 6, 7)))
sstPredictKey$layer<-1:nrow(sstPredictKey)

# Determine parameters needed for energy calculations
speciesNo<-6
paramNo<-16
modelParams<-tibble(species=rep(c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot"), paramNo))
modelParams$parameter<-rep(c("L1", "Th1", "Th2", "L1_colony", "dist_colony", "pLand", "c",
"RMR", "c1", "c2", "c3", "c4", "c5", "TC", "Beta_active", "Beta_rest"), each=6)
modelParams$values<-list(
240, 810, 90, 134, 88, 88, # Species-specific flight bout duration (minutes)
c(seq(0.95, 0.98, 0.01)), c(seq(0.95, 0.98, 0.01)), c(seq(0.84, 0.93, 0.01)), c(seq(0.84, 0.93, 0.01)), c(seq(0.84, 0.93, 0.01)), c(seq(0.84, 0.93, 0.01)), # Th1 % wet threshold for differenciating between behaviors
0, 0, 0, 0, 0, 0, # Th2 #% wet for determing dry vs. intermediate
c(0, 240), c(0, 810), c(0, 90), c(0, 134), c(0, 88), c(0, 88), # L1_colony: duration of dry bouts at start of tLand that can re-allocated to flight
500, 500, 500, 500, 500, 500,  # Distance to colony (km) below which it is considered possible to be on land
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, # pland: probability of dry being land or something else
c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), # coefficient for adjusting for leg-tucking
c(seq(1.60, 1.68, 0.01)), c(seq(0.96, 1.04, 0.01)), c(0), c(0), c(0), c(0), # RMR is just for a few of the species...BlKi / NoFu (Resting metabolic rate)
c(seq(2, 5.7, 0.1)), c(2.2), c(seq(106, 176, 1)), c(seq(5.4, 11.4, 0.1)), c(seq(106, 176, 1)), c(seq(106, 176, 1)), # c1 is the cost of Flight
c(seq(0.3, 3.8, 0.1)), c(seq(0.3, 3.8, 0.1)), c(0), c(0), c(0), c(0), # c2 is the cost of foraging...
c(seq(0.1, 1.1, 0.1)), c(0.8), c(seq(2.7, 15.3, 0.1)), c(seq(2.7, 15.3, 0.1)), c(seq(2.7, 15.3, 0.1)), c(seq(2.7, 15.3, 0.1)), # c3 is cost of being on land...
c(seq(0.1, 2.8, 0.1)), c(2), c(seq(2.7, 15.3, 0.1)), c(seq(2.7, 15.3, 0.1)), c(seq(2.7, 15.3, 0.1)), c(seq(2.7, 15.3, 0.1)), # c4 is the cost of resting on the water...
c(0), c(0), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), # c5 is the cost of being active on water when thermoneutral...
c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), 
c(seq(0.0997-0.1*0.0997, 0.0997 + 0.1*0.0997, 0.01)), c(seq(0.07-0.1*0.07, 0.07 + 0.1*0.07, 0.01)), # TC is thermal conductivity (code is TC - 1.96*TCError)
c(seq(2.75-0.1*2.75, 2.75 + 0.1*2.75, 0.01)), c(seq(2.75-0.1*2.75, 2.75 + 0.1*2.75, 0.01)),
c(seq(2.75-0.1*2.75, 2.75 + 0.1*2.75, 0.01)), c(seq(2.75-0.1*2.75, 2.75 + 0.1*2.75, 0.01)),
c(0), c(0), c(seq(118-118*0.1, 118 + 118*0.1, 1)),  # Intercepts of resting metabolic rate at 0°C during different activities (active)
c(seq(118-118*0.1, 118 + 118*0.1, 1)), c(seq(118-118*0.1, 118 + 118*0.1, 1)), c(seq(118-118*0.1, 118 + 118*0.1, 1)),
c(seq(1.87-1.87*0.1, 1.87 + 1.87*0.1, 1)), c(seq(1.34-1.34*0.1, 1.34+1.34*0.1, 1)),# # Intercepts of resting metabolic rate at 0°C during different activities (rest)
c(seq(72.2-72.2*0.1, 72.2 + 72.2*0.1, 0.1)), c(seq(72.2-72.2*0.1, 72.2 + 72.2*0.1, 0.1)),
c(seq(72.2-72.2*0.1, 72.2 + 72.2*0.1, 0.1)), c(seq(72.2-72.2*0.1, 72.2 + 72.2*0.1, 0.1)))

### Step 2: Calculate energy & map it many times ####

# Determine number of model colonies #
modelColony_no<-nrow(metaAll)

# Make a list to save output file in
allResults<-list()

for (p in 1:nrow(metaAll)) {

# Determine name of model colony for later
colonyName<-metaAll[p,]

# Determine number of months to loop through
Months<-unique(actBudgets$month)

# Make an an empty list to save the results in
PopEnergyFinal_All<-list()

for (i in 1:length(Months)) {

# Subset to month i
monthSub<-Months[i]

# Prepare corresponding SST surface
sstKeySub<-subset(sstPredictKey, month == monthSub)
sstSub <- raster(sstMonthly_predict[1], varname="tos", band=sstKeySub$layer) # SST
lon <- ncvar_get(sstMaps_predict, "longitude")
lat <- ncvar_get(sstMaps_predict, "latitude")
sstMean<-raster::calc(sstSub, fun = mean, na.rm = TRUE)

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

print(paste0("Mapping energy for month ", i, " : iteration ", j))

### We generate activity budget from distribution ###

print("Randomising individual activity budget...")

# To do this, we select a behaviour value at random from the CIs and adjust the last one accordingly
speciesSub<-actBudgets$species[1]
if (speciesSub %in% c("Black-legged kittiwake", "Northern fulmar")) {
behaviors<-c("tFlight", "tForage", "tRestWater", "tLand") 
actMonth$tActive<-0
} else {
behaviors<-c("tFlight", "tActive", "tRestWater", "tLand") 
actMonth$tForage<-0
}
behaviors_random<-sample(behaviors, 4, replace=FALSE) # Choose at random the order in which these will be determined

for (b in behaviors_random) {
  col_lower <- paste0(b, "_lower")
  col_upper <- paste0(b, "_upper")
  
  actMonth[[b]] <- runif(
    nrow(actMonth),
    min = actMonth[[col_lower]],
    max = actMonth[[col_upper]]
  )
}

# Change the value of the last behavior so they add up to something that is correct
meanBeh<-actMonth %>%
dplyr::select(month, tForage, tFlight, tActive, tRestWater, tLand, totHrs)

# make sure nothing is negative
meanBeh$tForage<-ifelse(meanBeh$tForage<0, 0, meanBeh$tForage)
meanBeh$tFlight<-ifelse(meanBeh$tFlight<0, 0, meanBeh$tFlight)
meanBeh$tActive<-ifelse(meanBeh$tActive<0, 0, meanBeh$tActive)
meanBeh$tLand<-ifelse(meanBeh$tLand<0, 0, meanBeh$tLand)
meanBeh$tRestWater<-ifelse(meanBeh$tRestWater<0, 0, meanBeh$tRestWater)

meanBeh[[behaviors_random[4]]]<-meanBeh$totHrs-rowSums(meanBeh[, behaviors_random[1:3]]) # We have a problem... and that is that some columns are negative meaning we need to remove hours from another behaviour

# Figure out if there are any negative rows
negativeRows<-subset(meanBeh, meanBeh[[behaviors_random[4]]]<0)

if (nrow(negativeRows)>0) {

rowChange<-3
rowOld<-4

repeat {

meanBeh[[behaviors_random[rowChange]]]<-ifelse(meanBeh[[behaviors_random[rowOld]]]<0, meanBeh[[behaviors_random[rowChange]]] + meanBeh[[behaviors_random[rowOld]]], meanBeh[[behaviors_random[rowChange]]])
meanBeh[[behaviors_random[rowOld]]]<-ifelse(meanBeh[[behaviors_random[rowOld]]]<0, 0, meanBeh[[behaviors_random[rowOld]]])

negativeRows<-subset(meanBeh, meanBeh[[behaviors_random[rowChange]]]<0)

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
totHrs<-meanBeh$totHrs[1]
meanBeh$totHrs<-rowSums(meanBeh[,2:6])
meanBeh$match<-ifelse(all.equal(meanBeh$totHrs, meanBeh$totHrs, tolerance=0.2), 1, 0)
nomatch<-subset(meanBeh, match==0)
if(nrow(nomatch)>0) stop(print("error: randomization problem"))

# Extract mean activity values at pop-level & make sure they match the total hours in a given month
meanPop<-meanBeh %>%
  group_by(month) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# Calculate new total hours
meanPop$totHrs_new<-rowSums(meanPop[,2:6])

# make sure these behaviors sum to correct number (otherwise there are incomplete months #
sumBehNoMatch<-subset(meanPop, !all.equal(totHrs, totHrs_new, tolerance=0.2))
if (nrow(sumBehNoMatch) >0) stop(print("Error: activity hours weird"))
  
### Now we randomize the energy expenditure parameters ###
print("Randomizing energy parameters...")

actRes<-meanPop # Allocate activity Sum to ActRes because that fits with previous code
actRes$species<-speciesSub
actRes$colony<-actMonth$colony[1]
actRes$session_id<-1
actRes$rep<-j
actRes$individ_id<-"Pop"
modelParamsSub<-subset(modelParams, species==actRes$species[1]) # Subset parameter data frame to species of interest

if (actRes$species[1] %in% c("Northern fulmar", "Black-legged kittiwake")) {

actRes$c2<-mean(sample(subset(modelParamsSub, parameter=="c2")$values[[1]], 1), sample(c(seq(2, 5.7, 0.1)) ,1))
#actRes$c5<-0

} else {

actRes$c2<-sample(subset(modelParamsSub, parameter=="c2")$values[[1]], 1) # choose at random from uniform distribution
#actRes$c5<-mean(sample(subset(modelParamsSub, parameter=="c5")$values[[1]][[1]], 1), sample(subset(modelParamsSub, parameter=="c5")$values[[1]][[2]], 1)) # mean of two numbers

}

actRes$RMR<-sample(subset(modelParamsSub, parameter=="RMR")$values[[1]], 1) # choose at random from uniform distribution
actRes$c1<-sample(subset(modelParamsSub, parameter=="c1")$values[[1]], 1) # choose at random from uniform distribution
actRes$c3<-sample(subset(modelParamsSub, parameter=="c3")$values[[1]], 1) # choose at random from uniform distribution
actRes$c4<-sample(subset(modelParamsSub, parameter=="c4")$values[[1]], 1) # Choose at random from uniform distribution
actRes$TC<-sample(subset(modelParamsSub, parameter=="TC")$values[[1]], 1) # Choose at random from uniform distribution
actRes$Beta_active<-sample(subset(modelParamsSub, parameter=="Beta_active")$values[[1]], 1) # Choose at random from uniform distribution
actRes$Beta_rest<-sample(subset(modelParamsSub, parameter=="Beta_rest")$values[[1]], 1) # Choose at random from uniform distribution

### Calculate energetics (monthly) ####
print("Calculating energetics")

actRes <- actRes %>%
  rename_with(~ paste0(.x, "_month"), .cols = c(tForage, tRestWater, tActive, tFlight, tLand))

energyMap<-calculateEnergetics(species=actRes$species[1], data=actRes, colonySub=actRes$colony[1], sstVals=sstMean, type="daily", type2="map")
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
dplyr::select(x, y, species, month, colony, weight, sst, NoBirds_mean, randomBirds, DEEkJ) %>%
#dplyr::mutate(energyPopkJ=randomBirds*DEEkJ)
dplyr::mutate(energyPopkJ=NoBirds_mean*DEEkJ) # Multiply by mean numbers of birds so I can compare to other method... 

IndEnergy<-NoBirds %>%
dplyr::left_join(energyMap, by=c("x", "y")) %>%
dplyr::select(x, y, species, month, colony, weight, sst, DEEkJ) 

} else {

# This is a base layer
PopEnergy2<-NoBirds %>%
dplyr::left_join(energyMap, by=c("x", "y")) %>%
dplyr::select(x, y, species, colony, weight, sst, NoBirds_mean, DEEkJ) %>%
#dplyr::mutate(energyPopkJ=randomBirds*DEEkJ)
dplyr::mutate(energyPopkJ=NoBirds_mean*DEEkJ)
colname <- paste0("energyPopkJ_", j)   # creates "col_2025"
PopEnergy[[colname]] <- PopEnergy2$energyPopkJ

IndEnergy2<-NoBirds %>%
dplyr::left_join(energyMap, by=c("x", "y")) %>%
dplyr::select(x, y, species, colony, weight, sst, DEEkJ) 
colname <- paste0("DEEkJ_", j)   # creates "col_2025"
IndEnergy[[colname]] <- IndEnergy2$DEEkJ

}

}

# Summarize result and make one raster layer with mean energy, lower, upper # (multplied by number of birds) per month 
PopEnergyFinal<-PopEnergy %>%
dplyr::select(x, y, species, month, colony, weight, sst, NoBirds_mean)

PopEnergyFinal$energyPopkJ_mean<-rowMeans(PopEnergy[,11:ncol(PopEnergy)])
PopEnergyFinal$energyPopkJ_sd <- apply(PopEnergy[, 11:ncol(PopEnergy)], 1, sd)
PopEnergyFinal$energyPopkJ_se<-PopEnergyFinal$energyPopkJ_sd/sqrt(overall.iterations)

# Summarize the same for energy without multiplying by numbers of birds...
IndEnergyFinal<-IndEnergy %>%
dplyr::select(x, y, species, month, colony, weight, sst)

IndEnergyFinal$DEEkJ_mean<-rowMeans(IndEnergy[,8:ncol(IndEnergy)])
IndEnergyFinal$DEEkJ_sd <- apply(IndEnergy[, 8:ncol(IndEnergy)], 1, sd)

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

write.csv(PopEnergyFinal_All, file=paste0("tmp5/", speciesSave, "_", popSave, "_energyMap_v1.csv"), row.names=FALSE)

# Save metadata #
results<-data.frame(species=speciesSub, colony=PopEnergyFinal_All$colony[1], modelcolony=PopEnergyFinal_All$modelcolony[1])
allResults<-rbind(allResults, results)

}

print("Saving results...")

# Save result
output_file1 <- args[2]

# Save metadata as main output file
write.csv(allResults, file=output_file1, row.names=FALSE)