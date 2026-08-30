## Here we just try & visualize whether we will need to spatially map behaviours ##
## Before calculating energetics ##

library(dplyr)
library(ncdf4)
library(raster)
library(sp)
library(sf)
library(lubridate)
library(data.table)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggplot2)
library(gridExtra)
library(biscale)
library(cowplot)
library(classInt)
library(spdep)

#### Step 0: setting up basic conditions ####

# Set-up number of iterations...
overall.iterations<-50 # how many times this is calculated per individual
print(paste0("Mapping activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1] # This will read in a species name
print(input_file1)
#input_file1<-"Commonguillemot"

# Now we determine the list of activity/energy data to loop through
allFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp4/", full.names=TRUE)
speciesFiles<-allFiles[grep(input_file1, allFiles)] # Subset to species-specific ones
actFiles<-speciesFiles[grep("activityMap.csv", speciesFiles)] # extract activity files
energyFiles<-speciesFiles[grep("energy", speciesFiles)] # extract energy files

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# open colony coordinates to put on maps
colony.summary<-readRDS("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")
colonyCoords<-colony.summary %>%
dplyr::ungroup() %>%
dplyr::group_by(colony) %>%
dplyr::slice(1) %>%
dplyr::select(colony, col_lat, col_lon)
coordinates(colonyCoords)<-~col_lon + col_lat
proj4string(colonyCoords)<-projection_84
colonyCoordsTrans<-data.frame(spTransform(colonyCoords, projection_NA))

# Determine model colonies for activity plotting #
speciesMatch<-data.frame(speciesLatin=c("Alle_alle", "Fratercula_arctica", "Fulmarus_glacialis", "Rissa_tridactyla", "Uria_lomvia", "Uria_aalge"), species=c("Littleauk", "Atlanticpuffin", "Northernfulmar",
"Blackleggedkittiwake", "Brunnichsguillemot", "Commonguillemot"))
speciesSub<-subset(speciesMatch, species==input_file1)
speciesSelect<-paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/data/popdata_raw/SEATRACK_Abundance_Model_", speciesSub$speciesLatin, "_Ver_3_1.nc")
nc<-nc_open(speciesSelect) # nc<-nc_open("data/popdata_raw/SEATRACK_Abundance_Model_Rissa_tridactyla_Ver_3_1.nc")
modelcolonies <- (ncvar_get(nc,"SmcolName")) # Find list of populations
colonies <- (ncvar_get(nc,"colonyName")) # Find list of ALL colonies
colcode <- (ncvar_get(nc,"SmcolCode")) # Colony code
meta <- data.frame(modelcolonies, colonies, colcode) # Make some metadata so easier to understand raster structure

#### Step 0: Figure out largest extent with birds for comparisons ###

# This step is made so that I can subset my map to a smaller one, accross which I can make comparison's of moran's I and other things #

print("Step 0: figuring out a map accross which to conduct comparisons...")

# Here I estimate the min & max coordinates 
extentAll<-mapExtent(actFiles)

# Open up main dataset and provide an index #
actSub<-fread(actFiles[1])
actSub<-subset(actSub, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))
actSubIndex<-actSub %>%
ungroup() %>%
dplyr::group_by(month) %>%
arrange(x, y) %>%
dplyr::mutate(index=row_number()) 

# Attach extent all file & create an index for x & y coordinates so I can easily filter future maps

# 1. Unique sorted x and y lists
actSubIndex <- actSubIndex %>%
  dplyr::mutate(
    x_clean = round(x, 0),   # round to nearest meter
    y_clean = round(y, 0)
  )

x_index <- actSubIndex %>%
ungroup() %>%
  dplyr::distinct(x_clean) %>%
  arrange(x_clean) %>%
  dplyr::mutate(indexX = row_number())

y_index <- actSubIndex %>%
ungroup() %>%
  distinct(y_clean) %>%
  arrange(y_clean) %>%
  mutate(indexY = row_number())

# 2. Join back
actSubIndex2 <- actSubIndex %>%
  left_join(x_index, by = "x_clean") %>%
  left_join(y_index, by = "y_clean") %>%
  dplyr::left_join(extentAll, by=c("index"))

# Determine min/max x y coords
extents<-subset(actSubIndex2, !is.na(type))
minY<-min(extents$indexY)
minX<-min(extents$indexX)
maxY<-max(extents$indexY)
maxX<-max(extents$indexX)

#### Step 1: Create species-weighted mean for activity (add SD later) ####

modelColonies<-list()
colonyLox<-list()
pairwisepop_res<-list()

print("Starting for loop")

for (i in 1:length(actFiles)) {

print(paste0("Mapping file", i, "/", length(actFiles), "..."))

#### Open file i ###
actSub<-read.csv(actFiles[i])
actSub<-subset(actSub, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))
energySub<-read.csv(energyFiles[i], nrow=10)

# Calculate month-specific weights
actSubWeights<-actSub %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(totalBirds) %>%
dplyr::mutate(weight=1/totalBirds)

# We number the rows so we can join data frames later ( doesn't work so well with floating numbers)
actSub<-actSub %>%
  ungroup() %>%
  arrange(month, x, y) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(index=row_number()) %>%
  ungroup() 

# We also number the x/y coordinates #  
actSubIndex <- actSub %>%
  dplyr::mutate(
    x_clean = round(x, 0),   # round to nearest meter
    y_clean = round(y, 0)
  )

x_index <- actSubIndex %>%
ungroup() %>%
  dplyr::distinct(x_clean) %>%
  arrange(x_clean) %>%
  dplyr::mutate(indexX = row_number())

y_index <- actSubIndex %>%
ungroup() %>%
  distinct(y_clean) %>%
  arrange(y_clean) %>%
  mutate(indexY = row_number())

# 2. Join back
actSubIndex2 <- actSubIndex %>%
  left_join(x_index, by = "x_clean") %>%
  left_join(y_index, by = "y_clean")

# Determine the model colony
metaSub<-subset(meta, colonies==actSub$colony[1])
modelColSub<-metaSub$modelcolonies[1]
print(modelColSub)

# Does the list already have this model colony?
if (modelColSub %in% modelColonies) {
print("Next")
next
}

# Can we match it with colonyCoords?
colonySub<-subset(colonyCoordsTrans, colony==actSub$colony_og[1])

if(nrow(colonySub<1)) {
colonyCoordsTrans$colony<-gsub("Bjørnøya", "Bjornoya", colonyCoordsTrans$colony)
colonySub<-subset(colonyCoordsTrans, colony==actSub$colony_og[1])
}

# Otherwise add the modelcolony name to the list
modelColonies<-append(modelColonies, modelColSub)

# Plot population-specific behaviors first
#mapFlight<-mapActivity_pop(actSub, colonyCoordsTrans, "Flight", energySub$species[1])
#mapRest<-mapActivity_pop(actSub, colonyCoordsTrans, "RestWater", energySub$species[1])
#if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
#mapForage<-mapActivity_pop(actSub, colonyCoordsTrans, "Forage", energySub$species[1])
#} else {
#mapActive<-mapActivity_pop(actSub, colonyCoordsTrans, "Active", energySub$species[1])
#}

# Susbet dataset to a minimum extent where there are birds accross all populations
actBirds<-subset(actSubIndex2, indexX <= maxX & indexX >= minX & indexY <= maxY & indexY >= minY) # Subset to squares where this data

# Conduct pair-wise correlation calculations #
pairwisePop_Flight<-popcorr(actBirds, actFiles, meta, "Flight")
pairwisePop_RestWater<-popcorr(actBirds, actFiles, meta, "RestWater")
pairwisePop_Land<-popcorr(actBirds, actFiles, meta, "Land")
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
pairwisePop_Forage<-popcorr(actBirds, actFiles, meta, "Forage")
pairwisePop_all<-rbind(pairwisePop_Flight, pairwisePop_RestWater, pairwisePop_Land, pairwisePop_Forage)
} else {
pairwisePop_Active<-popcorr(actBirds, actFiles, meta, "Active")
pairwisePop_all<-rbind(pairwisePop_Flight, pairwisePop_RestWater, pairwisePop_Land, pairwisePop_Active)
}

# Calculate Moran's i
actBirdsMoran_flight<-computeMoran_pop(actBirds, "Flight")
actBirdsMoran_rest<-computeMoran_pop(actBirds, "RestWater")
actBirdsMoran_land<-computeMoran_pop(actBirds, "Land")
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
actBirdsMoran_forage<-computeMoran_pop(actBirds, "Forage")
actBirdsMoran_active<-NA
} else {
actBirdsMoran_active<-computeMoran_pop(actBirds, "Active")
actBirdsMoran_forage<-NA
}

# Calculate correlation between mean & SD & add values back to main data frame
actBirdsCor<-subset(actSubIndex2, NoBirds >1)
actBirdsCorr<-corrMonth(actBirdsCor) # Calculate correlations for different behaviours
actSub<-actSub %>%
dplyr::left_join(actBirdsCorr, by=c("month")) %>%
dplyr::left_join(actBirdsMoran_flight, by=c("month")) %>%
dplyr::left_join(actBirdsMoran_rest, by=c("month")) %>%
dplyr::left_join(actBirdsMoran_land, by=c("month"))

if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
actSub<-actSub %>%
dplyr::left_join(actBirdsMoran_forage, by=c("month")) %>%
dplyr::mutate(Moran_Active=NA)
} else {
actSub<-actSub %>%
dplyr::left_join(actBirdsMoran_active, by=c("month")) %>%
dplyr::mutate(Moran_Forage=NA)
}

# Add min & max of these correlations to get an idea of data spread
actSub$moranFlight_min<-actSub$Moran_Flight
actSub$moranFlight_max<-actSub$Moran_Flight
actSub$moranRestWater_min<-actSub$Moran_RestWater
actSub$moranRestWater_max<-actSub$Moran_RestWater
actSub$moranActive_min<-actSub$Moran_Active
actSub$moranActive_max<-actSub$Moran_Active
actSub$moranForage_min<-actSub$Moran_Forage
actSub$moranForage_max<-actSub$Moran_Forage
actSub$moranLand_min<-actSub$Moran_Land
actSub$moranLand_max<-actSub$Moran_Land

actSub$corFlight_min<-actSub$corFlight
actSub$corFlight_max<-actSub$corFlight
actSub$corRestWater_min<-actSub$corRestWater
actSub$corRestWater_max<-actSub$corRestWater
actSub$corActive_min<-actSub$corActive
actSub$corActive_max<-actSub$corActive
actSub$corForage_min<-actSub$corForage
actSub$corForage_max<-actSub$corForage
actSub$corLand_min<-actSub$corLand
actSub$corLand_max<-actSub$corLand

# We will multiply time spent in different behaviours by 1/no birds
actSub[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]<-actSub[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]/actSub$totalBirds
actSub[c("timeFlight_sd", "timeForage_sd", "timeActive_sd", "timeLand_sd", "timeRestWater_sd")]<-actSub[c("timeFlight_sd", "timeForage_sd", "timeActive_sd", "timeLand_sd", "timeRestWater_sd")]/actSub$totalBirds
actSub[c("corFlight", "corForage", "corActive", "corLand", "corRestWater")]<-actSub[c("corFlight", "corForage", "corActive", "corLand", "corRestWater")]/actSub$totalBirds
actSub[c("Moran_Flight", "Moran_Forage", "Moran_Active", "Moran_Land", "Moran_RestWater")]<-actSub[c("Moran_Flight", "Moran_Forage", "Moran_Active", "Moran_Land", "Moran_RestWater")]/actSub$totalBirds
actSub$weight<-ifelse(actSub$NoBirds>0, 1/(actSub$totalBirds), 0) # According to birds which are in a given square
actSub$weightTotal<-1/(actSub$totalBirds) # all weights
#actSub[is.na(actSub)] <- 0
actSub$species<-input_file1
actSub$colonies<-ifelse(actSub$NoBirds>0, 1, 0)

if (i ==1) {

# Make a master dataset that we will join to
actSubTot<-actSub
actSubTot<-actSubTot %>%
dplyr::select(-c(weight, totalBirds)) %>%
dplyr::left_join(actSubWeights, by=c("month")) %>%
dplyr::select(species, month, index, NoBirds, timeFlight, timeRestWater, timeLand, timeForage, timeActive, timeFlight_sd, timeRestWater_sd, timeLand_sd, timeForage_sd, timeActive_sd, 
corFlight, corRestWater, corLand, corForage, corActive, corFlight_min, corFlight_max, corRestWater_min, corRestWater_max, corLand_min, corLand_max, corForage_min, corForage_max, corActive_min, corActive_max, 
Moran_Flight, Moran_RestWater, Moran_Land, Moran_Forage, Moran_Active, moranFlight_min, moranFlight_max, moranRestWater_min, moranRestWater_max, moranLand_min, moranLand_max, moranForage_min, moranForage_max, 
moranActive_min, moranActive_max, timeTotal, colonies, weight, totalBirds)
actSubTot$weight<-ifelse(actSubTot$NoBirds>0, 1/(actSubTot$totalBirds), 0)
actSubTot$weightTotal<-1/actSubTot$totalBirds

} else {

# We sum everything together (which will be averaged by total number of birds later... #
actSub<-actSub %>%
dplyr::select(species, month, index, NoBirds, timeFlight, timeRestWater, timeLand, timeForage, timeActive, timeFlight_sd, timeRestWater_sd, timeLand_sd, timeForage_sd, timeActive_sd, 
corFlight, corRestWater, corLand, corForage, corActive, corFlight_min, corFlight_max, corRestWater_min, corRestWater_max, corLand_min, corLand_max, corForage_min, corForage_max, corActive_min, corActive_max, 
Moran_Flight, Moran_RestWater, Moran_Land, Moran_Forage, Moran_Active, moranFlight_min, moranFlight_max, moranRestWater_min, moranRestWater_max, moranLand_min, moranLand_max, moranForage_min, moranForage_max, 
moranActive_min, moranActive_max, timeTotal, colonies) %>%
dplyr::left_join(actSubWeights, by=c("month"))%>%
dplyr::select(-n)
actSub$weight<-ifelse(actSub$NoBirds>0, 1/(actSub$totalBirds), 0)
actSub$weightTotal<-1/actSub$totalBirds

# join to main dataset & summarize findings
actSubTot<-rbind(actSubTot, actSub)
actSubTot<-actSubTot %>%
ungroup () %>%
dplyr::group_by(species, month, index) %>%
dplyr::summarise(NoBirds=sum(NoBirds, na.rm=TRUE), timeFlight=sum(timeFlight, na.rm=TRUE), 
timeRestWater=sum(timeRestWater, na.rm=TRUE), timeLand=sum(timeLand, na.rm=TRUE), timeForage=sum(timeForage, na.rm=TRUE), timeActive=sum(timeActive, na.rm=TRUE), 
timeFlight_sd=sum(timeFlight_sd, na.rm=TRUE), 
timeRestWater_sd=sum(timeRestWater_sd, na.rm=TRUE), timeLand_sd=sum(timeLand_sd, na.rm=TRUE), timeForage_sd=sum(timeForage_sd, na.rm=TRUE), timeActive_sd=sum(timeActive_sd, na.rm=TRUE),
corFlight=sum(corFlight), corRestWater=sum(corRestWater), corForage=sum(corForage), corActive=sum(corActive), corLand=sum(corLand), 
corFlight_min=min(corFlight_min), corFlight_max=max(corFlight_max), corActive_min=min(corActive_min), corActive_max=max(corActive_max),
corRestWater_min=min(corRestWater_min), corRestWater_max=max(corRestWater_max), corForage_min=min(corForage_min), corForage_max=max(corForage_max),
corLand_min=min(corLand_min), corLand_max=max(corLand_max), 
Moran_Flight=sum(Moran_Flight), Moran_RestWater=sum(Moran_RestWater), Moran_Forage=sum(Moran_Forage), Moran_Active=sum(Moran_Active), Moran_Land=sum(Moran_Land), 
moranFlight_min=min(moranFlight_min), moranFlight_max=max(moranFlight_max), moranActive_min=min(moranActive_min), moranActive_max=max(moranActive_max),
moranRestWater_min=min(moranRestWater_min), moranRestWater_max=max(moranRestWater_max), moranForage_min=min(moranForage_min), moranForage_max=max(moranForage_max),
moranLand_min=min(moranLand_min), moranLand_max=max(moranLand_max), 
timeTotal=sum(timeTotal, na.rm=TRUE), colonies=sum(colonies, na.rm=TRUE), 
weight=sum(weight, na.rm=TRUE), weightTotal=sum(weightTotal, na.rm=TRUE)) 

}

pairwisepop_res<-rbind(pairwisepop_res, pairwisePop_all)
colonyLox<-rbind(colonyLox, colonySub)

}

# Estimate species-weighted mean for the next calculation & mean within-population SD in time spent in different behaviours
actSubTot<-actSubTot %>%
ungroup() 
actSubTot[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]<-actSubTot[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]/actSubTot$weight #species, weighted-mean
actSubTot[c("timeFlight_sd", "timeForage_sd", "timeActive_sd", "timeLand_sd", "timeRestWater_sd")]<-actSubTot[c("timeFlight_sd", "timeForage_sd", "timeActive_sd", "timeLand_sd", "timeRestWater_sd")]/actSubTot$weight # within-population mean SD (weighted too)
actSubTot[c("corFlight", "corForage", "corActive", "corLand", "corRestWater")]<-actSubTot[c("corFlight", "corForage", "corActive", "corLand", "corRestWater")]/actSubTot$weightTotal # mean pop correlation
actSubTot[c("Moran_Flight", "Moran_Forage", "Moran_Active", "Moran_Land", "Moran_RestWater")]<-actSubTot[c("Moran_Flight", "Moran_Forage", "Moran_Active", "Moran_Land", "Moran_RestWater")]/actSubTot$weightTotal # mean pop moran's I
actSubTot[is.na(actSubTot)] <- 0 # Remove any NAs

# Prepare coordinates so we can join them back again
actSub<-read.csv(actFiles[i]) # Re-determine an index for the coordinates (this is to solve floating 0 problems)
actSub<-actSub %>%
  ungroup() %>%
  arrange(month, x, y) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(index=row_number())
actSubcoords<-actSub %>%
ungroup() %>%
dplyr::filter(month==1) %>%
dplyr::select(index, x, y)

# Rename metrics so it's obvious what they are & join coordinates for plotting
actSubTotFinal<-actSubTot %>%
ungroup() %>%
rename(timeFlightMean=timeFlight, timeForageMean=timeForage, timeActiveMean=timeActive, timeLandMean=timeLand, timeRestWaterMean=timeRestWater, timeTotalMean=timeTotal) %>%
rename(timeFlightSD_pop=timeFlight_sd, timeForageSD_pop=timeForage_sd, timeActiveSD_pop=timeActive_sd, timeLandSD_pop=timeLand_sd, timeRestWaterSD_pop=timeRestWater_sd,
corFlight_pop=corFlight, corActive_pop=corActive, corRestWater_pop=corRestWater, corLand_pop=corLand, corForage_pop=corForage) %>%
dplyr::select(-c(NoBirds, weight)) %>%
dplyr::left_join(actSubcoords, by=c("index"))

### Step 2: Do this again to estimate species-weighted SD ###

modelColonies2<-list()

for (i in 1:length(actFiles)) {

print(paste0("Mapping file SD", i, "/", length(actFiles), "..."))

#### Open file i ###
actSub2<-read.csv(actFiles[i])
actSub2<-subset(actSub2, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))
energySub2<-read.csv(energyFiles[i], nrow=10)

# Calculate month-specific weights
actSubWeights2<-actSub2 %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(totalBirds) %>%
dplyr::mutate(weight=1/totalBirds)

# We number the rows so we can join data frames later ( doesn't work so well with floating numbers)
actSub2<-actSub2 %>%
  ungroup() %>%
  arrange(month, x, y) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(index=row_number())

# Determine the model colony
metaSub<-subset(meta, colonies==actSub2$colony[1])
modelColSub2<-metaSub$modelcolonies[1]
print(modelColSub2)

# Does the list already have this model colony?
if (modelColSub2 %in% modelColonies2) {
print("Next")
next
}

# Otherwise add the modelcolony name to the list
modelColonies2<-append(modelColonies2, modelColSub2)

# We will multiply time spent in different behaviours by 1/no birds
#actSub2$month<-factor(actSub2$month)
actSub2$species<-input_file1
#actSub2$weight<-ifelse(actSub2$NoBirds>0, 1/actSub2$totalBirds[1], 0)
actSub2<-actSub2 %>%
dplyr::select(-c(x, y)) %>%
dplyr::left_join(actSubTotFinal, by=c("species", "month", "index")) %>%
dplyr::select(-totalBirds) %>%
dplyr::left_join(actSubWeights2, by=c("month"))
actSub2$weight<-ifelse(actSub2$NoBirds>0, 1/(actSub2$totalBirds), 0)
actSub2$timeFlight<-(actSub2$timeFlight - actSub2$timeFlightMean)^2*actSub2$weight
actSub2$timeLand<-(actSub2$timeLand - actSub2$timeLandMean)^2*actSub2$weight
actSub2$timeForage<-(actSub2$timeForage - actSub2$timeForageMean)^2*actSub2$weight
actSub2$timeActive<-(actSub2$timeActive - actSub2$timeActiveMean)^2*actSub2$weight
actSub2$timeRestWater<-(actSub2$timeRestWater - actSub2$timeRestWaterMean)^2*actSub2$weight
actSub2$timeTotal<-(actSub2$timeTotal - actSub2$timeTotalMean)^2*actSub2$weight
actSub2$colonies<-ifelse(actSub2$NoBirds>0, 1, 0)

#sub<-subset(actSub2, month==4 & index==1456)

#if (sub$NoBirds > 0) {
#print("Success")
#break
#}

if (i ==1) {

# Make a master dataset that we will join to
actSubTot2<-actSub2

} else {

# We sum everything together (which will be averaged by total number of birds later... #
actSub2<-actSub2 %>%
dplyr::select(species, month, index, NoBirds, timeFlight, timeRestWater, timeLand, timeForage, timeActive, timeTotal, weight, colonies, totalBirds)
actSub2$weight<-ifelse(actSub2$NoBirds>0, 1/(actSub2$totalBirds), 0)
actSubTot2<-rbind(actSubTot2, actSub2)
actSubTot2<-actSubTot2 %>%
ungroup() %>%
dplyr::group_by(species, month, index) %>%
dplyr::summarise(NoBirds=sum(NoBirds, na.rm=TRUE), timeFlight=sum(timeFlight, na.rm=TRUE), timeRestWater=sum(timeRestWater, na.rm=TRUE), 
timeLand=sum(timeLand, na.rm=TRUE), timeForage=sum(timeForage, na.rm=TRUE), timeActive=sum(timeActive, na.rm=TRUE), timeTotal=sum(timeTotal, na.rm=TRUE), weight=sum(weight, na.rm=TRUE), colonies=sum(colonies, na.rm=TRUE))

}

}

# Estimate species-weighted SD for the next calculation
actSubTot2<-actSubTot2 %>%
ungroup() 
actSubTot2[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]<-sqrt(actSubTot2[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]/actSubTot2$weight)
actSubTot2[is.na(actSubTot2)] <- 0
actSubTotFinal2<-actSubTot2 %>%
#dplyr::select(-c(timeActive_sd, timeActive_sd, timeRestWater_sd, timeLand_sd, timeActive_sd)) %>%
rename(timeFlight_sd=timeFlight, timeForage_sd=timeForage, timeActive_sd=timeActive, timeLand_sd=timeLand, timeRestWater_sd=timeRestWater, timeTotal_sd=timeTotal) %>%
dplyr::left_join(actSubcoords, by=c("index"))

# Conduct species-level correlations & plot results (species vs. population)

pop_means<-actSubTotFinal %>%
ungroup() %>%
dplyr::select(month, index, timeFlightMean, timeRestWaterMean, timeActiveMean, timeForageMean, timeLandMean)

actSubTotFinal2<-actSubTotFinal2 %>%
dplyr::left_join(pop_means, by=c("month", "index"))

actBirds_pop<-subset(actSubTotFinal2, colonies>1) # Subset to squares where this data
#indices<-data.frame(index=unique(actBirds$index))
#actBirds_pop<-actSubTotFinal2 %>%
#dplyr::inner_join(indices, by=c("index"))
actBirdsCorrSpecies<-corrMonth_species(actBirds_pop) # Calculate correlations for different behaviours
actSubTotFinal2<-actSubTotFinal2 %>%
dplyr::left_join(actBirdsCorrSpecies, by=c("month")) # Join back to main dataset

pop1<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corFlight_pop, corFlight_min, corFlight_max) %>%
dplyr::rename(corMean=corFlight_pop, corMin=corFlight_min, corMax=corFlight_max) %>%
dplyr::mutate(behaviour="Flight", type="Population")

pop2<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corActive_pop, corActive_min, corActive_max) %>%
dplyr::rename(corMean=corActive_pop, corMin=corActive_min, corMax=corActive_max) %>%
dplyr::mutate(behaviour="Active", type="Population")

pop3<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corRestWater_pop, corRestWater_min, corRestWater_max) %>%
dplyr::rename(corMean=corRestWater_pop, corMin=corRestWater_min, corMax=corRestWater_max) %>%
dplyr::mutate(behaviour="RestWater", type="Population")

pop4<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corForage_pop, corForage_min, corForage_max) %>%
dplyr::rename(corMean=corForage_pop, corMin=corForage_min, corMax=corForage_max) %>%
dplyr::mutate(behaviour="Forage", type="Population")

pop5<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corLand_pop, corLand_min, corLand_max) %>%
dplyr::rename(corMean=corLand_pop, corMin=corLand_min, corMax=corLand_max) %>%
dplyr::mutate(behaviour="Land", type="Population")

popAll_corr<-rbind(pop1, pop2, pop3, pop4, pop5)

species1<-actSubTotFinal2 %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corFlight) %>%
dplyr::rename(corMean=corFlight) %>%
dplyr::mutate(corMin=corMean, corMax=corMean) %>%
dplyr::mutate(behaviour="Flight", type="species")

species2<-actSubTotFinal2 %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corActive) %>%
dplyr::rename(corMean=corActive) %>%
dplyr::mutate(corMin=corMean, corMax=corMean) %>%
dplyr::mutate(behaviour="Active", type="species")

species3<-actSubTotFinal2 %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corRestWater) %>%
dplyr::rename(corMean=corRestWater) %>%
dplyr::mutate(corMin=corMean, corMax=corMean) %>%
dplyr::mutate(behaviour="RestWater", type="species")

species4<-actSubTotFinal2 %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corForage) %>%
dplyr::rename(corMean=corForage) %>%
dplyr::mutate(corMin=corMean, corMax=corMean) %>%
dplyr::mutate(behaviour="Forage", type="species")

species5<-actSubTotFinal2 %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(corLand) %>%
dplyr::rename(corMean=corLand) %>%
dplyr::mutate(corMin=corMean, corMax=corMean) %>%
dplyr::mutate(behaviour="Land", type="species")

speciesAll_corr<-rbind(species1, species2, species3, species4, species5)

# Join both types of data #

speciesPop_metrics<-rbind(speciesAll_corr, popAll_corr)

# Make the plot comparing populations & species #

corrPop_vs_species<-speciesPop_metrics %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))) %>%
ggplot(aes(x=factor(month), y=corMean)) +
geom_pointrange(aes(ymin=corMin, ymax=corMax, color=type, group=type), position=position_dodge2(width=0.5))  +
geom_hline(yintercept=-0.6, linetype="dashed") +
geom_hline(yintercept=0.6, linetype="dashed") +
geom_hline(yintercept=0, linetype="dashed", color="darkred") +
facet_wrap(~ behaviour) +
#ylim(-1, 1) +
theme_bw() +
xlab("Month") +
ylab("Correlation mean vs. sd") +
ggtitle("Are there shared hotspots? (neg corr)")

# Attach population SD to this data frame so we can compute a comparison metrics
pop_sd<-actSubTotFinal %>%
ungroup() %>%
dplyr::select(month, index, timeFlightSD_pop, timeRestWaterSD_pop, timeActiveSD_pop, timeForageSD_pop, timeLandSD_pop)

actSubTotFinal3<-actSubTotFinal2 %>%
dplyr::left_join(pop_sd, by=c("month", "index"))

actSubTotFinal3$timeFlight_sd_diff<-actSubTotFinal3$timeFlight_sd/actSubTotFinal3$timeFlightSD_pop
actSubTotFinal3$timeFlight_sd_diff<-ifelse(actSubTotFinal3$timeFlight_sd_diff %in% c("Inf"), 0, actSubTotFinal3$timeFlight_sd_diff)
actSubTotFinal3$timeFlight_sd_diff<-ifelse(actSubTotFinal3$timeFlight_sd_diff %in% c("-Inf"), 0, actSubTotFinal3$timeFlight_sd_diff)

actSubTotFinal3$timeForage_sd_diff<-actSubTotFinal3$timeForage_sd/actSubTotFinal3$timeForageSD_pop
actSubTotFinal3$timeForage_sd_diff<-ifelse(actSubTotFinal3$timeForage_sd_diff %in% c("Inf"), 0, actSubTotFinal3$timeForage_sd_diff)
actSubTotFinal3$timeForage_sd_diff<-ifelse(actSubTotFinal3$timeForage_sd_diff %in% c("-Inf"), 0, actSubTotFinal3$timeForage_sd_diff)

actSubTotFinal3$timeRestWater_sd_diff<-actSubTotFinal3$timeRestWater_sd/actSubTotFinal3$timeRestWaterSD_pop
actSubTotFinal3$timeRestWater_sd_diff<-ifelse(actSubTotFinal3$timeRestWater_sd_diff %in% c("Inf"), 0, actSubTotFinal3$timeRestWater_sd_diff)
actSubTotFinal3$timeRestWater_sd_diff<-ifelse(actSubTotFinal3$timeRestWater_sd_diff %in% c("-Inf"), 0, actSubTotFinal3$timeRestWater_sd_diff)

actSubTotFinal3$timeActive_sd_diff<-actSubTotFinal3$timeActive_sd/actSubTotFinal3$timeActiveSD_pop
actSubTotFinal3$timeActive_sd_diff<-ifelse(actSubTotFinal3$timeActive_sd_diff %in% c("Inf"), 0, actSubTotFinal3$timeActive_sd_diff)
actSubTotFinal3$timeActive_sd_diff<-ifelse(actSubTotFinal3$timeActive_sd_diff %in% c("-Inf"), 0, actSubTotFinal3$timeActive_sd_diff)

actSubTotFinal3$timeLand_sd_diff<-actSubTotFinal3$timeLand_sd/actSubTotFinal3$timeLandSD_pop
actSubTotFinal3$timeLand_sd_diff<-ifelse(actSubTotFinal3$timeLand_sd_diff %in% c("Inf"), 0, actSubTotFinal3$timeLand_sd_diff)
actSubTotFinal3$timeLand_sd_diff<-ifelse(actSubTotFinal3$timeLand_sd_diff %in% c("-Inf"), 0, actSubTotFinal3$timeLand_sd_diff)

actSubTotFinal3[is.na(actSubTotFinal3)] <- 0 # remove any NAs for plotting

# Now we make plots showing species-weighted mean, difference between within-pop SD & between-pop SD & the same but shown using boxplots (for every behaviour)

# Make this plot so I can make extents for my maps #
actMultiPops<-subset(actSubTotFinal, colonies>0)
actMultiPops2<-subset(actSubTotFinal2, colonies>0)

# make maps for every behaviour
mapFlightSpecies<-mapActivity_species(actSubTotFinal, actSubTotFinal3, actMultiPops, "Flight", energySub$species[1])
mapRestSpecies<-mapActivity_species(actSubTotFinal, actSubTotFinal3, actMultiPops, "RestWater", energySub$species[1])
if(energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
mapForageSpecies<-mapActivity_species(actSubTotFinal, actSubTotFinal3, actMultiPops, "Forage", energySub$species[1])
} else {
mapActiveSpecies<-mapActivity_species(actSubTotFinal, actSubTotFinal3, actMultiPops, "Active", energySub$species[1])
}
#mapLandSpecies<-mapActivity_species(actSubTotFinal, actSubTotFinal3, actMultiPops, "Land", energySub$species[1])

### R_ratio distribution (between/within pop SD) ###

R_summary1 <- actSubTotFinal3 %>%
dplyr::mutate(behaviour="Flight") %>%
ungroup() %>%
dplyr::filter(colonies>1) %>%
dplyr::mutate(betweenvswithin=ifelse(timeFlight_sd_diff > 1, 1, 0)) %>%
group_by(behaviour, month) %>%
dplyr::mutate(cells=n_distinct(index)) %>%
  summarize(
    median_R = median(timeFlight_sd_diff, na.rm=TRUE),
    pct_R_gt_1 = sum(betweenvswithin, na.rm=TRUE)/mean(cells),
    .groups="drop"
  )
  
R_summary2 <- actSubTotFinal3 %>%
dplyr::mutate(behaviour="RestWater") %>%
ungroup() %>%
dplyr::filter(colonies>1) %>%
dplyr::mutate(betweenvswithin=ifelse(timeRestWater_sd_diff > 1, 1, 0)) %>%
group_by(behaviour, month) %>%
dplyr::mutate(cells=n_distinct(index)) %>%
  summarize(
    median_R = median(timeRestWater_sd_diff, na.rm=TRUE),
    pct_R_gt_1 = sum(betweenvswithin, na.rm=TRUE)/mean(cells),
    .groups="drop"
  )
  
 R_summary3 <- actSubTotFinal3 %>%
dplyr::mutate(behaviour="Land") %>%
ungroup() %>%
dplyr::filter(colonies>1) %>%
dplyr::mutate(betweenvswithin=ifelse(timeLand_sd_diff > 1, 1, 0)) %>%
group_by(behaviour, month) %>%
dplyr::mutate(cells=n_distinct(index)) %>%
  summarize(
    median_R = median(timeLand_sd_diff, na.rm=TRUE),
    pct_R_gt_1 = sum(betweenvswithin, na.rm=TRUE)/mean(cells),
    .groups="drop"
  )
  
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {

R_summary4 <- actSubTotFinal3 %>%
dplyr::mutate(behaviour="Forage") %>%
ungroup() %>%
dplyr::filter(colonies>1) %>%
dplyr::mutate(betweenvswithin=ifelse(timeForage_sd_diff > 1, 1, 0)) %>%
group_by(behaviour, month) %>%
dplyr::mutate(cells=n_distinct(index)) %>%
  summarize(
    median_R = median(timeForage_sd_diff, na.rm=TRUE),
    pct_R_gt_1 = sum(betweenvswithin, na.rm=TRUE)/mean(cells),
    .groups="drop")

} else {
  
  R_summary4 <- actSubTotFinal3 %>%
dplyr::mutate(behaviour="Active") %>%
ungroup() %>%
dplyr::filter(colonies>1) %>%
dplyr::mutate(betweenvswithin=ifelse(timeActive_sd_diff > 1, 1, 0)) %>%
group_by(behaviour, month) %>%
dplyr::mutate(cells=n_distinct(index)) %>%
  summarize(
    median_R = median(timeActive_sd_diff, na.rm=TRUE),
    pct_R_gt_1 = sum(betweenvswithin, na.rm=TRUE)/mean(cells),
    .groups="drop"
  )
  
  }

R_summary<-rbind(R_summary1, R_summary2, R_summary3, R_summary4)

between_vs_within_plot<-R_summary %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))) %>%
ggplot(aes(x=month, y=pct_R_gt_1)) +
geom_point(aes(color=behaviour), position=position_dodge2(width=0.5)) +
ylim(0, 1) +
ylab("Prop cells SD1>SD2") +
theme_bw() +
ggtitle("Between vs. within pop SD")

#### Compute Moran's I: ####

# Subset datasets to same extent as pop-level moran's I

indices<-data.frame(index=unique(actBirds$index))

actSubTotFinal_ext<-actSubTotFinal %>%
dplyr::inner_join(indices, by=c("index"))

actSubTotFinal3_ext3<-actSubTotFinal3 %>%
dplyr::inner_join(indices, by=c("index"))

moranFlight<-computeMoran(actSubTotFinal_ext, actSubTotFinal3_ext3, "Flight") 
moranRest<-computeMoran(actSubTotFinal_ext, actSubTotFinal3_ext3, "RestWater") 
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
moranForage<-computeMoran(actSubTotFinal_ext, actSubTotFinal3_ext3, "Forage")
} else {
moranActive<-computeMoran(actSubTotFinal_ext, actSubTotFinal3_ext3, "Active")}
moranLand<-computeMoran(actSubTotFinal_ext, actSubTotFinal3_ext3, "Land")

# Make a species data frame with all behaviours for plotting #
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
moranSpecies<-rbind(moranFlight, moranRest, moranForage,  moranLand)
} else {
moranSpecies<-rbind(moranFlight, moranRest, moranActive,  moranLand)
}
moranSpeciesFinal<-moranSpecies %>%
dplyr::select(month, moranI_speciesMean, behaviour) %>%
dplyr::rename(moranI=moranI_speciesMean) %>%
dplyr::mutate(moranMin=moranI, moranMax=moranI) %>%
dplyr::mutate(type="Species") %>%
dplyr::mutate(behaviour=ifelse(behaviour=="Rest", "RestWater", behaviour))

# Make a population-specific data frame for plotting
pop1_mor<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(Moran_Flight, moranFlight_min, moranFlight_max) %>%
dplyr::rename(moranI=Moran_Flight, moranMin=moranFlight_min, moranMax=moranFlight_max) %>%
dplyr::mutate(behaviour="Flight", type="Population")

pop2_mor<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(Moran_RestWater, moranRestWater_min, moranRestWater_max) %>%
dplyr::rename(moranI=Moran_RestWater, moranMin=moranRestWater_min, moranMax=moranRestWater_max) %>%
dplyr::mutate(behaviour="RestWater", type="Population")

pop3_mor<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(Moran_Land, moranLand_min, moranLand_max) %>%
dplyr::rename(moranI=Moran_Land, moranMin=moranLand_min, moranMax=moranLand_max) %>%
dplyr::mutate(behaviour="Land", type="Population")

if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {

pop4_mor<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(Moran_Forage, moranForage_min, moranForage_max) %>%
dplyr::rename(moranI=Moran_Forage, moranMin=moranForage_min, moranMax=moranForage_max) %>%
dplyr::mutate(behaviour="Forage", type="Population")

} else {

pop4_mor<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::count(Moran_Active, moranActive_min, moranActive_max) %>%
dplyr::rename(moranI=Moran_Active, moranMin=moranActive_min, moranMax=moranActive_max) %>%
dplyr::mutate(behaviour="Active", type="Population")

}

moranPopFinal<-rbind(pop1_mor, pop2_mor, pop3_mor, pop4_mor)

# plot these together #
moranAll<-rbind(moranPopFinal, moranSpeciesFinal)

moranPop_vs_species<-moranAll %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))) %>%
ggplot(aes(x=factor(month), y=moranI)) +
geom_pointrange(aes(ymin=moranMin, ymax=moranMax, color=type, group=type), position=position_dodge2(width=0.5))  +
#geom_hline(yintercept=-0.6, linetype="dashed") +
#geom_hline(yintercept=0.6, linetype="dashed") +
geom_hline(yintercept=0, linetype="dashed", color="darkred") +
facet_wrap(~ behaviour) +
#ylim(-1, 1) +
theme_bw() +
xlab("Month") +
ylab("Moran's I (mean, min, max)") +
ggtitle("Are there hotspots?")

### Now I also make a corrPop to pop figrue ###
pairwiseRes<-pairwisepop_res %>%
ungroup() %>%
dplyr::group_by(behaviour, month) %>%
dplyr::summarise(meancor=mean(cor, na.rm=TRUE), sdcor=sd(cor, na.rm=TRUE), mincor=min(cor, na.rm=TRUE), maxcor=max(cor, na.rm=TRUE)) %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3, 4)))

popcor<-ggplot() +
geom_pointrange(data=pairwiseRes, aes(x=factor(month), y=meancor, ymin=meancor-sdcor, ymax=meancor + sdcor, color=behaviour), position=position_dodge2(width=0.6)) +
ylim(-1.4, 1.4) +
theme_bw() +
ylab("Pair-wise corr (mean/SD)") +
geom_hline(yintercept=0, linetype="dashed") +
geom_hline(yintercept=-0.6, linetype="dashed", color="darkred") +
geom_hline(yintercept=0.6, linetype="dashed", color="darkred") +
ggtitle("Pop-pop correlation") +
xlab("Month")

### Save three metrics ### : Moran's I, correlation between mean & SD + proportion of squares where between R is higher than within R

pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/activity/", energySub$species[1], "_summaryMetrics.pdf"), width=10, height=5)
plot(moranPop_vs_species)
plot(corrPop_vs_species) 
grid.arrange(between_vs_within_plot, popcor, nrow=2)
dev.off()

# Save outputs

output_file1 <- args[2]
write.csv(actSubTotFinal, file=output_file1, row.names=FALSE)

output_file2 <- args[3]
write.csv(actSubTotFinal2, file=output_file2, row.names=FALSE)


