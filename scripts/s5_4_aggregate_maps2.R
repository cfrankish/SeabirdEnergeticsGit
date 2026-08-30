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
#input_file1<-"Atlanticpuffin"

# Now we determine the list of activity/energy data to loop through
allFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp4/", full.names=TRUE)
speciesFiles<-allFiles[grep(input_file1, allFiles)] # Subset to species-specific ones
#print("species")
actFiles<-speciesFiles[grep("activityMap_yearly.csv", speciesFiles)] # extract activity files
#print("act")
energyFiles<-speciesFiles[grep("energy", speciesFiles)] # extract energy files
#print("energy")

# Create a world map for plotting
#print("Prepping map info")
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# open colony coordinates to put on maps
#print("Getting colony coords")
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

#### Step 1: Create species-weighted mean for activity (add SD later) ####

modelColonies<-list()
colonyLox<-list()

print("Starting for loop")

for (i in 1:length(actFiles)) {

print(paste0("Mapping file", i, "/", length(actFiles), "..."))

#### Open file i ###
actSub<-read.csv(actFiles[i])
#actSub<-subset(actSub, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))
energySub<-read.csv(energyFiles[i], nrow=10)

# Calculate month-specific weights
actSubWeights<-actSub %>%
ungroup() %>%
#dplyr::group_by(month) %>%
dplyr::count(totalBirds) %>%
dplyr::mutate(weight=1/totalBirds)

# We number the rows so we can join data frames later ( doesn't work so well with floating numbers)
actSub<-actSub %>%
  ungroup() %>%
  arrange(x, y) %>%
  #dplyr::group_by(month) %>%
  dplyr::mutate(index=row_number())

# Determine the model colony
metaSub<-subset(meta, colonies==actSub$colony[1])
modelColSub<-metaSub$modelcolonies[1]
print(modelColSub)

# Does the list already have this model colony?
if (modelColSub %in% modelColonies) {
print("Next")
next
}

# We plot population-level activity budgets here #
#actSub$month<-factor(actSub$month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))

# Can we match it with colonyCoords?
colonySub<-subset(colonyCoordsTrans, colony==actSub$colony_og[1])

if(nrow(colonySub<1)) {
colonyCoordsTrans$colony<-gsub("Bjørnøya", "Bjornoya", colonyCoordsTrans$colony)
colonySub<-subset(colonyCoordsTrans, colony==actSub$colony_og[1])
}

# Otherwise add the modelcolony name to the list
modelColonies<-append(modelColonies, modelColSub)

# We will multiply time spent in different behaviours by 1/no birds
actSub[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]<-actSub[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]/actSub$totalBirds
actSub$weight<-ifelse(actSub$NoBirds>0, 1/(actSub$totalBirds), 0)
actSub[is.na(actSub)] <- 0
actSub$species<-input_file1
actSub$colonies<-ifelse(actSub$NoBirds>0, 1, 0)

if (i ==1) {

# Make a master dataset that we will join to
actSub<-actSub %>%
dplyr::select(species, index, NoBirds, timeFlight, timeRestWater, timeLand, timeForage, timeActive, timeTotal, colonies, totalBirds, weight) 
actSubTot<-actSub
actSubTot<-actSubTot %>%
dplyr::select(-c(weight, totalBirds)) 
actSubTot$totalBirds<-actSubWeights$totalBirds
actSubTot$weight<-actSubWeights$weight
#actSubTot$timeFlight<-ifelse(actSubTot$NoBirds==0, NA, actSubTot$timeFlight)
#actSubTot$timeForage<-ifelse(actSubTot$NoBirds==0, NA, actSubTot$timeForage)
#actSubTot$timeLand<-ifelse(actSubTot$NoBirds==0, NA, actSubTot$timeLand)
#actSubTot$timeRestWater<-ifelse(actSubTot$NoBirds==0, NA, actSubTot$timeRestWater)
#actSubTot$timeActive<-ifelse(actSubTot$NoBirds==0, NA, actSubTot$timeActive)
#actSubTot$timeLand<-ifelse(actSubTot$NoBirds==0, NA, actSubTot$timeLand)
#actSubTot$timeTotal<-ifelse(actSubTot$NoBirds==0, NA, actSubTot$timeTotal)
actSubTot$weight<-ifelse(actSubTot$NoBirds==0, 0, actSubTot$weight)

} else {

# We sum everything together (which will be averaged by total number of birds later... #
actSub<-actSub %>%
dplyr::select(species, index, NoBirds, timeFlight, timeRestWater, timeLand, timeForage, timeActive, timeTotal, colonies, totalBirds) 
actSub$totalBirds<-actSubWeights$totalBirds
actSub$weight<-actSubWeights$weight
#actSub$timeFlight<-ifelse(actSub$NoBirds==0, NA, actSub$timeFlight)
#actSub$timeForage<-ifelse(actSub$NoBirds==0, NA, actSub$timeForage)
#actSub$timeLand<-ifelse(actSub$NoBirds==0, NA, actSub$timeLand)
#actSub$timeRestWater<-ifelse(actSub$NoBirds==0, NA, actSub$timeRestWater)
#actSub$timeActive<-ifelse(actSub$NoBirds==0, NA, actSub$timeActive)
#actSub$timeLand<-ifelse(actSub$NoBirds==0, NA, actSub$timeLand)
#actSub$timeTotal<-ifelse(actSub$NoBirds==0, NA, actSub$timeTotal)
actSub$weight<-ifelse(actSub$NoBirds==0, 0, actSub$weight)
actSubTot<-rbind(actSubTot, actSub)
actSubTot<-actSubTot %>%
ungroup () %>%
dplyr::group_by(species, index) %>%
dplyr::summarise(NoBirds=sum(NoBirds, na.rm=TRUE), timeFlight=sum(timeFlight, na.rm=TRUE), 
timeRestWater=sum(timeRestWater, na.rm=TRUE), timeLand=sum(timeLand, na.rm=TRUE), timeForage=sum(timeForage, na.rm=TRUE), timeActive=sum(timeActive, na.rm=TRUE), timeTotal=sum(timeTotal, na.rm=TRUE), colonies=sum(colonies, na.rm=TRUE),
weight=sum(weight, na.rm=TRUE)) 

}

colonyLox<-rbind(colonyLox, colonySub)

}

# Estimate species-weighted mean for the next calculation
actSubTot<-actSubTot %>%
ungroup() 
actSubTot[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]<-actSubTot[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater", "timeTotal")]/actSubTot$weight
actSubTot[is.na(actSubTot)] <- 0
actSub<-read.csv(actFiles[i])
actSub<-actSub %>%
  ungroup() %>%
  arrange(x, y) %>%
  #dplyr::group_by(month) %>%
  dplyr::mutate(index=row_number())
actSubcoords<-actSub %>%
ungroup() %>%
#dplyr::filter(month==1) %>%
dplyr::select(index, x, y)
actSubTotFinal<-actSubTot %>%
ungroup() %>%
rename(timeFlightMean=timeFlight, timeForageMean=timeForage, timeActiveMean=timeActive, timeLandMean=timeLand, timeRestWaterMean=timeRestWater, timeTotalMean=timeTotal) %>%
dplyr::select(-c(NoBirds, weight)) %>%
dplyr::left_join(actSubcoords, by=c("index"))

### Step 2: Do this again to estimate species-weighted SD ###

modelColonies2<-list()

for (i in 1:length(actFiles)) {

print(paste0("Mapping file SD", i, "/", length(actFiles), "..."))

#### Open file i ###
actSub2<-read.csv(actFiles[i])
#actSub2<-subset(actSub2, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))
energySub2<-read.csv(energyFiles[i], nrow=10)

# Calculate month-specific weights
actSubWeights2<-actSub2 %>%
ungroup() %>%
#dplyr::group_by(month) %>%
dplyr::count(totalBirds) %>%
dplyr::mutate(weight=1/totalBirds)

# We number the rows so we can join data frames later ( doesn't work so well with floating numbers)
actSub2<-actSub2 %>%
  ungroup() %>%
  arrange( x, y) %>%
 # dplyr::group_by(month) %>%
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
dplyr::left_join(actSubTotFinal, by=c("species", "index")) %>%
dplyr::select(species, index, NoBirds, timeFlight, timeFlightMean, timeRestWater, timeRestWaterMean, timeLand, timeLandMean, timeForage, timeForageMean, timeActive, timeActiveMean, timeTotal, timeTotalMean, colonies)
actSub2$weight<-actSubWeights$weight
actSub2$totalBirds<-actSubWeights$totalBirds
actSub2$weight<-ifelse(actSub2$NoBirds==0, 0, actSub2$weight)
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
dplyr::select(species, index, NoBirds, timeFlight, timeFlightMean, timeRestWater, timeRestWaterMean, timeLand, timeLandMean, timeForage, timeForageMean, timeActive, timeActiveMean, timeTotal, timeTotalMean, colonies)
actSub2$weight<-actSubWeights$weight
actSub2$totalBirds<-actSubWeights$totalBirds
actSub2$weight<-ifelse(actSub2$NoBirds==0, 0, actSub2$weight)
actSubTot2<-rbind(actSubTot2, actSub2)
actSubTot2<-actSubTot2 %>%
ungroup() %>%
dplyr::group_by(species, index) %>%
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
#dplyr::select(-c(timeFlight_sd, timeActive_sd, timeRestWater_sd, timeLand_sd, timeActive_sd)) %>%
rename(timeFlight_sd=timeFlight, timeForage_sd=timeForage, timeActive_sd=timeActive, timeLand_sd=timeLand, timeRestWater_sd=timeRestWater, timeTotal_sd=timeTotal) %>%
dplyr::left_join(actSubcoords, by=c("index"))

#actSubTotFinal$month<-factor(actSubTotFinal$month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))
#actSubTotFinal2$month<-factor(actSubTotFinal2$month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))

# Scale #
actMultiPops<-subset(actSubTotFinal, colonies>0)
actMultiPops2<-subset(actSubTotFinal2, colonies>0)
minFlight1<-min(actMultiPops$timeFlightMean)
minFlight2<-min(actMultiPops2$timeFlight_sd)
maxFlight1<-max(actMultiPops$timeFlightMean)
maxFlight2<-max(actMultiPops2$timeFlight_sd)
minFlight<-ifelse(minFlight1<minFlight2, minFlight2, minFlight2)
maxFlight<-ifelse(maxFlight1>maxFlight2, maxFlight1, maxFlight2)

# make map 1
flight1<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeFlightMean)) +
  #geom_text(data=actSubTotFinal, aes(x=x, y=y, label=colonies)) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent flight (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  scale_fill_gradientn('Prop time flight', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxFlight)) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_flightMean_winter.pdf"))
#plot(flight1)
#dev.off()

# SD #

flight2<-ggplot() +
  #geom_tile(data=filter(actSubTotFinal, timeFlightMean>0), aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeFlight_sd)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent flight (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeFlightMean), max(actSubTotFinal$timeFlightMean))) +
  scale_fill_gradientn('Prop time flight (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxFlight)) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_flightSD_winter.pdf"))
#plot(flight2)
#dev.off()

# make map 1

minRest1<-min(actMultiPops$timeRestWaterMean)
minRest2<-min(actMultiPops2$timeRestWater_sd)
maxRest1<-max(actMultiPops$timeRestWaterMean)
maxRest2<-max(actMultiPops2$timeRestWater_sd)
minRest<-ifelse(minRest1<minRest2, minRest2, minRest2)
maxRest<-ifelse(maxRest1>maxRest2, maxRest1, maxRest2)

actSubTotFinal$timeRestWaterMean2<-ifelse(actSubTotFinal$colonies<2, 0, actSubTotFinal$timeRestWaterMean)

RestWater1<-ggplot() +
  #geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeRestWaterMean)) +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeRestWaterMean)) +
  #scale_fill_gradientn('Time spent rest (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  scale_fill_gradientn('Prop time rest ', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxRest)) +
  #scale_fill_gradientn('Time spent RestWater (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
   coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_RestWaterMean_winter.pdf"))
#plot(RestWater1)
#dev.off()

# SD #

RestWater2<-ggplot() +
  #geom_tile(data=filter(actSubTotFinal, timeRestWaterMean>0), aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeRestWater_sd)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  scale_fill_gradientn('Prop time rest (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxRest)) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
   #scale_fill_gradientn('Time spent rest (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeRestWaterMean), max(actSubTotFinal$timeRestWaterMean))) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_RestWaterSD_winter.pdf"))
#plot(RestWater2)
#dev.off()

# make map 1

minLand1<-min(actMultiPops$timeLandMean)
minLand2<-min(actMultiPops2$timeLand_sd)
maxLand1<-max(actMultiPops$timeLandMean)
maxLand2<-max(actMultiPops2$timeLand_sd)
minLand<-ifelse(minLand1<minLand2, minLand2, minLand2)
maxLand<-ifelse(maxLand1>maxLand2, maxLand1, maxLand2)

Land1<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeLandMean)) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Land (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  scale_fill_gradientn('Prop time land', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minLand, maxLand)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_LandMean_winter.pdf"))
#plot(Land1)
#dev.off()

# Mean behaviour #

Land2<-ggplot() +
  #geom_tile(data=filter(actSubTotFinal, timeLandMean>0), aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeLand_sd)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Land (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeLandMean), max(actSubTotFinal$timeLandMean))) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  scale_fill_gradientn('Prop time land (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minLand, maxLand)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
 geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_LandSD_winter.pdf"))
#plot(Land2)
#dev.off()

if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {

# make map 1

minForage1<-min(actMultiPops$timeForageMean)
minForage2<-min(actMultiPops2$timeForage_sd)
maxForage1<-max(actMultiPops$timeForageMean)
maxForage2<-max(actMultiPops2$timeForage_sd)
minForage<-ifelse(minForage1<minForage2, minForage2, minForage2)
maxForage<-ifelse(maxForage1>maxForage2, maxForage1, maxForage2)

Forage1<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeForageMean)) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Forage (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
   scale_fill_gradientn('Prop time forage', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxForage)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_ForageMean_winter.pdf"))
#plot(Forage1)
#dev.off()

# Mean behaviour #

Forage2<-ggplot() +
  #geom_tile(data=filter(actSubTotFinal, timeForageMean>0), aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeForage_sd)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Forage (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeForageMean), max(actSubTotFinal$timeForageMean))) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
   scale_fill_gradientn('Prop time forage (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxForage)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_ForageSD_winter.pdf"))
#plot(Forage2)
#dev.off()

} else {

# make map 1

minActive1<-min(actMultiPops$timeActiveMean)
minActive2<-min(actMultiPops2$timeActive_sd)
maxActive1<-max(actMultiPops$timeActiveMean)
maxActive2<-max(actMultiPops2$timeActive_sd)
minActive<-ifelse(minActive1<minActive2, minActive2, minActive2)
maxActive<-ifelse(maxActive1>maxActive2, maxActive1, maxActive2)

Active1<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeActiveMean)) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Active (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  scale_fill_gradientn('Prop time active', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxActive)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_ActiveMean_winter.pdf"))
#plot(Active1)
#dev.off()

# Mean behaviour #

Active2<-ggplot() +
  #geom_tile(data=actSubTotFinal, aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeActive_sd)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Active (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeActiveMean), max(actSubTotFinal$timeActiveMean))) +
  scale_fill_gradientn('Prop time active (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxActive)) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_ActiveSD_winter.pdf"))
#plot(Active2)
#dev.off()

}

# Plot No of birds #

noBirds<-ggplot() +
  #geom_tile(data=actSubTotFinal, aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=NoBirds)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Active (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeActiveMean), max(actSubTotFinal$timeActiveMean))) +
  scale_fill_gradientn('No Birds', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_NoBirds_winter.pdf"))
#plot(noBirds)
#dev.off()

# Plot No of pops #

noPops<-ggplot() +
  #geom_tile(data=actSubTotFinal, aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=colonies)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Active (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeActiveMean), max(actSubTotFinal$timeActiveMean))) +
  scale_fill_gradientn('No pops', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_NoPops_winter.pdf"))
#plot(noPops)
#dev.off()

# Plot time spent in cells overall #

maxTime1<-max(actMultiPops$timeTotalMean)
maxTime2<-max(actMultiPops2$timeTotal_sd)
maxTime<-ifelse(maxTime1>maxTime2, maxTime1, maxTime2)

timeTotal1<-ggplot() +
  #geom_tile(data=actSubTotFinal, aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeTotalMean)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Active (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeActiveMean), max(actSubTotFinal$timeActiveMean))) +
  scale_fill_gradientn('Total time', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, maxTime)) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
 # facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_propTimeMean_winter.pdf"))
#plot(timeTotal1)
#dev.off()

timeTotal2<-ggplot() +
  #geom_tile(data=actSubTotFinal, aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeTotal_sd)) +
  #geom_text(data=actSubTotFinal2, aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  #scale_fill_gradientn('Time spent Active (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeActiveMean), max(actSubTotFinal$timeActiveMean))) +
  scale_fill_gradientn('Total time (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"),limits=c(0, maxTime)) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonyLox, aes(x=coords.x1, y=coords.x2), fill="yellow",cex=1, shape=21) +
  labs(colour="") +
  #facet_wrap(~ factor(month, levels = c(9, 10, 11, 12, 1, 2, 3, 4))) +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

#pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_propTimeSD_winter.pdf"))
#plot(timeTotal2)
#dev.off()

if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {

pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_allMean.pdf"),width=10)
grid.arrange(flight1, RestWater1, Forage1, timeTotal1, noPops, noBirds, nrow=2)
dev.off()

pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_allSD.pdf"),width=10)
grid.arrange(flight2, RestWater2, Forage2, timeTotal1, noPops, noBirds, nrow=2)
dev.off()


} else {

pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_allMean.pdf"), width=10)
grid.arrange(flight1, RestWater1, Active1, timeTotal1, noPops, noBirds, nrow=2)
dev.off()

pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_allSD.pdf"), width = 10)
grid.arrange(flight2, RestWater2, Active2, timeTotal2, noPops, noBirds, nrow = 2)
dev.off()


}


# Save outputs

output_file1 <- args[2]
write.csv(actSubTotFinal, file=output_file1, row.names=FALSE)

output_file2 <- args[3]
write.csv(actSubTotFinal2, file=output_file2, row.names=FALSE)


