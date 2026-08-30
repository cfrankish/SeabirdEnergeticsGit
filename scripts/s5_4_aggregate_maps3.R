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
#input_file1<-"Littleauk"

# Now we determine the list of activity/energy data to loop through
allFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp4/", full.names=TRUE)
speciesFiles<-allFiles[grep(input_file1, allFiles)] # Subset to species-specific ones
#print("species")
actFiles<-speciesFiles[grep("activity", speciesFiles)] # extract activity files
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

print("Starting for loop")

for (i in 1:length(actFiles)) {

print(paste0("Mapping file", i, "/", length(actFiles), "..."))

#### Open file i ###
actSub<-read.csv(actFiles[i])
actSub<-subset(actSub, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))
energySub<-read.csv(energyFiles[i], nrow=10)

# We number the rows so we can join data frames later ( doesn't work so well with floating numbers)
actSub<-actSub %>%
  ungroup() %>%
  arrange(month, x, y) %>%
  dplyr::group_by(month) %>%
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
actSub$month<-factor(actSub$month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))

# Can we match it with colonyCoords?
colonySub<-subset(colonyCoordsTrans, colony==energySub$colony[1])

if(nrow(colonySub<1)) {
colonyCoordsTrans$colony<-gsub("Bjørnøya", "Bjornoya", colonyCoordsTrans$colony)
colonySub<-subset(colonyCoordsTrans, colony==energySub$colony[1])
}

# Otherwise add the modelcolony name to the list
modelColonies<-append(modelColonies, modelColSub)

# We will multiply time spent in different behaviours by 1/no birds
actSub$species<-input_file1
actSub$colonies<-ifelse(actSub$NoBirds>0, 1, 0)

if (i ==1) {

# Make a master dataset that we will join to
actSubTot<-actSub %>%
dplyr::ungroup() %>%
dplyr::group_by(species, index) %>%
dplyr::summarise(NoBirds=sum(NoBirds, na.rm=TRUE), timeFlight=sum(timeFlight, na.rm=TRUE), 
timeRestWater=sum(timeRestWater, na.rm=TRUE), timeLand=sum(timeLand, na.rm=TRUE), timeForage=sum(timeForage, na.rm=TRUE), timeActive=sum(timeActive, na.rm=TRUE), colonies=max(colonies, na.rm=TRUE))
actSubTot[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater")]<-actSubTot[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater")]/actSubTot$NoBirds
actSubTot$weight<-ifelse(actSubTot$NoBirds>0, 1/actSubTot$NoBirds, 0)

} else {

# We sum everything together (which will be averaged by total number of birds later... #
actSub2<-actSub %>%
dplyr::ungroup() %>%
dplyr::group_by(species, index) %>%
dplyr::summarise(NoBirds=sum(NoBirds, na.rm=TRUE), timeFlight=sum(timeFlight, na.rm=TRUE), 
timeRestWater=sum(timeRestWater, na.rm=TRUE), timeLand=sum(timeLand, na.rm=TRUE), timeForage=sum(timeForage, na.rm=TRUE), timeActive=sum(timeActive, na.rm=TRUE), colonies=max(colonies, na.rm=TRUE))

actSub2[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater")]<-actSub2[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater")]/actSub2$NoBirds
actSub2$weight<-ifelse(actSub2$NoBirds>0, 1/actSub2$NoBirds, 0)

actSubTot<-rbind(actSubTot, actSub2)
actSubTot<-actSubTot %>%
ungroup () %>%
dplyr::group_by(species, index) %>%
dplyr::summarise(NoBirds=sum(NoBirds, na.rm=TRUE), timeFlight=sum(timeFlight, na.rm=TRUE), 
timeRestWater=sum(timeRestWater, na.rm=TRUE), timeLand=sum(timeLand, na.rm=TRUE), timeForage=sum(timeForage, na.rm=TRUE), timeActive=sum(timeActive, na.rm=TRUE), colonies=sum(colonies, na.rm=TRUE), weight=sum(weight, na.rm=TRUE))

}

}

# Estimate species-weighted mean for the next calculation
actSubTot<-actSubTot %>%
ungroup() 
actSubTot[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater")]<-actSubTot[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater")]/actSubTot$weight
actSubTot[is.na(actSubTot)] <- 0
actSub<-read.csv(actFiles[i])
actSub<-actSub %>%
  ungroup() %>%
  arrange(month, x, y) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(index=row_number())
actSubcoords<-actSub %>%
ungroup() %>%
dplyr::filter(month==1) %>%
dplyr::select(index, x, y)
actSubTotFinal<-actSubTot %>%
ungroup() %>%
rename(NoBirdsTot=NoBirds) %>%
rename(timeFlightMean=timeFlight, timeForageMean=timeForage, timeActiveMean=timeActive, timeLandMean=timeLand, timeRestWaterMean=timeRestWater) %>%
dplyr::left_join(actSubcoords, by=c("index"))

### Step 2: Do this again to estimate species-weighted SD ###

modelColonies2<-list()

for (i in 1:length(actFiles)) {

print(paste0("Mapping file SD", i, "/", length(actFiles), "..."))

#### Open file i ###
actSub2<-read.csv(actFiles[i])
actSub2<-subset(actSub2, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))
energySub2<-read.csv(energyFiles[i], nrow=10)

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
actSub2$month<-factor(actSub2$month)
actSub2$species<-input_file1
#actSub2$weight<-ifelse(actSub2$NoBirds>0, 1/actSub2$NoBirds, 0)
actSub3<-actSub2%>%
dplyr::ungroup() %>%
dplyr::group_by(species, index) %>%
dplyr::summarise(NoBirds=sum(NoBirds, na.rm=TRUE), timeFlight=sum(timeFlight, na.rm=TRUE), 
timeRestWater=sum(timeRestWater, na.rm=TRUE), timeLand=sum(timeLand, na.rm=TRUE), timeForage=sum(timeForage, na.rm=TRUE), timeActive=sum(timeActive, na.rm=TRUE), colonies=max(colonies, na.rm=TRUE))
actSub3$weight<-ifelse(actSub3$NoBirds>0, 1/actSub3$NoBirds, 0)
actSub3<-actSub3 %>%
#dplyr::select(-c(x, y)) %>%
dplyr::left_join(actSubTotFinal, by=c("species", "index"))
actSub3$timeFlight<-(actSub3$timeFlight - actSub3$timeFlightMean)^2*actSub3$weight.x
actSub3$timeLand<-(actSub3$timeLand - actSub3$timeLandMean)^2*actSub3$weight.x
actSub3$timeForage<-(actSub3$timeForage - actSub3$timeForageMean)^2*actSub3$weight.x
actSub3$timeActive<-(actSub3$timeActive - actSub3$timeActiveMean)^2*actSub3$weight.x
actSub3$timeRestWater<-(actSub3$timeRestWater - actSub3$timeRestWaterMean)^2*actSub3$weight.x
actSub3$colonies<-ifelse(actSub3$NoBirds>0, 1, 0)

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
actSub3<-actSub3 %>%
dplyr::select(species, index, NoBirds, timeFlight, timeRestWater, timeLand, timeForage, timeActive, colonies, weight.x)
actSubTot2<-rbind(actSubTot2, actSub3)
actSubTot2<-actSubTot2 %>%
ungroup() %>%
dplyr::group_by(species, index) %>%
dplyr::summarise(NoBirds=sum(NoBirds, na.rm=TRUE), timeFlight=sum(timeFlight, na.rm=TRUE), timeRestWater=sum(timeRestWater, na.rm=TRUE), 
timeLand=sum(timeLand, na.rm=TRUE), timeForage=sum(timeForage, na.rm=TRUE), timeActive=sum(timeActive, na.rm=TRUE), colonies=sum(colonies, na.rm=TRUE), weight=sum(weight.x, na.rm=TRUE))

}

}

# Estimate species-weighted SD for the next calculation
actSubTot2<-actSubTot2 %>%
ungroup() 
actSubTot2[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater")]<-sqrt(actSubTot2[c("timeFlight", "timeForage", "timeActive", "timeLand", "timeRestWater")]/actSubTot2$weight)
actSubTot2[is.na(actSubTot2)] <- 0
actSubTotFinal2<-actSubTot2 %>%
#dplyr::select(-c(timeFlight_sd, timeActive_sd, timeRestWater_sd, timeLand_sd, timeActive_sd)) %>%
rename(timeFlight_sd=timeFlight, timeForage_sd=timeForage, timeActive_sd=timeActive, timeLand_sd=timeLand, timeRestWater_sd=timeRestWater) %>%
dplyr::left_join(actSubcoords, by=c("index")) %>%
rename(NoBirdsTot=NoBirds)

#actSubTotFinal$month<-factor(actSubTotFinal$month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))
#actSubTotFinal2$month<-factor(actSubTotFinal2$month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))

# Scale : find maximum & minimum values so that the fill scale is the same on the mean & SD plots #
actSubTotFinal3 <- actSubTotFinal2 %>%
  mutate(across(where(is.numeric), ~ replace(., is.infinite(.), 0)))
actMultiPops<-subset(actSubTotFinal, colonies>0)
actMultiPops2<-subset(actSubTotFinal3, colonies>0)
minFlight1<-min(actMultiPops$timeFlightMean)
minFlight2<-min(actMultiPops2$timeFlight_sd)
maxFlight1<-max(actMultiPops$timeFlightMean)
maxFlight2<-max(actMultiPops2$timeFlight_sd)
minFlight<-ifelse(minFlight1<minFlight2, minFlight2, minFlight2)
maxFlight<-ifelse(maxFlight1>maxFlight2, maxFlight1, maxFlight2)

# make map 1
flight1<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeFlightMean), color="black") +
  #geom_text(data=subset(actSubTotFinal, NoBirdsTot>1& timeFlightMean >0), aes(x=x, y=y, label=NoBirdsTot), cex=1.5) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time flight (hrs)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #scale_fill_viridis_c(limits=c(minFlight, maxFlight), option="plasma") +
  #+ scale_fill_viridis_c(option = "plasma") 
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")

# SD #

flight2<-ggplot() +
  #geom_tile(data=filter(actSubTotFinal, timeFlightMean>0), aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeFlight_sd), color="black") +
  #geom_text(data=subset(actSubTotFinal2, NoBirdsTot>1& timeFlight_sd >0), aes(x=x, y=y, label=colonies), cex=1) +
  #scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent flight (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #scale_fill_gradientn('Time spent flight (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(minFlight, maxFlight)) +
  #scale_fill_viridis_c(limits=c(minFlight, maxFlight)) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")

# make map 1

minRest1<-min(actMultiPops$timeRestWaterMean)
minRest2<-min(actMultiPops2$timeRestWater_sd)
maxRest1<-max(actMultiPops$timeRestWaterMean)
maxRest2<-max(actMultiPops2$timeRestWater_sd)
minRest<-ifelse(minRest1<minRest2, minRest2, minRest2)
maxRest<-ifelse(maxRest1>maxRest2, maxRest1, maxRest2)

RestWater1<-ggplot() +
  #geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeRestWaterMean)) +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeRestWaterMean), color="black") +
  #geom_text(data=filter(actSubTotFinal, NoBirdsTot>1& timeRestWaterMean >0), aes(x=x, y=y, label=NoBirdsTot), cex=1.5) +
  scale_fill_gradientn('Time rest (hrs)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #scale_fill_viridis_c(limits=c(minRest, maxRest)) +
  #scale_fill_gradientn('Time spent RestWater (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
   coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")


# SD #

RestWater2<-ggplot() +
  #geom_tile(data=filter(actSubTotFinal, timeRestWaterMean>0), aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeRestWater_sd), color="black") +
  #geom_text(data=filter(actSubTotFinal2, timeRestWater_sd>0 & colonies >1), aes(x=x, y=y, label=colonies), cex=1) +
  scale_fill_gradientn('Time spent rest (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
   #scale_fill_gradientn('Time spent rest (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(min(actSubTotFinal$timeRestWaterMean), max(actSubTotFinal$timeRestWaterMean))) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")

# make map 1

if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {

# make map 1

minForage1<-min(actMultiPops$timeForageMean)
minForage2<-min(actMultiPops2$timeForage_sd)
maxForage1<-max(actMultiPops$timeForageMean)
maxForage2<-max(actMultiPops2$timeForage_sd)
minForage<-ifelse(minForage1<minForage2, minForage2, minForage2)
maxForage<-ifelse(maxForage1>maxForage2, maxForage1, maxForage2)

Forage1<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeForageMean), color="black") +
  #geom_text(data=filter(actSubTotFinal, NoBirdsTot>1 & timeForageMean>0), aes(x=x, y=y, label=NoBirdsTot), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent Forage (hours)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  #scale_fill_viridis_c(limits=c(minForage, maxForage)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")


# Mean behaviour #

Forage2<-ggplot() +
  #geom_tile(data=filter(actSubTotFinal, timeForageMean>0), aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeForage_sd), color="black") +
  #geom_text(data=filter(actSubTotFinal2, timeForage_sd>0 & NoBirdsTot>1), aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent Forage (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minForage, maxForage)) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  #scale_fill_viridis_c(limits=c(minForage, maxForage)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")


} else {

# make map 1

minActive1<-min(actMultiPops$timeActiveMean)
minActive2<-min(actMultiPops2$timeActive_sd)
maxActive1<-max(actMultiPops$timeActiveMean)
maxActive2<-max(actMultiPops2$timeActive_sd)
minActive<-ifelse(minActive1<minActive2, minActive2, minActive2)
maxActive<-ifelse(maxActive1>maxActive2, maxActive1, maxActive2)

Active1<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=timeActiveMean), color="black") +
  #geom_text(data=filter(actSubTotFinal, NoBirdsTot>1 & timeActiveMean>0), aes(x=x, y=y, label=NoBirdsTot), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time Active (hrs)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  #scale_fill_viridis_c(limits=c(minActive, maxActive)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")

# Mean behaviour #

Active2<-ggplot() +
  #geom_tile(data=filter(actSubTotFinal, timeActiveMean>0), aes(x=x, y=y)) +
  geom_tile(data=actSubTotFinal2, aes(x=x, y=y, fill=timeActive_sd), color="black") +
  #geom_text(data=filter(actSubTotFinal2, timeActive_sd>0 & NoBirdsTot>1), aes(x=x, y=y, label=colonies), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent Active (SD)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #scale_fill_viridis_c(limits=c(minActive, maxActive)) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")
  
  }

### Add total time spent in each square ###

actSubTotSum<-actSubTotFinal %>%
ungroup() %>%
dplyr::group_by(species, x, y) %>%
dplyr::summarise(totTime=sum(timeFlightMean, timeActiveMean, timeRestWaterMean, timeLandMean, timeForageMean))

totTime<-ggplot() +
  geom_tile(data=actSubTotSum, aes(x=x, y=y, fill=totTime), color="black") +
  #geom_text(data=filter(actSubTotFinal, NoBirdsTot>1 & timeActiveMean>0), aes(x=x, y=y, label=NoBirdsTot), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Total time (hrs)', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  #scale_fill_viridis_c(limits=c(minActive, maxActive)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")


### Number of birds per square ###

birdNos<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=NoBirdsTot), color="black") +
  #geom_text(data=filter(actSubTotFinal, NoBirdsTot>1 & timeActiveMean>0), aes(x=x, y=y, label=NoBirdsTot), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Number of birds', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  #scale_fill_viridis_c(limits=c(minActive, maxActive)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1]))+
  theme(legend.position="bottom")

### Total colonies per square ###

colonyNos<-ggplot() +
  geom_tile(data=actSubTotFinal, aes(x=x, y=y, fill=colonies), color="black") +
  #geom_text(data=filter(actSubTotFinal, NoBirdsTot>1 & timeActiveMean>0), aes(x=x, y=y, label=NoBirdsTot), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Number of pops', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  #scale_fill_viridis_c(limits=c(minActive, maxActive)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(actMultiPops$x)-400000, max(actMultiPops$x) + 400000), ylim=c(min(actMultiPops$y)-400000, max(actMultiPops$y) + 400000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(energySub$species[1])) +
  theme(legend.position="bottom")

# Save plots

if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {

pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_allMean.pdf"),width=10)
grid.arrange(flight1, RestWater1, Forage1, totTime, colonyNos, birdNos, nrow=2)
dev.off()

pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_allSD.pdf"),width=10)
grid.arrange(flight2, RestWater2, Forage2, totTime, colonyNos, birdNos, nrow=2)
dev.off()


} else {

pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_allMean.pdf"), width=10)
grid.arrange(flight1, RestWater1, Active1, totTime, colonyNos, birdNos, nrow=2)
dev.off()

pdf(paste0("results/figures/speciesMaps/activity/", energySub$species[1], "_allSD.pdf"), width=10)
grid.arrange(flight2, RestWater2, Active2, totTime, colonyNos, birdNos, nrow=2)
dev.off()


}


}

# Save outputs

output_file1 <- args[2]
write.csv(actSubTotFinal, file=output_file1, row.names=FALSE)

output_file2 <- args[3]
write.csv(actSubTotFinal2, file=output_file2, row.names=FALSE)


