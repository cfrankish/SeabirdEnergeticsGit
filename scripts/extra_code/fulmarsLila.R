## Here we map energetics using static activity (just changes by month but spatial SST ##

library(dplyr)
library(data.table)
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
library(terra)
#library(biscale)
#library(cowplot)
#library(classInt)
#library(spdep)

#### Step 0: setting up basic conditions ####

# Set-up number of iterations...
print(paste0("Mapping energy..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1] # This will read in a species name
print(input_file1)
#input_file1<-"Commonguillemot"

# Now we determine the list of activity/energy data to loop through
allFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp5/", full.names=TRUE)
speciesFiles<-allFiles[grep(input_file1, allFiles)] # Subset to species-specific ones
#print("species")
#actFiles<-speciesFiles[grep("activityMap.csv", speciesFiles)] # extract activity files
#print("act")
energyFiles<-speciesFiles[grep("energyMap_v1", speciesFiles)] # extract energy files
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
country <- (ncvar_get(nc,"colonyCountry")) # Colony code
colLon<-(ncvar_get(nc,"colonyLongitude"))
colLat<-(ncvar_get(nc,"colonyLatitude"))
meta <- data.frame(modelcolonies, colonies, colcode, country, colLon, colLat) # Make some metadata so easier to understand raster structure

#### Step 1: Add energy data together (initial visualization ####

# Figure out if it spends any time in NOrwegian EEZ
eez <- st_read("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/data/eez_norway/eez.shp")
eez_vector<-vect(eez)

# make empty list to save results in
energySpeciesMonth<-list() # Spatial maps
energyMonthCI_all<-list() # Temporal sums (for adding between populations)
colonyCoords<-list()

# Now we loop through the energy files #

for (j in 1:length(energyFiles)) {

print(paste0("File ", j, "/", length(energyFiles)))

# open energy file j 
energy_j<-fread(energyFiles[j])

# Subset to month i
energy_j_i<-subset(energy_j, month %in% c(11, 12, 1, 2))

# Sum accross study period for simplicity
energy_j_i_sum<-energy_j_i %>%
dplyr::group_by(species, modelcolony, colony, x, y) %>%
dplyr::summarise(energyPopkJ_mean=sum(energyPopkJ_mean, na.rm=TRUE), totBirds=sum(NoBirds_mean, na.rm=TRUE))

# Figure out if it's Norwegian or not #
metaCountry<-subset(meta, colonies==energy_j_i_sum$colony[1])

rasterdf<-energy_j_i_sum %>%
ungroup() %>%
dplyr::select(x, y, energyPopkJ_mean) %>%
rename(z=energyPopkJ_mean) # turn energy df into raster
energyRast<-rast(rasterdf)
birdsNorway<-sum(terra::extract(energyRast, vect(eez)))

# Essentially if it has no birds in norwegian EEZ and it's not from NOrway then next

if (birdsNorway==0 & !metaCountry$country[1] %in% c("Norway")) {
next}

#Save colony coordinates #
colonyCoords<-rbind(colonyCoords, metaCountry)

# Add an index for joining because I'm unsure whether floating points join well
energy_j_i_sum<-energy_j_i_sum %>%
ungroup() %>%
dplyr::arrange(x, y) %>%
dplyr::mutate(index=row_number()) %>%
dplyr::select(-c(x, y))

# Sum monthly energy many times #
energyTest<-subset(energy_j_i, !is.na(energyPopkJ_mean))
energyEstimatePop<-energyEstimate(energyTest, 50)
energyMonthCI_all<-rbind(energyMonthCI_all, energyEstimatePop)

# Make a master data frame for joining the others on
if (j==1) {
energySpecies<-energy_j_i_sum 
} else {
energy_j_i_sum$colony<-as.character(energy_j_i_sum$colony)
energy_j_i_sum$colony<-"All"
energySpecies<-rbind(energySpecies, energy_j_i_sum)
energySpecies<-energySpecies %>%
ungroup() %>%
dplyr::mutate(colony="All") %>%
dplyr::group_by(index, species, colony) %>%
dplyr::summarise(energyPopkJ_mean=sum(energyPopkJ_mean, na.rm=TRUE), NoBirds_mean=sum(totBirds, na.rm=TRUE))

}

}

# Now we re-add the x & y coordinates
energy_j<-fread(energyFiles[1]) # Open random file

# Subset to month i
energy_j_i<-subset(energy_j, month==9)

# Add an index for joining because I'm unsure whether floating points join well
energy_j_i<-energy_j_i %>%
ungroup() %>%
dplyr::arrange(x, y) %>%
dplyr::mutate(index=row_number()) %>%
dplyr::select(x, y, index)

# Join to main data frame
energySpeciesCoords<-energySpecies %>%
dplyr::left_join(energy_j_i, by=c("index"))

### Step 2: plot outputs ####

pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/energy/", energySpeciesMonth$species[1], "_energy_v1.pdf"), width=10, height=7)

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

EnergyExpenditure<-ggplot() +
 geom_tile(data=energySpeciesCoords, aes(x=x, y=y, fill=energyPopkJ_mean/1e06)) +
 scale_fill_gradientn(colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
 geom_sf(data=coast) +
 coord_sf(crs=projection_84, xlim=c(min(energySpeciesCoords$x), max(energySpeciesCoords$x)), ylim=c(min(energySpeciesCoords$y), max(energySpeciesCoords$y))) +
 theme_minimal() +
 geom_point(data=colonyCoords, aes(x=colLon, y=colLat, color="Populations"), size=0.5, shape=19) +
 scale_color_manual(values=c("yellow")) +  
 ggtitle("Fulmar energy expenditure (Nov-Feb)") +
 labs(fill="Energy expenditure (kJx1e06)") +
 theme(legend.position="bottom") +
 labs(color="", shape="") +
 xlab("") +
 ylab("")
 
pdf("./results/figures/forLila/populationEnergy.pdf")
plot(EnergyExpenditure)
dev.off()
 
### Now we extract variation in activity budgets & energy expenditure for fulmars from Norway ###

average_activity<-read.csv("./results/tables/main/table2_budgets_species.csv") # species level
average_activity_colony<-read.csv("./results/tables/main/table3_budgets_population.csv") # colony-level

# Make some kind of summary #

average_activity_all_colony<-average_activity_colony %>%
  ungroup() %>%
  dplyr::filter(species=="Northern fulmar" & colony %in% c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Karmøy")) %>%
  dplyr::group_by(species, colony, weekNo) %>%
  dplyr::summarise(reps=n_distinct(rep), Flight=mean(meanFlight), sdFlight=sd(meanFlight), seFlight=sdFlight/sqrt(reps),
                   Forage=mean(meanForage), sdForage=sd(meanForage), seForage=sdForage/sqrt(reps),
                   Land=mean(meanLand), sdLand=sd(meanLand), seLand=sdLand/sqrt(reps),
                   Rest=mean(meanRest), sdRest=sd(meanRest), seRest=sdRest/sqrt(reps),
                   Active=mean(meanActive), sdActive=sd(meanActive), seActive=sdActive/sqrt(reps),
                   SST=mean(meansst), sdSST=sd(meansst), seSST=sdSST/sqrt(reps), SST_max=max(maxsst), SST_min=min(minsst),
                   energy=mean(meanDEE), sdenergy=sd(meanDEE), seenergy=sdenergy/sqrt(reps)) 

# Here is the species average that we will also add #

# Make intermediary data frame
average_activity_all<-average_activity %>%
  ungroup() %>%
  dplyr::group_by(species, weekNo) %>%
  dplyr::summarise(reps=n_distinct(rep), Flight=mean(meanFlight), sdFlight=sd(meanFlight), seFlight=sdFlight/sqrt(reps),
                   Forage=mean(meanForage), sdForage=sd(meanForage), seForage=sdForage/sqrt(reps),
                   Land=mean(meanLand), sdLand=sd(meanLand), seLand=sdLand/sqrt(reps),
                   Rest=mean(meanRest), sdRest=sd(meanRest), seRest=sdRest/sqrt(reps),
                   Active=mean(meanActive), sdActive=sd(meanActive), seActive=sdActive/sqrt(reps), meanBirds=mean(birdsTot),
                   SST=mean(meansst), sdSST=sd(meansst), seSST=sdSST/sqrt(reps), SST_max=max(maxsst), SST_min=min(minsst),
                   energy=mean(meanDEE), sdenergy=sd(meanDEE), seenergy=sdenergy/sqrt(reps))

average_activity_all_species<-average_activity_all %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::filter(species=="Northern fulmar")												 

# Determine study period #
startDate = "2021-09-15"
endDate = "2022-04-15"

# Create list of weeks to roll through
dates<-data.frame(dateKeep=seq(as.Date(startDate), as.Date(endDate), 1))
dates$doy<-1:nrow(dates)
dates$month<-as.numeric(substr(dates$date, 6, 7))
dates$day<-as.numeric(substr(dates$date, 9, 10))

# Add week number for summarizing information
dates_weekly<-dates %>%
  dplyr::mutate(weekNo=ceiling(doy/7)) %>%
  dplyr::group_by(weekNo) %>%
  dplyr::mutate(days=n_distinct(dateKeep)) %>%
  dplyr::filter(days==7) %>%
  dplyr::select(-days)
  
# Activity plot #

startMonth<-dates_weekly %>%
  dplyr::filter(day==1)

FigureS9<-average_activity_all_colony %>%
  ggplot() +
  geom_pointrange(aes(x=weekNo, y=Flight, colour="Flight",  ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight), cex=0.1, alpha=0.05) +
  geom_line(aes(colour="Flight", x=weekNo, y=Flight, group=colony), alpha=0.05) +
  geom_ribbon(aes(x=weekNo, y=Flight, ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight, fill="Flight", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2) +
  labs(color="", fill="") +
  #geom_pointrange(data=filter(average_activity_all_colony, Active>0) ,aes(x=weekNo, y=Active, colour="Active", ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive), cex=0.1, alpha=0.05) +
  #geom_line(data=filter(average_activity_all_colony, Active>0) ,aes(colour="Active", x=weekNo, y=Active, group=colony), alpha=0.05) +
  #geom_ribbon(data=filter(average_activity_all_colony ,Active>0) ,aes(x=weekNo, y=Active, ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive, fill="Active", group=colony), alpha=0.05) +
  geom_pointrange(data=filter(average_activity_all_colony, Forage>0) ,aes(x=weekNo, y=Forage, colour="Forage", ymin=Forage-1.96*seForage, ymax=Forage + 1.96*seForage, group=colony), cex=0.1, alpha=0.05) +
  geom_line(data=filter(average_activity_all_colony, Forage>0) ,aes(colour="Forage", x=weekNo, y=Forage, group=colony), alpha=0.05) +
  geom_ribbon(data=filter(average_activity_all_colony, Forage>0) ,aes(x=weekNo, y=Forage, ymin=Forage-1.96*seForage, ymax=Forage + 1.96*seForage, fill="Forage", group=colony), alpha=0.05) +
  geom_pointrange(aes(x=weekNo, y=Rest, colour="Water", ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest, group=colony), cex=0.1, alpha=0.05) +
  geom_line(aes(colour="Water", x=weekNo, y=Rest, group=colony), alpha=0.05) +
  geom_ribbon(aes(x=weekNo, y=Rest, ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest, fill="Water", group=colony), alpha=0.05) +
  geom_pointrange(aes(x=weekNo, y=Land, colour="Land", ymin=Land-1.96*seLand, ymax=Land + 1.96*seLand, group=colony), alpha=0.05, cex=0.1) +
  geom_line(aes(colour="Land", x=weekNo, y=Land, group=colony), alpha=0.05) +
  geom_ribbon(aes(x=weekNo, y=Land, ymin=Land-1.96*seLand, ymax=Land + 1.96*seLand, fill="Land", group=colony), alpha=0.05) +
  theme_bw() +
  geom_pointrange(data=average_activity_all_species, aes(x=weekNo, y=Land, colour="Land", ymin=Land-1.96*seLand, ymax=Land + 1.96*seLand), cex=0.1) +
  geom_line(data=average_activity_all_species,aes(colour="Land", x=weekNo, y=Land)) +
  geom_ribbon(data=average_activity_all_species,aes(x=weekNo, y=Land, ymin=Land-1.96*seLand, ymax=Land + 1.96*seLand, fill="Land"), alpha=0.2) +
  geom_pointrange(data=filter(average_activity_all_species, Forage>0), aes(x=weekNo, y=Forage, colour="Forage", ymin=Forage-1.96*seForage, ymax=Forage + 1.96*seForage), cex=0.1) +
  geom_line(data=filter(average_activity_all_species, Forage>0),aes(colour="Forage", x=weekNo, y=Forage)) +
  geom_ribbon(data=filter(average_activity_all_species, Forage>0),aes(x=weekNo, y=Forage, ymin=Forage-1.96*seForage, ymax=Forage + 1.96*seForage, fill="Forage"), alpha=0.2) +
  #geom_pointrange(data=filter(average_activity_all_species, Active>0), aes(x=weekNo, y=Active, colour="Active", ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive), cex=0.1) +
  #geom_line(data=filter(average_activity_all_species, Active>0),aes(colour="Active", x=weekNo, y=Active)) +
  geom_ribbon(data=filter(average_activity_all_species, Active>0),aes(x=weekNo, y=Active, ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive, fill="Active"), alpha=0.2) +
  geom_pointrange(data=filter(average_activity_all_species, Rest>0), aes(x=weekNo, y=Rest, colour="Water", ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest), cex=0.1) +
  geom_line(data=filter(average_activity_all_species, Rest>0),aes(colour="Water", x=weekNo, y=Rest)) +
  geom_ribbon(data=filter(average_activity_all_species, Rest>0),aes(x=weekNo, y=Rest, ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest, fill="Water"), alpha=0.2) +
  geom_pointrange(data=filter(average_activity_all_species, Flight>0), aes(x=weekNo, y=Flight, colour="Flight", ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight), cex=0.1) +
  geom_line(data=filter(average_activity_all_species, Flight>0),aes(colour="Flight", x=weekNo, y=Flight)) +
  geom_ribbon(data=filter(average_activity_all_species, Flight>0),aes(x=weekNo, y=Flight, ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight, fill="Flight"), alpha=0.2) +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Proportion of day spent in behaviour") +
  scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))+
  scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") 

pdf("./results/figures/forLila/activity.pdf", width=6, height=5)
grid.arrange(FigureS9)
dev.off()

# Add maximum & min countours #

Contours<-average_activity_all_colony %>%
  ungroup() %>%
  dplyr::mutate(upper=energy + 1.96*seenergy, lower=energy - 1.96*seenergy) %>%
  dplyr::group_by(species, weekNo) %>%
  dplyr::summarise(minEnergy=min(lower), maxEnergy=max(upper))

# And minimum contour #

FigureS12<-average_activity_all_colony %>%
  ggplot(aes(x=weekNo, y=energy)) +
  geom_pointrange(aes(colour="DEE", x=weekNo, y=energy, ymin=energy-1.96*seenergy, ymax=energy + 1.96*seenergy), cex=0.1, alpha=0.05) +
  geom_line(aes(colour="DEE", x=weekNo, y=energy, group=colony), alpha=0.05) +
  geom_ribbon(aes(x=weekNo, y=energy, ymin=energy-1.96*seenergy, ymax=energy + 1.96*seenergy, fill="DEE", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free_y") +
  labs(color="Behaviour", fill="Behaviour") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  #geom_pointrange(aes(x=weekNo, y=SST, ymin=(SST-seSST*1.96), ymax=(SST+seSST*1.96), colour="scaled SST", group=colony), cex=0.1, alpha=0.05) +
  #geom_line(aes(colour="scaled SST", x=weekNo, y=SST, group=colony), alpha=0.05) +
  #geom_ribbon(aes(x=weekNo, y=(SST), ymin=(SST-1.96*seSST), ymax=(SST + 1.96*seSST), fill="scaled SST", group=colony), alpha=0.05) +
  #scale_y_continuous(sec.axis = sec_axis(~.1, name="SST (degrees C)")) +
  ylab("DEE (kJ.day-1)") +
  labs(colour="", fill="") +
  scale_color_manual(values=c("#0072b2", "#E25822"))+
  scale_fill_manual(values=c( "#0072b2", "#E25822")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  #ggtitle("B) Energy budgets") +
  #geom_pointrange(data=average_activity_all_species,aes(x=weekNo, y=SST, ymin=(SST-seSST*1.96), ymax=(SST+seSST*1.96), colour="scaled SST"), cex=0.1, alpha=0.5) +
  #geom_line(data=average_activity_all_species,aes(colour="scaled SST", x=weekNo, y=SST), alpha=0.5) +
  #geom_ribbon(data=average_activity_all_species,aes(x=weekNo, y=(SST), ymin=(SST-1.96*seSST), ymax=(SST + 1.96*seSST), fill="scaled SST"), alpha=0.2) +
  geom_pointrange(data=average_activity_all_species, aes(colour="DEE", x=weekNo, y=energy, ymin=energy-1.96*seenergy, ymax=energy + 1.96*seenergy), cex=0.1) +
  geom_line(data=average_activity_all_species, aes(colour="DEE", x=weekNo, y=energy)) +
  geom_ribbon(data=average_activity_all_species, aes(x=weekNo, y=energy, ymin=energy-1.96*seenergy, ymax=energy + 1.96*seenergy, fill="DEE"), alpha=0.2) +
  geom_line(data=Contours, aes(colour="DEE", x=weekNo, y=maxEnergy), linetype="dashed", alpha=0.8) +
  geom_line(data=Contours, aes(colour="DEE", x=weekNo, y=minEnergy), linetype="dashed", alpha=0.8) 

pdf("./results/figures/forLila/energy.pdf", width=6, height=5)
plot(FigureS12)
dev.off()




