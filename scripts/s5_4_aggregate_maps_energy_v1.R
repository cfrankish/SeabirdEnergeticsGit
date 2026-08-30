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
library(biscale)
library(cowplot)
library(classInt)
library(spdep)

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
allFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp4/", full.names=TRUE)
speciesFiles<-allFiles[grep(input_file1, allFiles)] # Subset to species-specific ones
#print("species")
#actFiles<-speciesFiles[grep("activityMap.csv", speciesFiles)] # extract activity files
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

#### Step 1: Add energy data together (initial visualization ####

# make empty list to save results in
energySpeciesMonth<-list()

# Determine unique months I would like to loop through #
months<-c(9, 10, 11, 12, 1, 2, 3, 4)

for (i in 1:length(months)) {

print(paste0("Summing energy accross pops for month ", i, "...")) 

# Now we loop through the energy files #

for (j in 1:length(energyFiles)) {

print(paste0("File ", j, "/", length(energyFiles)))

# open energy file j 
energy_j<-fread(energyFiles[j])

# Subset to month i
energy_j_i<-subset(energy_j, month==months[i])

# Add an index for joining because I'm unsure whether floating points join well
energy_j_i<-energy_j_i %>%
ungroup() %>%
dplyr::arrange(x, y) %>%
dplyr::mutate(index=row_number()) %>%
dplyr::select(-c(x, y))

# Make a master data frame for joining the others on
if (j==1) {
energySpecies<-energy_j_i 
} else {
energySpecies<-rbind(energySpecies, energy_j_i)
energySpecies<-energySpecies %>%
ungroup() %>%
dplyr::mutate(colony="All") %>%
dplyr::group_by(index, species, colony, month, weight, sst) %>%
dplyr::summarise(energyPopkJ_mean=sum(energyPopkJ_mean), energyPopkJ_sd=mean(energyPopkJ_sd), energyPopkJ_se=mean(energyPopkJ_se), NoBirds_mean=sum(NoBirds_mean))

}

}

# Now we re-add the x & y coordinates
energy_j<-fread(energyFiles[1]) # Open random file

# Subset to month i
energy_j_i<-subset(energy_j, month==months[i])

# Add an index for joining because I'm unsure whether floating points join well
energy_j_i<-energy_j_i %>%
ungroup() %>%
dplyr::arrange(x, y) %>%
dplyr::mutate(index=row_number()) %>%
dplyr::select(x, y, index)

# Join to main data frame
energySpeciesCoords<-energySpecies %>%
dplyr::left_join(energy_j_i, by=c("index"))

# Save monthly value
energySpeciesMonth<-rbind(energySpeciesMonth, energySpeciesCoords)

}

### Step 2: plot outputs ####

print("Step 2: plotting outputs...")

scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# So I will make seperate datasets for energy expenditure, no of birds and the difference between these two

energySpecies1<-energySpeciesMonth %>%
dplyr::ungroup() %>%
dplyr::select(species, x, y, month, energyPopkJ_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric="Energy_expenditure") %>%
dplyr::select(-energyPopkJ_mean)

energySpecies2<-energySpeciesMonth %>%
dplyr::ungroup() %>%
dplyr::select(species, x, y, month, NoBirds_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(NoBirds_mean)) %>%
dplyr::mutate(metric="No_Birds") %>%
dplyr::select(-NoBirds_mean)

energySpeciesAll<-rbind(energySpecies1, energySpecies2)

energySpecies3<-energySpeciesMonth %>%
dplyr::ungroup() %>%
dplyr::select(species, x, y, month, energyPopkJ_mean, NoBirds_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled1=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric_scaled2=scale01(NoBirds_mean)) %>%
dplyr::mutate(metric_diff=metric_scaled1 - metric_scaled2)%>%
dplyr::mutate(metric="Energy_expenditure - No_Birds") 

# Now we plot them side by side by month

# Plot these together #

pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/energy/", energySpeciesMonth$species[1], "_energy_v1.pdf"), width=10, height=7)

monthsLoop<-unique(energySpeciesMonth$month)

for (m in 1:length(monthsLoop)) {

plot1_energy<-ggplot() +
geom_tile(data=filter(energySpeciesAll, metric=="Energy_expenditure" & month==monthsLoop[m]), aes(x=x, y=y, fill=metric_scaled)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, 1)) +
facet_wrap(~metric) +
theme_bw()

plot1_density<-ggplot() +
geom_tile(data=filter(energySpeciesAll, metric=="No_Birds"& month==monthsLoop[m]), aes(x=x, y=y, fill=metric_scaled)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, 1)) +
facet_wrap(~metric) +
theme_bw()

# Now I make a plot showing the difference between the two

plot2_diff<-ggplot() +
geom_tile(data=filter(energySpecies3, month==monthsLoop[m]), aes(x=x, y=y, fill=metric_diff)) +
scale_fill_gradientn(
  colours = c(
    "#364B9A", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#A50026"
  ),
  limits = c(min(energySpecies3$metric_diff, na.rm=TRUE), max(energySpecies3$metric_diff, na.rm=TRUE)),
  values = scales::rescale(c(min(energySpecies3$metric_diff, na.rm=TRUE), 0, max(energySpecies3$metric_diff, na.rm=TRUE)))
) +
facet_wrap(~metric) +
theme_bw()

# Finally we make an SST plot? #

energySpeciesMonth$metric<-"SST"

plot3_sst<-ggplot() +
geom_tile(data=filter(energySpeciesMonth, month==monthsLoop[m]), aes(x=x, y=y, fill=sst)) +
scale_fill_gradientn(
  colours = c(
    "#364B9A", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#A50026"
  ),
  limits = c(min(energySpeciesMonth$sst, na.rm=TRUE), max(energySpeciesMonth$sst, na.rm=TRUE)),
  values = scales::rescale(c(min(energySpeciesMonth$sst, na.rm=TRUE), 5, max(energySpeciesMonth$sst, na.rm=TRUE)))
) +
facet_wrap(~metric) +
theme_bw()

grid.arrange(plot1_energy, plot1_density, plot2_diff, plot3_sst, nrow=2)

}

dev.off()

# Save outputs

print("Saving outputs...")

output_file1 <- args[2]
write.csv(energySpeciesMonth, file=output_file1, row.names=FALSE)



