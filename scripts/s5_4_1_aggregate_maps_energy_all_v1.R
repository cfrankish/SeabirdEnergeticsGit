## Here we aggregate energy maps to make estimates accross all species ##

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
#input_file1 <- args[1] # This will read in a species name
#print(input_file1)
#input_file1<-"Commonguillemot"

# Now we determine the list of activity/energy data to loop through
allFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/", full.names=TRUE)
energyFiles<-allFiles[grep("energyMap_monthly_map_v1", allFiles)] # extract energy files

# Create a world map for plotting
print("Prepping map info")
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#### Step 0.2 prep some diet stuff for aggregation component ####

# Main functional groups, AE and energy content in Kj.g of wet weight #
seabirdDiet<-data.frame(functionalGroups=c("Forage fish", "Macrozooplankton", "Herbivorous zooplankton", "Other"))
seabirdDiet$AE<-0.75 # Assimilation efficiency -> assume same for all seabird species & prey groups for now
seabirdDiet$energyDensity<-c(4.5, 4, 3.1, 2) # (kJ per g wet mass) #
# Energy density for forage fish based on this paper for Capelin: https://d1wqtxts1xzle7.cloudfront.net/46124708/Assimilation_efficiency_of_adult_Kittiwa20160601-4336-10nyk6f-libre.pdf?1464771742=&response-content-disposition=inline%3B+filename%3DAssimilation_efficiency_of_adult_Kittiwa.pdf&Expires=1776686451&Signature=FiV260FL8sA6XfT611m3V764hctijC01w2o6Cx9yHcrlMSq-hw3PRJQT-duIKO3~V9lh63q3HFgDZO2kv2tTPPH54yr9~xc7cSAmUGZTisvbtVaO5i5Dd3UUGIHqehy~dG4AyzRfAQYGZNFpN~s2ZgO--emdvI~GM88W3jbDydiJf-2EEbz3bcBav70q9l~4AdDjWEXb-2UfqHWx9-tVWGDvJTtJd4hx6hpYOtO1nDklqsIXrmcqDydleX1lijaaFGHDoUJENqx3U12SVPfzjL4mF5-9IwM9wDc7A4GMRiFe2UDN0vQHQX0rP4u97wOOYWq-IyvckPnI49hu9pVL5w__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA
# Energy density for herb. zoo. based on this paper for Calanus: https://www.sciencedirect.com/science/article/pii/S0079661116302105#s0010
# Energy density for macro. zoo. based on this paper for amphipods: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190345

# Determine proportion of each class in diet #
seabirdSpecies<-data.frame(species=rep(c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot"), each=4))
seabirdSpecies$functionalGroups<-rep(c("Forage fish", "Macrozooplankton", "Herbivorous zooplankton", "Other"), 6)
seabirdSpecies$propGroup<-c(0.94, 0.06, 0, 0, # BLK https://www.int-res.com/articles/meps2007/349/m349p269.pdf (I made an average of the Other category)
0.93, 0, 0, 0.07, # NoFu # file:///C:/Users/caitlin.frankish/Downloads/s002270050613.pdf
0.66, 0.02, 0, 0.41, # AtPu : https://onlinelibrary.wiley.com/doi/pdf/10.1111/ibi.12272
0.12, 0.87, 0.01, 0,# LiAu: https://link.springer.com/article/10.1007/s00300-013-1379-4/tables/2
1, 0, 0, 0,# CoGu # https://www.tandfonline.com/doi/epdf/10.1080/17451000802279636?needAccess=true
0.39, 0.55, 0, 0.06) #BrGu: https://cdnsciencepub.com/doi/epdf/10.1139/cjz-2021-0120 # But there are 6% other species so I will jus attribute this to the two classes  

# Attach energy density values #
seabirdSpeciesDiet<-seabirdSpecies %>%
dplyr::left_join(seabirdDiet, by=c("functionalGroups"))

#### Step 1: Loop through species files, open & add ###

for (i in 1:length(energyFiles)) {

# Print status #
print(paste0("Opening file", i, "/", length(energyFiles)))

# open file i
energyFile_i<-fread(energyFiles[i])

# Find out area size of every cell
energySpatialSub<-subset(energyFile_i, month==9) # subset to a random month
energyRast<-energySpatialSub %>%
dplyr::select(x, y, energyPopkJ_mean) %>%
rename(z=energyPopkJ_mean)
energyRast<-rasterFromXYZ(energyRast) # Turn into a raster
area_raster <- area(energyRast) # calculate area of every cell
area_raster_df<-as.data.frame(area_raster, xy=TRUE) # Turn into a data frame
colnames(area_raster_df)<-c("x", "y", "areaKm2") # Change col names

# Join this to main dataset? 
energySpatialArea<-energyFile_i %>%
dplyr::left_join(area_raster_df, by=c("x", "y"))

# Subset to species of interest
seabirdSpeciesDiet_sub<-subset(seabirdSpeciesDiet, species==energyFile_i$species[1])

# Determine monthly prey requirements #
fishInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Forage fish")
energySpatialArea$ForageFish_g<-((energySpatialArea$energyPopkJ_mean/(fishInfo$AE[1]*fishInfo$energyDensity))*fishInfo$propGroup)/energySpatialArea$areaKm2

herbZooInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Herbivorous zooplankton")
energySpatialArea$HerbZoo_g<-((energySpatialArea$energyPopkJ_mean/(herbZooInfo$AE[1]*herbZooInfo$energyDensity))*herbZooInfo$propGroup)/energySpatialArea$areaKm2

macroZooInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Macrozooplankton")
energySpatialArea$MacroZoo_g<-((energySpatialArea$energyPopkJ_mean/(macroZooInfo$AE[1]*macroZooInfo$energyDensity))*macroZooInfo$propGroup)/energySpatialArea$areaKm2

# if it's the first file, we create a master file #
if (i==1) {

energyAll<-energySpatialArea %>%
dplyr::select(month, index, energyPopkJ_mean, NoBirds_mean, ForageFish_g, HerbZoo_g, MacroZoo_g, areaKm2) %>%
dplyr::mutate(speciesNo=ifelse(NoBirds_mean>0, 1, 0))

} else {

# Otherwise they are summed together #

energySpatialArea<-energySpatialArea %>%
dplyr::select(month, index, energyPopkJ_mean, NoBirds_mean, ForageFish_g, HerbZoo_g, MacroZoo_g, areaKm2)%>%
dplyr::mutate(speciesNo=ifelse(NoBirds_mean>0, 1, 0))

energyAll<-rbind(energyAll, energySpatialArea)

energyAll<-energyAll %>%
ungroup() %>%
dplyr::group_by(month, index) %>%
dplyr::summarise(energyPopkJ_mean=sum(energyPopkJ_mean, na.rm=TRUE), NoBirds_mean=sum(NoBirds_mean, na.rm=TRUE), ForageFish_g=sum(ForageFish_g, na.rm=TRUE), HerbZoo_g=sum(HerbZoo_g, na.rm=TRUE), MacroZoo_g=sum(MacroZoo_g, na.rm=TRUE),
speciesNo=sum(speciesNo, na.rm=TRUE), areaKm2=mean(areaKm2, na.rm=TRUE))

}

}

### Step 2: plot outputs ####

print("Step 2: plotting outputs...")

scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Open up random file to get x & y coordinates 
energyFile_i<-fread(energyFiles[1])

coords<-energyFile_i %>%
ungroup() %>%
dplyr::group_by(index) %>%
dplyr::slice(1) %>%
dplyr::select(index, x, y)

energySave<-energyAll %>%
dplyr::ungroup() %>%
dplyr::left_join(coords, by=c("index"))

# So I will make seperate datasets for energy expenditure, no of birds and the difference between these two

energySpecies1<-energyAll %>%
dplyr::ungroup() %>%
dplyr::left_join(coords, by=c("index")) %>%
dplyr::select(x, y, month, energyPopkJ_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric_scaled=ifelse(energyPopkJ_mean==0, NA, metric_scaled)) %>%
dplyr::mutate(metric="Energy_expenditure") %>%
dplyr::select(-energyPopkJ_mean)

energySpecies2<-energyAll %>%
dplyr::ungroup() %>%
dplyr::left_join(coords, by=c("index")) %>%
dplyr::select(x, y, month, NoBirds_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(NoBirds_mean)) %>%
dplyr::mutate(metric_scaled=ifelse(NoBirds_mean==0, NA, metric_scaled)) %>%
dplyr::mutate(metric="No_Birds") %>%
dplyr::select(-NoBirds_mean)

energySpeciesAll<-rbind(energySpecies1, energySpecies2)

energySpecies3<-energyAll %>%
dplyr::ungroup() %>%
dplyr::left_join(coords, by=c("index")) %>%
dplyr::select(x, y, month, energyPopkJ_mean, NoBirds_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled1=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric_scaled2=scale01(NoBirds_mean)) %>%
dplyr::mutate(metric_diff=metric_scaled1 - metric_scaled2)%>%
dplyr::mutate(metric_diff=ifelse(energyPopkJ_mean==0, NA, metric_diff)) %>%
dplyr::mutate(metric="Energy_expenditure - No_Birds") 

# Now we plot them side by side by month

# Plot these together #

pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/energy/allBirds_energy_v1.pdf"), width=10, height=7)

print("Saving PDF")

monthsLoop<-c(9, 10, 11, 12, 3, 4)

for (m in 1:length(monthsLoop)) {

print(paste0("Month", m))

plot1_energy<-ggplot() +
geom_tile(data=filter(energySpeciesAll, metric=="Energy_expenditure" & month==monthsLoop[m] & metric_scaled>0), aes(x=x, y=y, fill=metric_scaled)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, 1), na.value = "grey80") +
facet_wrap(~metric) +
theme_bw() +
ggtitle(paste0("Month ", monthsLoop[m]))

plot1_density<-ggplot() +
geom_tile(data=filter(energySpeciesAll, metric=="No_Birds"& month==monthsLoop[m] & metric_scaled >0), aes(x=x, y=y, fill=metric_scaled)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, 1)) +
facet_wrap(~metric) +
theme_bw()

# Now I make a plot showing the difference between the two

plot2_diff<-ggplot() +
geom_tile(data=filter(energySpecies3, month==monthsLoop[m] & metric_scaled1>0), aes(x=x, y=y, fill=metric_diff)) +
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

grid.arrange(plot1_energy, plot1_density, plot2_diff, nrow=2)

}

dev.off()

# Save outputs

print("Saving outputs...")

output_file1 <- args[1]
write.csv(energySave, file=output_file1, row.names=FALSE)