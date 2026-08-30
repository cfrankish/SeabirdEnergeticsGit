## Here we add species-specific diets & energy content so that we can relate prey requirements to energy demand ##

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
input_file1 <- args[1] # This will read in a species name
print(input_file1)
#input_file1<-"Commonguillemot"

# Now we determine the list of activity/energy data to loop through
allFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/", full.names=TRUE)
energyFiles<-allFiles[grep("energyMap_monthly_map_v1", allFiles)] # extract energy files
energyFilesSpecies<-energyFiles[grep(input_file1, energyFiles)] # Subset to species of interest

#### Step 1: determine possible diet for our species #### (start by simple)

# Main functional groups, AE and energy content in Kj.g of wet weight #
seabirdDiet<-data.frame(functionalGroups=c("Forage fish", "Macrozooplankton", "Herbivorous zooplankton", "Other"))
seabirdDiet$AE<-0.75 # Assimilation efficiency -> assume same for all seabird species & prey groups for now
seabirdDiet$energyDensity<-c(4.5, 4, 3.1, 2) # (kJ per g wet mass) #
# Energy density for forage fish based on this paper for Capelin: https://d1wqtxts1xzle7.cloudfront.net/46124708/Assimilation_efficiency_of_adult_Kittiwa20160601-4336-10nyk6f-libre.pdf?1464771742=&response-content-disposition=inline%3B+filename%3DAssimilation_efficiency_of_adult_Kittiwa.pdf&Expires=1776686451&Signature=FiV260FL8sA6XfT611m3V764hctijC01w2o6Cx9yHcrlMSq-hw3PRJQT-duIKO3~V9lh63q3HFgDZO2kv2tTPPH54yr9~xc7cSAmUGZTisvbtVaO5i5Dd3UUGIHqehy~dG4AyzRfAQYGZNFpN~s2ZgO--emdvI~GM88W3jbDydiJf-2EEbz3bcBav70q9l~4AdDjWEXb-2UfqHWx9-tVWGDvJTtJd4hx6hpYOtO1nDklqsIXrmcqDydleX1lijaaFGHDoUJENqx3U12SVPfzjL4mF5-9IwM9wDc7A4GMRiFe2UDN0vQHQX0rP4u97wOOYWq-IyvckPnI49hu9pVL5w__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA
# Energy density for herb. zoo. based on this paper for Calanus: https://www.sciencedirect.com/science/article/pii/S0079661116302105#s0010
# Energy density for macro. zoo. based on this paper for amphipods: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190345

# Determine proportion of each class in diet #
seabirdSpecies<-data.frame(species=rep(c("Blackleggedkittiwake", "Northernfulmar", "Atlanticpuffin", "Littleauk", "Commonguillemot", "Brunnichsguillemot"), each=4))
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

# Subset to species of interest
seabirdSpeciesDiet_sub<-subset(seabirdSpeciesDiet, species==input_file1)

#### Step 2: Determine area of cells in raster ####

# Figure out resolution of every 0.25*0.25 grid #

energySpatial<-read.csv(energyFilesSpecies[1]) # open up specific raster

energySpatialSub<-subset(energySpatial, month==9) # subset to a random month

energyRast<-energySpatialSub %>%
dplyr::select(x, y, energyPopkJ_mean) %>%
rename(z=energyPopkJ_mean)
energyRast<-rasterFromXYZ(energyRast) # Turn into a raster

area_raster <- area(energyRast) # calculate area of every cell

area_raster_df<-as.data.frame(area_raster, xy=TRUE) # Turn into a data frame

colnames(area_raster_df)<-c("x", "y", "areaKm2") # Change col names

# Join this to main dataset? 
energySpatialArea<-energySpatial %>%
dplyr::left_join(area_raster_df, by=c("x", "y"))

#### Step 3: Open up energy demand file & add a new column for estimating mass of diet items ###

# Determine monthly prey requirements #
fishInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Forage fish")
energySpatialArea$ForageFish_g<-((energySpatialArea$energyPopkJ_mean/(fishInfo$AE[1]*fishInfo$energyDensity))*fishInfo$propGroup)/energySpatialArea$areaKm2

herbZooInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Herbivorous zooplankton")
energySpatialArea$HerbZoo_g<-((energySpatialArea$energyPopkJ_mean/(herbZooInfo$AE[1]*herbZooInfo$energyDensity))*herbZooInfo$propGroup)/energySpatialArea$areaKm2

macroZooInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Macrozooplankton")
energySpatialArea$MacroZoo_g<-((energySpatialArea$energyPopkJ_mean/(macroZooInfo$AE[1]*macroZooInfo$energyDensity))*macroZooInfo$propGroup)/energySpatialArea$areaKm2





