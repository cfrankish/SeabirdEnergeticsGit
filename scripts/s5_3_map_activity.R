library(dplyr)
library(ncdf4)
library(raster)
library(sp)
library(sf)
library(lubridate)
library(data.table)
library(tidyr)
library(ggplot2)
library(rnaturalearthdata)
library(rnaturalearth)

#### Step 0: setting up basic conditions ####

# Set-up number of iterations...
overall.iterations<-10 # how many times this is calculated per individual
print(paste0("Determining activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1] # This will read in a population-specific activity budget  

# open input file 1
actBudgets<-read.csv(input_file1) # actBudgets<-read.csv("tmp3/Commonguillemot_lundy_actbudgets_daily.csv")

### Step 1: Make world grid to loop through ###

print("Making world grid...")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Load world polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

# Reproject world to metric CRS
world_proj <- st_transform(world, crs = projection_NA)

# Define extent in projected units (meters) - here: North Atlantic area
bbox <- st_bbox(world_proj)
extent_proj <- extent(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])

# Create raster template in metric CRS with 10 km resolution
res_m <- 200000  # resolution in meters
r <- raster(extent_proj, res = res_m, crs = projection_NA)

# Rasterize
world_raster <- rasterize(world_proj, r, field = 1, background = NA)

# Land = 1, Sea = 10000
grid<-as.data.frame(world_raster, xy=TRUE)

# Set-up extent 
#fake.rast<-raster(xmn=(-180), xmx=180, ymn=(-60), ymx=90, res=5)
#values(fake.rast)<-0
#proj4string(fake.rast)<-projection_84
#WorldProject<-projectRaster(fake.rast, crs=projection_NA)
#grid<-as.data.frame(WorldProject, xy=TRUE)
grid$loxCol1<-0	
grid$NoBirds<-0
grid$timeFlight<-0 # time spent in flight in hours at colony 1
grid$timeForage<-0 # time spent foraging in hours at colony 1
grid$timeActive<-0 # time spent foraging in hours at colony 1
grid$timeLand<-0 # time spent foraging in hours at colony 2
grid$timeRestWater<-0 # time spent foraging in hours at colony 2
grid$timeFlight_sd<-0 # time spent in flight in hours at colony 1
grid$timeForage_sd<-0 # time spent foraging in hours at colony 1
grid$timeActive_sd<-0 # time spent foraging in hours at colony 1
grid$timeLand_sd<-0 # time spent foraging in hours at colony 2
grid$timeRestWater_sd<-0 # time spent foraging in hours at colony 2
grid$timeTotal<-0 # total time spent
grid$timeTotal_sd<-0 # total time spent SD
grid<-grid %>%
arrange(x, y) 

### Step 2: Calculate energy & map it many times ####

# Determine number of months to loop through
Months<-unique(actBudgets$month)

# Make an an empty list to save the results in
PopActivityFinal_All<-list()

for (i in 1:length(Months)) {

print(paste0("Saving for month", i))

# Subset to month i
monthSub<-Months[i]

# Cut actBudgets
actMonth<-subset(actBudgets, month==monthSub)

# Calculate total time in different behaviours per month
totalTime<-actMonth %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id, month) %>%
dplyr::summarise(totFlight=sum(tFlightMean), totActive=sum(tActiveMean), totRest=sum(tRestWaterMean), totLand=sum(tLandMean), totForage=sum(tForageMean), totHrs=sum(totFlight, totActive, totRest, totLand, totForage))

# Project coordinates
coordinates(actMonth)<-~meanLon + mean.lat
proj4string(actMonth)<-projection_84
actMonthTrans<-data.frame(spTransform(actMonth, projection_NA))

# Determine total birds accross month
totBirds<-n_distinct(actMonth$individ_id)

# Rename our grid
gridMonth<-grid

for (m in 1:nrow(grid)) {

# Subset grid x
gridSub<-gridMonth[m,]

resx<-res(world_raster)[1]
resy<-res(world_raster)[2]

# Subset coordinates which fit #
loxSub1<-subset(actMonthTrans,   coords.x1 > gridSub$x & coords.x1< gridSub$x + resx & coords.x2 > gridSub$y & coords.x2 < gridSub$y + resy)

if (nrow(loxSub1)>0) {

# Calculate time spent in flight
timeFlight<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tFlightMean), totTimeAll=sum(tFlightMean, tActiveMean, tRestWaterMean, tLandMean, tForageMean)) %>%
#dplyr::full_join(actMonthTrans, by=c("species", "colony", "individ_id")) %>%
#dplyr::left_join(totalTime, by=c("species", "colony", "individ_id")) %>%
#replace_na(list(totTime=0)) %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent foraging
timeForage<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tForageMean), totTimeAll=sum(tFlightMean, tActiveMean, tRestWaterMean, tLandMean, tForageMean)) %>%
#dplyr::full_join(actMonthTrans, by=c("species", "colony", "individ_id")) %>%
#dplyr::left_join(totalTime, by=c("species", "colony", "individ_id")) %>%
#replace_na(list(totTime=0)) %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeLand<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tLandMean), totTimeAll=sum(tFlightMean, tActiveMean, tRestWaterMean, tLandMean, tForageMean)) %>%
#dplyr::full_join(actMonthTrans, by=c("species", "colony", "individ_id")) %>%
#dplyr::left_join(totalTime, by=c("species", "colony", "individ_id")) %>%
#replace_na(list(totTime=0)) %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeActive<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tActiveMean), totTimeAll=sum(tFlightMean, tActiveMean, tRestWaterMean, tLandMean, tForageMean)) %>%
#dplyr::full_join(actMonthTrans, by=c("species", "colony", "individ_id")) %>%
#dplyr::left_join(totalTime, by=c("species", "colony", "individ_id")) %>%
#replace_na(list(totTime=0)) %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeRest<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tRestWaterMean), totTimeAll=sum(tFlightMean, tActiveMean, tRestWaterMean, tLandMean, tForageMean)) %>%
#dplyr::full_join(actMonthTrans, by=c("species", "colony", "individ_id")) %>%
#dplyr::left_join(totalTime, by=c("species", "colony", "individ_id")) %>%
#replace_na(list(totTime=0)) %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate total time
timeTot<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tFlightMean, tRestWaterMean, tForageMean, tActiveMean, tLandMean)) %>%
#dplyr::full_join(actMonthTrans, by=c("species", "colony", "individ_id")) 
dplyr::full_join(totalTime, by=c("species", "colony", "individ_id")) %>%
replace_na(list(totTime=0)) %>%
dplyr::mutate(propTime=totTime/totHrs)

#nas<-subset(timeLand, is.na(propTime))
#if(nrow(nas)>0) {break}

# Save results in grid for plotting
gridMonth$timeFlight[m]<-mean(timeFlight$propTime)
gridMonth$timeFlight_sd[m]<-sd(timeFlight$propTime)
gridMonth$timeForage[m]<-mean(timeForage$propTime)
gridMonth$timeForage_sd[m]<-sd(timeForage$propTime)
gridMonth$timeActive[m]<-mean(timeActive$propTime)
gridMonth$timeActive_sd[m]<-sd(timeActive$propTime)
gridMonth$timeLand[m]<-mean(timeLand$propTime)
gridMonth$timeLand_sd[m]<-sd(timeLand$propTime)
gridMonth$timeRestWater[m]<-mean(timeRest$propTime)
gridMonth$timeRestWater_sd[m]<-sd(timeRest$propTime)
gridMonth$timeTotal[m]<-mean(timeTot$propTime)
gridMonth$timeTotal_sd[m]<-sd(timeTot$propTime)
gridMonth$loxCol1[m]<-nrow(loxSub1)
gridMonth$NoBirds[m]<-n_distinct(loxSub1$individ_id)

# Replace NAs with 0
gridMonth[is.na(gridMonth)] <- 0

}

}

# Add number of locations
gridMonth$month<-Months[i]
gridMonth$totalBirds<-totBirds

# Save results
PopActivityFinal_All<-rbind(PopActivityFinal_All, gridMonth)

}

print("Saving results...")

PopActivityFinal_All$colony<-actMonthTrans$colony[1]
PopActivityFinal_All$colony_og<-actMonth$og_colony[1]

# Save result
output_file1 <- args[2]

# Save metadata as main output file
write.csv(PopActivityFinal_All, file=output_file1, row.names=FALSE)