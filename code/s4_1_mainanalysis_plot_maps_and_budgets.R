# This script does the following: makes maps for Figures 1, S1 & S2, and generates species and population-level activity & energy budgets #
### Input files is the id catalogue (table1_idcatalogue.csv) and start & end dates of the study period (which can be changed manually in the workflow file) ###
### The script maps where the study populations are (Figures 1, S1 & S2) ###
### It also outputs some figures showing variation in actvity budgets, sst & energy over time for my own visualisation ### 
### The data used to produce these are saved in the following output files  : 
#output_files1 = f"./results/tables/main/table2_budgets_species.csv" # Table with species-specific activity and energy budgets
#output_files2 = f"./results/tables/main/table3_budgets_population.csv" # Table with population-specific activity and energy budgets
#output_files3 = f"./results/tables/main/table4_budgets_individual.csv" # Table with population-specific activity and energy budgets

library(ggplot2)
library(dplyr)
library(fields)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(fasterize)
library(gdistance)
library(sf)
library(sp)
library(tidyr)
library(suncalc)
library(lubridate)
library(data.table)
library(terra)
library(ncdf4)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

# Set up minimum sample size & number of iterations

# minSampleSize<-1

minSampleSize<-5
reps<-50

### Step 1: Load all individual Ids from result files ###

# This is so we can loop through them later #

print("Step 1: Load id catalogue")
input_file <- args[1]
energyAll <- read.csv(input_file)

# Code below for testing purposes only
#ids<-unique(energyAll$individ_id)
#idsSub<-sample(ids, 100, replace=F)
#energyAll<-subset(energyAll, individ_id %in% c(idsSub))

### Step 2: Make a map of study site ###

print("Step 2: Figure 1 - make a map of study site")

#### First we start by making an sst map ####

# Make a list of where SST data is hidden #
print("Assembling SST layer...")
sstfiles<-list.files("./data/sst/", full.names=TRUE)
sstfiles_df<-as.data.frame(sstfiles)
colnames(sstfiles_df)<-c("FileName")
sstfiles_df$Month<-substr(sstfiles_df$FileName, 28, 29)
sstfiles_df$Month<-gsub("\\.", "", sstfiles_df$Month)

# Subset to winter months
sstfiles_df_winter<-subset(sstfiles_df, Month %in% c(9, 10, 11, 12, 1, 2, 3, 4))

# Open these one by one, stack & average
for (i in 1:nrow(sstfiles_df_winter)) {
  
  print(paste0("Averaging layer", i, "/", nrow(sstfiles_df_winter)))  
  
  # Open original raster
  if (i ==1) {  
    rastSub<-raster(sstfiles_df_winter$FileName[i]) 
    
  } else {
    rastSub2<-raster(sstfiles_df_winter$FileName[i])
    rastSub<-stack(rastSub2, rastSub)
    rastSub<-calc(rastSub, mean)
    
  }
  
}

# Convert from kelvin to degrees
values(rastSub)<-values(rastSub)-273.15

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Set CRS: Orthographic projection centered on the North Pole
north_pole_proj <- "+proj=ortho +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +no_defs"

# Set-up extent 
fake.rast<-crop(rastSub, extent(-200, 200, -60, 90), res=0.5)
SSTProject<-projectRaster(fake.rast, crs=projection_NA)
# Fill missing values using nearest neighbor interpolation (or other methods)
raster_data_filled <- focal(SSTProject, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)

# Turn into a data frame
SSTProject_df<-as.data.frame(raster_data_filled, xy=TRUE)

#### Then we make a polygon covering the location of the birds ####

# Make polygon around locations of birds #
print("Making polygon around locations...")

# open up locations of individual species #
idsPlot<-unique(energyAll$individ_id)

# Determine where locations are stored
allResults<-list.files("./tmp/", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# make a list to save locations in
allLocations<-list()

for (i in 1:length(idsPlot)) {
  
  print(paste0("Extracting lox from bird", i, "/", length(idsPlot)))
  
  # Subset to individual i
  birdID<-idsPlot[i]
  birdID<-gsub("-", "_", birdID)
  birdID<-gsub("ø", "o", birdID)
  birdSub<-fread(energyRes_day[grepl(birdID, energyRes_day)])
  
  # Susbet to specific months of interest
  birdSub_studyperiod<-subset(birdSub, Month %in% c(9, 10, 11, 12, 1, 2, 3, 4))
  
  # Determine number of years
  birdSub_studyperiod$year<-as.numeric(substr(birdSub_studyperiod$date, 1, 4))
  birdSub_studyperiod$day<-as.numeric(substr(birdSub_studyperiod$date, 9, 10))
  
  # Filter to one rep per year-month combination
  birdSub_studyperiod_lox<-birdSub_studyperiod %>%
    ungroup() %>%
    dplyr::group_by(year, Month, day) %>%
    dplyr::slice_sample(n=1) %>%
    dplyr::select(species, colony, individ_id, year, Month, day, mean.lon, mean.lat) %>%
    ungroup() %>%
    dplyr::group_by(species, colony, individ_id, Month) %>%
    dplyr::summarise(meanLon=mean(mean.lon), meanLat=mean(mean.lat))
  
  # Skip to next if no data
  
  if (nrow(birdSub_studyperiod_lox) < 1) {
    next
  }
  
  # Otherwise save in list
  allLocations<-rbind(allLocations, birdSub_studyperiod_lox)
  
}

# project coordinates & make a polygon around them #
coordinates(allLocations)<-~meanLon + meanLat
proj4string(allLocations)<-projection_84
allLocations_trans<-data.frame(spTransform(allLocations, north_pole_proj))
allLocations_polygon<-allLocations_trans %>%
  st_as_sf(coords = c("coords.x1", "coords.x2")) %>%
  group_by(species) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_convex_hull() %>%
  st_set_crs(north_pole_proj)
allLocations_polygon$species<-factor(allLocations_polygon$species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot"))

# Extract coordinates so I can define the max x and y limits of my plot
extentLox<-as.data.frame(st_coordinates(allLocations_polygon))

#### Now we make our map with labels ####

# Open country layers
print("Making map S1...")
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Add colony locations
colonies<-energyAll %>%
  dplyr::group_by(species) %>%
  dplyr::count(colony)

# open colony information & join so we have longitude & latitude & project into correct format
colony.summary<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")

speciesMatch<-data.frame(speciesLatin=c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle"), 
                         species=c("Brünnich's guillemot", "Black-legged kittiwake", "Common guillemot", "Atlantic puffin", "Northern fulmar", "Litte auk"))

colonyMatch<-colony.summary %>%
  dplyr::filter(species %in% c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle")) %>%
  dplyr::group_by(species, colony) %>%
  dplyr::slice(1) %>%
  dplyr::select(species, colony) %>%
  rename(speciesLatin=species) %>%
  dplyr::left_join(speciesMatch, by=c("speciesLatin"))

colonyNames<-colony.summary %>%
  dplyr::filter(species %in% c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle")) %>%
  ungroup() %>%
  dplyr::group_by(colony) %>%
  dplyr::slice(1) %>%
  dplyr::select(colony, col_lon, col_lat) %>%
  arrange(desc(col_lat)) %>%
  ungroup() %>%
  dplyr::mutate(country=c("RU", "NO", "NO", "NO", "GR", "RU", "NO", "GR", "GR", "RU", "NO", "CA", "GR", "CA", "NO", "NO", "GR",
                          "RU", "NO", "GR", "RU", "RU", "NO", "RU", "NO", "CA", "IC", "IC", "IC", "IC", "IC", "IC", 
                          "GR", "NO", "IC", "IC", "GR", "IC", "IC", "NO", "IC", "CA", "CA", "NO", "FA", "GR", "UK", "NO", "UK", "UK", "DK",
                          "UK", "UK", "UK", "IR", "IR", "CA", "IR", "IR", "IR", "UK", "CA", "CA")) %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(colonyNo=row_number()) %>%
  dplyr::mutate(colonyName=paste0(country, colonyNo)) # Here i make a special naming system with country first and colony second. It is ordered from North to South

colonies_lox<-colonyMatch %>%
  dplyr::left_join(colonyNames, by=c("colony")) %>%
  ungroup() %>%
  arrange(colonyName) %>%
  dplyr::group_by(colonyName) %>%
  dplyr::slice(1)

# Determine amount to jitter by to get all coordinates & labels on the same plot
amountJitter<-data.frame(colonyName=c(unique(colonies_lox$colonyName)), xAmount=c(1), yAmount=c(150000))
amountJitter$yAmount[4]<-200000 # CA4
amountJitter$xAmount[5]<-200000 # CA5
amountJitter$yAmount[9]<-(-200000) # DK1
amountJitter$xAmount[10]<-400000 # FA1
amountJitter$yAmount[10]<--200000 # FA1
amountJitter$yAmount[11]<-150000 # GR1
amountJitter$yAmount[12]<-(-180000) # GR2
amountJitter$xAmount[12]<-200000 # GR2
amountJitter$xAmount[13]<-(-200000) # GR3
amountJitter$yAmount[13]<-(-200000) # GR3
amountJitter$xAmount[14]<-(200000) # GR4
amountJitter$yAmount[14]<-(-180000) # GR4
amountJitter$yAmount[15]<-(400000) # GR5
amountJitter$xAmount[18]<-(200000) # GR8
amountJitter$yAmount[18]<-(-180000) # GR8
amountJitter$xAmount[21]<-(-290000) # IC10
amountJitter$yAmount[21]<-(-50000) # IC10
amountJitter$xAmount[22]<-(290000) # IC11
amountJitter$yAmount[22]<-(-200000) # IC11
amountJitter$xAmount[23]<-(-180000) # IC2
amountJitter$xAmount[24]<-(180000) # IC2
amountJitter$yAmount[25]<-(-800000) # IC4
amountJitter$yAmount[26]<-(-90000) # IC5
amountJitter$xAmount[26]<-(180000) # IC5
amountJitter$yAmount[27]<-(-30000) # IC6
amountJitter$xAmount[27]<-(-180000) # IC6
amountJitter$yAmount[28]<-600000 # IC7
amountJitter$xAmount[29]<-600000 # IC8
amountJitter$xAmount[30]<-600000 # IC9
amountJitter$yAmount[30]<--200000 # IC9
amountJitter$xAmount[31]<-(-1500000) # IR1
amountJitter$xAmount[32]<-(-600000) # IR2
amountJitter$xAmount[33]<-(-400000) # IR3
amountJitter$yAmount[33]<-(250000) # IR3
amountJitter$yAmount[34]<-(-250000) # IR4
amountJitter$xAmount[35]<-(-400000) # IR5
amountJitter$yAmount[36]<-(250000) # NO1
amountJitter$yAmount[37]<-(250000) # N10
amountJitter$xAmount[37]<-(-250000) # N10
amountJitter$xAmount[38]<-(450000) # NO11
amountJitter$xAmount[39]<-(-250000) # NO12
amountJitter$xAmount[40]<-(450000) # NO13
amountJitter$xAmount[42]<-(400000) # N15
amountJitter$xAmount[43]<-(-250000) # NO2
amountJitter$xAmount[44]<-(550000) # NO3
amountJitter$yAmount[44]<-(250000) # NO3
amountJitter$xAmount[45]<-(250000) # NO4
amountJitter$yAmount[48]<-(400000) # NO8
amountJitter$yAmount[49]<-(-200000) # NO9
amountJitter$xAmount[49]<-(-500000) # NO9
amountJitter$yAmount[50]<-(600000) # NO9
amountJitter$xAmount[52]<-(-250000) # RU3
amountJitter$yAmount[52]<-(-150000) # RU3
amountJitter$xAmount[53]<-(-250000) # RU3
amountJitter$yAmount[54]<-(650000) # RU5
amountJitter$yAmount[55]<-(-200000) # RU5
amountJitter$yAmount[56]<-(-600000) # RU6
amountJitter$xAmount[57]<-(300000) # UK1
amountJitter$yAmount[57]<-(-200000) # UK1
amountJitter$xAmount[58]<-(500000) # UK2
amountJitter$yAmount[58]<-(-600000) # UK2
amountJitter$yAmount[59]<-(300000) # UK3
amountJitter$xAmount[60]<-(-150000) # UK4
amountJitter$yAmount[61]<-(-1000000) # UK5
amountJitter$xAmount[61]<-(1000000) # UK5
amountJitter$yAmount[62]<-(-1200000) # UK6
amountJitter$yAmount[63]<-(-250000) # UK7
amountJitter$xAmount[63]<-(-550000) # UK7

coordinates(colonies_lox)<-~col_lon + col_lat
proj4string(colonies_lox)<-projection_84
colonies_lox_trans<-data.frame(spTransform(colonies_lox, projection_NA))
colonies_lox_trans<-colonies_lox_trans %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::left_join(amountJitter, by=c("colonyName")) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("GR3"), coords.x1-70000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("GR2"), coords.x1 + 30000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("IC5"), coords.x1 + 30000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("IC8"), coords.x1 - 100000, coords.x1)) 

# Make a map to save

FigureS1<-ggplot() +
  geom_tile(data=SSTProject_df, aes(x=x, y=y, fill=layer)) +
  scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(-2.5, 26)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_segment(data=colonies_lox_trans, aes(x=coords.x1, xend=coords.x1 + xAmount, y=coords.x2, yend=coords.x2 + yAmount, group=colonyName), color="darkgrey") + 
  geom_label(data=colonies_lox_trans, aes(x=coords.x1 + xAmount, y=coords.x2 + yAmount, label=colonyName), cex=2.5) +
  geom_point(data=colonies_lox_trans, aes(x=coords.x1, y=coords.x2),  cex=1.5, fill="yellow", shape=21) + 
  coord_sf(crs=projection_NA, xlim=c(-4095718.87 - 10000, 1963494.28 + 90000), ylim=c(-1025375.3 - 1600000, 2837995.2  + 220000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  theme(legend.position="bottom")

pdf("./results/figures/supplementary/FigureS1.pdf")
plot(FigureS1)
dev.off()

# Now we make map # 1 

print("Making map #1...")

# Set CRS: Orthographic projection centered on the North Pole
north_pole_proj <- "+proj=ortho +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +no_defs"

# Convert to a line object so I can change thickness of outline...
polygon_lines <- st_cast(allLocations_polygon, "MULTILINESTRING")

# Project SST into north pole projection
SSTProject2<-projectRaster(fake.rast, crs=north_pole_proj)
# Fill missing values using nearest neighbor interpolation (or other methods)
raster_data_filled <- focal(SSTProject2, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
# Turn into a data frame
SSTProject_df<-as.data.frame(raster_data_filled, xy=TRUE)

# project colony coordinates into north pole projection
colonies_lox_trans<-data.frame(spTransform(colonies_lox, north_pole_proj))
colonies_lox_trans<-colonies_lox_trans %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::left_join(amountJitter, by=c("colonyName")) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("GR3"), coords.x1-70000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("GR2"), coords.x1 + 30000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("IC5"), coords.x1 + 30000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("IC8"), coords.x1 - 100000, coords.x1)) 

Figure1<-ggplot() +
  geom_raster(data=subset(SSTProject_df, !is.na(layer)), aes(x=x, y=y, fill=layer), alpha=0.4, color=NA) +
  scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(-2.5, 33)) +
  geom_sf(data=polygon_lines, aes(color=species), fill=NA, linewidth=1) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=colonies_lox_trans, aes(x=coords.x1, y=coords.x2),  cex=2, fill="yellow", shape=21) + 
  scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs = north_pole_proj, datum=NA, clip="on") +
  xlab("") +
  ylab("") +
  labs(colour="Species", tag="A) Study system") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color="none")

pdf("./results/figures/main/Figure1.pdf", width = 8, height=12)
plot(Figure1)
dev.off()

##### Now we make individual maps with kernels ###### (collage of distirbution)

Figure1_v2<-ggplot() +
  geom_tile(data=SSTProject_df, aes(x=x, y=y, fill=layer)) +
  scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(-2.5, 26)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  #geom_segment(data=colonies_lox_trans, aes(x=coords.x1, xend=coords.x1 + xAmount, y=coords.x2, yend=coords.x2 + yAmount, group=colonyName), color="darkgrey") + 
  #geom_label(data=colonies_lox_trans, aes(x=coords.x1 + xAmount, y=coords.x2 + yAmount, label=colonyName), cex=2.5) +
  geom_point(data=colonies_lox_trans, aes(x=coords.x1, y=coords.x2),  cex=3, fill="yellow", shape=21) + 
  coord_sf(crs=north_pole_proj, xlim=c(min(colonies_lox_trans$coords.x1) - 10000, max(colonies_lox_trans$coords.x1) + 10000), 
  ylim=c(min(colonies_lox_trans$coords.x2) - 10000, max(colonies_lox_trans$coords.x2) + 10000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  theme(legend.position="bottom")
  
# Now for every species I need to make kernels & a seperate map #

speciesList<-unique(allLocations_trans$species)

colonyLoxes<-list()
speciesRasters<-list()

# Make a grid to plot accross
cellsize<-200000
coordinates(allLocations_trans)<-~coords.x1 + coords.x2
proj4string(allLocations_trans)<-north_pole_proj
null.grid <- expand.grid(x=seq(min(coordinates(allLocations_trans)[,1])-100000, max(coordinates(allLocations_trans)[,1])+100000, by = cellsize),
                         y=seq(min(coordinates(allLocations_trans)[,2])-100000, max(coordinates(allLocations_trans)[,2])+100000, by = cellsize))
coordinates(null.grid) <- ~x+y
gridded(null.grid) <- TRUE

# Make the kernels
h<-200000

for (i in 1:length(speciesList)) {

# Subset to species i
allLocations_trans<-as.data.frame(allLocations_trans)
speciesSub<-subset(allLocations_trans, species==speciesList[i])
speciesSub<-speciesSub %>%
dplyr::group_by(individ_id) %>%
dplyr::mutate(locs=n_distinct(Month)) %>%
dplyr::filter(locs>5) %>%
dplyr::select(individ_id, coords.x1, coords.x2)
coordinates(speciesSub)<-~coords.x1 + coords.x2
proj4string(speciesSub)<-north_pole_proj


# Create weights to weigh the kernels by
weights<-as.data.frame(allLocations_trans) %>%
dplyr::filter(species==speciesList[i])%>%
dplyr::group_by(colony) %>%
dplyr::summarise(birds=n_distinct(individ_id))

# Colonies
colonies<-unique(weights$colony)

# list for saving #
colonyUDs<-list()

for (k in 1:length(colonies)) {

print(paste0("Species", i, " Colony", k, "/", length(colonies)))

birdsSub<-subset(allLocations_trans, colony==colonies[k] & species==speciesList[i])
birdsSubIds<-unique(birdsSub$individ_id)
colonySub<-subset(allLocations_trans, individ_id %in% birdsSubIds)
colonySub<-colonySub %>%
dplyr::group_by(individ_id) %>%
dplyr::mutate(locs=n_distinct(Month)) %>%
dplyr::filter(locs==8) %>%
dplyr::select(individ_id, coords.x1, coords.x2)
coordinates(colonySub)<-~coords.x1 + coords.x2
proj4string(colonySub)<-north_pole_proj
kernels<-adehabitatHR::kernelUD(colonySub, grid = null.grid, h = h)

# Subset weights
weightsSub<-subset(weights, colony==colonies[k])

# Merge rasters 
uds<-list()

for (j in 1:length(kernels)) {

print(j)

# Extract individual rasters & merge (here i give colonies equal weighting)
rast1<-raster(kernels[[j]])
values(rast1)<-values(rast1)*(1/weightsSub$birds)

# Save result
uds<-append(uds, rast1)

}

# Merge the rasters
alluds<-stack(uds)
alludsMean<-mean(alluds)

# Save results
colonyUDs<-append(colonyUDs, alludsMean)

}

# Merge the rasters
colStack<-stack(colonyUDs)
colMean<-mean(colStack)


# Transform into a data frame
colMeanDf<-as.data.frame(colMean, xy=TRUE)

# Order colMeanDf
colMeanDf_order<-colMeanDf %>%
arrange(desc(layer))

twentyfive<-quantile(colMeanDf$layer, 0.75)
fifty<-quantile(colMeanDf$layer, 0.5)
ninetyfive<-quantile(colMeanDf$layer, 0.05)
colMeanDf$quantile<-ifelse(colMeanDf$layer>ninetyfive, 95, "Other")
colMeanDf$quantile<-ifelse(colMeanDf$layer>fifty, 50, colMeanDf$quantile)
#colMeanDf$quantile<-ifelse(colMeanDf$layer>twentyfive, 25, colMeanDf$quantile)

# Save result
colMeanDf$species<-speciesList[i]
speciesRasters<-rbind(speciesRasters, colMeanDf)

# Subset colonies & Save
colonyLoxSub<-subset(colonies_lox_trans, colony %in% c(colonies))
colonyLoxSub<-colonyLoxSub %>%
dplyr::select(colony, coords.x1, coords.x2) %>%
dplyr::mutate(species=speciesList[i])
colonyLoxes<-rbind(colonyLoxes, colonyLoxSub)

}

speciesRasters<-speciesRasters %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))
												 
colonyLoxes<-colonyLoxes %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
      "Little auk", "Common guillemot", "Brünnich's guillemot")))

minmax_scale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    return(rep(0, length(x)))  # handle case where all values are identical
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

speciesScale<-speciesRasters %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::mutate(layerScale=minmax_scale(layer))

FigureS2<-ggplot() +
  #geom_tile(data=filter(speciesRasters,quantile==50), aes(x=x, y=y, fill=species)) +
  geom_contour_filled(data=filter(speciesScale, layerScale>0), aes(x=x, y=y, z=layerScale), limits=c(0.01, 1)) +
  #geom_density_2d(data=filter(speciesRasters, layer>0), aes(x=x, y=y, fill=layer)) +
  #scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  #geom_sf(data=world, color = "#E5E5E5", fill = "white") +
  geom_sf(data=coast) +
  labs(fill="Kernel UD") +
  #geom_segment(data=colonies_lox_trans, aes(x=coords.x1, xend=coords.x1 + xAmount, y=coords.x2, yend=coords.x2 + yAmount, group=colonyName), color="darkgrey") + 
  #geom_label(data=colonies_lox_trans, aes(x=coords.x1 + xAmount, y=coords.x2 + yAmount, label=colonyName), cex=2.5) +
  geom_point(data=colonyLoxes, aes(x=coords.x1, y=coords.x2),  cex=2, fill="yellow", shape=21) + 
  #scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=north_pole_proj, xlim=c(min(allLocations_trans$coords.x1) - 10000, max(allLocations_trans$coords.x1) + 10000), 
  ylim=c(min(allLocations_trans$coords.x2) - 10000, max(allLocations_trans$coords.x2) + 10000)) +
  scale_y_continuous(breaks=c(0, 30)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  theme(legend.position="bottom") +
  facet_wrap(~species) +
  theme(plot.background = element_rect(fill="white"))
  
  pdf("./results/figures/supplementary/FigureS2.pdf")
  plot(FigureS2)
  dev.off()


# Now we add any extra one with currents & sea names for better describing what's going on #

print("Making map S2...")

currents<-st_read("./data/currents/Major_Ocean_Currents.shp")
currents_trans<-st_transform(currents, crs = projection_NA)

colonies_lox_trans<-data.frame(spTransform(colonies_lox, projection_NA))
colonies_lox_trans<-colonies_lox_trans %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::left_join(amountJitter, by=c("colonyName")) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("GR3"), coords.x1-70000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("GR2"), coords.x1 + 30000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("IC5"), coords.x1 + 30000, coords.x1)) %>%
  dplyr::mutate(coords.x1=ifelse(colonyName %in% c("IC8"), coords.x1 - 100000, coords.x1)) 

FigureS25<-ggplot() +
  geom_sf(data=currents_trans, aes(fill=TEMP), alpha=0.2) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_text(aes(label="Arctic Ocean", x=1963494.28 - 1700000, y=2837995.2 + 90000), color="darkblue", size=3 ,fontface = "italic")+
  geom_text(aes(label="Barents Sea", x=1963494.28 -700000, y=2837995.2 - 900000), color="darkblue", size=2 ,fontface = "italic", angle=45)+
  geom_text(aes(label="Greenland Sea", x=1963494.28 -1800000, y=2837995.2 - 800000), color="darkblue", size=2 ,fontface = "italic", angle=80)+
  geom_text(aes(label="Norwegian Sea", x=1963494.28 -1800000, y=2837995.2 - 2000000), color="darkblue", size=2 ,fontface = "italic", angle=80)+
  geom_text(aes(label="North Sea", x=1963494.28 -1200000, y=2837995.2 - 3300000), color="darkblue", size=2 ,fontface = "italic")+
  geom_text(aes(label="Labrador Sea", x=1963494.28 -4300000, y=2837995.2 - 3000000), color="darkblue", size=2 ,fontface = "italic")+
  geom_text(aes(label="Baffin Bay", x=1963494.28 -3800000, y=2837995.2 - 1300000), color="darkblue", size=2 ,fontface = "italic", angle=65)+
  scale_fill_manual(values=c("darkblue", "red")) +
  geom_point(data=colonies_lox_trans, aes(x=coords.x1, y=coords.x2),  cex=3, fill="yellow", shape=21) + 
  coord_sf(crs=projection_NA, xlim=c(-4095718.87 - 10000, 1963494.28 + 90000), ylim=c(-1025375.3 - 1600000, 2837995.2  + 220000)) +
  geom_sf_text(data=currents_trans, aes(label=NAME), size=2) +
  xlab("") +
  ylab("") +
  labs(colour="", fill="Type of current") +
  guides(
    color = guide_legend(position = "bottom"),  # Move color legend to bottom
    fill = guide_legend(position = "bottom")    # Keep shape legend on the right
  ) 

# This map is not appearing in the MS currently #

#pdf("./results/figures/supplementary/FigureS25.pdf")
#plot(FigureS25)
#dev.off()

# make a blank one so I can draw myself

FigureS1_part3<-ggplot() +
  #geom_sf(data=currents_trans, aes(fill=TEMP), alpha=0.2) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_text(aes(label="Arctic Ocean", x=1963494.28 - 1700000, y=2837995.2 + 90000), color="darkblue", size=3 ,fontface = "italic")+
  geom_text(aes(label="Barents Sea", x=1963494.28 -700000, y=2837995.2 - 900000), color="darkblue", size=2 ,fontface = "italic", angle=45)+
  geom_text(aes(label="Greenland Sea", x=1963494.28 -1800000, y=2837995.2 - 800000), color="darkblue", size=2 ,fontface = "italic", angle=80)+
  geom_text(aes(label="Norwegian Sea", x=1963494.28 -1800000, y=2837995.2 - 2000000), color="darkblue", size=2 ,fontface = "italic", angle=80)+
  geom_text(aes(label="North Sea", x=1963494.28 -1200000, y=2837995.2 - 3300000), color="darkblue", size=2 ,fontface = "italic")+
  geom_text(aes(label="Labrador Sea", x=1963494.28 -4300000, y=2837995.2 - 3000000), color="darkblue", size=2 ,fontface = "italic")+
  geom_text(aes(label="Baffin Bay", x=1963494.28 -3800000, y=2837995.2 - 1300000), color="darkblue", size=2 ,fontface = "italic", angle=65)+
  scale_fill_manual(values=c("darkblue", "red")) +
  geom_point(data=colonies_lox_trans, aes(x=coords.x1, y=coords.x2),  cex=2, fill="yellow", shape=21) + 
  coord_sf(crs=projection_NA, xlim=c(-4095718.87 - 10000, 1963494.28 + 90000), ylim=c(-1025375.3 - 1600000, 2837995.2  + 220000)) +
  #geom_sf_text(data=currents_trans, aes(label=NAME), size=2) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  guides(
    color = guide_legend(position = "bottom"),  # Move color legend to bottom
    fill = guide_legend(position = "bottom")    # Keep shape legend on the right
  ) 

#pdf("./results/figures/supplementary/FigureS25_blank.pdf")
#plot(FigureS1_part3)
#dev.off()

### Step 4: Species-specific average activity & possible energy ####

# Here we make general maps showing temporal variation in activity budgets, SST & energy (but mainly for checking everything looks good! #

print("Step 4: calculating species-specific budgets...")

# Determine where daily files are
allResults<-list.files("./tmp/", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# Determine study period
startDate<-args[2] # Read-in start of study period
endDate<-args[3] # Read-in end date of study period

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
 
day1<-as.Date(min(dates_weekly$dateKeep)-1)
day2<-as.Date(max(dates_weekly$dateKeep)+1)
date1<-data.frame(dateKeep=day1, doy=0, month=as.numeric(substr(day1, 6, 7)), day=as.numeric(substr(day1, 9, 10)), weekNo=0)
date2<-data.frame(dateKeep=day2, doy=0, month=as.numeric(substr(day2, 6, 7)), day=as.numeric(substr(day2, 9, 10)), weekNo=0)
dates_weekly2<-rbind(date1, dates_weekly, date2)

# Within sessions, subset to 'trackyears' if possible
day1_month<-as.numeric(substr(day1, 6, 7))
day2_month<-as.numeric(substr(day2, 6, 7))

# Determine number of individuals
ids<-unique(energyAll$individ_id)

# Define number of species
speciesList<-unique(energyAll$species)

# Define coefficients to divide energy by
speciesWeightDivider<-data.frame(species=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                           "Little auk", "Common guillemot", "Brünnich's guillemot"), divider=c(0.717, 0.765, 0.689, 0.689, 0.689, 0.689),
                                 newWeight=c(365, 651, 370, 149, 989, 989))

# 365 & 561 are from Gabrielsen's paper
# 989 is average weight of BrGu from our dataset
# 370 is from annette's apper
# 151 is from st-marie paper

# Make list to save results in
average_activity<-list()

# Loop through species

for (j in 1:length(speciesList)) {
  
  print(paste0("Step 4: Species ", j, "/", length(speciesList)))
  
  speciesSub<-speciesList[j]
  
  # Assemble birds 
  
  # Find all relevant ids
  ids<-subset(energyAll, species==speciesSub)
  ids<-unique(ids$individ_id)
  
  # Create list to save results
  speciesDaily<-list()
  
  # Update message
  print("Assembling birds...")
  
  for (k in 1:length(ids)) {
    
    # Print update message
    print(paste0("Step4: Assembling... Species ", j, " Bird ", k, "/", length(ids)))  
    
    # Open id k
    birdSub<-ids[k]
	birdSub<-gsub("-", "_", birdSub)
	birdSub<-gsub("ø", "o", birdSub)
    
    # Find csv file
    print("Opening file...")
    birdSub<-fread(energyRes_day[grep(birdSub, energyRes_day)])
    
	# Determine month and year
	birdSub$day<-as.numeric(substr(birdSub$date, 9, 10))
	birdSub$month<-as.numeric(substr(birdSub$date, 6, 7))
	birdSub$year<-as.numeric(substr(birdSub$date, 1, 4))
	birdSub$track_year<-ifelse(birdSub$month<day1_month, paste0(birdSub$year-1, "_", substr(birdSub$year, 3, 4)), paste0(birdSub$year, "_", substr(birdSub$year + 1, 3, 4)))
    birdSub$session_year<-paste0(birdSub$session_id, "_", birdSub$track_year)

# Determine session_year combos with enough data	
dates_sessions<-birdSub%>%
dplyr::mutate(date=substr(date, 1, 10)) %>%
dplyr::group_by(rep, session_id, track_year) %>%
dplyr::count(date) %>%
dplyr::mutate(month=as.numeric(substr(date, 6, 7))) %>%
dplyr::mutate(day=as.numeric(substr(date, 9, 10))) %>%
dplyr::inner_join(dates_weekly2, by=c("month", "day")) %>%
ungroup() %>%
dplyr::group_by(rep, session_id, track_year) %>%
dplyr::summarise(daysTot=n_distinct(date)) %>%
dplyr::filter(daysTot==nrow(dates_weekly2)) %>%
dplyr::mutate(session_year=paste0(session_id, "_", track_year))

if (nrow(dates_sessions)<1) {
next}
    
# Choose 50 reps & session_years at random 
repsSelect<-sample(1:max(dates_sessions$rep), reps, replace=FALSE)
sessionSelect<-sample(unique(dates_sessions$session_year), reps, replace=TRUE)
randomSelect<-data.frame(rep=repsSelect, session_year=sessionSelect)
    
birdSub_random<-birdSub %>%
   dplyr::inner_join(randomSelect, by=c("rep", "session_year")) %>%
   ungroup() %>%
   dplyr::group_by(individ_id) %>%
   dplyr::mutate(rep=as.numeric(factor(rep)))
    
# Check if a column with a specific name exists
    
if ("tActive" %in% colnames(birdSub) == FALSE) {
      
birdSub_random$tActive<-0
    }
    
weighted_activity<-birdSub_random %>%
    dplyr::left_join(dates_weekly2, by=c("day", "month")) %>%
	dplyr::filter(weekNo>0) %>%
    dplyr::left_join(speciesWeightDivider, by=c("species")) %>%
    dplyr::mutate(DEE=DEEkJ/weight^divider) %>%
    dplyr::mutate(DEEkJ_2=DEE*newWeight^divider) %>%
    dplyr::group_by(rep, species, colony, weekNo, individ_id) %>%
    dplyr::summarise(Flight_tot=sum(tFlight)/168, Rest_tot=sum(tRestWater)/168, Active_tot=sum(tActive)/168, Land_tot=sum(tLand)/168,
                       Forage_tot=sum(tForage)/168, sstMean=mean(sst_random), weight=mean(weight), energyTot=sum(DEEkJ_2)/7, energyTot_weekly=sum(DEEkJ)) %>%
    ungroup() %>%
    dplyr::group_by(rep, species, colony, individ_id) %>%
    dplyr::mutate(meanDEE=mean(energyTot_weekly), deviance_weekly=abs(energyTot_weekly-meanDEE)/meanDEE, totdeviance_weekly=sum(deviance_weekly)) %>%
    ungroup() 
    
    # Save results
    speciesDaily<-rbind(speciesDaily, weighted_activity)
    
    
  }
  
  # create weights related to colony size to conduct a weighted average of activity budgets
  colonySampleSize<-speciesDaily %>%
    dplyr::group_by(rep, colony) %>%
    dplyr::summarise(birds=n_distinct(individ_id)) %>%
    ungroup() %>%
    dplyr::mutate(weights=1/birds) 
  
  # Check if a column with a specific name exists
  
  if ("tActive" %in% colnames(speciesDaily) == FALSE) {
    
    speciesDaily$tActive<-0
  }
  
  weighted_activity2<-speciesDaily %>%
    ungroup() %>%
    dplyr::group_by(rep, species, colony, weekNo) %>%
    dplyr::summarise(weightedFlight=mean(Flight_tot), weightedRest=mean(Rest_tot),  
                     weightedActive=mean(Active_tot), weightedLand=mean(Land_tot), 
                     weightedForage=mean(Forage_tot), weightedsst=mean(sstMean), maxsst=max(sstMean), minsst=min(sstMean), energy=mean(energyTot), deviance_weekly=mean(totdeviance_weekly), birds_sample=n_distinct(individ_id)) %>%
    ungroup() %>%
    dplyr::left_join(colonySampleSize, by=c("rep", "colony")) %>%
    dplyr::filter(birds>=minSampleSize) %>%
    dplyr::group_by(rep, species, weekNo) %>%
    dplyr::summarise(sumWeights=sum(weights), meanFlight=sum(weightedFlight*weights)/sumWeights,
                     meanRest=sum(weightedRest*weights)/sumWeights, meanActive=sum(weightedActive*weights)/sumWeights, meanLand=sum(weightedLand*weights)/sumWeights, 
                     meanForage=sum(weightedForage*weights)/sumWeights, meansst=sum(weightedsst*weights)/sumWeights,
                     minsst=min(minsst), maxsst=max(maxsst), meanDEE=sum(energy*weights)/sumWeights, meanDeviance=sum(deviance_weekly*weights)/sumWeights, birdsTot=sum(birds_sample))
  
  average_activity<-rbind(average_activity, weighted_activity2)
  
} 

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

### Step 5: the same as above but for colonies ####

print("Step 5: estimating population-specific parameters...")

# Run from 1st of September through to end of April # 

# Determine number of individuals
ids<-unique(energyAll$individ_id)

# Define number of species
speciesList<-unique(energyAll$species)

# Define coefficients to divide energy by
speciesWeightDivider<-data.frame(species=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                           "Little auk", "Common guillemot", "Brünnich's guillemot"), divider=c(0.717, 0.765, 0.689, 0.689, 0.689, 0.689),
                                 newWeight=c(365, 651, 370, 151, 989, 989))

# Make list to save results in
average_activity_colony<-list()
average_activity_colony_individual<-list()

# Loop through species

for (j in 1:length(speciesList)) {
  
  print(paste0("Step 5: Species ", j, "/", length(speciesList)))
  
  speciesSub<-speciesList[j]
  
  # Assemble birds 
  
  # Find all relevant ids
  ids<-subset(energyAll, species==speciesSub)
  ids<-unique(ids$individ_id)
  
  # Create list to save results
  speciesDaily<-list()
  
  # Update message
  print("Assembling birds...")
  
  for (k in 1:length(ids)) {
    
# Print update message
print(paste0("Step5: Assembling... Species ", j, " Bird ", k, "/", length(ids)))  
    
# Open id k
    birdSub<-ids[k]
	birdSub<-gsub("-", "_", birdSub)
	birdSub<-gsub("ø", "o", birdSub)
    
    # Find csv file
    print("Opening file...")
    birdSub<-fread(energyRes_day[grep(birdSub, energyRes_day)])
    
# Determine month and year
birdSub$day<-as.numeric(substr(birdSub$date, 9, 10))
birdSub$month<-as.numeric(substr(birdSub$date, 6, 7))
birdSub$year<-as.numeric(substr(birdSub$date, 1, 4))
birdSub$track_year<-ifelse(birdSub$month<day1_month, paste0(birdSub$year-1, "_", substr(birdSub$year, 3, 4)), paste0(birdSub$year, "_", substr(birdSub$year + 1, 3, 4)))
birdSub$session_year<-paste0(birdSub$session_id, "_", birdSub$track_year)

# Determine session_year combos with enough data	
dates_sessions<-birdSub%>%
dplyr::mutate(date=substr(date, 1, 10)) %>%
dplyr::group_by(rep, session_id, track_year) %>%
dplyr::count(date) %>%
dplyr::mutate(month=as.numeric(substr(date, 6, 7))) %>%
dplyr::mutate(day=as.numeric(substr(date, 9, 10))) %>%
dplyr::inner_join(dates_weekly2, by=c("month", "day")) %>%
ungroup() %>%
dplyr::group_by(rep, session_id, track_year) %>%
dplyr::summarise(daysTot=n_distinct(date)) %>%
dplyr::filter(daysTot==nrow(dates_weekly2)) %>%
dplyr::mutate(session_year=paste0(session_id, "_", track_year))

if (nrow(dates_sessions)<1) {
next}
    
# Choose 50 reps & session_years at random 
repsSelect<-sample(1:max(dates_sessions$rep), reps, replace=FALSE)
sessionSelect<-sample(unique(dates_sessions$session_year), reps, replace=TRUE)
randomSelect<-data.frame(rep=repsSelect, session_year=sessionSelect)
    
birdSub_random<-birdSub %>%
   dplyr::inner_join(randomSelect, by=c("rep", "session_year")) %>%
   ungroup() %>%
   dplyr::group_by(individ_id) %>%
   dplyr::mutate(rep=as.numeric(factor(rep)))
    
# Check if a column with a specific name exists
    
if ("tActive" %in% colnames(birdSub) == FALSE) {
      
birdSub_random$tActive<-0
    }
    
    weighted_activity<-birdSub_random %>%
      ungroup() %>%
      dplyr::filter(!is.na(DEEkJ)) %>%
      dplyr::left_join(dates_weekly2, by=c("day", "month")) %>%
	  dplyr::filter(weekNo>0) %>%
      dplyr::left_join(speciesWeightDivider, by=c("species")) %>%
      dplyr::mutate(DEE=DEEkJ/weight^divider) %>%
      dplyr::mutate(DEEkJ_2=DEE*newWeight^divider) %>%
      dplyr::group_by(rep, species, colony, weekNo, individ_id) %>%
      dplyr::summarise(Flight_tot=sum(tFlight)/168, Rest_tot=sum(tRestWater)/168, Active_tot=sum(tActive)/168, Land_tot=sum(tLand)/168,
                       Forage_tot=sum(tForage)/168, sstMean=mean(sst_random), weight=mean(weight), energyTot=sum(DEEkJ_2)/7) %>%
      ungroup()
    
    # Save results
    speciesDaily<-rbind(speciesDaily, weighted_activity)
    
    
  }
  
  # create weights related to colony size to conduct a weighted average of activity budgets
  colonySampleSize<-speciesDaily %>%
    dplyr::group_by(rep, colony) %>%
    dplyr::summarise(birds=n_distinct(individ_id)) %>%
    ungroup() %>%
    dplyr::mutate(weights=1/birds) 
  
  # Check if a column with a specific name exists
  
  if ("tActive" %in% colnames(speciesDaily) == FALSE) {
    
    speciesDaily$tActive<-0
  }
  
  weighted_activity2<-speciesDaily %>%
    dplyr::left_join(colonySampleSize, by=c("rep", "colony")) %>%
    dplyr::filter(birds>=minSampleSize) %>%
    dplyr::group_by(rep, species, colony, weekNo) %>%
    dplyr::summarise(meanFlight=mean(Flight_tot), meanRest=mean(Rest_tot),  
                     meanActive=mean(Active_tot), meanLand=mean(Land_tot), 
                     meanForage=mean(Forage_tot), meansst=mean(sstMean), maxsst=max(sstMean), minsst=min(sstMean), meanDEE=mean(energyTot), birds_sample=n_distinct(individ_id)) %>%
    ungroup() 
  
  average_activity_colony<-rbind(average_activity_colony, weighted_activity2)
  
  weighted_activity3<-speciesDaily %>%
    dplyr::left_join(colonySampleSize, by=c("rep", "colony")) %>%
    dplyr::filter(birds>=minSampleSize) %>%
    dplyr::group_by(species, colony, individ_id, weekNo) %>%
    dplyr::summarise(meanFlight=mean(Flight_tot), meanRest=mean(Rest_tot),  
                     meanActive=mean(Active_tot), meanLand=mean(Land_tot), 
                     meanForage=mean(Forage_tot), meansst=mean(sstMean), maxsst=max(sstMean), minsst=min(sstMean), meanDEE=mean(energyTot), birds_sample=n_distinct(individ_id)) %>%
    ungroup() 
  
  average_activity_colony_individual<-rbind(average_activity_colony_individual, weighted_activity3)
  
} 
# Create an intermediate data frame before plotting

average_activity_all_colony<-average_activity_colony %>%
  ungroup() %>%
  dplyr::group_by(species, colony, weekNo) %>%
  dplyr::summarise(reps=n_distinct(rep), Flight=mean(meanFlight), sdFlight=sd(meanFlight), seFlight=sdFlight/sqrt(reps),
                   Forage=mean(meanForage), sdForage=sd(meanForage), seForage=sdForage/sqrt(reps),
                   Land=mean(meanLand), sdLand=sd(meanLand), seLand=sdLand/sqrt(reps),
                   Rest=mean(meanRest), sdRest=sd(meanRest), seRest=sdRest/sqrt(reps),
                   Active=mean(meanActive), sdActive=sd(meanActive), seActive=sdActive/sqrt(reps),
                   SST=mean(meansst), sdSST=sd(meansst), seSST=sdSST/sqrt(reps), SST_max=max(maxsst), SST_min=min(minsst),
                   energy=mean(meanDEE), sdenergy=sd(meanDEE), seenergy=sdenergy/sqrt(reps)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 

# Here is the species average that we will also add #

average_activity_all_species<-average_activity_all %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 

average_activity_individual_all<-average_activity_colony_individual %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 

# Figure S9

print("Making supplementary plots - activity")

startMonth<-dates_weekly %>%
  dplyr::filter(day==1)

FigureS9<-average_activity_all_colony %>%
  ggplot() +
  geom_pointrange(aes(x=weekNo, y=Flight, colour="Flight",  ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight), cex=0.1, alpha=0.05) +
  geom_line(aes(colour="Flight", x=weekNo, y=Flight, group=colony), alpha=0.05) +
  geom_ribbon(aes(x=weekNo, y=Flight, ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight, fill="Flight", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2) +
  labs(color="", fill="") +
  geom_pointrange(data=filter(average_activity_all_colony, Active>0) ,aes(x=weekNo, y=Active, colour="Active", ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive), cex=0.1, alpha=0.05) +
  geom_line(data=filter(average_activity_all_colony, Active>0) ,aes(colour="Active", x=weekNo, y=Active, group=colony), alpha=0.05) +
  geom_ribbon(data=filter(average_activity_all_colony ,Active>0) ,aes(x=weekNo, y=Active, ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive, fill="Active", group=colony), alpha=0.05) +
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
  geom_pointrange(data=filter(average_activity_all_species, Active>0), aes(x=weekNo, y=Active, colour="Active", ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive), cex=0.1) +
  geom_line(data=filter(average_activity_all_species, Active>0),aes(colour="Active", x=weekNo, y=Active)) +
  geom_ribbon(data=filter(average_activity_all_species, Active>0),aes(x=weekNo, y=Active, ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive, fill="Active"), alpha=0.2) +
  geom_pointrange(data=filter(average_activity_all_species, Rest>0), aes(x=weekNo, y=Rest, colour="Water", ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest), cex=0.1) +
  geom_line(data=filter(average_activity_all_species, Rest>0),aes(colour="Water", x=weekNo, y=Rest)) +
  geom_ribbon(data=filter(average_activity_all_species, Rest>0),aes(x=weekNo, y=Rest, ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest, fill="Water"), alpha=0.2) +
  geom_pointrange(data=filter(average_activity_all_species, Flight>0), aes(x=weekNo, y=Flight, colour="Flight", ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight), cex=0.1) +
  geom_line(data=filter(average_activity_all_species, Flight>0),aes(colour="Flight", x=weekNo, y=Flight)) +
  geom_ribbon(data=filter(average_activity_all_species, Flight>0),aes(x=weekNo, y=Flight, ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight, fill="Flight"), alpha=0.2) +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Proportion of day spent in behaviour") +
  scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))+
  scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") 

pdf("./results/figures/supplementary/activity.pdf", width=10, height=7)
grid.arrange(FigureS9)
dev.off()

# Supplementary plot showing time in flight

print("Making supplementary plots - flight")

FigureS10<-average_activity_all_colony %>%
  ggplot() +
  #geom_image(data=birdsimage, aes(x=weekNo, y=maxDEE, image=image), size=0.26, alpha=1) +
  geom_pointrange(aes(x=weekNo, y=Flight, colour="Flight",  ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight), cex=0.1, alpha=0.05) +
  geom_line(aes(colour="Flight", x=weekNo, y=Flight, group=colony), alpha=0.05) +
  geom_ribbon(aes(x=weekNo, y=Flight, ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight, fill="Flight", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free") +
  labs(color="", fill="") +
  # geom_pointrange(data=filter(average_activity_all_colony, Active>0) ,aes(x=weekNo, y=Active, colour="Active", ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive), cex=0.1, alpha=0.05) +
  #geom_line(data=filter(average_activity_all_colony, Active>0) ,aes(colour="Active", x=weekNo, y=Active, group=colony), alpha=0.05) +
  #geom_ribbon(data=filter(average_activity_all_colony ,Active>0) ,aes(x=weekNo, y=Active, ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive, fill="Active", group=colony), alpha=0.05) +
  #geom_pointrange(data=filter(average_activity_all_colony, Forage>0) ,aes(x=weekNo, y=Forage, colour="Forage", ymin=Forage-1.96*seForage, ymax=Forage + 1.96*seForage, group=colony), cex=0.1, alpha=0.05) +
  #geom_line(data=filter(average_activity_all_colony, Forage>0) ,aes(colour="Forage", x=weekNo, y=Forage, group=colony), alpha=0.05) +
  #geom_ribbon(data=filter(average_activity_all_colony, Forage>0) ,aes(x=weekNo, y=Forage, ymin=Forage-1.96*seForage, ymax=Forage + 1.96*seForage, fill="Forage", group=colony), alpha=0.05) +
  #geom_pointrange(aes(x=weekNo, y=Rest, colour="Water", ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest, group=colony), cex=0.1, alpha=0.05) +
  #geom_line(aes(colour="Water", x=weekNo, y=Rest, group=colony), alpha=0.05) +
  #geom_ribbon(aes(x=weekNo, y=Rest, ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest, fill="Water", group=colony), alpha=0.05) +
  #geom_pointrange(aes(x=weekNo, y=Land, colour="Land", ymin=Land-1.96*seLand, ymax=Land + 1.96*seLand, group=colony), alpha=0.05, cex=0.1) +
  #geom_line(aes(colour="Land", x=weekNo, y=Land, group=colony), alpha=0.05) +
  #geom_ribbon(aes(x=weekNo, y=Land, ymin=Land-1.96*seLand, ymax=Land + 1.96*seLand, fill="Land", group=colony), alpha=0.05) +
  theme_bw() +
  ### Now for species-specific trends
  #geom_pointrange(data=average_activity_all_species, aes(x=weekNo, y=Land, colour="Land", ymin=Land-1.96*seLand, ymax=Land + 1.96*seLand), cex=0.1) +
  #geom_line(data=average_activity_all_species,aes(colour="Land", x=weekNo, y=Land)) +
  #geom_ribbon(data=average_activity_all_species,aes(x=weekNo, y=Land, ymin=Land-1.96*seLand, ymax=Land + 1.96*seLand, fill="Land"), alpha=0.2) +
  #geom_pointrange(data=filter(average_activity_all_species, Forage>0), aes(x=weekNo, y=Forage, colour="Forage", ymin=Forage-1.96*seForage, ymax=Forage + 1.96*seForage), cex=0.1) +
  #geom_line(data=filter(average_activity_all_species, Forage>0),aes(colour="Forage", x=weekNo, y=Forage)) +
  #geom_ribbon(data=filter(average_activity_all_species, Forage>0),aes(x=weekNo, y=Forage, ymin=Forage-1.96*seForage, ymax=Forage + 1.96*seForage, fill="Forage"), alpha=0.2) +
  #geom_pointrange(data=filter(average_activity_all_species, Active>0), aes(x=weekNo, y=Active, colour="Active", ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive), cex=0.1) +
  #geom_line(data=filter(average_activity_all_species, Active>0),aes(colour="Active", x=weekNo, y=Active)) +
  #geom_ribbon(data=filter(average_activity_all_species, Active>0),aes(x=weekNo, y=Active, ymin=Active-1.96*seActive, ymax=Active + 1.96*seActive, fill="Active"), alpha=0.2) +
  #geom_pointrange(data=filter(average_activity_all_species, Rest>0), aes(x=weekNo, y=Rest, colour="Water", ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest), cex=0.1) +
  #geom_line(data=filter(average_activity_all_species, Rest>0),aes(colour="Water", x=weekNo, y=Rest)) +
  #geom_ribbon(data=filter(average_activity_all_species, Rest>0),aes(x=weekNo, y=Rest, ymin=Rest-1.96*seRest, ymax=Rest + 1.96*seRest, fill="Water"), alpha=0.2) +
  geom_pointrange(data=filter(average_activity_all_species, Flight>0), aes(x=weekNo, y=Flight, colour="Flight", ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight), cex=0.1) +
  geom_line(data=filter(average_activity_all_species, Flight>0),aes(colour="Flight", x=weekNo, y=Flight)) +
  geom_ribbon(data=filter(average_activity_all_species, Flight>0),aes(x=weekNo, y=Flight, ymin=Flight-1.96*seFlight, ymax=Flight + 1.96*seFlight, fill="Flight"), alpha=0.2) +
  
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Proportion of day spent in behaviour") +
  scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))+
  scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") 

pdf("./results/figures/supplementary/flight.pdf", width=10, height=7)
plot(FigureS10)
dev.off()

# Making supplementary plots - SST

print("Making plot S11")

average_activity_all_colony2<-average_activity_colony %>%
  ungroup() %>%
  dplyr::group_by(species, colony, rep) %>%
  #dplyr::mutate(SST_scale=scale_to_range(meansst, old_min=min(meansst), old_max=max(meansst), new_min=min(meanDEE), new_max=max(meanDEE))) %>%
  dplyr::group_by(species, colony, weekNo) %>%
  dplyr::summarise(reps=n_distinct(rep), Flight=mean(meanFlight), sdFlight=sd(meanFlight), seFlight=sdFlight/sqrt(reps),
                   Forage=mean(meanForage), sdForage=sd(meanForage), seForage=sdForage/sqrt(reps),
                   Land=mean(meanLand), sdLand=sd(meanLand), seLand=sdLand/sqrt(reps),
                   Rest=mean(meanRest), sdRest=sd(meanRest), seRest=sdRest/sqrt(reps),
                   Active=mean(meanActive), sdActive=sd(meanActive), seActive=sdActive/sqrt(reps),
                   SST=mean(meansst), sdSST=sd(meansst), seSST=sdSST/sqrt(reps), SST_max=max(maxsst), SST_min=min(minsst),
                   energy=mean(meanDEE), sdenergy=sd(meanDEE), seenergy=sdenergy/sqrt(reps)) %>%
  #filter(row_number() %% 15 == 0) %>%
  #dplyr::filter(meanBirds>10 & reps==10) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 

average_activity_all2<-average_activity %>%
  ungroup() %>%
  dplyr::group_by(species, rep) %>%
  #dplyr::mutate(SST_scale=scale_to_range(meansst, old_min=min(meansst), old_max=max(meansst), new_min=min(meanDEE), new_max=max(meanDEE))) %>%
  ungroup() %>%
  dplyr::group_by(species, weekNo) %>%
  dplyr::summarise(reps=n_distinct(rep), Flight=mean(meanFlight), sdFlight=sd(meanFlight), seFlight=sdFlight/sqrt(reps),
                   Forage=mean(meanForage), sdForage=sd(meanForage), seForage=sdForage/sqrt(reps),
                   Land=mean(meanLand), sdLand=sd(meanLand), seLand=sdLand/sqrt(reps),
                   Rest=mean(meanRest), sdRest=sd(meanRest), seRest=sdRest/sqrt(reps),
                   Active=mean(meanActive), sdActive=sd(meanActive), seActive=sdActive/sqrt(reps), meanBirds=mean(birdsTot),
                   SST=mean(meansst), sdSST=sd(meansst), seSST=sdSST/sqrt(reps), SST_max=max(maxsst), SST_min=min(minsst),
                   energy=mean(meanDEE), sdenergy=sd(meanDEE), seenergy=sdenergy/sqrt(reps)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

FigureS11<-average_activity_all_colony2 %>%
  ggplot(aes(x=weekNo, y=SST)) +
  geom_pointrange(aes(colour="SST", x=weekNo, y=SST, ymin=SST-1.96*seSST, ymax=SST + 1.96*seSST), cex=0.1, alpha=0.05) +
  geom_line(aes(colour="SST", x=weekNo, y=SST, group=colony), alpha=0.05) +
  geom_ribbon(aes(x=weekNo, y=SST, ymin=SST-1.96*seSST, ymax=SST + 1.96*seSST, fill="SST", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free_y") +
  labs(color="Behaviour", fill="Behaviour") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  #geom_pointrange(aes(x=weekNo, y=SST, ymin=(SST-seSST*1.96), ymax=(SST+seSST*1.96), colour="scaled SST", group=colony), cex=0.1, alpha=0.05) +
  #geom_line(aes(colour="scaled SST", x=weekNo, y=SST, group=colony), alpha=0.05) +
  #geom_ribbon(aes(x=weekNo, y=(SST), ymin=(SST-1.96*seSST), ymax=(SST + 1.96*seSST), fill="scaled SST", group=colony), alpha=0.05) +
  #scale_y_continuous(sec.axis = sec_axis(~.1, name="SST (degrees C)")) +
  ylab("SST (°C)") +
  labs(colour="", fill="") +
  scale_color_manual(values=c( "#E25822"))+
  scale_fill_manual(values=c(  "#E25822")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  #ggtitle("B) Energy budgets") +
  geom_pointrange(data=average_activity_all2,aes(x=weekNo, y=SST, ymin=(SST-seSST*1.96), ymax=(SST+seSST*1.96), colour="SST"), cex=0.1, alpha=0.5) +
  geom_line(data=average_activity_all2,aes(colour="SST", x=weekNo, y=SST), alpha=0.5) +
  geom_ribbon(data=average_activity_all2,aes(x=weekNo, y=(SST), ymin=(SST-1.96*seSST), ymax=(SST + 1.96*seSST), fill="SST"), alpha=0.2) 
#geom_pointrange(data=average_activity_all_species, aes(colour="DEE", x=weekNo, y=energy, ymin=energy-1.96*seenergy, ymax=energy + 1.96*seenergy), cex=0.1) +
#geom_line(data=average_activity_all_species, aes(colour="DEE", x=weekNo, y=energy)) +
#geom_ribbon(data=average_activity_all_species, aes(x=weekNo, y=energy, ymin=energy-1.96*seenergy, ymax=energy + 1.96*seenergy, fill="DEE"), alpha=0.2) 


pdf("./results/figures/supplementary/sst.pdf", width=10, height=7)
plot(FigureS11)
dev.off()

# Making supplementary plots - energy-related variation

print("Making figure 12...")

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
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
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

pdf("./results/figures/supplementary/energy.pdf", width=10, height=7)
plot(FigureS12)
dev.off()

# print
print(paste0("Saving files... output file1 size: ", nrow(energyAll), " rows"))
print(paste0("Saving files... output file 2 size: ", nrow(average_activity), " rows"))
print(paste0("Saving files... output file 2 size: ", nrow(average_activity_colony), " rows"))

# Save output files
print("Saving output files")

# Number 2
output_file1 <- args[4]
print("Saving output file 1")
write.csv(average_activity, file = output_file1, row.names = FALSE) # species-specific budgets

# Number 3
output_file2 <- args[5]
print("Saving output file 2")
write.csv(average_activity_colony, file = output_file2, row.names = FALSE) # Population-specific budgets

# Number 3
output_file3 <- args[6]
print("Saving output file 3")
write.csv(average_activity_colony_individual, file = output_file3, row.names = FALSE) # individual-specific budgets

print("DONE")

