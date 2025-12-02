# This script mainly calculates migratory distance #
# As well as weekly variation in energy expenditure used to make Figure 3B #
# Input files are the id catalogue & start/end date of the study peirod #
# Output files are : f"./results/tables/main/table5_migratory_distance.csv" -> which is migratory distance calculated for all individuals
# as well as species (table7_species_mean_deviance.csv) & population-level (table8_population_mean_deviance.csv) weekly deviance in energy expenditure for Figure 3B 
# It also outputs all supplementary figures showing weekly deviation in different behaviours & SST # (Figures S15-S20)


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
#library(lme4)
#library(MuMIn)
library(mgcv)
library(ggplot2)
library(igraph)

args <- commandArgs(trailingOnly = TRUE) # This allows R to read in arguments written in the workflow file

### Step 0: set-up sample size & iteration number ###

print("Step 0: setting up initial parameters")

# Set up minimum sample size & number of iterations
minSampleSize<-5
print(paste0("min sample size per colony is: ", minSampleSize))
reps<-50
print(paste0("min iteration number is: ", reps))

### Step 1: read in id catalogue ###

# This is so we can loop through them later #

print("Step 1: Load id catalogue")
input_file <- args[1]
energyAll <- read.csv(input_file)
#energyAll<-energyAll %>%
#dplyr::group_by(rep, species, colony) %>%
#dplyr::slice_sample(n=5) %>%
#dplyr::mutate(birds=n_distinct(individ_id)) %>%
#dplyr::filter(birds==5)

### Step 2: Calculate migratory distance ###

print("Step 2: calculate migratory distance")

# First we create a raster that will allow us to calculate distances around land 
#world<-ne_countries(scale = "medium", returnclass = "sf") # open world map
#res <- 0.5  # Determine resolution in degrees (is quite coarse to speed things up)
#extent_raster <- extent(-180, 180, -90, 90)  # Full extent for global raster, including the poles
#r <- raster(extent_raster, resolution = res) # Create a raster with the new extent
#world_raster <- rasterize(world, r, field = "labelrank", fun = "last") # rasterize
#world_raster[!is.na(world_raster)] <- 1 # Change na to land
#world_raster[is.na(world_raster)] <- 10000 # Change sea to high number (this means they can still cross land but not preferable
#projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61" # Projection for North Atlantic
#projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # Flat projection
#crs(world_raster)<-projection_84
#NAExtent<-extent(-180, 180, 10, 90) # Crop a little smaller
#naRaster<-crop(world_raster, NAExtent) # Crop to extent round colony to reduce calculation time
#proj4string(naRaster)<-projection_84
#template<-raster::projectExtent(naRaster, projection_NA)
#rasterTrans<-raster::projectRaster(naRaster, template) # project to metric system so I can calculate distances

# Load world polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define projections
projection_NA <- "+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

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
world_raster[!is.na(world_raster)] <- 1
world_raster[is.na(world_raster)] <- 1000000000
rasterTrans<-world_raster

# Now we open up location of breeding colonies so we have starting coordinates
colony.summary.irma<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")
colony.summary.match<-colony.summary.irma %>%
  dplyr::select(colony, col_lon, col_lat) %>%
  dplyr::group_by(colony) %>%
  dplyr::slice(1)

# To do this we will loop through every species/individual
speciesList<-unique(energyAll$species)

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

# make list to save results
allLocations<-list()
migratoryDistance<-list()

for (i in 1:length(speciesList)) {
  
  # Subset species i
  speciesSub<-speciesList[i]
  
  # subset weekly energy file
  speciesSub_ids<-subset(energyAll, species==speciesSub)
  
  # Go through individuals one by one
  birdIds<-unique(speciesSub_ids$individ_id)
  
  # Make list to save results
  birdResAll<-list()
  birdResAll_locations<-list()
  
  for (j in 1:length(birdIds)) {
    
    print(paste0("Step2: calculating for species ", i, ": Bird ", j))  
    
    # Subset bird j
    birdSub<-birdIds[j]
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
    dplyr::group_by(session_id, track_year) %>%
    dplyr::count(date) %>%
    dplyr::mutate(month=as.numeric(substr(date, 6, 7))) %>%
    dplyr::mutate(day=as.numeric(substr(date, 9, 10))) %>%
    dplyr::inner_join(dates_weekly2, by=c("month", "day")) %>%
    ungroup() %>%
    dplyr::group_by(session_id, track_year) %>%
    dplyr::summarise(daysTot=n_distinct(date)) %>%
    dplyr::filter(daysTot==nrow(dates_weekly2)) %>%
    dplyr::mutate(session_year=paste0(session_id, "_", track_year))
    
	if (nrow(dates_sessions) < 1) {
      next}
	
    # Subset to one location per date to reduce file size
    birdcsv_reducted<-birdSub %>%
      ungroup() %>%
      dplyr::group_by(species, colony, individ_id, date) %>%
      dplyr::slice(1) %>%
      dplyr::select(species, colony, individ_id, date, mean.lon, mean.lat, day, month, year, track_year, session_year) %>%
	  dplyr::filter(session_year %in% c(dates_sessions$session_year)) %>%
      dplyr::inner_join(dates_weekly2, by=c("day", "month")) %>%
	  dplyr::filter(weekNo>0) %>%
	  dplyr::select(-c("dateKeep", "doy"))
	  
    print("Calculating distance...")
    
    # transform colony locations
    col_locations<-birdcsv_reducted %>%
      ungroup() %>%
      dplyr::left_join(colony.summary.match, by=c("colony")) %>%
      dplyr::slice(1)
    coordinates(col_locations)<-~col_lon + col_lat
    proj4string(col_locations)<-projection_84
    point1Trans<-spTransform(col_locations, projection_NA)
    
    # Transform it's coordinates
    coordinates(birdcsv_reducted)<-~mean.lon + mean.lat
    proj4string(birdcsv_reducted)<-projection_84
    point2Trans<-spTransform(birdcsv_reducted, projection_NA)
    point2Trans.df<-as.data.frame(point2Trans)
    point1Trans.df<-as.data.frame(point1Trans)
    minLon<-ifelse(min(point2Trans.df$coords.x1)-800000<point1Trans.df$coords.x1 - 800000, min(point2Trans.df$coords.x1)-800000, point1Trans.df$coords.x1-800000)
    maxLon<-ifelse(max(point2Trans.df$coords.x1)+800000>point1Trans.df$coords.x1 + 800000, max(point2Trans.df$coords.x1)+800000, point1Trans.df$coords.x1+800000)
    minLat<-ifelse(min(point2Trans.df$coords.x2)-800000<point1Trans.df$coords.x2 - 800000, min(point2Trans.df$coords.x2)-800000, point1Trans.df$coords.x2-800000)
    maxLat<-ifelse(max(point2Trans.df$coords.x2)+800000>point1Trans.df$coords.x2+800000, max(point2Trans.df$coords.x2)+800000, point1Trans.df$coords.x2+800000)
    cropextent<-extent(minLon, maxLon, minLat, maxLat)
    
    # Crop map quickly
    rasterCrop<-crop(rasterTrans, cropextent)
    transitionRaster <- transition(rasterCrop, mean, directions = 16) # create transition layer
    tr <- geoCorrection(transitionRaster, "c")
    
    # Calculate shortest distance to every point around land
    distance<- gdistance::shortestPath(tr, coordinates(point1Trans), coordinates(point2Trans), output ="SpatialLines")
    distancesf<-st_as_sf(distance)
    lengthKm<-data.frame(distanceKm=st_length(distancesf)/1000)
    
	#pdf("example.pdf")
    #plot(raster(tr))
    #plot(point1Trans, add=T)
    #plot(distance[10], add=T)
	#dev.off()
    
    print("Attaching dataset...")
    
	# Save location of birds for next step
	birdResLox<-data.frame(birdcsv_reducted) %>%
      dplyr::bind_cols(lengthKm)
	
    # Attach results
    birdRes<-data.frame(birdcsv_reducted) %>%
      dplyr::bind_cols(lengthKm) %>%
      ungroup() %>%
      dplyr::group_by(species, colony, individ_id, session_year, track_year) %>%
      dplyr::summarise(MigratoryDistKm=max(distanceKm), startDate=min(date), endDate=max(date))
    
    # save results
    birdResAll<-rbind(birdResAll, birdRes)
	#birdResAll_locations<-rbind(birdResAll_locations, birdResLox)
    
    
  }
  
  migratoryDistance<-rbind(migratoryDistance, birdResAll)
  #allLocations<-rbind(allLocations, birdResAll_locations)
  
}

output_file1 <- args[4]
print("Saving output file 1")
write.csv(migratoryDistance, file = output_file1, row.names = FALSE) # Migratory characterisitcs including negative SST anomaly

### Step 3: Deviance - species-mean  ####

print("Step 3: Calculating deviance for all species")

# Determine where daily files are
allResults<-list.files("./tmp/", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# We will use this to turn DEE into kj.g
species<-data.frame(species=c("Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot", "Little auk", "Atlantic puffin"))
species$allometryCoef<-c(0.765, 0.717, 0.689, 0.689, 0.689, 0.689)

# Define number of species
speciesList<-unique(migratoryDistance$species)

# Make list to save results in
migratoryStrategyRange_weekly<-list()
migratoryStrategyRange_species_weekly<-list()

for (j in 1:length(speciesList)) {
  
  print(paste0("Calculating deviance for species ", j ))
  
  # Find all relevant ids
  speciesSub<-subset(migratoryDistance, species==speciesList[j])
  ids<-unique(speciesSub$individ_id)
  
  # Create list to save results
  speciesDaily<-list()
  
  # Update message
  print("Assembling birds...")
  
  for (k in 1:length(ids)) {
    
    # Print update message
    print(paste0(" Step 5: Species ", j, " Bird ", k))    
    
    # Open id k
	birdID<-ids[k]
	birdID<-gsub("-", "_", birdID)
	birdID<-gsub("ø", "o", birdID)
    birdSub<-fread(energyRes_day[grepl(birdID, energyRes_day)])
	birdSub$year<-as.numeric(substr(birdSub$date, 1, 4))
	birdSub$month<-as.numeric(substr(birdSub$date, 6, 7))
	birdSub$day<-as.numeric(substr(birdSub$date, 9, 10))
	birdSub$track_year<-ifelse(birdSub$month<day1_month, paste0(birdSub$year-1, "_", substr(birdSub$year, 3, 4)), paste0(birdSub$year, "_", substr(birdSub$year + 1, 3, 4)))
    birdSub$session_year<-paste0(birdSub$session_id, "_", birdSub$track_year)
	birdSub<-subset(birdSub, Duration > 100)
	
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
  dplyr::select(-doy) %>%
   dplyr::inner_join(randomSelect, by=c("rep", "session_year")) %>%
   dplyr::inner_join(dates_weekly2, by=c("month", "day")) %>%
   dplyr::filter(weekNo >0) %>%
   ungroup() 
    
    # Subset migratory characteristics to choose the same information
    migrSub<-subset(migratoryDistance, individ_id %in% c(birdSub$individ_id[1]))
    
    # If no match then we go next
    if (nrow(migrSub)<1) {
      next}
    
    #migrSubReps<-unique(migrSub$rep)
    migrSubInfo<-migrSub %>%
      dplyr::ungroup() %>%
      #dplyr::group_by(rep, individ_id) %>%
	  dplyr::group_by(individ_id) %>%
	  dplyr::select(session_year, individ_id, MigratoryDistKm)
      #dplyr::select(rep, track_year, individ_id, MigratoryDistKm, negDevianceSST, meanSST, minSST, maxSST, devianceDEE_nb, negDevianceSST_nb, DEE_cov_nb)
    
    birdSub_random2<-birdSub_random %>%
      #dplyr::filter(rep %in% migrSubReps) %>%
      #dplyr::inner_join(migrSubInfo, by=c("rep", "individ_id", "track_year")) %>%
	  dplyr::inner_join(migrSubInfo, by=c("individ_id", "session_year")) 
    
    # Make sure there are 50 reps
    repsTot<-n_distinct(birdSub_random2$rep)
    
    if (repsTot < reps) {
      stop(print("Error: wrong number of reps..."))
    }
	
	# make sure only one track year per rep
	trackyears<-birdSub_random2 %>%
	ungroup() %>%
	dplyr::group_by(rep) %>%
	dplyr::summarise(years=n_distinct(track_year)) %>%
	dplyr::filter(years>1)
	
	if (nrow(trackyears)>0) {
      stop(print("Error: wrong number of track years..."))
    }
    
    # Save results
    speciesDaily<-rbind(speciesDaily, birdSub_random2)
    
    
  }
  
  
  # Define doy in a standardized way that starts on 1st of July through to end of June (?) 
  speciesDaily$month<-as.numeric(substr(speciesDaily$date, 6, 7))
  speciesDaily$day<-as.numeric(substr(speciesDaily$date, 9, 10))
  
  # Calculate colony weights to attach
  colonySampleSize<-speciesDaily %>%
    dplyr::group_by(colony) %>%
    dplyr::summarise(birds=n_distinct(individ_id)) %>%
    ungroup() %>%
    dplyr::mutate(weights=1/birds) 
  
  # # Make species-specific file for looping through with stats
  
  DevianceCalculation_base<-speciesDaily %>%
    ungroup() %>%
    dplyr::left_join(species, by=c("species")) %>%
    dplyr::mutate(DEEg=DEEkJ/(weight^allometryCoef)) %>%
    #dplyr::inner_join(dates_weekly, by=c("month", "day")) %>%
    dplyr::group_by(rep, species, colony, individ_id, weekNo, track_year) %>%
    dplyr::mutate(days=n_distinct(date)) %>%
    dplyr::filter(days==7) %>%
    dplyr::summarise(weeklyDEE=sum(DEEg), MigratoryDistKm=mean(MigratoryDistKm), weeklySST=mean(sst_random)) %>%
    ungroup() %>%
    dplyr::group_by(individ_id) %>%
    dplyr::mutate(rep=as.numeric(as.factor(rep))) %>%
    dplyr::group_by(rep, species, colony, individ_id) %>%
    dplyr::mutate(meanWeeklyDEE=mean(weeklyDEE), meanWeeklySST=mean(weeklySST), deviationDEE=abs(weeklyDEE-meanWeeklyDEE)/meanWeeklyDEE, sdWeeklyDEE=sd(weeklyDEE), sdWeeklySST=sd(weeklySST)) %>%
    dplyr::group_by(species, rep, colony, individ_id) %>%
    dplyr::summarise(devianceWeekly=sum(deviationDEE),
                     MigratoryDistKm=mean(MigratoryDistKm), covWeekly=(first(sdWeeklyDEE)/first(meanWeeklyDEE))*100, covSSTWeekly=(first(sdWeeklySST)/first(meanWeeklySST))*100) %>%
    dplyr::left_join(colonySampleSize, by=c("colony")) %>%
    ungroup() %>%
    dplyr::group_by(individ_id) %>%
    dplyr::mutate(repsTot=n_distinct(rep))
  
  # Make sure there are still enough reps
  
  minReps<-min(DevianceCalculation_base$repsTot)
  
  if (minReps < reps) {
    stop(print("Error: not enough reps..."))
  }
  
  saveRDS(DevianceCalculation_base, file=paste0("./results/tables/main/deviance_", speciesList[j], ".rds"))
  
  # Save species-specific result
  migratoryStrategyRange_species<-DevianceCalculation_base %>%
    ungroup() %>%
    dplyr::group_by(species, colony, individ_id) %>%
    dplyr::mutate(repsTot=n_distinct(rep)) %>%
    ungroup() %>%
    dplyr::group_by(rep, species, colony) %>%
    dplyr::reframe(averageDeviance=mean(devianceWeekly), averageCOV=mean(covWeekly), sampleSize=n_distinct(individ_id), weight=1/sampleSize, birds=n_distinct(individ_id)) %>%
    ungroup() %>%
    dplyr::group_by(rep, species) %>%
    dplyr::reframe(weightedDeviance=sum(averageDeviance*weight)/sum(weight), weightedCov=sum(averageCOV*weight)/sum(weight))%>%
    ungroup() %>%
    dplyr::group_by(species) %>%
    dplyr::reframe(mean_weightedDeviance=mean(weightedDeviance), mean_weightedCov=mean(weightedCov))%>%
    dplyr::mutate(scale="weekly") %>%
    dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))
  
  # Make a colony-specific results
  
  migratoryStrategyRange<-DevianceCalculation_base %>%
    ungroup() %>%
    dplyr::group_by(species, colony, individ_id) %>%
    dplyr::summarise(meanDeviance=mean(devianceWeekly), meanCov=mean(covWeekly)) %>%
    ungroup() %>%
    dplyr::group_by(species, colony) %>%
    dplyr::summarise(meanDeviance2=mean(meanDeviance), minDeviance=min(meanDeviance), maxDeviance=max(meanDeviance), 
	meanCov2=mean(meanCov), minCov=min(meanCov), maxCov=max(meanCov),birds=n_distinct(individ_id)) %>%
    dplyr::filter(birds>=minSampleSize) %>%
    dplyr::mutate(scale="weekly") %>%
    dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))
  
  # Make a weekly trend for plotting purposes
  
  if ("tActive" %in% colnames(speciesDaily) == FALSE) {
    
    speciesDaily$tActive<-0
  }
  
  if ("DEEkJ_active" %in% colnames(speciesDaily) == FALSE) {
    
    speciesDaily$DEEkJ_active<-0
  }
  
  weeklyDEE<-speciesDaily %>%
    ungroup() %>%
    dplyr::left_join(species, by=c("species")) %>%
    dplyr::mutate(DEEg=DEEkJ/(weight^allometryCoef)) %>%
	dplyr::mutate(DEEg_col=DEEkJ_col/(weight^allometryCoef)) %>%
    dplyr::mutate(DEEg_flight=DEEkJ_flight/(weight^allometryCoef)) %>%
    dplyr::mutate(DEEg_forage=DEEkJ_forage/(weight^allometryCoef)) %>%
    dplyr::mutate(DEEg_active=DEEkJ_active/(weight^allometryCoef)) %>%
    dplyr::mutate(DEEg_land=DEEkJ_restland/(weight^allometryCoef)) %>%
    dplyr::mutate(DEEg_water=DEEkJ_rest/(weight^allometryCoef)) %>%
    #dplyr::select(-doy) %>%
    #dplyr::inner_join(dates_weekly, by=c("month", "day")) %>%
    dplyr::group_by(rep, species, colony, individ_id, track_year, weekNo) %>%
    dplyr::mutate(days=n_distinct(date)) %>%
    dplyr::filter(days==7) %>%
    dplyr::summarise(weeklyDEE=sum(DEEg), weeklyDEE_col=sum(DEEg_col), meanSST=mean(sst_random), meanSST_colony=mean(sst_random_colony), propFlight=sum(tFlight)/168, propActive=sum(tActive)/168, propRest=sum(tRestWater)/168, propForage=sum(tForage)/168,
                     propLand=sum(tLand)/168, totFlightCost=sum(DEEg_flight), totActiveCost=sum(DEEg_active), totRestCost=sum(DEEg_water), totLandCost=sum(DEEg_land), totForageCost=sum(DEEkJ_forage)) %>%
    ungroup() %>%
    dplyr::mutate(propDEE_flight=totFlightCost/weeklyDEE, propDEE_forage=totForageCost/weeklyDEE, propDEE_active=totActiveCost/weeklyDEE, propDEE_water=totRestCost/weeklyDEE, propDEE_land=totLandCost/weeklyDEE) %>%
    dplyr::group_by(individ_id) %>%
    dplyr::mutate(rep=as.numeric(as.factor(rep))) %>%
    ungroup() %>%
    dplyr::group_by(rep, species, colony, individ_id, track_year) %>%
    dplyr::mutate(meanWeeklyDEE=mean(weeklyDEE), sdWeeklyDEE=sd(weeklyDEE), devianceDEE=(weeklyDEE-meanWeeklyDEE)/meanWeeklyDEE, devianceDEE_abs=abs(weeklyDEE-meanWeeklyDEE)/meanWeeklyDEE, totDeviance=sum(devianceDEE_abs),
	covDEE=(first(sdWeeklyDEE)/first(meanWeeklyDEE))*100) %>%
    ungroup() %>%
    dplyr::group_by(rep, species, colony) %>%
    dplyr::mutate(birds=n_distinct(individ_id)) 
  
  repsCheck<-weeklyDEE %>%
    ungroup() %>%
    dplyr::group_by(individ_id) %>%
    dplyr::mutate(repsTot=n_distinct(rep)) 
  
  if (min(repsCheck$repsTot)<reps) {
    stop(print("Error: not enough reps.."))
  }			  
  
  saveRDS(weeklyDEE, file=paste0("./results/tables/main/weeklydeviance_", speciesList[j], ".rds"))
  
  
  # Make species-specific file for looping through with stats
  
  remove(weeklyDEE)
  remove(DevianceCalculation_base)
  
  # Save results 
  migratoryStrategyRange_species_weekly<-rbind(migratoryStrategyRange_species_weekly, migratoryStrategyRange_species)
  migratoryStrategyRange_weekly<-rbind(migratoryStrategyRange_weekly, migratoryStrategyRange)
  
}

### Step 4: make some final plots ####

print("Step 4: making some plots...") 

# Now we make some plots showing weekly deviance in energy & other behaviours #

# First we open the appropriate files
allResults_deviation<-list.files("./results/tables/main/", full.names=TRUE)
deviation_weekly<-allResults_deviation[grepl("/weeklydeviance_", allResults_deviation)]

# make an empty list to save data in
devianceMean_weekly<-list()

for (m in 1:length(deviation_weekly)) {
  
  deviance_sub<-readRDS(deviation_weekly[m])
  
  devianceMean_weekly<-rbind(devianceMean_weekly, deviance_sub)
  
}

# Make a species-mean (weighted for varying sample size):
print("Calculating species weighted averages by week...")

devianceSpeciesRes<-list()

for (p in 1:reps) {

print(paste0("Rep ", p))

devianceSpeciesRep<-devianceMean_weekly %>%
  ungroup() %>%
  dplyr::filter(rep==p) %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::mutate(birds=n_distinct(individ_id)) %>%
  dplyr::filter(birds >=minSampleSize) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, individ_id) %>%
  dplyr::mutate(meanPropFlight=mean(propFlight), meanPropActive=mean(propActive), meanPropRest=mean(propRest), meanPropLand=mean(propLand), meanPropForage=mean(propForage), annualSST=mean(meanSST)) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, individ_id, weekNo) %>%
  dplyr::mutate(devianceFlight=(propFlight - meanPropFlight)/meanPropFlight, devianceActive=(propActive - meanPropActive)/meanPropActive, devianceRest=(propRest - meanPropRest)/meanPropRest, devianceLand=(propLand - meanPropLand)/meanPropLand, devianceForage=(propForage - meanPropForage)/meanPropForage, devianceSST=(meanSST - annualSST)) %>%
  dplyr::group_by(species, colony, individ_id) %>%
  dplyr::mutate(repsTot=n_distinct(rep)) %>%
  ungroup() %>%
  replace_na(list(devianceFlight=0, devianceActive=0, devianceLand=0)) %>%
  dplyr::group_by(rep, species, colony, weekNo) %>%
  dplyr::reframe(averageDeviance=mean(devianceDEE), averageDevianceFlight=mean(devianceFlight), averageDevianceActive=mean(devianceActive, na.rm=TRUE), 
  averageDevianceRest=mean(devianceRest), averageDevianceLand=mean(devianceLand), averageDevianceForage=mean(devianceForage), averageDevianceSST=mean(devianceSST), sampleSize=n_distinct(individ_id), weight=1/sampleSize, birds=n_distinct(individ_id)) %>%
  ungroup() %>%
  #dplyr::filter(birds>=minSampleSize) %>%
  dplyr::group_by(rep, species, weekNo) %>%
  dplyr::reframe(weightedDeviance=sum(averageDeviance*weight)/sum(weight), weightedDevianceFlight=sum(averageDevianceFlight*weight)/sum(weight), weightedDevianceActive=sum(averageDevianceActive*weight)/sum(weight),
  weightedDevianceRest=sum(averageDevianceRest*weight)/sum(weight), weightedDevianceLand=sum(averageDevianceLand*weight)/sum(weight), weightedDevianceForage=sum(averageDevianceForage*weight)/sum(weight),
  weightedDevianceSST=sum(averageDevianceSST*weight)/sum(weight))%>%
  ungroup()
  
  devianceSpeciesRes<-rbind(devianceSpeciesRes, devianceSpeciesRep)
  
  }

devianceSpecies<-devianceSpeciesRes %>%
  dplyr::group_by(species, weekNo) %>%
  dplyr::reframe(mean_weightedDeviance=mean(weightedDeviance), sd_weightedDeviance=sd(weightedDeviance), se=sd_weightedDeviance/sqrt(reps),
  mean_weightedDevianceFlight=mean(weightedDevianceFlight), sd_weightedDevianceFlight=sd(weightedDevianceFlight), seFlight=sd_weightedDevianceFlight/sqrt(reps),
  mean_weightedDevianceActive=mean(weightedDevianceActive), sd_weightedDevianceActive=sd(weightedDevianceActive), seActive=sd_weightedDevianceActive/sqrt(reps),
  mean_weightedDevianceRest=mean(weightedDevianceRest), sd_weightedDevianceRest=sd(weightedDevianceRest), seRest=sd_weightedDevianceRest/sqrt(reps),
  mean_weightedDevianceLand=mean(weightedDevianceLand), sd_weightedDevianceLand=sd(weightedDevianceLand), seLand=sd_weightedDevianceLand/sqrt(reps),
  mean_weightedDevianceForage=mean(weightedDevianceForage), sd_weightedDevianceForage=sd(weightedDevianceForage), seForage=sd_weightedDevianceForage/sqrt(reps),
  mean_weightedDevianceSST=mean(weightedDevianceSST), sd_weightedDevianceSST=sd(weightedDevianceSST), seSST=sd_weightedDevianceSST/sqrt(reps))%>%
  dplyr::mutate(scale="weekly") %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

print("Calculating colony weighted averages...")

devianceColonyRes<-list()

for (q in 1:reps) {

print(paste0("Rep ", q))

devianceColonyRep<-devianceMean_weekly %>%
  ungroup() %>%
  dplyr::filter(rep==q) %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::mutate(birds=n_distinct(individ_id)) %>%
  dplyr::filter(birds >=minSampleSize) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, individ_id) %>%
  dplyr::mutate(meanPropFlight=mean(propFlight), meanPropActive=mean(propActive), meanPropRest=mean(propRest), meanPropLand=mean(propLand), meanPropForage=mean(propForage), annualSST=mean(meanSST)) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, individ_id, weekNo) %>%
  dplyr::mutate(devianceFlight=(propFlight - meanPropFlight)/meanPropFlight, devianceActive=(propActive - meanPropActive)/meanPropActive, devianceRest=(propRest - meanPropRest)/meanPropRest, devianceLand=(propLand - meanPropLand)/meanPropLand, devianceForage=(propForage - meanPropForage)/meanPropForage, devianceSST=(meanSST - annualSST)) %>%
  dplyr::group_by(species, colony, individ_id) %>%
  dplyr::mutate(repsTot=n_distinct(rep)) %>%
  ungroup() %>%
  replace_na(list(devianceFlight=0, devianceActive=0, devianceLand=0)) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, weekNo) %>%
  dplyr::reframe(averageDeviance=mean(devianceDEE), averageDevianceFlight=mean(devianceFlight), averageDevianceActive=mean(devianceActive), 
  averageDevianceRest=mean(devianceRest), averageDevianceLand=mean(devianceLand), averageDevianceForage=mean(devianceForage), averageDevianceSST=mean(devianceSST), sampleSize=n_distinct(individ_id), weight=1/sampleSize, birds=n_distinct(individ_id)) %>%
  ungroup()  
  
  devianceColonyRes<-rbind(devianceColonyRes, devianceColonyRep)
  
  }


deviancePopulation<-devianceColonyRes %>%
  ungroup() %>%
  dplyr::group_by(species, colony, weekNo) %>%
  dplyr::summarise(mean_weightedDeviance=mean(averageDeviance), sd_weightedDeviance=sd(averageDeviance), se=sd_weightedDeviance/sqrt(reps),
  mean_weightedDevianceFlight=mean(averageDevianceFlight), sd_weightedDevianceFlight=sd(averageDevianceFlight), seFlight=sd_weightedDevianceFlight/sqrt(reps),
  mean_weightedDevianceActive=mean(averageDevianceActive), sd_weightedDevianceActive=sd(averageDevianceActive), seActive=sd_weightedDevianceActive/sqrt(reps),
  mean_weightedDevianceRest=mean(averageDevianceRest), sd_weightedDevianceRest=sd(averageDevianceRest), seRest=sd_weightedDevianceRest/sqrt(reps),
  mean_weightedDevianceLand=mean(averageDevianceLand), sd_weightedDevianceLand=sd(averageDevianceLand), seLand=sd_weightedDevianceLand/sqrt(reps),
  mean_weightedDevianceForage=mean(averageDevianceForage), sd_weightedDevianceForage=sd(averageDevianceForage), seForage=sd_weightedDevianceForage/sqrt(reps),
  mean_weightedDevianceSST=mean(averageDevianceSST), sd_weightedDevianceSST=sd(averageDevianceSST), seSST=sd_weightedDevianceSST/sqrt(reps)) %>%
  dplyr::mutate(scale="weekly") %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))


# Now we make weekly plots for different metrics

colony.summary<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")

speciesMatch<-data.frame(speciesLatin=c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle"), 
                         species=c("Brünnich's guillemot", "Black-legged kittiwake", "Common guillemot", "Atlantic puffin", "Northern fulmar", "Little auk"))

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


startMonth<-dates_weekly %>%
  dplyr::filter(day==1)
  
# Here are some dates that I will add #

# Arrival & departure dates - IRMA data
SEATRACK_BreedingDates_20241120<-readRDS("./data/positionsIRMA/SEATRACK_BreedingDates_20241120.rds")

ColonyAttendance<-SEATRACK_BreedingDates_20241120 %>%
  ungroup() %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(doyArrive=ceiling(as.numeric(difftime(start.date.mean, as.Date(paste0(substr(start.date.mean, 1, 4), "-01-01")))))) %>%
  dplyr::mutate(doyDepart=ceiling(as.numeric(difftime(end.date.mean, as.Date(paste0(substr(end.date.mean, 1, 4), "-01-01")))))) %>%
  dplyr::summarise(Depart_early=min(doyDepart), Depart_late=max(doyDepart), Arrive_early=min(doyArrive), Arrive_late=max(doyArrive)) %>%
  rename(speciesLatin=species) %>%
  dplyr::inner_join(speciesMatch, by=c("speciesLatin")) %>%
  dplyr::mutate(Depart_early=ifelse(Depart_early<=244, (0-(244-Depart_early))/7 + 1, (Depart_early-244)/7 + 1)) %>%
  dplyr::mutate(Depart_late=ifelse(Depart_late<=244, (0-(244-Depart_late))/7 + 1, (Depart_late-244)/7 + 1)) %>%
  dplyr::mutate(Arrive_early2=(365-244 + Arrive_early)/7+1) %>%
  dplyr::mutate(Arrive_late2=(365-244 + Arrive_late)/7+1) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

# We will change guillemot dates to those from Merkel et al. 2019 (28 Jan - 18 April)
Arrive_early<-as.Date("2020-01-28")
Arrive_early_doy<-ceiling(as.numeric(difftime(Arrive_early, as.Date(paste0(substr(Arrive_early, 1, 4), "-01-01")))))
Arrive_early_doy2<-(365-244 + Arrive_early_doy)/7+1
Arrive_late<-as.Date("2020-04-18")
Arrive_late_doy<-ceiling(as.numeric(difftime(Arrive_late, as.Date(paste0(substr(Arrive_late, 1, 4), "-01-01")))))
Arrive_late_doy2<-(365-244 + Arrive_late_doy)/7+1

ColonyAttendance$Arrive_early2[5]<-Arrive_early_doy2
ColonyAttendance$Arrive_early2[6]<-Arrive_early_doy2
#ColonyAttendance$Arrive_late2[5]<-Arrive_late_doy2
#ColonyAttendance$Arrive_late2[6]<-Arrive_late_doy2

# We will change the fulmars so that the first arrival is end of December as according to https://onlinelibrary.wiley.com/doi/full/10.1111/ibi.12714
Arrive_early_nofu<-as.Date("2020-12-15")
Arrive_early_doy_nofu<-ceiling(as.numeric(difftime(Arrive_early_nofu, as.Date(paste0(substr(Arrive_early_nofu, 1, 4), "-01-01")))))
Arrive_early_doy2_nofu<-(Arrive_early_doy_nofu - 244)/7+1

ColonyAttendance$Arrive_early2[3]<-Arrive_early_doy2_nofu

# Moulting for fulmars: grissot 

# Earliest date: 6th July and latest end date was 23 December 
startDate<-as.Date("2024-07-06")
endDate<-as.Date("2024-12-23")
moltFulmar<-data.frame(start=ceiling(as.numeric(difftime(startDate, as.Date(paste0(substr(startDate, 1, 4), "-01-01"))))),
end=ceiling(as.numeric(difftime(endDate, as.Date(paste0(substr(endDate, 1, 4), "-01-01"))))),
species="Northern fulmar")
moltFulmar$start<-ifelse(moltFulmar$start<=244, (0-(244-moltFulmar$start))/7 + 1, (moltFulmar$start-244)/7 + 1)
moltFulmar$end<-ifelse(moltFulmar$end<=244, (0-(244-moltFulmar$end))/7 + 1, (moltFulmar$end-244)/7 + 1)
moltFulmar<-moltFulmar %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))
												 
# Visiting fulmars #

visitLand<-subset(devianceSpecies, mean_weightedDevianceLand>0.2)

# Moulting for Puffins (Darby) - based on Table 2

# Moult # 1
startDate_moult1<-data.frame(Start=as.Date("2024-09-01"))
startDate_moult1$doyStart<-ceiling(as.numeric(difftime(startDate_moult1$Start, as.Date(paste0(substr(startDate_moult1$Start, 1, 4), "-01-01")))))
lateStarts<-data.frame(Start=c("2024-09-16", "2024-09-09", "2024-09-19", "2024-09-21"), duration=c(45, 63, 35, 46))
lateStarts$doyStart<-ceiling(as.numeric(difftime(lateStarts$Start, as.Date(paste0(substr(lateStarts$Start, 1, 4), "-01-01")))))
lateStarts$doyEnd<-lateStarts$doyStart + lateStarts$duration
molt1<-data.frame(Start=startDate_moult1$doyStart, End=max(lateStarts$doyEnd))
molt1$Start2<-ifelse(molt1$Start<=244, (0-(244-molt1$Start))/7 + 1, (molt1$Start-244)/7 + 1)
molt1$End2<-ifelse(molt1$End<=244, (0-(244-molt1$End))/7 + 1, (molt1$End-244)/7 + 1)
molt1$species<-"Atlantic puffin"
molt1<-molt1 %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

# Moult # 2
startDate_moult2<-data.frame(Start=as.Date("2024-02-09"))
startDate_moult2$doyStart<-ceiling(as.numeric(difftime(startDate_moult2$Start, as.Date(paste0(substr(startDate_moult2$Start, 1, 4), "-01-01")))))
lateStarts<-data.frame(Start=c("2024-02-20", "2024-02-14"), duration=c(32, 37))
lateStarts$doyStart<-ceiling(as.numeric(difftime(lateStarts$Start, as.Date(paste0(substr(lateStarts$Start, 1, 4), "-01-01")))))
lateStarts$doyEnd<-lateStarts$doyStart + lateStarts$duration
molt2<-data.frame(Start=startDate_moult2$doyStart, End=max(lateStarts$doyEnd))
molt2$Start2<-(365-244 + molt2$Start)/7+1
molt2$End2<-(365-244 + molt2$End)/7+1
molt2$species<-"Atlantic puffin"
molt2<-molt2 %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

# Moulting for guillemots?

# Lila: Synchrony paper 
moltGuillemots<-data.frame(Start=as.Date(c("2019-08-16")), End=as.Date(c("2019-09-15")))
moltGuillemots$doyStart<-ceiling(as.numeric(difftime(moltGuillemots$Start, as.Date(paste0(substr(moltGuillemots$Start, 1, 4), "-01-01")))))
moltGuillemots$doyEnd<-ceiling(as.numeric(difftime(moltGuillemots$End, as.Date(paste0(substr(moltGuillemots$End, 1, 4), "-01-01")))))

moltGuillemots$Start2<-ifelse(moltGuillemots$doyStart<=244, (0-(244-moltGuillemots$doyStart))/7 + 1, (moltGuillemots$doyStart-244)/7 + 1)
moltGuillemots$End2<-ifelse(moltGuillemots$doyEnd<=244, (0-(244-moltGuillemots$doyEnd))/7 + 1, (moltGuillemots$doyEnd-244)/7 + 1)
moltGuillemots$species<-"Common guillemot"
moltGuillemots<-moltGuillemots %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))
# Assume the same for Brunnichs at the moment												 
moltGuillemots2<-moltGuillemots
moltGuillemots2$species<-"Brünnich's guillemot"
moltGuillemots2<-moltGuillemots2 %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

# Mouting for little auks: https://link.springer.com/article/10.1007/s00300-011-1064-4

# Lila: Synchrony paper 
moltLittleAuk<-data.frame(Start=as.Date(c("2019-08-15")), End=as.Date(c("2019-09-15")))
moltLittleAuk$doyStart<-ceiling(as.numeric(difftime(moltLittleAuk$Start, as.Date(paste0(substr(moltLittleAuk$Start, 1, 4), "-01-01")))))
moltLittleAuk$doyEnd<-ceiling(as.numeric(difftime(moltLittleAuk$End, as.Date(paste0(substr(moltLittleAuk$End, 1, 4), "-01-01")))))

moltLittleAuk$Start2<-ifelse(moltLittleAuk$doyStart<=244, (0-(244-moltLittleAuk$doyStart))/7 + 1, (moltLittleAuk$doyStart-244)/7 + 1)
moltLittleAuk$End2<-ifelse(moltLittleAuk$doyEnd<=244, (0-(244-moltLittleAuk$doyEnd))/7 + 1, (moltLittleAuk$doyEnd-244)/7 + 1)
moltLittleAuk$species<-"Little auk"
moltLittleAuk<-moltLittleAuk %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))
												 
# Moulting for black-legged kittiwakes: https://www.int-res.com/articles/meps2011/435/m435p251.pdf + https://www.sciencedirect.com/science/article/pii/S0048969724060510#bb0145
moltBLK<-data.frame(Start=as.Date(c("2019-05-01")), End=as.Date(c("2019-12-31")))
moltBLK$doyStart<-ceiling(as.numeric(difftime(moltBLK$Start, as.Date(paste0(substr(moltBLK$Start, 1, 4), "-01-01")))))
moltBLK$doyEnd<-ceiling(as.numeric(difftime(moltBLK$End, as.Date(paste0(substr(moltBLK$End, 1, 4), "-01-01")))))

moltBLK$Start2<-ifelse(moltBLK$doyStart<=244, (0-(244-moltBLK$doyStart))/7 + 1, (moltBLK$doyStart-244)/7 + 1)
moltBLK$End2<-ifelse(moltBLK$doyEnd<=244, (0-(244-moltBLK$doyEnd))/7 + 1, (moltBLK$doyEnd-244)/7 + 1)
moltBLK$species<-"Black-legged kittiwake"
moltBLK<-moltBLK %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

# Determine max values for plotting text #

maxVals<-devianceSpecies %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::mutate(upper=mean_weightedDeviance + 1.96*se) %>%
arrange(desc(upper)) %>%
dplyr::slice(1) %>%
dplyr::select(species, upper) %>%
dplyr::mutate(height=upper + 0.7*upper, height2=upper + 0.2*upper)

ColonyAttendance2<-ColonyAttendance %>%
dplyr::inner_join(maxVals, by=c("species"))


Figure3B<-deviancePopulation %>%
  ggplot(aes(x=weekNo, y=mean_weightedDeviance)) +
  #geom_pointrange(aes(colour="Deviance", x=weekNo, y=mean_weightedDeviance, ymin=mean_weightedDeviance-1.96*se, ymax=mean_weightedDeviance + 1.96*se), cex=0.1, alpha=0.05) +
  geom_line(aes(colour="Population", x=weekNo, y=mean_weightedDeviance, group=colony), alpha=0.1) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=weekNo, y=mean_weightedDeviance, ymin=mean_weightedDeviance-1.96*se, ymax=mean_weightedDeviance + 1.96*se, fill="Population", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free_y") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Deviance in weekly energy expenditure") +
  scale_color_manual(values=c("#0072b2", "#e6550d", "#bcbddc"))+
  scale_fill_manual(values=c( "#0072b2", "#e6550d", "#bcbddc")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  ggtitle("B)") +
  geom_line(data=devianceSpecies, aes(colour="Species", x=weekNo, y=mean_weightedDeviance)) +
  geom_ribbon(data=devianceSpecies, aes(x=weekNo, y=mean_weightedDeviance, ymin=mean_weightedDeviance-1.96*se, ymax=mean_weightedDeviance + 1.96*se, fill="Species"), alpha=0.2) +
  geom_segment(data=ColonyAttendance2, aes(x=Depart_early, xend=Depart_late, y=height2), linetype="dashed") +
  geom_text(data=ColonyAttendance2, aes(x=Depart_early + 3, y=height, label="Departure"), size=2) +
  geom_segment(data=ColonyAttendance2, aes(x=Arrive_early2, xend=Arrive_late2, y=height2), linetype="dashed") +
  geom_text(data=ColonyAttendance2, aes(x=Arrive_early2 + 3, y=height, label="Return"), size=2) +
  geom_segment(data=moltFulmar, aes(x=start, xend=end, y=0.3), linetype="dashed", color="red") +
  geom_text(data=moltFulmar, aes(x=start + 3, y=0.35, label="Moult"), size=2, color="red") +
  geom_segment(data=molt1, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=molt1, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") +
  geom_segment(data=molt2, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=molt2, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots2, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots2, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") +
  geom_segment(data=moltLittleAuk, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=moltLittleAuk, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") +
  geom_segment(data=moltBLK, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=moltBLK, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") +
  labs(color="", fill="")

pdf("./results/figures/main/Figure3B.pdf", width=9, height=6)
grid.arrange(Figure3B)
dev.off()

# Flight #

FigureS16<-deviancePopulation %>%
  ggplot(aes(x=weekNo, y=mean_weightedDevianceFlight)) +
  geom_line(aes(colour="Population", x=weekNo, y=mean_weightedDevianceFlight, group=colony), alpha=0.1) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=weekNo, y=mean_weightedDevianceFlight, ymin=mean_weightedDevianceFlight-1.96*seFlight, ymax=mean_weightedDevianceFlight + 1.96*seFlight, fill="Population", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free_y") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Deviance in weekly time spent in flight") +
  scale_color_manual(values=c("#0072b2", "#e6550d", "#bcbddc"))+
  scale_fill_manual(values=c( "#0072b2", "#e6550d", "#bcbddc")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  geom_line(data=devianceSpecies, aes(colour="Species", x=weekNo, y=mean_weightedDevianceFlight)) +
  geom_ribbon(data=devianceSpecies, aes(x=weekNo, y=mean_weightedDevianceFlight, ymin=mean_weightedDevianceFlight-1.96*seFlight, ymax=mean_weightedDevianceFlight + 1.96*seFlight, fill="Species"), alpha=0.2) +
  labs(color="", fill="") +
  geom_segment(data=ColonyAttendance, aes(x=Depart_early, xend=Depart_late, y=1), linetype="dashed") +
  geom_text(data=ColonyAttendance, aes(x=Depart_early + 3, y=1.2, label="Departure"), size=2) +
  geom_segment(data=ColonyAttendance, aes(x=Arrive_early2, xend=Arrive_late2, y=1), linetype="dashed") +
  geom_text(data=ColonyAttendance, aes(x=Arrive_early2 + 3, y=1.2, label="Return"), size=2) +
  geom_segment(data=moltFulmar, aes(x=start, xend=end, y=2), linetype="dashed", color="red") +
  geom_text(data=moltFulmar, aes(x=start + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=molt1, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=molt1, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=molt2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=molt2, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots2, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltLittleAuk, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltLittleAuk, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltBLK, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=moltBLK, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") 

pdf("./results/figures/supplementary/FigureS16.pdf", width=9, height=6)
grid.arrange(FigureS16)
dev.off()

# Active #

FigureS17<-deviancePopulation %>%
  dplyr::filter(!species %in% c("Northern fulmar", "Black-legged kittiwake")) %>%
  ggplot(aes(x=weekNo, y=mean_weightedDevianceActive)) +
  #geom_pointrange(aes(colour="Deviance", x=weekNo, y=mean_weightedDeviance, ymin=mean_weightedDeviance-1.96*se, ymax=mean_weightedDeviance + 1.96*se), cex=0.1, alpha=0.05) +
  geom_line(aes(colour="Population", x=weekNo, y=mean_weightedDevianceActive, group=colony), alpha=0.1) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=weekNo, y=mean_weightedDevianceActive, ymin=mean_weightedDevianceActive-1.96*seActive, ymax=mean_weightedDevianceActive + 1.96*seActive, fill="Population", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free_y") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Deviance in weekly time spent active") +
  scale_color_manual(values=c("#0072b2", "#e6550d", "#bcbddc"))+
  scale_fill_manual(values=c( "#0072b2", "#e6550d", "#bcbddc")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  geom_line(data=subset(devianceSpecies, !species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(colour="Species", x=weekNo, y=mean_weightedDevianceActive)) +
  geom_ribbon(data=subset(devianceSpecies, !species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=weekNo, y=mean_weightedDevianceActive, ymin=mean_weightedDevianceActive-1.96*seActive, ymax=mean_weightedDevianceActive + 1.96*seActive, fill="Species"), alpha=0.2) +
  labs(color="", fill="") +
  geom_segment(data=subset(ColonyAttendance, !species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=Depart_early, xend=Depart_late, y=1), linetype="dashed") +
  geom_text(data=subset(ColonyAttendance, !species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=Depart_early + 3, y=1.2, label="Departure"), size=2) +
  geom_segment(data=subset(ColonyAttendance, !species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=Arrive_early2, xend=Arrive_late2, y=1), linetype="dashed") +
  geom_text(data=subset(ColonyAttendance, !species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=Arrive_early2 + 3, y=1.2, label="Return"), size=2) +
  #geom_segment(data=moltFulmar, aes(x=start, xend=end, y=2), linetype="dashed", color="red") +
  #geom_text(data=moltFulmar, aes(x=start + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=molt1, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=molt1, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=molt2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=molt2, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots2, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltLittleAuk, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltLittleAuk, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") 
  #geom_segment(data=moltBLK, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  #geom_text(data=moltBLK, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") 

pdf("./results/figures/supplementary/FigureS17.pdf", width=9, height=6)
grid.arrange(FigureS17)
dev.off()

# Rest #

FigureS18<-deviancePopulation %>%
  ggplot(aes(x=weekNo, y=mean_weightedDevianceRest)) +
  geom_line(aes(colour="Population", x=weekNo, y=mean_weightedDevianceRest, group=colony), alpha=0.1) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=weekNo, y=mean_weightedDevianceRest, ymin=mean_weightedDevianceRest-1.96*seRest, ymax=mean_weightedDevianceRest + 1.96*seRest, fill="Population", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free_y") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Deviance in weekly time spent resting (on water)") +
  scale_color_manual(values=c("#0072b2", "#e6550d", "#bcbddc"))+
  scale_fill_manual(values=c( "#0072b2", "#e6550d", "#bcbddc")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  geom_line(data=devianceSpecies, aes(colour="Species", x=weekNo, y=mean_weightedDevianceRest)) +
  geom_ribbon(data=devianceSpecies, aes(x=weekNo, y=mean_weightedDevianceRest, ymin=mean_weightedDevianceRest-1.96*seRest, ymax=mean_weightedDevianceRest + 1.96*seRest, fill="Species"), alpha=0.2) +
geom_segment(data=ColonyAttendance, aes(x=Depart_early, xend=Depart_late, y=1), linetype="dashed") +
  geom_text(data=ColonyAttendance, aes(x=Depart_early + 3, y=1.2, label="Departure"), size=2) +
  geom_segment(data=ColonyAttendance, aes(x=Arrive_early2, xend=Arrive_late2, y=1), linetype="dashed") +
  geom_text(data=ColonyAttendance, aes(x=Arrive_early2 + 3, y=1.2, label="Return"), size=2) +
  geom_segment(data=moltFulmar, aes(x=start, xend=end, y=2), linetype="dashed", color="red") +
  geom_text(data=moltFulmar, aes(x=start + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=molt1, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=molt1, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=molt2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=molt2, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots2, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltLittleAuk, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=moltLittleAuk, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltBLK, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=moltBLK, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") +
  labs(color="", fill="")

pdf("./results/figures/supplementary/FigureS18.pdf", width=9, height=6)
grid.arrange(FigureS18)
dev.off()

# Land #

FigureS19<-deviancePopulation %>%
  ggplot(aes(x=weekNo, y=mean_weightedDevianceLand)) +
  geom_line(aes(colour="Population", x=weekNo, y=mean_weightedDevianceLand, group=colony), alpha=0.1) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=weekNo, y=mean_weightedDevianceLand, ymin=mean_weightedDevianceLand-1.96*seLand, ymax=mean_weightedDevianceLand + 1.96*seLand, fill="Population", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free_y") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Deviance in weekly time spent on land") +
  scale_color_manual(values=c("#0072b2", "#e6550d", "#bcbddc"))+
  scale_fill_manual(values=c( "#0072b2", "#e6550d", "#bcbddc")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  geom_line(data=devianceSpecies, aes(colour="Species", x=weekNo, y=mean_weightedDevianceLand)) +
  geom_ribbon(data=devianceSpecies, aes(x=weekNo, y=mean_weightedDevianceLand, ymin=mean_weightedDevianceLand-1.96*seLand, ymax=mean_weightedDevianceLand + 1.96*seLand, fill="Species"), alpha=0.2) +
  geom_segment(data=ColonyAttendance, aes(x=Depart_early, xend=Depart_late, y=1), linetype="dashed") +
  geom_text(data=ColonyAttendance, aes(x=Depart_early + 3, y=1.9, label="Departure"), size=2) +
  geom_segment(data=ColonyAttendance, aes(x=Arrive_early2, xend=Arrive_late2, y=1), linetype="dashed") +
  geom_text(data=ColonyAttendance, aes(x=Arrive_early2 + 3, y=1.9, label="Return"), size=2) +
  geom_segment(data=moltFulmar, aes(x=start, xend=end, y=3), linetype="dashed", color="red") +
  geom_text(data=moltFulmar, aes(x=start + 3, y=3.5, label="Moult"), size=2, color="red") +
  geom_segment(data=molt1, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=molt1, aes(x=Start2 + 3, y=2.8, label="Moult"), size=2, color="red") +
  geom_segment(data=molt2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  geom_text(data=molt2, aes(x=Start2 + 3, y=2.8, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots, aes(x=Start2, xend=End2, y=3), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots, aes(x=Start2 + 3, y=3.5, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots2, aes(x=Start2, xend=End2, y=3), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots2, aes(x=Start2 + 3, y=3.5, label="Moult"), size=2, color="red") +
  geom_segment(data=moltLittleAuk, aes(x=Start2, xend=End2, y=3), linetype="dashed", color="red") +
  geom_text(data=moltLittleAuk, aes(x=Start2 + 3, y=3.5, label="Moult"), size=2, color="red") +
  geom_segment(data=moltBLK, aes(x=Start2, xend=End2, y=0.4), linetype="dashed", color="red") +
  geom_text(data=moltBLK, aes(x=Start2 + 3, y=1, label="Moult"), size=2, color="red") +
  labs(color="", fill="")

pdf("./results/figures/supplementary/FigureS19.pdf", width=9, height=6)
grid.arrange(FigureS19)
dev.off()

# Forage #

FigureS20<-deviancePopulation %>%
filter(species %in% c("Northern fulmar", "Black-legged kittiwake")) %>%
  ggplot(aes(x=weekNo, y=mean_weightedDevianceForage)) +
  geom_line(aes(colour="Population", x=weekNo, y=mean_weightedDevianceForage, group=colony), alpha=0.1) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=weekNo, y=mean_weightedDevianceForage, ymin=mean_weightedDevianceForage-1.96*seForage, ymax=mean_weightedDevianceForage + 1.96*seForage, fill="Population", group=colony), alpha=0.05) +
  facet_wrap(~species, scales="free_y") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Deviance in weekly time spent foraging") +
  scale_color_manual(values=c("#0072b2", "#e6550d", "#bcbddc"))+
  scale_fill_manual(values=c( "#0072b2", "#e6550d", "#bcbddc")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  geom_line(data=subset(devianceSpecies, species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(colour="Species", x=weekNo, y=mean_weightedDevianceForage)) +
  geom_ribbon(data=subset(devianceSpecies, species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=weekNo, y=mean_weightedDevianceForage, ymin=mean_weightedDevianceForage-1.96*se, ymax=mean_weightedDevianceForage + 1.96*se, fill="Species"), alpha=0.2) +
  geom_segment(data=subset(ColonyAttendance, species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=Depart_early, xend=Depart_late, y=1), linetype="dashed") +
  geom_text(data=subset(ColonyAttendance, species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=Depart_early + 3, y=1.2, label="Departure"), size=2) +
  geom_segment(data=subset(ColonyAttendance, species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=Arrive_early2, xend=Arrive_late2, y=1), linetype="dashed") +
  geom_text(data=subset(ColonyAttendance, species %in% c("Northern fulmar", "Black-legged kittiwake")), aes(x=Arrive_early2 + 3, y=1.2, label="Return"), size=2) +
  geom_segment(data=moltFulmar, aes(x=start, xend=end, y=2), linetype="dashed", color="red") +
  geom_text(data=moltFulmar, aes(x=start + 3, y=2.2, label="Moult"), size=2, color="red") +
  #geom_segment(data=molt1, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  #geom_text(data=molt1, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  #geom_segment(data=molt2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  #geom_text(data=molt2, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  #geom_segment(data=moltGuillemots, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  #geom_text(data=moltGuillemots, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  #geom_segment(data=moltGuillemots2, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  #geom_text(data=moltGuillemots2, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  #geom_segment(data=moltLittleAuk, aes(x=Start2, xend=End2, y=2), linetype="dashed", color="red") +
  #geom_text(data=moltLittleAuk, aes(x=Start2 + 3, y=2.2, label="Moult"), size=2, color="red") +
  labs(color="", fill="") +
  geom_segment(data=moltBLK, aes(x=Start2, xend=End2, y=0.3), linetype="dashed", color="red") +
  geom_text(data=moltBLK, aes(x=Start2 + 3, y=0.35, label="Moult"), size=2, color="red") 

pdf("./results/figures/supplementary/FigureS20.pdf", width=6, height=4)
grid.arrange(FigureS20)
dev.off()

# SST #

maxVals<-devianceSpecies %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::mutate(upper=mean_weightedDevianceSST + 1.96*seSST) %>%
arrange(desc(upper)) %>%
dplyr::slice(1) %>%
dplyr::select(species, upper) %>%
dplyr::mutate(height=upper + 0.7*upper, height2=upper + 0.2*upper, heightM=upper + 1.6*upper, heightM2=upper + 2.4*upper)

ColonyAttendance2<-ColonyAttendance %>%
dplyr::inner_join(maxVals, by=c("species"))

moltFulmar2<-moltFulmar %>%
dplyr::inner_join(maxVals, by=c("species"))

molt1_2<-molt1 %>%
dplyr::inner_join(maxVals, by=c("species"))

molt2_2<-molt2 %>%
dplyr::inner_join(maxVals, by=c("species"))

moltGuillemots_2<-moltGuillemots %>%
dplyr::inner_join(maxVals, by=c("species"))

moltGuillemots2_2<-moltGuillemots2 %>%
dplyr::inner_join(maxVals, by=c("species"))

moltBLK2<-moltBLK %>%
dplyr::inner_join(maxVals, by=c("species"))

moltLittleAuk2<-moltLittleAuk %>%
dplyr::inner_join(maxVals, by=c("species"))

FigureS15<-deviancePopulation %>%
  ggplot(aes(x=weekNo, y=mean_weightedDevianceSST)) +
  geom_line(aes(colour="Population", x=weekNo, y=mean_weightedDevianceSST, group=colony), alpha=0.1) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=weekNo, y=mean_weightedDevianceSST, ymin=mean_weightedDevianceSST-1.96*seSST, ymax=mean_weightedDevianceSST + 1.96*seSST, fill="Population", group=colony), alpha=0.05) +
  facet_wrap(~species, nrow=2, scales="free_y") +
  theme_bw() +
  scale_x_continuous(breaks=startMonth$weekNo, labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")) +
  xlab("") +
  ylab("Deviance in weekly SST") +
  scale_color_manual(values=c("#0072b2", "#e6550d", "#bcbddc"))+
  scale_fill_manual(values=c( "#0072b2", "#e6550d", "#bcbddc")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position ="bottom") +
  geom_line(data=devianceSpecies, aes(colour="Species", x=weekNo, y=mean_weightedDevianceSST)) +
  geom_ribbon(data=devianceSpecies, aes(x=weekNo, y=mean_weightedDevianceSST, ymin=mean_weightedDevianceSST-1.96*seSST, ymax=mean_weightedDevianceSST + 1.96*seSST, fill="Species"), alpha=0.2) +
  labs(color="", fill="") +
  geom_segment(data=ColonyAttendance2, aes(x=Depart_early, xend=Depart_late, y=height2), linetype="dashed") +
  geom_text(data=ColonyAttendance2, aes(x=Depart_early + 3, y=height, label="Departure"), size=2) +
  geom_segment(data=ColonyAttendance2, aes(x=Arrive_early2, xend=Arrive_late2, y=height2), linetype="dashed") +
  geom_text(data=ColonyAttendance2, aes(x=Arrive_early2 + 3, y=height, label="Return"), size=2) +
  geom_segment(data=moltFulmar2, aes(x=start, xend=end, y=heightM), linetype="dashed", color="red") +
  geom_text(data=moltFulmar2, aes(x=start + 3, y=heightM2, label="Moult"), size=2, color="red") +
  geom_segment(data=molt1_2, aes(x=Start2, xend=End2, y=heightM), linetype="dashed", color="red") +
  geom_text(data=molt1_2, aes(x=Start2 + 3, y=heightM2, label="Moult"), size=2, color="red") +
  geom_segment(data=molt2_2, aes(x=Start2, xend=End2, y=heightM), linetype="dashed", color="red") +
  geom_text(data=molt2_2, aes(x=Start2 + 3, y=heightM2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots_2, aes(x=Start2, xend=End2, y=heightM), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots_2, aes(x=Start2 + 3, y=heightM2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltGuillemots2_2, aes(x=Start2, xend=End2, y=heightM), linetype="dashed", color="red") +
  geom_text(data=moltGuillemots2_2, aes(x=Start2 + 3, y=heightM2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltLittleAuk2, aes(x=Start2, xend=End2, y=heightM), linetype="dashed", color="red") +
  geom_text(data=moltLittleAuk2, aes(x=Start2 + 3, y=heightM2, label="Moult"), size=2, color="red") +
  geom_segment(data=moltBLK2, aes(x=Start2, xend=End2, y=heightM), linetype="dashed", color="red") +
  geom_text(data=moltBLK2, aes(x=Start2 + 3, y=heightM2, label="Moult"), size=2, color="red") 

pdf("./results/figures/supplementary/FigureS15.pdf", width=9, height=6)
grid.arrange(FigureS15)
dev.off()

# Save output files
print("Saving output files...")

# Number 3
output_file3 <- args[5]
print("Saving output file 3")
write.csv(migratoryStrategyRange_species_weekly, file = output_file3, row.names = FALSE) # Species-mean deviance

# Number 4
output_file4 <- args[6]
print("Saving output file 4")
write.csv(migratoryStrategyRange_weekly, file = output_file4, row.names = FALSE) # Population-mean deviance





