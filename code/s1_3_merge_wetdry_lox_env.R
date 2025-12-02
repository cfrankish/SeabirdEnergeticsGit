### This script merged location, environmental & wet-dry data. It is parallelized by species ###
### The data is time-matched with IRMA data (by ID & session ) - a few more tracks might thus disapear ###
### It is location & time-matched with SST & ice data using rasters from copernicus ###
### Each ten-minute wet-dry data is then assigned to day/night/twilight ###
### Files are then split into individual files for the next step
### Input will be a list of species to iterate through (maybe based existence of activity data? ###
### Main output files will be individual csv files stored in 'birddata_ind/speciesname/id.csv" #### 
### A list of unique ids is also stored in 'birddata_raw/speciesname_IDs.txt"

# load functions
library(ggplot2)
library(dplyr)
library(fields)
#library(rnaturalearth)
#library(rnaturalearthdata)
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
library(seatrackR)

### Step 0: Determine species of interest ####

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_species <- args[1] # Read-in species i
print(paste0("Processing species: ", input_species))

# Match species name with how it is stored in SEATRACK database
possibleSpecies<-data.frame(input_species_names=c("Littleauk", "Northernfulmar", "Atlanticpuffin", "Blackleggedkittiwake", "Commonguillemot", "Brunnichsguillemot"), speciesNormal=c("Little auk", "Northern fulmar",
"Atlantic puffin", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot"))
speciesSub<-subset(possibleSpecies, input_species_names==input_species)
speciesName<-speciesSub$speciesNormal[1] # for now

# Make a new species name that will be used for saving files so there are no special characters
speciesNameSave<-speciesName
speciesNameSave<-gsub(" ", "", speciesNameSave)
speciesNameSave<-gsub("ü", "u", speciesNameSave)
speciesNameSave<-gsub("'", "", speciesNameSave)
speciesNameSave<-gsub("-", "", speciesNameSave)

# Source all necessary functions
source("./scripts/functions.R")

### Step 1: Open up correct wet-dry data ###

print("Step 1: opening up wet-dry file..")

# determine where act data is
list.act.data<-list.files("./data/wetdry_raw/", full.names=TRUE)
act.data<-list.act.data[grepl("csv", list.act.data)] # Subset to the big csv files
act.data.species<-act.data[grepl(input_species, act.data)] # Subset to species of interest

# Open species-specific immersion data
act.data.sub<-fread(act.data.species[1])

### Step 2: Merge wet-dry with positional data by session_id & extract environmental data ###

print("Step 2: merge wet-dry & env data...")

# Define some projections for later
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Determine location IRMA data (for locations)
irma.files<-list.files("./data/positionsIRMA/", full.names=TRUE)
irma.files<-irma.files[grepl("IRMAlocs", irma.files)]
irma.files.df<-data.frame(irma.files)
colnames(irma.files.df)<-c("FileName")
irma.files.df$species<-c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot")

# Determine location of environmental data 
lox.ice<-list.files("./data/ice/", full.names=TRUE)
lox.sst<-list.files("./data/sst/", full.names=TRUE)
lox.distance<-list.files("./data/distCoast/", full.names=TRUE) # I don't think I am using this anymore but just in case

# Summarise existing colonies
colonies<-unique(act.data.sub$colony)

# Open species-specific IRMA data (all individuals)
speciesSub<-act.data.sub$species[1]
irma.file.sub<-subset(irma.files.df, species==speciesSub)
irmaFile<-readRDS(irma.file.sub$FileName)

# open up colony locations
#load("./data/raw/IRMA/Seatrack_export_20220307_summaryTables.R")
colInfo<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")
colInfo<-colInfo %>%
dplyr::group_by(species, colony) %>%
dplyr::count(col_lon, col_lat)

# Subset to species of interest
speciesMatch<-data.frame(species=c("Fratercula_arctica","Uria_lomvia","Uria_aalge","Alle_alle",
           "Fulmarus_glacialis","Rissa_tridactyla"), speciesOG=c("Atlantic puffin", "Brünnich's guillemot", "Common guillemot", "Little auk", "Northern fulmar", "Black-legged kittiwake"))
colInfoSub<-colInfo %>%
dplyr::left_join(speciesMatch, by=c("species")) %>%
dplyr::filter(speciesOG==speciesName)

# Run loop through every colony & save in empty data frame below
infocol<-list()

for (j in 1:length(colonies)) {

#for (j in 1:2) {
 
# Print status update to make it easier to de-bug 
print(paste("Step 2: ", speciesName, ": Colony ", j, "/", length(colonies), sep=" ")) 
  
# Subset to colony j
dataSpeciesCol<-subset(act.data.sub, colony==colonies[j]) 

# Set up list of unique sessions
birdsIDs<-unique(dataSpeciesCol$individ_id)
  
# Run loop through list of names 
infoid<-list()

for (k in 1:length(birdsIDs)) {

#for (k in 1:2) {
  
print(paste("Colony ", j, "/", length(colonies), ", Bird", k, "/", length(birdsIDs), sep=" "))   
  
# Subset individual i
idSub<-birdsIDs[k] 

# Here i make some modifications so that id in my dataset matches the IRMA dataset
idIRMA<-subset(dataSpeciesCol, individ_id==idSub)
idIRMA<-idIRMA$individ_id_IRMA[1]

# Subset immersion & gls file 
dataSpeciesIdSub<-subset(dataSpeciesCol, individ_id==idSub)
dataSpeciesIdSub$date_time<-as.character(dataSpeciesIdSub$date_time)
dataSpeciesIdSub$nchar<-nchar(dataSpeciesIdSub$date_time)
dataSpeciesIdSub$date_time<-ifelse(dataSpeciesIdSub$nchar < 11, paste0(dataSpeciesIdSub$date_time, " 00:00:00"), dataSpeciesIdSub$date_time)
dataSpeciesIdSub$date_time<-as.POSIXct(dataSpeciesIdSub$date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="GMT")

# Subset IRMA file
irmaSub<-subset(irmaFile, individ_id==idIRMA)

# If no match then jump to next...
if (nrow(irmaSub)<1) {
  next
}

# Create empty list for sessions
sessions<-unique(dataSpeciesIdSub$session_id)
lox_allsessions<-list()

for (l in 1:length(sessions)) {
  
print(paste0("Colony ", j, "/", length(colonies), ", Bird", k, "/", length(birdsIDs), ", Sessions", l, "/", length(sessions), sep=" "))   
  
# Subset to specific session (l)
sessionSub<-subset(dataSpeciesIdSub, session_id %in% sessions[l])
date1<-as.POSIXct(min(sessionSub$date_time), tz="UTC") # determine start of time-series
date2<-as.POSIXct(max(sessionSub$date_time), tz="UTC") # determine end  of time-series

# We subset the IRMA file to start & end date of the immersion dataset with some buffer time before and after
irmaSub$timestamp<-as.POSIXct(irmaSub$timestamp, tz="UTC")
irmaSessionSub<-subset(irmaSub, session_id %in% sessions[l])

if (nrow(irmaSessionSub) <1) {

print("No matching IRMA data..")

next}

irmaSub_shorten<-subset(irmaSessionSub, timestamp >= date1 -48*60*60 & timestamp <= date2 + 48*60*60) # Allow a day either side... otherwise June disappears from most datasets

# If no match then next
if (nrow(irmaSub_shorten)<1) {
print("No data")  
next
}

# Now we attach colony location info
colonySub<-subset(colInfoSub, colony==dataSpeciesIdSub$colony[1])

# error message so I can catch mismatches between colony names
if(nrow(colonySub)<1) {
  next
stop(print("No colony match"))
  
}

# We also shorten the activity data so that it matches the length of the IRMA data 
sessionSub$date_time<-as.POSIXct(sessionSub$date_time, tz="UTC")
dataSpeciesIdSub_shorten<-subset(sessionSub, date_time >= min(irmaSub_shorten$timestamp) & date_time <=max(irmaSub_shorten$timestamp))

# If no match then next
if (nrow(dataSpeciesIdSub_shorten)<1) {
print("No data")
next  
}

# Now I extract data for the positions
irmaSub_shorten$month<-as.numeric(substr(irmaSub_shorten$timestamp, 6, 7))
irmaSub_shorten$year<-as.numeric(substr(irmaSub_shorten$timestamp, 1, 4))
irmaSub_shorten$year_month<-paste(irmaSub_shorten$year, "Month", irmaSub_shorten$month, sep="_")

# Determine unique set of month-year for extracting sst & ice data
dataSpeciesIdSub_shorten$year<-substr(dataSpeciesIdSub_shorten$date_time, 1, 4)
dataSpeciesIdSub_shorten$month<-as.numeric(substr(dataSpeciesIdSub_shorten$date_time, 6, 7))
dataSpeciesIdSub_shorten$year_month<-paste(dataSpeciesIdSub_shorten$year, "Month", dataSpeciesIdSub_shorten$month, sep="_")
monthyears<-unique(irmaSub_shorten$year_month)

# Create empty list to save results in
lox_allmonths<-list()

print("Extracting env data...")

for (m in 1:length(monthyears)) {
  
print(paste("MonthYear", m, "/", length(monthyears)))  
  
# Subset specific month-year
monthyearSub<-subset(irmaSub_shorten, year_month %in% monthyears[m])  

# Find relevant distance to coast
distCoast<-rast("./data/distCoast/distCoast_projected.tiff")

# Find relevant SST info
sstRast<-rast(lox.sst[grepl(paste0(monthyearSub$year_month[1], ".nc"), lox.sst)])
values(sstRast)<-values(sstRast) - 273 # to change from kelvin to degrees C
crs(sstRast)<-projection_84
sstProj<-terra::project(sstRast, distCoast)

# Find relevant sea ice info
iceRast<-rast(lox.ice[grepl(paste0(monthyearSub$year_month[1], ".nc"), lox.ice)])
#values(iceRast)<-values(iceRast) - 273 # to change from kelvin to degrees C
crs(iceRast)<-projection_84
iceProj<-project(iceRast, distCoast)

# Project Bird coordinates
coordinates(monthyearSub)<-~lon + lat
projection(monthyearSub)<-projection_84
monthProj<-spTransform(monthyearSub, projection_NA)
lox_sf<-st_as_sf(monthProj)

# Project colony coordinates (so we have an idea of sst if the bird hadn't "left" it's colony
coordinates(colonySub)<-~col_lon + col_lat
projection(colonySub)<-projection_84
colProj<-spTransform(colonySub, projection_NA)
lox_sf_col<-st_as_sf(colProj)

# Create extraction buffer
buffered_lox <- st_buffer(lox_sf, dist = 200000) # buffer of 200 km radius
buffered_lox_col <- st_buffer(lox_sf_col, dist = 200000) # buffer of 200 km radius

# Extract mean + sd (locations)
loxdf<-as.data.frame(monthyearSub)

# Around bird positions
sstVals<-terra::extract(sstProj, buffered_lox)
sstVals_mean<-sstVals %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(sst_mean=mean(sst, na.rm=TRUE), sst_sd=sd(sst, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(sst_mean, sst_sd)
 
# Around colony locations 
sstVals_col<-terra::extract(sstProj, buffered_lox_col)
sstVals_mean_col<-sstVals_col %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(sst_mean_col=mean(sst, na.rm=TRUE), sst_sd_col=sd(sst, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(sst_mean_col, sst_sd_col)  

# Attach results  
loxdf$sst_lox_mean<-sstVals_mean$sst_mean
loxdf$sst_lox_sd<-sstVals_mean$sst_sd

loxdf$sst_col_mean<-sstVals_mean_col$sst_mean_col
loxdf$sst_col_sd<-sstVals_mean_col$sst_sd_col

iceVals<-terra::extract(iceProj, buffered_lox)
iceVals_mean<-iceVals %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(ice_mean=mean(siconc, na.rm=TRUE), ice_sd=sd(siconc, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(ice_mean, ice_sd)
loxdf$ice_mean<-iceVals_mean$ice_mean
loxdf$ice_sd<-iceVals_mean$ice_sd

distVals<-terra::extract(distCoast, buffered_lox)
distVals_mean<-distVals %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(dist_mean=mean(distCoast_projected, na.rm=TRUE), dist_sd=sd(distCoast_projected, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(dist_mean, dist_sd)
loxdf$distCoastKm_mean<-distVals_mean$dist_mean
loxdf$distCoastKm_sd<-distVals_mean$dist_sd

# Check no Nas in all this
na1<-subset(loxdf, is.na(sst_lox_mean))
na2<-subset(loxdf, is.na(ice_mean))

if(nrow(na1)>0) {
  stop(print("NAs in sst layer"))
}

if(nrow(na2)>0) {
  stop(print("NAs in ice layer"))
}

# Add colony locations
colonySub<-subset(colInfoSub, colony==dataSpeciesIdSub$colony[1])
colonySub<-colonySub %>%
dplyr::select(col_lon, col_lat)
loxdf_col<-bind_cols(loxdf, colonySub)

# Add distance to colony
loxdf_col<-loxdf_col %>%
  dplyr::select(individ_id:colony, col_lon, col_lat)
distanceCol<-rdist.earth(x1=as.matrix(loxdf_col[,6:7]), x2=as.matrix(loxdf_col[,21:22]), miles=FALSE) 
loxdf_col$distColonyKm<-diag(distanceCol)

# save to ID list
lox_allmonths<-rbind(lox_allmonths, loxdf_col)
  
} # end of month loop

# Attach to main dataset
lox_allmonths$date_time<-lox_allmonths$timestamp

# Figure out if any missing days. I am annotating 'parts' if there are missing locations in the IRMA dataset
lox_allmonths$date<-as.Date(substr(lox_allmonths$date_time, 1, 10))

# Sub-select to days which exist in the dataset first?
datesMatch<-unique(lox_allmonths$date)

# Divide dataset into sections based on gaps in data (partNo)
lox_allmonths2<-lox_allmonths %>%
ungroup() %>%
dplyr::mutate(daydiff=as.numeric(difftime(lead(date_time), date_time, unit=c("days")))) %>%
replace_na(list(daydiff=0.5)) %>%
dplyr::mutate(part=0) %>%
dplyr::mutate(part=ifelse(daydiff>2, 1, 0)) %>%
dplyr::mutate(partNo=cumsum(part)) %>%
dplyr::select(-c("session_id", "part", "daydiff", "timestamp", "date")) 

print(paste0("Session has ", n_distinct(lox_allmonths2$partNo), " sections"))

data_join<-dataSpeciesIdSub_shorten %>%
  ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::filter(date %in% c(datesMatch)) %>% # filter days without location data out
  dplyr::select(-c(date)) %>%
  dplyr::select(-c(month, year, colony)) %>%
  dplyr::full_join(lox_allmonths2, by=c("species", "date_time", "individ_id", "year_month")) %>%
  arrange(date_time) %>%
  fill(loc_type:partNo, .direction=c("down")) %>%
  dplyr::filter(!is.na(max.cond)) %>%
  dplyr::select(-c(V1, nchar)) %>%
  dplyr::filter(!is.na(lon))

# Check for nas in the colony names
nas<-subset(data_join, is.na(colony))

if(nrow(nas)>0) {
  stop(print("Colony NAs"))
}

# Assign day/night
lox_allmonths_day<-assignTimePeriod(data_join)

lox_allsessions<-rbind(lox_allsessions, lox_allmonths_day)

  
} # End of track year loop

if (length(lox_allsessions)<1) {
  print("No data for this bird")
  next
}

print("Saving individual csv file...")

write.csv(lox_allsessions, file=paste0("./data/birddata_ind/", tolower(speciesNameSave), "/", lox_allsessions$individ_id[1], ".csv"))  

info<-lox_allsessions %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, colony, data_responsible) %>%
  dplyr::count(individ_id)

infoid<-rbind(infoid, info)
  
} # End of individual-specific loop

infocol<-rbind(infocol, infoid)
  
} # End of colony specific loop

# Save unique ids from file as a txt file

print("Saving unique ids...")

ids<-unique(infocol$individ_id)
output_file1 <- args[2]
write.table(ids, file=output_file1, sep = "\t", row.names = FALSE, quote = FALSE, col.names=FALSE)
