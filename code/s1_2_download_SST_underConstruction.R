### This script downloads SST data to match with positions of tracked birds. It is parallelized by species ###
### The data is then time-matched with IRMA data (by session ID & time) - so locations are added & wet/dry data with positions is also removed ###
### Time-steps are characterized by occuring during day, night or twilight ###
### Species-specific files will be split into individual bird files for paralellizing the next step ###
### Input will be a list of species to iterate through ###
### Output will be ... ###

# load functions
library(ggplot2)
library(dplyr)
library(fields)
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
library(ecmwfr)

### Step 1: Determine dates for which we need environmental data ####

print("Step 1: determining dates for which we need data ...")

# determine where act data is
list.act.data<-list.files("./data/wetdry_raw/", full.names=TRUE)
act.data<-list.act.data[grepl("csv", list.act.data)] # Subset to the big csv files

# make empty list to save result in
alldatesSum<-list()

for (i in 1:length(act.data)) {

print(paste0("Species ", i, "/", length(act.data))

# Open species-specific immersion data
act.data.sub<-fread(act.data[i])

# Summarize all existing dates
act.data.sub$date<-substr(act.data.sub$date_time, 1, 10) # make a date column
alldates<-act.data.sub %>%
  ungroup() %>%
  dplyr::count(date)
  
# Save result
alldatesSum<-rbind(alldatesSum, alldates)

}

# Summarize unique dates (by month-year combo)
alldatesSum$year<-as.numeric(substr(alldatesSum$date, 1, 4)) # Substr year
alldatesSum$month<-as.numeric(substr(alldatesSum$date, 6, 7)) # Subst month

alldates_unique<-alldatesSum %>%
  dplyr::group_by(year) %>%
  dplyr::count(month) 
  
#### Step 2: Download ice/sst rasters for these month-year combinations ####

print("Step 2: Downloading environmental rasters...")

options(keyring_backend = "file")
wf_set_key(user='cc75d4bb-c933-4045-a4df-7f72bcee08cd', key='01a598e6-074d-449d-afeb-4e5d5a6e0757', service='cds')

# Set your credentials directly
wf_set_key(
  user = "cc75d4bb-c933-4045-a4df-7f72bcee08cd",
  key = "01a598e6-074d-449d-afeb-4e5d5a6e0757",
  service = "cds"
)

#wf_set_key(key = "7bf5a340-e1fc-4ba2-a0e9-3aee477b897b", )
#wf_get_key()

# To get around the keychain thing, I have made a yml file in the folder with my keychain credentials in it #



for (i in 1:nrow(alldates_unique)) {
  
print(paste("Downloading date", i, "/", nrow(alldates_unique)))  
  
# Subset to date i
datesSub<-alldates_unique[i,]  

# See if they exist or not
sstTarget<-list.files("./data/sst/", full.names=TRUE)
sstTargetSub<-sstTarget[grepl(datesSub$year, sstTarget)] # subset to year x
sstTargetSub_month<-sstTargetSub[grepl(paste("Month", datesSub$month, sep="_"), sstTargetSub)]

if (length(sstTargetSub_month)>0) {
  next
}

# Determine location of where data stored
sstTarget<-"./data/sst/"
seaiceTarget<-"./data/ice/"

# This is an example of a request as converted from 
request1 <- list(
  dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
  product_type = "monthly_averaged_reanalysis",
  variable = "sea_surface_temperature",
  year = datesSub$year[1],
  month = datesSub$month[1],
  time = "00:00",
  data_format = "netcdf",
  download_format = "unarchived",
  area = c(90, -180, 0, 180),
  target = paste0("sst_", datesSub$year, "_Month_", datesSub$month, ".nc")
)

# This is an example of a request as converted from 
request2 <- list(
  dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
  product_type = "monthly_averaged_reanalysis",
  variable = "sea_ice_cover",
  year = datesSub$year[1],
  month = datesSub$month[1],
  time = "00:00",
  data_format = "netcdf",
  download_format = "unarchived",
  area = c(90, -180, 0, 180),
  target = paste0("ice_", datesSub$year, "_Month_", datesSub$month, ".nc")
)

# If you have stored your user login information
# in the keyring by calling cds_set_key you can
# call:
file1 <- wf_request(
  request  = request1,  # the request
  transfer = TRUE,     # download the file
  path     = sstTarget       # store data in current working directory
)

file2 <- wf_request(
  request  = request2,  # the request
  transfer = TRUE,     # download the file
  path     = seaiceTarget,
  retry = 10# store data in current working directory
)

}