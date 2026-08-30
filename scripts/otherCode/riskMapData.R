# Download other data #

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

# Create a sequence of dates in December for every from 2013-2023

alldates<-data.frame(dates=seq(as.POSIXct("2013-01-01 00:00:00", tz="GMT"), as.POSIXct("2023-12-31 00:00:00", tz="GMT"), 86400))


#### SST & ICE & AIR TEMP ####

wf_get_key()

alldates_unique<-alldates
alldates_unique$year<-as.numeric(substr(alldates_unique$dates, 1, 4))
alldates_unique$month<-as.numeric(substr(alldates_unique$dates, 6, 7))
alldates_unique$day<-as.numeric(substr(alldates_unique$dates, 9, 10))

for (i in 897:nrow(alldates_unique)) {
  
  print(paste("Downloading date", i, "/", nrow(alldates_unique)))  
  
  datesSub<-alldates_unique[i,] 
  
  # Determine location of where data stored
  sstTarget<-"./data/sstDaily/"
  #seaiceTarget<-"./data/ice/"
  #airtempTarget<-"./data/temp/" # air temperature
  
  # This is an example of a request as converted from 
  
  request1 <- list(
    dataset_short_name = "reanalysis-era5-single-levels",
    product_type = "reanalysis",
    variable = "sea_surface_temperature",
    year = datesSub$year[1],
    month = datesSub$month[1],
    day = datesSub$day[1],
    time = "00:00",
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(90, -180, 0, 180),
    target = paste0("sst_", datesSub$year, "_Month_", datesSub$month, "_Day_", datesSub$day, ".nc")
  )
  
  # If you have stored your user login information
  # in the keyring by calling cds_set_key you can
  # call:
  
  file1 <- wf_request(
    request  = request1,  # the request
    transfer = TRUE,     # download the file
    path     = sstTarget       # store data in current working directory
  )
  
  
  
}
