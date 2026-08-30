### This script will be used to determine the appropriate study period for our birds 

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

#### Step 1: Determine total possible dates that I will look through ####

# Determine study period
startDate<-"2021-07-01" # Read-in start of study period
endDate<-"2022-06-30" # Read-in end date of study period

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
 # dplyr::filter(days==7) %>%
  dplyr::select(-days)
 
day1<-as.Date(min(dates_weekly$dateKeep)-1)
day2<-as.Date(max(dates_weekly$dateKeep)+1)
date1<-data.frame(dateKeep=day1, doy=0, month=as.numeric(substr(day1, 6, 7)), day=as.numeric(substr(day1, 9, 10)), weekNo=0)
date2<-data.frame(dateKeep=day2, doy=0, month=as.numeric(substr(day2, 6, 7)), day=as.numeric(substr(day2, 9, 10)), weekNo=0)
dates_weekly2<-rbind(date1, dates_weekly, date2)

# Within sessions, subset to 'trackyears' if possible
day1_month<-as.numeric(substr(day1, 6, 7))
day2_month<-as.numeric(substr(day2, 6, 7))

#### Step 2: gather all individual energy files & extract start/end date ####

# Determine where daily files are
allResults<-list.files("./tmp/", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# make an empty list to save results in
startend_dates<-list()

for (i in 3669:length(energyRes_day)) {

print(paste0("Joining bird ", i, "/", length(energyRes_day)))

# Open ID 1 & subset to rep 1
energySub<-fread(energyRes_day[i]) # Open up energy file i
energyRep<-subset(energySub, rep==1) # Subset to one rep

# Add month & day
energyRep$day<-as.numeric(substr(energyRep$date, 9, 10))
energyRep$month<-as.numeric(substr(energyRep$date, 6, 7))
energyRep$year<-as.numeric(substr(energyRep$date, 1, 4))

# Figure out if there are many track years
energyRepYears<-energyRep %>%
ungroup() %>%
dplyr::mutate(track_year=ifelse(month<=day1_month, paste0(year-1, "_", substr(year, 3, 4)), paste0(year, "_", substr(year+1, 3, 4))))

# Remove all end Dates which have a doy that is lower than the start doy
energyTestStart<-energyRepYears %>%
ungroup() %>%
dplyr::left_join(dates_weekly, by=c("month", "day")) %>%
dplyr::select(species, colony, individ_id, session_id, track_year, date) %>%
dplyr::group_by(species, colony, individ_id, session_id, track_year) %>%
dplyr::slice(1) # Subset to first date

energyTestEnd<-energyRepYears %>%
dplyr::select(-doy) %>%
dplyr::left_join(dates_weekly, by=c("month", "day")) %>%
dplyr::group_by(species, colony, individ_id, session_id, track_year) %>%
dplyr::filter(doy>first(doy)) %>% # Make sure end date doy is not lower than start date doy
dplyr::select(species, colony, individ_id, session_id, track_year, date) %>%
dplyr::slice(n()) # Subset to last date

# Join to date dataset to get adjusted doy
energyTestStart$day<-as.numeric(substr(energyTestStart$date, 9, 10)) # Add day for joining
energyTestStart$month<-as.numeric(substr(energyTestStart$date, 6, 7)) # Add month for joining
energyDoyStart<-energyTestStart %>%
dplyr::left_join(dates_weekly, by=c("month", "day")) %>%
dplyr::mutate(start_doy=doy) %>%
dplyr::mutate(startDate=date) %>%
dplyr::select(-c(weekNo, day, month, doy, date, dateKeep))

energyTestEnd$day<-as.numeric(substr(energyTestEnd$date, 9, 10)) # Add day for joining
energyTestEnd$month<-as.numeric(substr(energyTestEnd$date, 6, 7)) # Add month for joining
energyDoyEnd<-energyTestEnd %>%
dplyr::left_join(dates_weekly, by=c("month", "day")) %>%
dplyr::mutate(end_doy=doy) %>%
dplyr::mutate(endDate=date)%>%
dplyr::select(-c(weekNo, day, month, doy, date, dateKeep))

# Save findings in a one-line dataset
energyDoyAll<-energyDoyStart %>%
dplyr::left_join(energyDoyEnd, by=c("species", "colony", "individ_id", "session_id", "track_year")) %>%
dplyr::filter(!is.na(end_doy)) %>%
dplyr::filter(!is.na(start_doy))

# Save results
startend_dates<-rbind(startend_dates, energyDoyAll)

}

### Step 3: summarize results ###

# If multiple rows per bird, collapse to min start and max end
bird_ranges <- startend_dates %>%
  group_by(individ_id) %>%
  summarise(start_doy = min(start_doy),
            end_doy   = max(end_doy),
            .groups="drop")

# Candidate windows
candidates <- expand.grid(
  study_start = unique(bird_ranges$start_doy),
  study_end   = unique(bird_ranges$end_doy)
) %>%
  filter(study_start < study_end)

# Count overlap (vectorized)
results <- candidates %>%
ungroup() %>%
  rowwise() %>%
  dplyr::mutate(
    n_birds = sum(bird_ranges$start_doy <= study_start &
                  bird_ranges$end_doy   >= study_end)
  ) %>%
  ungroup()

# Rank the options
allbest <- results %>%
  arrange(desc(n_birds), desc(study_end - study_start))

# Add actual dates back in
start<-subset(dates_weekly, doy %in% allbest$study_start)
start<-start %>%
dplyr::group_by(doy) %>%
dplyr::slice(1) %>%
dplyr::rename(study_start=doy, study_start_date=dateKeep) %>%
dplyr::select(study_start, study_start_date)
end<-subset(dates_weekly, doy %in% allbest$study_end)
end<-end %>%
dplyr::group_by(doy) %>%
dplyr::slice(1) %>%
dplyr::rename(study_end=doy, study_end_date=dateKeep) %>%
dplyr::select(study_end, study_end_date)

allbestfinal<-allbest %>%
dplyr::left_join(start, by=c("study_start")) %>%
dplyr::left_join(end, by=c("study_end")) %>%
dplyr::mutate(duration=study_end-study_start) %>%
arrange(desc(duration)) %>%
dplyr::mutate(propAll=n_birds/(length(energyRes_day)))
#dplyr::filter(duration > 200)

write.csv(allbestfinal, file="studyperiod.csv")
