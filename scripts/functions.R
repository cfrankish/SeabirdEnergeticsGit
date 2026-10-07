# Various functions used for estimating activity budgets and energy expenditure # 

#### Activity functions ####

# Function to calculate time spent in different activities (tFlight, tRestWater, tActive, tLand and tForage)
# Input data 1: species is common species name (can only be be 'Black-legged kittiwake', 'Northern fulmar', 'Common guillemot', 'Brünnich's guillemot', 'Atlantic puffin' or 'Little auk') 
# Input data 2: Data is immersion data standardized to 0-1 (with 1= 100% wet and 0 = 100% dry) for a given individual
# Input data 3: IrmaData is SEATRACK's irma dataset subset to the individual bird of interest
# The function essentially re-directs to species-specific activity budget functions (described below) and returns two objects:
# Object 1 is the dataset annotated with different behaviours for exploratory use
# Object 2 is time spent in activity per day

calculateTimeInActivity<-function(species, data, irmaData) {
 
# Print some of the model parameters for checking purposes
print(paste0("L1 = ", data$L1[1], " mins"))
print(paste0("Th1 = ", data$Th1[1], " % wet")) 
print(paste0("Th2 = ", data$Th2[1], " % wet")) 
print(paste0("L1_colony; range =  ", data$L1_colony_min[1], "-", data$L1_colony_max[1], " mins"))
print(paste0("dist_colony = ", data$dist_colony[1], " km"))
print(paste0("PLand = ", data$pLand_prob[1], " %")) 
print(paste0("c = ", data$c[1], " %")) 
 
# Determine number of sessions

sessionNo<-unique(data$session_id)

# For loop estimates activity budgets for separate sessions 

timeActivity_sessions<-list() # List to save results in

for (session in 1:length(sessionNo)) {

print(paste0("Activity for session ", session))

dataSub<-subset(data, session_id %in% sessionNo[session]) 
  
  if (species=="Black-legged kittiwake") {
    
    print(paste("Calculating time in activity..."))
    
    timeActivity<-calculateTimeInActivity_BLK(dataSub,  irmaData)
    
  }  
  
  if (species=="Northern fulmar") {

    timeActivity<-calculateTimeInActivity_NF(dataSub, irmaData)  

  } 
  
  if (species=="Common guillemot") {
   
    timeActivity<-calculateTimeInActivity_CoGu(dataSub,  irmaData)  
    
  } 
  
  if (species=="Brünnich's guillemot") {
    
    timeActivity<-calculateTimeInActivity_BrGu(dataSub, irmaData)  
    
  }   
  
  if (species=="Little auk") {
    
    timeActivity<-calculateTimeInActivity_LiA(dataSub,  irmaData)  
    
  }
  
  if (species=="Atlantic puffin") {
    
    timeActivity<-calculateTimeInActivity_AP(dataSub,  irmaData)  
    
  }  
 
if (session>1) {
 
timeActivity1<-timeActivity[[1]] # Open current dataset
timeActivity1$session_id<-sessionNo[session] 
timeActivity1_dataset<-timeActivity_sessions[[1]] # Open what is already inside
timeActivity1_all<-rbind(timeActivity1_dataset, timeActivity1)
timeActivity_sessions[[1]]<-timeActivity1_all

timeActivity2<-timeActivity[[2]] # Open current dataset
timeActivity2$session_id<-sessionNo[session] 
timeActivity2_dataset<-timeActivity_sessions[[2]] # Open what is already inside
timeActivity2_all<-rbind(timeActivity2_dataset, timeActivity2)
timeActivity_sessions[[2]]<-timeActivity2_all

} else {

timeActivity1<-timeActivity[[1]] # Open current dataset
timeActivity1$session_id<-sessionNo[session] 
timeActivity_sessions[[1]]<-timeActivity1

timeActivity2<-timeActivity[[2]] # Open current dataset
timeActivity2$session_id<-sessionNo[session] 
timeActivity_sessions[[2]]<-timeActivity2

}
 
} 
  
  return(timeActivity_sessions)  
  
}

##### Black-legged Kittiwake #####

calculateTimeInActivity_BLK<-function(data, irmaData){
  
# PURPOSE: Classify time-series observations into behavioural states and calculate daily activity budgets and flight-bout statistics.

# Broad workflow:
#   1. Classify observations as Forage, RestWater, or Dry.
#   2. Identify continuous bouts of Dry observations.
#   3. Reclassify Dry bouts as Flight or Land according to bout duration.
#   4. Potentially reallocate the beginning of some Land bouts to Flight.
#   5. Recalculate final Flight bouts and check their maximum duration.
#   6. Calculate daily and darkness-period activity summaries.
  
# RETURNS

# A list containing:
#   [[1]] actResults  - daily activity/environment summaries
#   [[2]] boutResults - darkness-period flight summaries

# 1: INITIAL ACTIVITY CLASSIFICATION 
    
# Start by by assigning each observation to three behaviours: RestWater, Forage and Dry according to Th1 and Th2
  
dataCalc<-data %>%
  rename(new_cond=new.cond) %>% # Standadize the conductvitiy variable name
  dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
  dplyr::mutate(Activity=ifelse(new_cond<Th1, "Forage", "RestWater")) %>%
  dplyr::mutate(Activity=ifelse(new_cond<=Th2, "Dry", Activity)) %>%
  dplyr::mutate(Activity=ifelse(Th2==0 & new_cond==0, "Dry", Activity)) %>%
  dplyr::mutate(MaxDistColKm=max(distColonyKm)) # Save maximum distance reached from the colony

# 2: IDENTIFY THE START OF EACH DRY BOUT

# Consecutive dry observations are grouped into numbered bouts based on the time lag between sequential dry observations.
# If the time lag is more than 10 minutes, then a new numbered group is created. 
# The dataset is then subset to the first reading of each group. 
  
FlightBouts<-dataCalc %>%
  dplyr::filter(Activity=="Dry") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  dplyr::mutate(date_time=as.character(date_time))
	
# Make sure no NAs in the date time 8this was happening sometimes with changing the class of the date.time columns
nas_date<-subset(FlightBouts, is.na(date_time))
 
 if(nrow(nas_date)>0) {stop(print("Error: nas in date time")) }
  
# 3: PROPAGATE BOUT NUMBERS & CALCULATE THE DURATION OF EACH 'DRY BOUT'

# Join the identified bout starts back onto the complete dataset, propagate
# the bout number through the corresponding observations, and calculate
# the duration of each Dry bout.

FlightBoutLengths<-dataCalc %>%
  dplyr::mutate(date_time=as.character(date_time)) %>%
  dplyr::left_join(FlightBouts, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Dry", flightLengthMins, 0)) %>%
  ungroup()
  
# 4. CLASSIFY DRY BOUTS AS FLIGHT OR LAND

# For kittiwakes this is based on a dry bout being longer than L1 regardless of time of day
# or geographical position as assume they can roost on different structures at-sea/on land

activityAdjust1<-FlightBoutLengths %>%
  ungroup() %>%
  arrange(date_time) %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(NewActivity=ifelse(Activity=="Dry" & flightLengthMins > L1, "Land", NA)) %>%
  dplyr::mutate(NewActivity=ifelse(Activity=="Dry" & flightLengthMins <= L1, "Flight", NewActivity)) %>%
  dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity))
  
# Check for remaining dry bouts & stop if there are some as error

dryBouts<-subset(activityAdjust1, Activity=="Dry")
  
if (nrow(dryBouts)>0) {
    stop(print("Error: remaining dry bouts"))
  }

# 5. REALLOCATE THE START OF LAND BOUTS TO FLIGHT

# Long Dry bouts were classified as Land above. This section 
# accounts for the possibility that some time immediately before a landing
# was actually spent flying.
  
# Set possible amount of  values for L1 colony (this is to allow for a slightly different analysis for the sensitivity part)

uniqueValues<-unique(c(data$L1_colony_max[1], data$L1_colony_min[1]))

if (length(uniqueValues)>1){

# Possibility # 1: analysis conducted in main text  
   
# The start of land bouts can be re-allocated to flight
  
activityAdjust2_reallocate<-activityAdjust1 %>%
  dplyr::select(-NewActivity) %>%
  ungroup() %>%
  dplyr::mutate(firstLand=ifelse(Activity=="Land" & !lag(Activity)=="Land", 1, 0)) %>%  # Determine whether it's the first ten-minutes of a 'Land' bout
  dplyr::mutate(LandBoutNo=cumsum(firstLand)) %>% # Number the land bouts so we can do calculations by land bout number later on
  dplyr::mutate(LandBoutNo=ifelse(Activity=="Land", LandBoutNo, NA)) %>% # this just turns the number of all non-land bouts to NA
  dplyr::group_by(LandBoutNo) %>%
  dplyr::mutate(DurationLandMins=ifelse(Activity=="Land", n_distinct(date_time)*10, NA)) %>% # Determine duration of every land bout
  dplyr::ungroup() %>%
  dplyr::mutate(PrevFlight=ifelse(firstLand==1 & lag(Activity)=="Flight", 1, 0)) %>% # Determine whether the previous 10-mins was flight or not (re-allocation only occurs if previous 10-mins was something else)
  dplyr::mutate(LagMins=ifelse(Activity=="Land" & !PrevFlight %in% c(1) & firstLand==1, sample(c(seq(data$L1_colony_min[1], data$L1_colony_max[1], 10)), replace=TRUE), 0)) %>% # Determine a random number of 10-minute bouts to re-allocate from land to flight for every bout No
  dplyr::mutate(LagMinsAdj=ifelse(LagMins>=DurationLandMins, DurationLandMins-10, LagMins)) %>% # Make sure they are not longer than the actual land bout duration
  dplyr::group_by(LandBoutNo)  %>%
  dplyr::mutate(LandBoutRow=ifelse(Activity=="Land", row_number()*10, 0)) %>% # Annotate increasing duration of an individual land bout in minutes 
  replace_na((list(LagMinsAdj=0))) %>%
  dplyr::mutate(NewActivity=ifelse(Activity=="Land" & first(LagMinsAdj)>0 & LandBoutRow<=first(LagMinsAdj) & !is.na(LandBoutRow) & first(PrevFlight) %in% c(0), "Flight", NA)) %>%  # Re-allocate 10-minute segments of a land bout to flight if they are less in duration to 'LagMinsAdj' and are not preceded by flight
  dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity)) 
	
	} else {
	
# Possibility # 2: Sensitivity analysis (basically a pre-determined number of 10-minutes are changed for every land bout)
	  
# The start of land bouts can be re-allocated to flight

activityAdjust2_reallocate<-activityAdjust1 %>%
  dplyr::select(-NewActivity) %>%
  ungroup() %>%
  dplyr::mutate(firstLand=ifelse(Activity=="Land" & !lag(Activity)=="Land", 1, 0)) %>% # Determine whether it's the first ten-minutes of a 'Land' bout
  dplyr::mutate(LandBoutNo=cumsum(firstLand)) %>% # Now i number the land bouts so I get do some calculations by bout No later
  dplyr::mutate(LandBoutNo=ifelse(Activity=="Land", LandBoutNo, NA)) %>% # this just turns the number of all non-land bouts to NA
  dplyr::group_by(LandBoutNo) %>%
  dplyr::mutate(DurationLandMins=ifelse(Activity=="Land", n_distinct(date_time)*10, NA)) %>% # Determine duration of evey land bout
  dplyr::ungroup() %>%
  dplyr::mutate(PrevFlight=ifelse(firstLand==1 & lag(Activity)=="Flight", 1, 0)) %>% # Here I determine whether the previous bout was flight or not
  dplyr::mutate(LagMins=ifelse(Activity=="Land" & !PrevFlight %in% c(1) & firstLand==1, data$L1_colony_min[1], 0)) %>% # Here I determine a random number of 10-minute bouts to re-allocate from land to flight 
  dplyr::mutate(LagMinsAdj=ifelse(LagMins>=DurationLandMins, DurationLandMins-10, LagMins)) %>% # and here I make sure they are not longer than the actual land bout
  dplyr::group_by(LandBoutNo) %>%
  dplyr::mutate(LandBoutRow=ifelse(Activity=="Land", row_number()*10, 0)) %>% # Here i make a crazy system to re-allocate a certain number of rows...
  replace_na((list(LagMinsAdj=0))) %>%
  dplyr::mutate(NewActivity=ifelse(Activity=="Land" & first(LagMinsAdj)>0 & LandBoutRow<=first(LagMinsAdj) & !is.na(LandBoutRow) & first(PrevFlight) %in% c(0), "Flight", NA)) %>%
  dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity))
	
	}
  
# 6. REBUILD FLIGHT BOUTS AFTER REALLOCATION

# Some observations previously classified as Land may now be Flight.
# Therefore Flight bouts and their durations need to be recalculated from
# scratch just to make sure they are not too long. 
  
FlightBouts_2<-activityAdjust2_reallocate %>%
  dplyr::filter(Activity=="Flight") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  dplyr::mutate(date_time=as.character(date_time))
	
# Make sure no NAs in the date time
nas_date<-subset(FlightBouts_2, is.na(date_time))
 
if(nrow(nas_date)>0) {stop(print("Error: nas in date time")) }
  
# 7. CALCULATE FINAL FLIGHT-BOUT DURATIONS

FlightBoutLengths_final<-activityAdjust2_reallocate %>%
 ungroup() %>%
 dplyr::select(-BoutNo) %>%
 dplyr::left_join(FlightBouts_2, by=c("date_time")) %>%
 dplyr::group_by(Activity) %>%
 fill(BoutNo, .direction=c("down")) %>%
 dplyr::ungroup() %>%
 dplyr::group_by(BoutNo, Activity) %>%
 dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
 dplyr::mutate(flightLengthMins=ifelse(Activity=="Flight", flightLengthMins, 0)) %>%
 ungroup()
	
# Check no residual error
maxFlight<-max(FlightBoutLengths_final$flightLengthMins)
	
if(maxFlight > data$L1[1]) {
	stop(print("Error: flight bouts too long"))
}

# 8. CALCULATE DARKNESS-PERIOD FLIGHT STATISTICS 

# Produce a daily summary focused on flight during darkness (this is used in the supplementary analysis).
# Twilight is first merged into Daylight, leaving two effective periods:
# Daylight and Darkness.
  
dataCalcDay_period<-FlightBoutLengths_final %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::mutate(Period=ifelse(Period %in% c("Daylight", "Twilight"), "Daylight", "Darkness")) %>%
  dplyr::group_by(date, Period) %>%
  dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
  ungroup() %>%
  dplyr::group_by(species, colony, session_id, date, Period, Activity, BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(!Activity %in% c("Flight"), 0, flightLengthMins))%>%
  ungroup() %>%
  dplyr::group_by(species, colony, session_id, date, Period) %>%
  dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0)) %>%
  dplyr::summarise(tForage=sum(Forage)*10/60, tFlight=sum(Flight)*10/60, tRestWater=sum(RestWater)*10/60, tLand=sum(Land)*10/60,
                   tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, maxFlightBoutsMins_dark=max(flightLengthMins)) %>%
  dplyr::mutate(propDay_forage=tForage/tDaylight, propflight_dark=tFlight/tDarkness, flightTimeMins_dark=tFlight*60)  %>%
  ungroup() %>%
  dplyr::filter(Period=="Darkness") %>%
  dplyr::select(date, propflight_dark, maxFlightBoutsMins_dark, flightTimeMins_dark)
  
# Attach max flight bout length

daily_max_flightBout<-FlightBoutLengths_final %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::filter(Activity=="Flight") %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(maxFlightBoutsMins=max(flightLengthMins))
  
# 10. CALCULATE DAILY ACTIVITY BUDGETS

# Summarise hours spent in each activity and light period, together with
# environmental and spatial variables.

dataCalcDay<-FlightBoutLengths_final %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>% # Determine date
  dplyr::group_by(date) %>%
  dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
  dplyr::filter(Duration>=1430) %>% # Only keep full days %>%
  dplyr::group_by(species, colony, session_id, date) %>%
  dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                  Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0), sstRandom=sst_random_start) %>%
  dplyr::mutate(sstRandom=ifelse(sstRandom < -1.9, -1.9, sstRandom)) %>%
  dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tLand=sum(Land)*10/60,
                     tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, MaxDistColKm=max(MaxDistColKm),
                     sst_random=mean(sstRandom), ice_random=mean(ice_random, na.rm=TRUE), air_random=mean(air_mean, na.rm=TRUE), immersionType=mean(max.cond), distColonyKm_mean=mean(distColonyKm, na.rm=TRUE), mean.lon=mean(lon), mean.lat=mean(lat)) %>%
  ungroup() %>%
  dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
  dplyr::mutate(dayLengthHrs=tDaylight) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(DurationTot=sum(tForage, tRestWater, tLand, tFlight)) %>%
  dplyr::left_join(daily_max_flightBout, by=c("date")) %>%
  dplyr::left_join(dataCalcDay_period, by=c("date"))
  
# 11. PREPARE OUTPUT

# Main daily activity/environment results.
actResults<-dataCalcDay

# Darkness-period Flight-bout results.
boutResults<-dataCalcDay_period

# Return both outputs as a two-element list.  
allResults<-list(actResults, boutResults)

return(allResults)
  
}

##### Northern fulmar #####

calculateTimeInActivity_NF<-function(data, irmaData){

# PURPOSE: Classify time-series observations into behavioural states and calculate daily activity budgets and flight-bout statistics.
  
# Broad workflow:
# 1. Classify observations as Forage, RestWater, or Dry.
# 2. Calculate whether birds could potentially reach land based on their current and next locations.
# 3. Identify continuous bouts of Dry observations and calculate their duration.
# 4. Reclassify Dry bouts as Flight or Land according to bout duration, geographical position and probability of landing.
# 5. Potentially reallocate the beginning of some Land bouts to Flight.
# 6. Recalculate final Flight bouts and reclassify unrealistically long flights where landing is possible.
# 7. Calculate daily and daylight/darkness-period activity summaries.
# 8. Extract individual Flight bouts for supplementary analyses.
  
# RETURNS
  
# A list containing:
#   [[1]] actResults  - daily activity/environment summaries
#   [[2]] boutResults - activity summaries separated into daylight and darkness
  
  
# 1: INITIAL ACTIVITY CLASSIFICATION  
  
# First standardise the date-time column and assign each observation to three preliminary
# behaviours: RestWater, Forage and Dry according to conductivity thresholds Th1 and Th2.  
    
   
data$date_time<-as.character(data$date_time)
  
# Here we assign three behaviours: flight, forage & rest   
dataCalc<-data %>%
  dplyr::ungroup() %>%
  rename(new_cond=new.cond) %>%
  dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
  dplyr::mutate(Activity=ifelse(new_cond<Th1, "Forage", "RestWater")) %>%
  dplyr::mutate(Activity=ifelse(new_cond<Th2, "Dry", Activity)) %>%
  dplyr::mutate(Activity=ifelse(Th2==0 & new_cond==0, "Dry", Activity)) %>%
  dplyr::mutate(MaxDistColKm=max(distColonyKm)) 

# 2: CALCULATE DISTANCE TO COLONY AT THE NEXT LOCATION

# Whether a dry bout could represent time on land partly depends on whether the bird is close
# enough to the colony. Here we determine the colony distance associated with the next location,
# as a bird may be approaching the colony even if its current location is further away.
  
distances_next<-dataCalc %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  dplyr::group_by(distColonyKm) %>%
  dplyr::slice(1) %>%
  arrange(date_time) %>%
  dplyr::select(date_time, distColonyKm) %>%
  ungroup() %>%
  dplyr::mutate(distColonyKm_next=lead(distColonyKm)) %>%
  dplyr::mutate(distColonyKm_next=ifelse(is.na(distColonyKm_next), distColonyKm, distColonyKm_next)) %>%
  dplyr::select(-distColonyKm)
  
# 3: IDENTIFY THE START OF EACH DRY BOUT

# Consecutive dry observations are grouped into numbered bouts based on the time lag between sequential dry observations.
# If the time lag is more than 10 minutes, then a new numbered group is created.
# The dataset is then subset to the first reading of each group.

FlightBouts<-dataCalc %>%
  dplyr::filter(Activity=="Dry") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  dplyr::mutate(date_time=as.character(date_time))
	
# Check no nas in date.time 
na_dates<-subset(FlightBouts, is.na(date_time))	

if (nrow(na_dates)>1) {
stop(print("Error: na in dates"))
}

# 4: PROPAGATE BOUT NUMBERS & CALCULATE THE DURATION OF EACH DRY BOUT

# Join the identified bout starts back onto the complete dataset, propagate
# the bout number through the corresponding observations, and calculate
# the duration of each Dry bout.
  
FlightBoutLengths<-dataCalc %>%
  dplyr::left_join(FlightBouts, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo, distColonyKm) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Dry", flightLengthMins, NA)) %>%
  ungroup()

# 4: PROPAGATE BOUT NUMBERS & CALCULATE THE DURATION OF EACH DRY BOUT

# Join the identified bout starts back onto the complete dataset, propagate
# the bout number through the corresponding observations, and calculate
# the duration of each Dry bout.

FlightBoutLengths<-dataCalc %>%
  #dplyr::mutate(date_time=as.character(date_time)) %>%
  #dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  dplyr::left_join(FlightBouts, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo, distColonyKm) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>% # Calculate bout duration assuming observations represent 10-minute intervals
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Dry", flightLengthMins, NA)) %>% # Only retain calculated duration for Dry observations
  ungroup()


# 5: CLASSIFY DRY BOUTS AS FLIGHT OR LAND

# Unlike the BLK function, classification here depends on both the duration of the dry bout
# and whether the bird could plausibly be on land or sea ice.
#
# PossLand is set to 1 when there is sea ice, when the bird is within dist_colony of the colony,
# or when its next location is within dist_colony.
#
# For short bouts where landing is possible, a probabilistic allocation determines whether
# the observation is classified as Land or Flight. This combines pLand_prob with ice concentration.
  
activityAdjust1<-FlightBoutLengths %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  dplyr::full_join(distances_next, by=c("date_time")) %>%
  arrange(date_time) %>%
  fill(distColonyKm_next, .direction=c("down")) %>%
  dplyr::group_by(BoutNo, distColonyKm) %>%
  # make sure ice_random is between 0 and 1
  dplyr::mutate(ice_random=ifelse(!is.na(BoutNo), first(ice_random), NA)) %>%
  dplyr::mutate(ice_random=ifelse(ice_random<0, 0, ice_random)) %>%
  dplyr::mutate(ice_random=ifelse(ice_random>1, 1, ice_random)) %>%
  # Determine whether is was possible for a bird to be on land or ice
  dplyr::mutate(PossLand=ifelse(ice_random >0 | distColonyKm <= dist_colony | distColonyKm_next <=dist_colony, 1, 0)) %>%
  dplyr::mutate(p=runif(n=n_distinct(BoutNo))) %>% # Determine random probability for each bout no (used for land)
  dplyr::mutate(p2=runif(n=n_distinct(BoutNo))) %>% # Determine random probability for each bout no (used for sea-ice)
  # Create a probability for being on land or ice based on this for every bout
  dplyr::mutate(pLand=ifelse(p <= pLand_prob | p2 <= ice_random, 1, 0)) %>%
  replace_na(list(flightLengthMins=0)) %>%
  dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins >= L1 & PossLand==0, "Flight", NA)) %>% # If no land in sight, then it has to be flight
  dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins >= L1 & PossLand==1, "Land", NewActivity)) %>% # otherwise it's land
  dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==0, "Flight", NewActivity)) %>% # If it's short & no land then it's flight
  dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm <= dist_colony & pLand==1, "Land", NewActivity)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
  dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm <= dist_colony & pLand==0 , "Flight", NewActivity)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
  dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm > dist_colony & pLand==1, "Land", NewActivity)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
  dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm > dist_colony & pLand==0 , "Flight", NewActivity)) %>%
  # Create new column where I count how many rows this concerns (for supplementary analysis)
  dplyr::mutate(RandomAllocation=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm <= dist_colony & pLand==1, 1, 0)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
  dplyr::mutate(RandomAllocation=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm <= dist_colony & pLand==0 , 1, RandomAllocation)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
  dplyr::mutate(RandomAllocation=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm > dist_colony & pLand==1, 1, RandomAllocation)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
  dplyr::mutate(RandomAllocation=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm > dist_colony & pLand==0 , 1, RandomAllocation)) %>%
  dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity)) %>%
  ungroup() %>%
  dplyr::group_by(Activity) %>%
  dplyr::mutate(maxLength=ifelse(Activity %in% c("Flight", "Land"), max(flightLengthMins, na.rm=TRUE), NA))
  
# Check for remaining dry bouts & stop if there are some as error
dryBouts<-subset(activityAdjust1, Activity=="Dry")
  
if (nrow(dryBouts)>0) {
  stop(print("Error: remaining dry bouts"))
}
  
# 6: REALLOCATE THE START OF LAND BOUTS TO FLIGHT
  
# Long Dry bouts may have been classified as Land above. As in the BLK function,
# this section accounts for the possibility that some time immediately before
# landing was actually spent flying.
  
# Determine whether L1_colony is being sampled across a range of values
# or is fixed to a single value for the sensitivity analysis
  
uniqueVals<-unique(c(data$L1_colony_min[1], data$L1_colony_max[1]))
  
if (length(uniqueVals)>1) {
  
# Possibility #1: analysis conducted in main text
  
# The start of Land bouts can be reallocated to Flight using a randomly selected
# duration between L1_colony_min and L1_colony_max.
	
activityAdjust2_reallocate<-activityAdjust1 %>%
  dplyr::select(-NewActivity) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(firstLand=ifelse(Activity=="Land" & !lag(Activity)=="Land", 1, 0)) %>% # Determine whether it's the first ten-minutes of a 'Land' bout
  replace_na(list(firstLand=0)) %>%
  dplyr::mutate(LandBoutNo=cumsum(firstLand)) %>% # Number the land bouts to do some calculations by bout No later
  dplyr::mutate(LandBoutNo=ifelse(Activity=="Land", LandBoutNo, NA)) %>% # Change the number of all non-land bouts to NA
  dplyr::group_by(LandBoutNo) %>%
  dplyr::mutate(DurationLandMins=ifelse(Activity=="Land", n_distinct(date_time)*10, NA)) %>% # Determine duration of evey land bout
  dplyr::ungroup() %>%
  dplyr::mutate(PrevFlight=ifelse(firstLand==1 & lag(Activity)=="Flight", 1, 0)) %>% # Determine whether the previous bout was flight or not
  dplyr::mutate(LagMins=ifelse(Activity=="Land" & !PrevFlight %in% c(1) & firstLand==1, sample(seq(data$L1_colony_min[1], data$L1_colony_max[1], 10), replace=TRUE), 0)) %>% # Determine a random number of 10-minute bouts to re-allocate from land to flight 
  dplyr::mutate(LagMinsAdj=ifelse(LagMins>=DurationLandMins, DurationLandMins-10, LagMins)) %>% # Make sure this number is not longer than the actual land bout duration
  dplyr::group_by(LandBoutNo) %>%
  dplyr::mutate(LandBoutRow=ifelse(Activity=="Land", row_number()*10, 0)) %>% # Annotate increasing duration within each Land bout
  replace_na((list(LagMinsAdj=0))) %>%
  dplyr::mutate(NewActivity=ifelse(Activity=="Land" & first(LagMinsAdj)>0 & LandBoutRow<=first(LagMinsAdj) & !is.na(LandBoutRow) & first(PrevFlight) %in% c(0), "Flight", NA)) %>% #  Reallocate the beginning of eligible Land bouts to Flight
  dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity))
	
	
	} else {
	
# Possibility #2: sensitivity analysis
	  
# Here a predetermined value of L1_colony_min is used rather than randomly
# selecting a value between L1_colony_min and L1_colony_max.
	
activityAdjust2_reallocate<-activityAdjust1 %>%
  dplyr::select(-NewActivity) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(firstLand=ifelse(Activity=="Land" & !lag(Activity)=="Land", 1, 0)) %>% # Determine whether it's the first ten-minutes of a 'Land' bout
  replace_na(list(firstLand=0)) %>%
  dplyr::mutate(LandBoutNo=cumsum(firstLand)) %>% # Now i number the land bouts so I get do some calculations by bout No later
  dplyr::mutate(LandBoutNo=ifelse(Activity=="Land", LandBoutNo, NA)) %>% # this just turns the number of all non-land bouts to NA
  dplyr::group_by(LandBoutNo) %>%
  dplyr::mutate(DurationLandMins=ifelse(Activity=="Land", n_distinct(date_time)*10, NA)) %>% # Determine duration of evey land bout
  dplyr::ungroup() %>%
  dplyr::mutate(PrevFlight=ifelse(firstLand==1 & lag(Activity)=="Flight", 1, 0)) %>% # Here I determine whether the previous bout was flight or not
  dplyr::mutate(LagMins=ifelse(Activity=="Land" & !PrevFlight %in% c(1) & firstLand==1, data$L1_colony_min[1], 0)) %>% # Here I determine a random number of 10-minute bouts to re-allocate from land to flight 
  dplyr::mutate(LagMinsAdj=ifelse(LagMins>=DurationLandMins, DurationLandMins-10, LagMins)) %>% # and here I make sure they are not longer than the actual land bout
  dplyr::group_by(LandBoutNo) %>%
  dplyr::mutate(LandBoutRow=ifelse(Activity=="Land", row_number()*10, 0)) %>% # Here i make a crazy system to re-allocate a certain number of rows...
  replace_na((list(LagMinsAdj=0))) %>%
  dplyr::mutate(NewActivity=ifelse(Activity=="Land" & first(LagMinsAdj)>0 & LandBoutRow<=first(LagMinsAdj) & !is.na(LandBoutRow) & first(PrevFlight) %in% c(0), "Flight", NA)) %>%
  dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity))
	
	
	}
	
# 7: REBUILD FLIGHT BOUTS AFTER REALLOCATION

# Some observations previously classified as Land may now be Flight.
# Therefore Flight bouts and their durations need to be recalculated from scratch.

FlightBouts_final<-activityAdjust2_reallocate %>%
  dplyr::filter(Activity=="Flight") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >11) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  dplyr::mutate(date_time=as.character(date_time))
	
# 8: CALCULATE FINAL FLIGHT-BOUT DURATIONS & RECLASSIFY LONG FLIGHTS

# Recalculate Flight-bout duration without splitting bouts according to distance from the colony.
# Following reallocation, any Flight bout longer than L1 is changed back to Land where landing
# is geographically possible.
	
FlightBoutLengths_final<-activityAdjust2_reallocate %>%
  dplyr::mutate(date_time=as.character(date_time)) %>%
  dplyr::select(-BoutNo) %>%
  dplyr::left_join(FlightBouts_final, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Flight", flightLengthMins, NA)) %>%
  ungroup() %>%
  dplyr::mutate(Activity=ifelse(Activity=="Flight" & PossLand==1 & flightLengthMins > L1, "Land", Activity)) %>%
  ungroup()
	
# 9: CALCULATE DARKNESS-PERIOD FLIGHT STATISTICS

# Produce a daily summary focused on Flight during darkness (used in supplementary analyses).
# Twilight is merged into Daylight, leaving two effective periods: Daylight and Darkness.

dataCalcDay_period<-FlightBoutLengths_final %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::mutate(Period=ifelse(Period %in% c("Daylight", "Twilight"), "Daylight", "Darkness")) %>%
  dplyr::group_by(date, Period) %>%
  dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
  ungroup() %>%
  dplyr::group_by(species, colony, session_id, date, Period, Activity, BoutNo, distColonyKm) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(!Activity %in% c("Flight"), 0, flightLengthMins))%>%
  ungroup() %>%
  dplyr::group_by(species, colony, session_id, date, Period) %>%
  dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0)) %>%
  dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tRestWater=sum(RestWater)*10/60, tLand=sum(Land)*10/60,
                   tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, maxFlightBoutsMins_dark=max(flightLengthMins)) %>%
  dplyr::mutate(propDay_forage=tForage/tDaylight, propflight_dark=tFlight/tDarkness, flightTimeMins_dark=tFlight*60)  %>%
  ungroup() %>%
  dplyr::filter(Period=="Darkness") %>%
  dplyr::select(date, propflight_dark, maxFlightBoutsMins_dark, flightTimeMins_dark)
	  
# 10: CALCULATE MAXIMUM DAILY FLIGHT-BOUT LENGTH

daily_max_flightBout<-FlightBoutLengths_final %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::filter(Activity=="Flight") %>%
  dplyr::group_by(date, BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  ungroup() %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(maxFlightBoutsMins=max(flightLengthMins))

# Calculate a corresponding darkness-specific flight metric

daily_max_flightBout_dark<- dataCalcDay_period%>%
  dplyr::ungroup() %>%
  #dplyr::filter(Activity=="Flight") %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(maxFlightBoutsMins_dark=max(flightTimeMins_dark)) %>%
  ungroup() 

# 11: CALCULATE DAILY ACTIVITY BUDGETS

# Summarise hours spent in each activity and light period, together with
# environmental and spatial variables. Only essentially complete days
# (at least 1430 minutes of observations) are retained.
	  
dataCalcDay<-FlightBoutLengths_final %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
  dplyr::filter(Duration>=1430) %>%  # Only keep full days %>%
  dplyr::group_by(species, colony, session_id, date) %>%
  dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0), sstRandom=sst_random_start) %>%
  dplyr::mutate(sstRandom=ifelse(sstRandom < -1.9, -1.9, sstRandom)) %>%
  dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tRestWater=sum(RestWater)*10/60, tLand=sum(Land)*10/60,
                   tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, MaxDistColKm=max(MaxDistColKm),
                   sst_random=mean(sstRandom), ice_random=mean(ice_random, na.rm=TRUE), air_random=mean(air_mean, na.rm=TRUE), immersionType=mean(max.cond), distColonyKm_mean=mean(distColonyKm, na.rm=TRUE), mean.lon=mean(lon),
                   mean.lat=mean(lat), boutsRandom=(sum(RandomAllocation)/144)) %>%
  ungroup() %>%
  dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
  dplyr::mutate(dayLengthHrs=tDaylight) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(DurationTot=sum(tForage, tRestWater, tLand, tFlight)) %>%
  dplyr::left_join(daily_max_flightBout, by=c("date")) %>%
  dplyr::left_join(dataCalcDay_period, by=c("date")) 
	
  
# 12: CALCULATE ACTIVITY BUDGETS SEPARATELY FOR DAYLIGHT AND DARKNESS

# Create a second summary dataset containing time spent in each behaviour
# separately during Daylight and Darkness. Twilight is included within Daylight.

dataCalcDay_period2<-FlightBoutLengths_final %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::mutate(Period=ifelse(Period %in% c("Daylight", "Twilight"), "Daylight", "Darkness")) %>%
  dplyr::group_by(date, Period) %>%
  dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
  dplyr::group_by(species, colony, session_id, date,  Period) %>%
  dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0)) %>%
  dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tRestWater=sum(RestWater)*10/60, tLand=sum(Land)*10/60,
                   tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10) 
  
# 13: EXTRACT INDIVIDUAL FLIGHT BOUTS FOR SUPPLEMENTARY ANALYSIS

# Reconstruct individual Flight bouts for fulmars so that the distribution
# of Flight-bout durations across all times of day can be saved and analysed.

finalFlightBoutNos<-FlightBoutLengths_final %>%
  dplyr::filter(Activity=="Flight") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
	dplyr::mutate(date_time=as.character(date_time))
	
# Join the new bout numbers back onto Flight observations and calculate
# the final duration of every individual Flight bout.

finalFlightBoutLengths<-FlightBoutLengths_final %>%
  dplyr::select(-BoutNo) %>%
  dplyr::filter(Activity=="Flight") %>%
  dplyr::ungroup() %>%
  dplyr::left_join(finalFlightBoutNos, by=c("date_time")) %>%
  dplyr::group_by(species, colony, individ_id) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10)
	
# Save results
#write.csv(finalFlightBoutLengths, file=paste0("./results/tables/supplementary/fulmarBoutLengths/fulmar_flightbouts_", finalFlightBoutLengths$individ_id[1], "_rep", i, ".csv"))
  
# 14: PREPARE OUTPUT

# Main daily activity/environment results.
actResults<-dataCalcDay

# Daylight/darkness-specific activity results.
boutResults<-dataCalcDay_period2

# Return both outputs as a two-element list.
allResults<-list(actResults, boutResults)

return(allResults)
  
}

##### Methods for activity budgets - auks ####

methodCaitlin<-function(data, irmaData, speciesLatin) {

# Broad workflow:
# 1. Initially classify observations as Dry or Other according to conductivity.
# 2. Identify continuous Dry bouts and calculate their duration.
# 3. During darkness, determine whether Dry observations represent Land or RestWater.
# 4. During daylight, allocate remaining Dry observations to Flight, Land or RestWater.
# 5. Reclassify remaining observations as Active or RestWater according to conductivity.
# 6. Reallocate the beginning of some Land bouts to Flight.
# 7. Repeatedly check Flight-bout duration and correct bouts exceeding the maximum threshold.
# 8. Calculate daily activity/environment summaries and daylight/darkness activity budgets.
  
# RETURNS
  
# A list containing:
#   [[1]] actResults  - daily activity/environment summaries
#   [[2]] boutResults - activity summaries separated into daylight and darkness
  
  
# 1: INITIAL DATA PREPARATION
  
# Calculate day of year, maximum distance reached from the colony, and arrange
# observations chronologically.  
    
dataCalc<-data %>%
  dplyr::ungroup() %>%
  dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
  dplyr::mutate(MaxDistColKm=max(distColonyKm)) %>%
  arrange(date) %>%
  ungroup()  
  
# 2: CALCULATE DISTANCE TO COLONY AT THE NEXT LOCATION

# Whether a Dry observation could represent time on land partly depends on whether
# the bird is currently close to the colony or approaching it. Here we determine
# the colony distance associated with the next available location.

distances_next<-dataCalc %>%
  dplyr::group_by(distColonyKm) %>%
  dplyr::slice(1) %>%
  arrange(date_time) %>%
  dplyr::select(date_time, distColonyKm) %>%
  ungroup() %>%
  dplyr::mutate(distColonyKm_next=lead(distColonyKm)) %>%
  dplyr::mutate(distColonyKm_next=ifelse(is.na(distColonyKm_next), distColonyKm, distColonyKm_next)) %>%
  dplyr::select(-distColonyKm) %>%
  ungroup() %>%
  dplyr::mutate(date_time=as.character(date_time))
  
# 3: INITIAL DRY/OTHER CLASSIFICATION

# Initially separate observations into Dry and Other according to the lower
# conductivity threshold (Th2). Dry observations are subsequently divided
# into Flight, Land or RestWater.

dataCalc$Activity<-ifelse(dataCalc$new.cond<=data$Th2[1], "Dry", "Other")
  
# Now we allocate dry bouts to TFlight, TLand or TRest

# 4: IDENTIFY DRY BOUTS & CALCULATE THEIR DURATION

# Consecutive Dry observations are grouped into numbered bouts based on the time
# lag between observations. A gap greater than 10 minutes starts a new bout.
  
FlightBouts<-dataCalc %>%
  dplyr::filter(Activity=="Dry") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, session_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  distinct() %>%
  dplyr::mutate(date_time=as.character(date_time))
	
# Make sure no NAs in date_time
nas_date<-subset(FlightBouts, is.na(date_time))
  
if (nrow(nas_date)>0) {stop(print("Error: nas in date_time")) }

# Join bout numbers back onto the complete dataset and calculate the duration
# of each Dry bout assuming observations represent 10-minute intervals.

FlightBoutLengths<-dataCalc %>%
  dplyr::mutate(date_time=as.character(date_time)) %>%
  dplyr::left_join(FlightBouts, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Dry", flightLengthMins, 0)) %>%
  ungroup()
  
# 5: CLASSIFY DRY OBSERVATIONS DURING DARKNESS

# During darkness, Dry observations are assumed to represent either Land or
# RestWater rather than Flight. Land is only possible when the bird is dry for
# the entire night and is geographically able to be on land or sea ice.

# Sea ice is considered a possible resting surface only for Little auks and
# Brünnich's guillemots.
  
# Being on sea-ice is only applicable to little auks & Brunnich's guillemots 

FlightBoutLengths$use_seaice<-ifelse(dataCalc$species[1] %in% c("Little auk", "Brünnich's guillemot"), 1, 0)
  
# First identify individual periods of darkness ("nights"). 

darkness_boutNo<-FlightBoutLengths %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  dplyr::filter(Period=="Darkness")  %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo_dark=row_number()) %>%
  dplyr::select(date_time, BoutNo_dark) %>%
	dplyr::mutate(date_time=as.character(date_time))
  
# Determine the proportion of each darkness period that the bird was Dry.
# propDarkDry == 1 identifies nights during which the bird was dry throughout.

darkness_propDry<-dataCalc %>%
  dplyr::mutate(date_time=as.character(date_time)) %>%
  dplyr::left_join(darkness_boutNo, by=c("date_time")) %>%
  dplyr::group_by(Period) %>%
  fill(BoutNo_dark, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Period, BoutNo_dark) %>%
  dplyr::mutate(totDry=sum(1-new.cond)*10, darknessLength=n_distinct(date_time)*10) %>%
  dplyr::mutate(propDarkDry=ifelse(Period=="Darkness", totDry/darknessLength, NA)) %>%
  ungroup() %>%
  dplyr::select(date_time, BoutNo_dark, propDarkDry)
  
# Allocate Dry observations during darkness to Land or RestWater.
#
# PossLand is initially based on the presence of sea ice, current distance to
# the colony, or distance to the colony at the next location. It is then set
# to zero unless the bird remained Dry for the entire darkness period.
#
# Where landing is possible, pLand probabilistically determines whether the
# bird is classified as Land or RestWater.

activityAdjust1<-FlightBoutLengths %>%
  ungroup() %>%
  dplyr::left_join(darkness_propDry, by=c("date_time")) %>% # Join info on proportion of night dry
  dplyr::full_join(distances_next, by=c("date_time")) %>% # Join info on distance to colony
  arrange(date_time) %>%
  fill(distColonyKm_next, .direction=c("down")) %>%
  dplyr::group_by(BoutNo, BoutNo_dark, distColonyKm) %>%
  # Determine distance to land & ice concentration
  dplyr::mutate(ice_random=ifelse(!is.na(BoutNo), first(ice_random), NA)) %>%
  dplyr::mutate(ice_random=ifelse(!is.na(BoutNo) & ice_random<0, 0, ice_random)) %>%
  dplyr::mutate(ice_random=ifelse(!is.na(BoutNo) & ice_random>1, 1, ice_random)) %>%
  dplyr::mutate(ice_random=ifelse(use_seaice==1, ice_random, 0)) %>% # Make this null if species does not stand on ice
  # Determine whether the bird could have been on land or not (PossLand)
  dplyr::mutate(PossLand=ifelse( ice_random >0 | distColonyKm <= dist_colony | distColonyKm_next <= dist_colony, 1, 0)) %>%
  # Change this to zero if prop night that is dry is not equal to 1
  dplyr::mutate(PossLand=ifelse(propDarkDry==1, PossLand, 0)) %>%
  replace_na(list(PossLand=0)) %>%
  # Generate random probabilities & use these to determine whether a bird could have been on land or not
  dplyr::mutate(p=runif(n=n_distinct(BoutNo))) %>%
  dplyr::mutate(p2=runif(n=n_distinct(BoutNo))) %>%
  dplyr::mutate(pLand=ifelse(p <= pLand_prob | p2 <= ice_random, 1, 0)) %>%
  dplyr::mutate(NewActivity=ifelse(Period=="Darkness" & Activity=="Dry" & PossLand==0, "RestWater", NA)) %>%
  dplyr::mutate(NewActivity=ifelse(Period=="Darkness" & Activity=="Dry" & PossLand==1 & pLand==1, "Land", NewActivity)) %>%
  dplyr::mutate(NewActivity=ifelse(Period=="Darkness" & Activity=="Dry" & PossLand==1 & pLand==0, "RestWater", NewActivity)) %>%
  dplyr::mutate(RandomAllocation=ifelse(Period=="Darkness" & Activity=="Dry" & PossLand==1 & pLand==1, 1, 0)) %>%
  dplyr::mutate(RandomAllocation=ifelse(Period=="Darkness" & Activity=="Dry" & PossLand==1 & pLand==0, 1, RandomAllocation)) %>%
  dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity))
  
# 6: CLASSIFY REMAINING DRY OBSERVATIONS DURING DAYLIGHT

# Dry observations that remain after the darkness classification are re-numbered
# and their bout durations recalculated. These observations can subsequently be
# assigned to Flight, Land or RestWater.
  
FlightBouts2<-activityAdjust1 %>%
  dplyr::filter(Activity=="Dry") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  dplyr::mutate(date_time=as.character(date_time))
  
# Recalculate the duration of remaining Dry bouts.

FlightBoutLengths2<-activityAdjust1 %>%
  ungroup() %>%
  dplyr::select(-BoutNo) %>%
  dplyr::left_join(FlightBouts2, by=c("date_time")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo, distColonyKm) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Dry", flightLengthMins, 0)) %>%
  ungroup()

# Allocate daylight Dry observations according to bout duration, whether landing
# is possible, and the probabilistic landing parameter.
#
# Long bouts (>= L1):
#   - No possible landing -> RestWater
#   - Possible landing + pLand == 1 -> Land
#   - Possible landing + pLand == 0 -> RestWater
#
# Short bouts (< L1):
#   - No possible landing -> randomly Flight or RestWater
#   - Possible landing + pLand == 1 -> Land
#   - Possible landing + pLand == 0 -> randomly Flight or RestWater
  
activityAdjust2<-FlightBoutLengths2 %>%
  dplyr::group_by(BoutNo, distColonyKm) %>%
  # Determine distance to land & ice concentration
  dplyr::mutate(ice_random=ifelse(!is.na(BoutNo), ice_mean, NA)) %>%
  dplyr::mutate(ice_random=ifelse(ice_random<0, 0, ice_random)) %>%
  dplyr::mutate(ice_random=ifelse(ice_random>1, 1, ice_random)) %>%
  dplyr::mutate(ice_random=ifelse(use_seaice==1, ice_random, 0)) %>% # Change this to zero if species is incorrect
  # Determine whether it was possible for the bird to be on land or not (PossLand)
  dplyr::mutate(PossLand=ifelse(ice_random >0 | distColonyKm <= dist_colony | distColonyKm_next <=dist_colony, 1, 0)) %>%
  # Generate random probabilities which will be used to determine whether bird is on land or not (pLand)
  dplyr::mutate(p=runif(n=n_distinct(BoutNo))) %>%
  dplyr::mutate(p2=runif(n=n_distinct(BoutNo))) %>%
  dplyr::mutate(pLand=ifelse(p <= pLand_prob | p2 <= ice_random, 1, 0)) %>%
  dplyr::mutate(NewActivity=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins >= L1 & PossLand==0, "RestWater", NA)) %>% # If no land in sight, then it has to be flight
  dplyr::mutate(NewActivity=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins >= L1 & PossLand==1 & pLand == 1, "Land", NewActivity)) %>% # otherwise it's land
  dplyr::mutate(NewActivity=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins >= L1 & PossLand==1 & pLand == 0, "RestWater", NewActivity)) %>% # otherwise it's land
  dplyr::mutate(NewActivity=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins < L1 & PossLand==0, sample(c("Flight", "RestWater"), 1), NewActivity)) %>%
  dplyr::mutate(NewActivity=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & pLand==1, "Land", NewActivity)) %>%
  dplyr::mutate(NewActivity=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & pLand ==0, sample(c("Flight", "RestWater"), 1), NewActivity)) %>%
  # Sum how often these random allocations occurs
  dplyr::mutate(RandomAllocation=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins >= L1 & PossLand==1 & pLand == 1, 1, RandomAllocation)) %>% # otherwise it's land
  dplyr::mutate(RandomAllocation=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins >= L1 & PossLand==1 & pLand == 0, 1, RandomAllocation)) %>% # otherwise it's land
  dplyr::mutate(RandomAllocation=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins < L1 & PossLand==0, 1, RandomAllocation)) %>%
  dplyr::mutate(RandomAllocation=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & pLand==1, 1, RandomAllocation)) %>%
  dplyr::mutate(RandomAllocation=ifelse(!Period=="Darkness" & Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & pLand ==0, 1, RandomAllocation)) %>%
  dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity)) %>%
  ungroup() %>%
  dplyr::group_by(Activity) %>%
  dplyr::mutate(maxLength=ifelse(Activity %in% c("Flight", "Land"), max(flightLengthMins, na.rm=TRUE), NA)) %>%
    ungroup()
  
# Check for remaining dry bouts & stop if there are some as error
dryBouts<-subset(activityAdjust2, Activity=="Dry")
  
if (nrow(dryBouts)>0) {
  stop(print("Error: Dry bouts remain!"))
}
  
# 7: CLASSIFY REMAINING OBSERVATIONS AS ACTIVE OR RESTWATER

# Observations that were initially classified as Other are now divided according
# to Th1. Conductivity values between zero and Th1 become RestWater, while values
# greater than or equal to Th1 become Active.

dry_final<-activityAdjust2
dry_final$Activity<-ifelse(dry_final$new.cond<data$Th1[1] & dry_final$new.cond>0, "RestWater",dry_final$Activity)
dry_final$Activity<-ifelse(dry_final$new.cond>=data$Th1[1], "Active", dry_final$Activity)

# 8: REALLOCATE THE START OF LAND BOUTS TO FLIGHT

# Identify individual Land bouts. The beginning of a Land bout can subsequently
# be reallocated to Flight if the bird was not already flying immediately before
# the Land bout.
  
dry_lengths<-dry_final %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  dplyr::filter(Activity=="Land") %>%
  dplyr::mutate(gap=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list(gap=11)) %>%
  dplyr::filter(gap>10) %>%
  dplyr::mutate(DryBoutNo=row_number()) %>%
  dplyr::select(date_time, DryBoutNo) %>%
  dplyr::mutate(date_time=as.character(date_time))
  
# Determine whether L1_colony is being sampled across a range of values
# or fixed to one value for the sensitivity analysis.
  
uniqueValues<-unique(c(data$L1_colony_min, data$L1_colony_max))
  
if (length(uniqueValues)>1) {

# Possibility #1: main analysis
  
# Randomly select the amount of time to reallocate between L1_colony_min and
# L1_colony_max. Reallocation only occurs where a Land bout was not preceded
# by Flight, and at least the final 10 minutes remain classified as Land.  
    
reAllocate<-dry_final %>%
  dplyr::ungroup() %>%
  #dplyr::select(-DryBoutNo) %>%
  dplyr::left_join(dry_lengths, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(DryBoutNo, .direction=c("down")) %>%
  ungroup() %>%
  dplyr::mutate(DryBoutNo=ifelse(Activity=="Land", DryBoutNo, NA)) %>%
  dplyr::group_by(DryBoutNo) %>%
  dplyr::mutate(DurationLandMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(DurationLandMins=ifelse(Activity=="Land", DurationLandMins, NA)) %>%
  dplyr::mutate(OnLand=ifelse(Activity=="Land", 1, 0)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(FirstLand=ifelse(OnLand==1 & lag(OnLand)==0, 1, 0)) %>%
  dplyr::mutate(PrevAct=ifelse(OnLand==1 & lag(OnLand)==0, lag(Activity), "Other")) %>%
  dplyr::mutate(LagMins=ifelse(Activity=="Land" & !PrevAct %in% c("Flight") & FirstLand==1, sample(seq(data$L1_colony_min[1], data$L1_colony_max[1], 10), replace=TRUE), 0)) %>%
  dplyr::mutate(LagMinsAdj=ifelse(LagMins>=DurationLandMins, DurationLandMins-10, LagMins)) %>%
  ungroup() %>%
  dplyr::group_by(DryBoutNo) %>%
  dplyr::mutate(DryBoutRow=ifelse(Activity=="Land", row_number()*10, 0)) %>%
  replace_na((list(LagMinsAdj=0, PrevAct="Other"))) %>%
  dplyr::mutate(Activity=ifelse(Activity=="Land" & first(LagMinsAdj)>0 & DryBoutRow<=first(LagMinsAdj) & !is.na(DryBoutRow), "Flight", Activity))
	
	} else {
	
# Possibility #2: sensitivity analysis

# Use the fixed L1_colony_min value rather than randomly selecting the amount
# of time reallocated from Land to Flight.
	
reAllocate<-dry_final %>%
  dplyr::ungroup() %>%
  #dplyr::select(-DryBoutNo) %>%
  dplyr::left_join(dry_lengths, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(DryBoutNo, .direction=c("down")) %>%
  ungroup() %>%
  dplyr::mutate(DryBoutNo=ifelse(Activity=="Land", DryBoutNo, NA)) %>%
  dplyr::group_by(DryBoutNo) %>%
  dplyr::mutate(DurationLandMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(DurationLandMins=ifelse(Activity=="Land", DurationLandMins, NA)) %>%
  dplyr::mutate(OnLand=ifelse(Activity=="Land", 1, 0)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(FirstLand=ifelse(OnLand==1 & lag(OnLand)==0, 1, 0)) %>%
  dplyr::mutate(PrevAct=ifelse(OnLand==1 & lag(OnLand)==0, lag(Activity), "Other")) %>%
  dplyr::mutate(LagMins=ifelse(Activity=="Land" & !PrevAct %in% c("Flight") & FirstLand==1, data$L1_colony_min[1], 0)) %>%
  dplyr::mutate(LagMinsAdj=ifelse(LagMins>=DurationLandMins, DurationLandMins-10, LagMins)) %>%
  ungroup() %>%
  dplyr::group_by(DryBoutNo) %>%
  dplyr::mutate(DryBoutRow=ifelse(Activity=="Land", row_number()*10, 0)) %>%
  replace_na((list(LagMinsAdj=0, PrevAct="Other"))) %>%
  dplyr::mutate(Activity=ifelse(Activity=="Land" & first(LagMinsAdj)>0 & DryBoutRow<=first(LagMinsAdj) & !is.na(DryBoutRow), "Flight", Activity))
	
	}
  
# If there are activities which are 'other' or NA then STOP
nas<-subset(reAllocate, is.na(Activity))
other<-subset(reAllocate, Activity %in% c("Other"))
  
if(nrow(nas)>0  | nrow(other)>0) {
    
    stop(print("Error: nas in activity"))  
    
} 
  
# 9: CHECK & CORRECT FLIGHT BOUTS EXCEEDING L1

# Reallocation can create Flight bouts that exceed the maximum permitted Flight
# duration (L1). The following sections progressively correct these bouts.
#
# First, overlong Flight bouts are changed to Land where landing is possible
# during daylight. Remaining overlong Flight bouts are then shortened by
# reclassifying part of the bout as RestWater.
  
# 9a: IDENTIFY CURRENT FLIGHT BOUTS
  
FlightBouts_lastcheck<-reAllocate %>%
  dplyr::filter(Activity=="Flight") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  distinct() %>%
  dplyr::mutate(date_time=as.character(date_time))
	
# Make sure no NAs in date_time
nas_date<-subset(FlightBouts_lastcheck, is.na(date_time))
  
if (nrow(nas_date)>0) {stop(print("Error: nas in date_time")) }
  
# Join bout numbers back onto the dataset and calculate Flight-bout duration.

FlightBoutLengths_lastcheck<-reAllocate %>%
  dplyr::select(-BoutNo) %>%
  #dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  dplyr::left_join(FlightBouts_lastcheck, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Flight", flightLengthMins, 0)) %>%
  ungroup()
	
# 9b: RECLASSIFY OVERLONG FLIGHTS AS LAND WHERE POSSIBLE

# During daylight, overlong Flight bouts are reclassified as Land where the
# bird could plausibly have landed.

FlightBoutLengths_lastcheck_reclassify<-FlightBoutLengths_lastcheck %>%
  ungroup() %>%
  dplyr::mutate(Activity=ifelse(Activity=="Flight" & flightLengthMins > L1 & PossLand==1 & Period=="Daylight", "Land", Activity)) 
  
# 9c: RECALCULATE FLIGHT BOUTS AFTER LAND RECLASSIFICATION
  
FlightBouts_lastcheck2<-FlightBoutLengths_lastcheck_reclassify %>%
  dplyr::filter(Activity=="Flight") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  distinct() %>%
  dplyr::mutate(date_time=as.character(date_time))
	
# Make sure no NAs in date_time
nas_date<-subset(FlightBouts_lastcheck2, is.na(date_time))
  
if (nrow(nas_date)>0) {stop(print("Error: nas in date_time")) }
  

FlightBoutLengths_lastcheck2<-FlightBoutLengths_lastcheck_reclassify %>%
  dplyr::select(-BoutNo) %>%
  #dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  dplyr::left_join(FlightBouts_lastcheck2, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Flight", flightLengthMins, 0)) %>%
  ungroup()
	
# 9d: SHORTEN REMAINING OVERLONG FLIGHTS USING RESTWATER

# For Flight bouts that still exceed L1, calculate how many 10-minute observations
# must be removed to bring the bout back to the maximum permitted duration.
# Those observations at the beginning of the bout are changed to RestWater.

FlightBoutLengths_lastcheck_reclassify2<-FlightBoutLengths_lastcheck2 %>%
  ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(index=row_number()) %>%
  dplyr::mutate(Activity_change=ifelse(Activity=="Flight" & flightLengthMins > L1, 1, 0)) %>% # Determine whether parts of a bout need to be re-assigned
  dplyr::mutate(Activity_change_duration=ifelse(Activity_change==1, (flightLengthMins-L1)/10, 0)) %>% # Calculate how many rows need to be changed
  dplyr::mutate(Activity=ifelse(Activity_change==1 & index <=ceiling(Activity_change_duration), "RestWater", Activity )) %>%
  dplyr::select(-c(Activity_change, Activity_change_duration)) 
  
# 9e: FINAL FLIGHT-BOUT CHECK

# Recalculate Flight bouts one final time and stop the function if any Flight
# bout still exceeds L1.

FlightBouts_lastcheck3<-FlightBoutLengths_lastcheck_reclassify2 %>%
  dplyr::filter(Activity=="Flight") %>%
  ungroup() %>%
  dplyr::mutate(date_characters=nchar(date_time)) %>%
  dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
  dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  arrange(individ_id, date_time) %>%
  dplyr::mutate(timediff=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
  replace_na(list("timediff"=0)) %>%
  dplyr::filter(timediff==0 | timediff >10) %>%
  dplyr::mutate(BoutNo=row_number()) %>%
  dplyr::select(date_time, BoutNo) %>%
  distinct() %>%
  dplyr::mutate(date_time=as.character(date_time))
	
# Make sure no NAs in date_time
nas_date<-subset(FlightBouts_lastcheck3, is.na(date_time))
  
if (nrow(nas_date)>0) {stop(print("Error: nas in date_time")) }
  
# Calculate final Flight-bout durations. 
FlightBoutLengths_lastcheck3<-FlightBoutLengths_lastcheck_reclassify2 %>%
  ungroup() %>%
  dplyr::select(-BoutNo) %>%
  #dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
  dplyr::left_join(FlightBouts_lastcheck3, by=c("date_time")) %>%
  dplyr::group_by(Activity) %>%
  fill(BoutNo, .direction=c("down")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
  dplyr::mutate(flightLengthMins=ifelse(Activity=="Flight", flightLengthMins, 0)) %>%
  ungroup()
	
# Determine max flight bout length
maxLength<-max(FlightBoutLengths_lastcheck3$flightLengthMins)
 
 if(maxLength > data$L1[1]) {
 stop(print("Error: flights still too long")) }
  
# 10: CALCULATE MAXIMUM DAILY FLIGHT-BOUT LENGTH

# Extract Flight observations and calculate the longest Flight bout occurring
# on each date. If there are no Flight observations, return an empty summary
# structure with maxFlightBoutsMins set to zero.

daily_max_flightBout_temp<-FlightBoutLengths_lastcheck3 %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::filter(Activity=="Flight")

if (nrow(daily_max_flightBout_temp) >0) {

daily_max_flightBout<-daily_max_flightBout_temp %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(maxFlightBoutsMins=max(flightLengthMins, na.rm=TRUE))
	
	} else {
	
daily_max_flightBout<-daily_max_flightBout_temp %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(maxFlightBoutsMins=0)
	
	}
  
# 11: CALCULATE DAILY ACTIVITY BUDGETS

# Summarise hours spent in each activity and light period together with
# environmental and spatial variables. Only essentially complete days
# (at least 1430 minutes of observations) are retained.
#
# RestWater time is subsequently adjusted using coefficient c. The additional
# RestWater time is subtracted from Active time, without allowing Active time
# to fall below zero.
  
dataCalcDay<-FlightBoutLengths_lastcheck3 %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
  dplyr::filter(Duration>=1430) %>% # This allows for full days including those which have a weird thing where they finish at midnight...
  dplyr::group_by(species, colony, session_id, date) %>%
  dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Active=ifelse(Activity=="Active", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0), sstRandom=sst_random_start) %>%
  dplyr::mutate(sstRandom=ifelse(sstRandom < -1.9, -1.9, sstRandom)) %>%
  dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater1=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tActive=sum(Active)*10/60, tLand=sum(Land)*10/60,
                     tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, MaxDistColKm=max(MaxDistColKm),
                     sst_random=mean(sstRandom), ice_random=mean(ice_random, na.rm=TRUE), air_random=mean(air_mean), distColonyKm_mean=mean(distColonyKm, na.rm=TRUE), mean.lon=mean(lon), mean.lat=mean(lat),
                   boutsRandom=(sum(RandomAllocation)/144), immersionType=mean(max.cond)) %>%
  ungroup() %>%
  dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(tRestWater2=tRestWater1*data$c[1], amountSub=tRestWater2-tRestWater1) %>%
  dplyr::mutate(tActive2=ifelse(amountSub > tActive, 0, tActive-amountSub), tRestWater2=ifelse(amountSub>tActive, tRestWater1 + tActive, tRestWater2)) %>%
  dplyr::select(-c(tActive, tRestWater1)) %>%
  dplyr::rename(tActive=tActive2, tRestWater=tRestWater2) %>%
  dplyr::mutate(tActive=ifelse(tActive<0, 0, tActive)) %>%
  dplyr::mutate(DurationTot=sum(tFlight, tForage, tActive, tRestWater, tLand)) %>%
  dplyr::mutate(dayLengthHrs=tDaylight) %>%
  dplyr::left_join(daily_max_flightBout, by=c("date")) 
  
# 12: CALCULATE ACTIVITY BUDGETS SEPARATELY FOR DAYLIGHT AND DARKNESS

# Create a second summary dataset primarily for verification, containing time
# spent in each behaviour separately during Daylight and Darkness.
# Twilight is included within Daylight for this summary.

dataCalcDay_period2<-reAllocate %>%
  dplyr::ungroup() %>%
  dplyr::mutate(date=substr(date_time, 1, 10)) %>%
  dplyr::mutate(Period=ifelse(Period %in% c("Daylight", "Twilight"), "Daylight", "Darkness")) %>%
  dplyr::group_by(species, colony, session_id, date, Period) %>%
  dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Active=ifelse(Activity=="Active", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                  Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0), sstRandom=sst_random_start) %>%
  dplyr::mutate(sstRandom=ifelse(sstRandom < -1.9, -1.9, sstRandom)) %>%
  dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater1=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tActive=sum(Active)*10/60, tLand=sum(Land)*10/60,
                     tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, MaxDistColKm=max(MaxDistColKm),
                     sst_random=mean(sstRandom), ice_random=mean(ice_random, na.rm=TRUE), distColonyKm_mean=mean(distColonyKm, na.rm=TRUE), mean.lon=mean(lon), mean.lat=mean(lat),
                     c=sample(seq(1.5, 3, 0.1), 1), immersionType=mean(max.cond), air_random=mean(air_mean)) %>%
  rename(tRestWater=tRestWater1)
  
# 13: PREPARE OUTPUT

# Main daily activity/environment results.
actResults<-dataCalcDay

# Daylight/darkness-specific activity results.
boutResults<-dataCalcDay_period2

# Return both outputs as a two-element list.
allResults<-list(actResults, boutResults)

return(allResults)
  
}

##### Atlantic puffin #####

calculateTimeInActivity_AP<-function(data, irmaData){
    
    speciesLatin<-"Fratercula_arctica"    
    actResults<-methodCaitlin(data, irmaData, speciesLatin) 
  
  return(actResults)
  
  
}

##### Little auk #####

calculateTimeInActivity_LiA<-function(data,  irmaData){
    
    speciesLatin<-"Alle_alle"    
    actResults<-methodCaitlin(data,  irmaData, speciesLatin) 
  
  return(actResults)
  
  
}

##### Common guillemot #####

calculateTimeInActivity_CoGu<-function(data,  irmaData){
  
    speciesLatin<-"Uria_aalge"    
    actResults<-methodCaitlin(data,  irmaData, speciesLatin) 
  
  return(actResults)
  
  
}

##### Br guillemot #####

calculateTimeInActivity_BrGu<-function(data, irmaData){
  
    speciesLatin<-"Uria_lomvia"    
    actResults<-methodCaitlin(data, irmaData, speciesLatin) 
  
  return(actResults)
  
}

##### Assign day.night #####

assignTimePeriod<-function(data) {
  
  data$date<-as.Date(data$date)
  timeDay<-getSunlightTimes(data=data, keep=c("nightEnd", "sunrise", "sunset", "night"))
  data2<-data
  data2$Period<-ifelse(data2$date_time >= data2$sunrise & data2$date_time < data2$sunset, "Daylight", NA)
  data2$Period<-ifelse(is.na(data2$nightEnd) & is.na(data2$night) & data2$date_time >= data2$sunset , "Twilight", data2$Period)
  data2$Period<-ifelse(is.na(data2$nightEnd) & is.na(data2$night) & data2$date_time < data2$sunrise , "Twilight", data2$Period)
  data2$Period<-ifelse(!is.na(data2$night) & data2$date_time >= data2$sunset & data2$date_time < data2$night , "Twilight", data2$Period)
  data2$Period<-ifelse(!is.na(data2$nightEnd)  & data2$date_time < data2$sunrise & data2$date_time >= data2$nightEnd , "Twilight", data2$Period)
  data2$Period<-ifelse(!is.na(data2$nightEnd) & !is.na(data2$night) & data2$date_time < data2$nightEnd , "Darkness", data2$Period)
  data2$Period<-ifelse(!is.na(data2$nightEnd) & !is.na(data2$night) & data2$date_time >= data2$night , "Darkness", data2$Period)
  data2$Period<-ifelse(!is.na(data2$nightEnd) & is.na(data2$night) & data2$date_time < data2$nightEnd  , "Twilight", data2$Period)
  data2$Period<-ifelse(!is.na(data2$nightEnd) & is.na(data2$night) & data2$date_time >= data2$sunset  , "Twilight", data2$Period)
  #data2$Period<-ifelse(!is.na(data2$night) &  data2$date_time >= data2$sunset  , "Twilight", data2$Period)
  data2$month<-as.numeric(substr(data2$date, 6, 7))
  data2$Period<-ifelse(is.na(data2$nightEnd) & is.na(data2$sunrise) & is.na(data2$sunset) & is.na(data2$night) & data2$month %in% c(4, 5, 6, 7, 8, 9), "Daylight", data2$Period) # midnight sun
  #data2$Period<-ifelse(is.na(data2$nightEnd) & is.na(data2$sunrise) & is.na(data2$sunset) & !is.na(data2$night) & data2$month %in% c(4, 5, 6, 7, 8, 9), "Daylight", data2$Period) # midnight sun
  data2$Period<-ifelse( is.na(data2$sunrise) & is.na(data2$sunset) & data2$month %in% c(10, 11, 12, 1, 2) & data2$date_time <data2$nightEnd, "Darkness", data2$Period) # Polar night
  data2$Period<-ifelse( is.na(data2$sunrise) & is.na(data2$sunset) & data2$month %in% c(10, 11, 12, 1, 2) & data2$date_time >data2$night, "Darkness", data2$Period) # Polar night
  data2$Period<-ifelse( is.na(data2$sunrise) & is.na(data2$sunset) & data2$month %in% c(10, 11, 12, 1, 2) & data2$date_time >=data2$nightEnd & data2$date_time < data2$night, "Twilight", data2$Period) # Polar night
  data2$Period<-ifelse( is.na(data2$sunrise) & is.na(data2$sunset) & is.na(data2$nightEnd) & data2$month %in% c(10, 11, 12, 1, 2, 3) & data2$date_time >=data2$nauticalDusk, "Darkness", data2$Period) # Polar night
  data2$Period<-ifelse( is.na(data2$sunrise) & is.na(data2$sunset) & is.na(data2$nightEnd) & data2$month %in% c(10, 11, 12, 1, 2, 3) & data2$date_time <data2$nauticalDusk, "Twilight", data2$Period) # Polar night
  #data2$Period<-ifelse( is.na(data2$sunrise) & is.na(data2$sunset) & !is.na(datanightEnd) & !is.na(night) & data2$date_time>= data2$nightEnd & data2$date_time < data2$night,  "Twilight", data2$Period) # Polar night
  #data2$Period<-ifelse( is.na(data2$sunrise) & is.na(data2$sunset) & data2$month %in% c(12) & data2$date_time >data2$night, "Darkness", data2$Period) #
  data2$Period<-ifelse(is.na(data2$sunrise) & is.na(data2$sunset) & !is.na(data2$nightEnd) & is.na(data2$night) & data2$date_time>= data2$nightEnd,  "Twilight", data2$Period) # Polar night
  data2$Period<-ifelse(is.na(data2$sunrise) & is.na(data2$sunset) & !is.na(data2$nightEnd) & is.na(data2$night) & data2$date_time< data2$nightEnd,  "Darkness", data2$Period)
  data2$Period<-ifelse(is.na(data2$sunrise) & is.na(data2$sunset) & is.na(data2$nightEnd) & !is.na(data2$night) & data2$month %in% c(11, 12, 1, 2),  "Darkness", data2$Period) # Polar night
  data2$Period<-ifelse(is.na(data2$sunrise) & is.na(data2$sunset) & is.na(data2$nightEnd) & is.na(data2$night) & data2$month %in% c(10, 11, 12, 1, 2),  "Darkness", data2$Period) # Polar night
  
  # Search for NA twilights & stop if there is as is an error
  naTwilights<-subset(data2, is.na(Period))
  if(nrow(naTwilights)>0) {
    print("NA Twilights")
    break}  
  
  return(data2)
  
}

#### Energetic functions ####

# Function to calculate energetics based on time spent in activity
# species is one of the following: Black-legged kittiwake, Northern fulmar, Atlantic puffin, Little auk, Common guillemot, Brünnich's guillemot
# Data is data frame containing activity budgets
# ColonySub is the colony
# sstVals is just the main data frame (this will be different in the following mapping project)
# Type is daily or monthly resolution (type="daily" or "monthly")
# Type 2 is using sst from pop maps or individual locations ("pop", "ind")
# Function re-directs to individual species functions but will return daily energy expenditure values

calculateEnergetics<-function(species, data, colonySub, sstVals, type, type2) {
  
# PURPOSE: Calculate energy expenditure from activity budgets by directing each
# species to the appropriate species-specific energetic function.
  
# Broad workflow:
#   1. Print the energetic parameter values used in the calculation.
#   2. Identify individual tracking sessions.
#   3. Loop through each session separately.
#   4. Assign species-specific body mass and other required constants.
#   5. Select the appropriate energetic function according to species and data source.
#   6. Add body mass and session ID to the results.
#   7. Combine results across sessions.
  
# INPUTS
#   species   - species being analysed
#   data      - data frame containing activity budgets and energetic parameters
#   colonySub - colony associated with the individual
#   sstVals   - data used to provide SST values for map-based calculations
#   type      - temporal resolution of the calculation (e.g. "daily")
#   type2     - source of environmental data ("ind" for individual locations,
#               "map" for population-level maps)
  
# RETURNS
# A data frame containing energy-expenditure estimates for all tracking sessions.
  
  
# 1: PRINT ENERGETIC PARAMETERS
  
# Print the parameter values selected for this iteration. These include
# activity-specific energetic costs and parameters describing thermoregulatory
# costs in air and water.
 
# Print randomized parameter values
print(paste0("RMR = ", data$RMR[1]))
print(paste0("c1 = ", data$c1[1]))
print(paste0("c2 = ", data$c2[1]))
print(paste0("c3 = ", data$c3[1]))
print(paste0("c4 = ", data$c4[1]))
#print(paste0("c5 = ", data$c5[1]))
print(paste0("TC_air = ", data$TC_air[1]))
print(paste0("TC_water = ", data$TC_water[1]))
print(paste0("Beta_active = ", data$Beta_active[1]))
print(paste0("Beta_rest = ", data$Beta_rest[1]))
print(paste0("LCT_water = ", data$LCT_water[1]))
print(paste0("LCT_air = ", data$LCT_air[1]))

# 2: IDENTIFY TRACKING SESSIONS

# Energetic calculations are carried out separately for each tracking session
# before being combined into a single output.

sessionNo<-unique(data$session_id)

# Create an object to store results across sessions.
energyAll<-list() # Make a list to save results in

# 3: LOOP THROUGH TRACKING SESSIONS

for (session in 1:length(sessionNo)) {

print(paste0("Calculating energy for session ", session))

# Subset activity data to the current tracking session.
dataSub<-subset(data, session_id %in% sessionNo[session]) 
 
if (species=="Black-legged kittiwake") {
 
  weightG<-392
    
  if (type =="daily" & type2=="ind") {
      
  energySpent<-calculateEnergetics_BLK_daily(dataSub, weightG)  
      
  } 
    
  if (type =="daily" & type2=="map") {
      
  energySpent<-calculateEnergetics_BLK_daily_map(data, weightG, sstVals)  
      
    } 
    
  }  
  
  if (species=="Northern fulmar") {
    
  weightG<-728
    
    if (type =="daily"& type2=="ind") {
      
  energySpent<-calculateEnergetics_NF_daily(dataSub, weightG)  
      
    } 
	
	if (type =="daily" & type2=="map") {
      
  energySpent<-calculateEnergetics_NF_daily_map(data, weightG, sstVals)  
      
    } 
    
  } 
  
  if (species=="Common guillemot") {
    
  CostDivider<-803 # g (mean weight of all BrG as I can't find the mass of birds in Kyle's paper)  
  weightG<-940
    
    if (type =="daily" & type2=="ind") {
      
  energySpent<-calculateEnergetics_CoGu_daily(dataSub,  CostDivider,  weightG)  
      
    } 
    
    if (type =="daily" & type2=="map") {
      
  energySpent<-calculateEnergetics_CoGu_daily_map(data, CostDivider, weightG, sstVals)  
      
    } 
    
  } 
  
  if (species=="Brünnich's guillemot") {
    
  CostDivider<-803 # g (mean weight of all BrG as I can't find the mass of birds in Kyle's paper)  
  weightG<-980
    
    if (type =="daily" & type2=="ind") {
      
  energySpent<-calculateEnergetics_BrGu_daily(dataSub,  CostDivider, weightG)  
      
    } 
    
    if (type =="daily" & type2=="map") {
      
  energySpent<-calculateEnergetics_BrGu_daily_map(data, CostDivider, weightG, sstVals)  
      
    } 
    
  }   
  
  if (species=="Little auk") {
    
  CostDivider<-803 # g (mean weight of all BrG as I can't find the mass of birds in Kyle's paper)  
  weightG<-149
    
    if (type =="daily" & type2=="ind") {
      
  energySpent<-calculateEnergetics_LiA_daily(dataSub,  CostDivider,  weightG)  
      
    } 
    
    if (type =="daily" & type2=="map") {
      
  energySpent<-calculateEnergetics_LiA_daily_map(data, CostDivider, weightG, sstVals)  
      
    } 
    
  }    
  
  if (species=="Atlantic puffin") {
    
  CostDivider<-803 # g (mean weight of all BrG as I can't find the mass of birds in Kyle's paper)  
  weightG<-395
    
    if (type =="daily" & type2=="ind") {
      
  energySpent<-calculateEnergetics_AP_daily(dataSub,  CostDivider, weightG)  
      
    } 
    
    if (type =="daily" & type2=="map") {
      
  energySpent<-calculateEnergetics_AP_daily_map(data, CostDivider, weightG,  sstVals)  
      
    } 
    
  } 

# 4: COMBINE SESSION RESULTS

# Add the species-specific body mass and original session identifier to the
# energetic estimates before combining results across tracking sessions.
  
  energySpent$weight<-weightG 
  energySpent$session_id<-sessionNo[session]
  energyAll<-rbind(energyAll, energySpent)
  
}

# 5: PREPARE OUTPUT
  
  return(energyAll)    
  
}

##### Kittiwake #####

calculateEnergetics_BLK_daily<-function(data, weightG) {
  
# PURPOSE: Calculate daily energy expenditure for Black-legged kittiwakes from
# activity budgets, body mass, and environmental temperature.
  
# Broad workflow:
# 1. Extract activity-specific energetic and thermoregulatory parameters.
# 2. Convert activity-cost coefficients to hourly energetic costs.
# 3. Scale energetic parameters from reference body masses to the focal body mass.
# 4. Define temperature-dependent energetic costs below lower critical temperatures.
# 5. Calculate energetic costs separately for each activity.
# 6. Sum activity-specific costs to estimate total daily energy expenditure.
  
# INPUTS
# data- data frame containing daily activity budgets, environmental
# conditions, and energetic parameters
# weightG - body mass used to scale energetic costs (g)
  
# RETURNS
# A data frame containing the original daily activity data plus activity-specific
# and total daily energy expenditure estimates (kJ).
  
# 1: EXTRACT & CONVERT ENERGETIC PARAMETERS
  
# Caloric conversion factor used to convert oxygen consumption to energetic
# expenditure (20.1 J per mL O2; Schmidt-Nielsen 1997).
cf<-20.1
  
# Extract activity-specific energetic coefficients. Coefficients derived from
# Tremblay et al. 2024 are converted from daily to hourly energetic costs.

# Rest on Water:
restCoef<-data$c4[1]  # kJ.g.day
restCoef<-restCoef/24 # kJ.g.hr
  
# Foraging: 
forageCoef<-data$c2[1]/24 # kJ.g.hr
  
# On land: 
landCoef<-data$c3[1]
landCoef<-landCoef/24 # kJ.g.hr
  
# Flight: 
flightCoef<-data$c1[1]
flightCoef<-flightCoef/24 # kJ.g.hr
  
# Extract parameters used to calculate temperature-dependent energetic costs.
TCCoef_water<-data$TC_water[1] # kJ g-1 hr-1 C-1
TCCoef_air<-data$TC_air[1] # kJ g-1 hr-1 C-1
  
# 2: SCALE ENERGETIC PARAMETERS TO BODY MASS

# Scale activity-specific energetic costs from the reference body mass of 450 g
# to the focal body mass using a mass-scaling exponent of 0.717.
flightConstantx<-((flightCoef*450)/450^0.717)*weightG^0.717
restConstant2x<-((restCoef*450)/450^0.717)*weightG^0.717
forageConstantx<-((forageCoef*450)/450^0.717)*weightG^0.717
landConstantx<-((landCoef*450)/450^0.717)*weightG^0.717

# Convert and scale the resting intercept and thermal-conductance parameters.
# These parameters use a reference body mass of 365 g.
TCx_water<-(((TCCoef_water*cf/1000)*365)/365^0.717)*weightG^0.717
TCx_air<-(((TCCoef_air*cf/1000)*365)/365^0.717)*weightG^0.717

# 3: DEFINE TEMPERATURE-DEPENDENT RESTING COSTS

# Below the lower critical temperature (LCT), energetic expenditure increases
# as environmental temperature decreases. Separate LCTs are used for birds
# resting on water and on land.
LCT_water<-data$LCT_water[1] # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1474-919X.2006.00618.x
LCT_air<-data$LCT_air[1] # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1474-919X.2006.00618.x

# Calculate intercepts so that temperature-dependent energetic costs below the
# LCT meet the thermoneutral activity cost at the corresponding LCT.
restConstant1x<-(LCT_water*TCx_water + restConstant2x)
beta_land<-landConstantx + LCT_air*TCx_air

# 4: CALCULATE DAILY ENERGY EXPENDITURE

# Calculate energetic expenditure separately for each activity. Resting costs
# on water and land increase below their respective lower critical temperatures,
# whereas flight and foraging costs are independent of temperature here.
  
energySub2<-data %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(DEEkJ_active=0, DEEkJ_active_col=0) %>%
  dplyr::mutate(DEEkJ_rest=ifelse(sst_random <=LCT_water, (restConstant1x - TCx_water*sst_random)*tRestWater, restConstant2x*tRestWater)) %>%
  dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <=LCT_water, (restConstant1x - TCx_water*sst_random_colony)*tRestWater, restConstant2x*tRestWater)) %>%
  dplyr::mutate(DEEkJ_flight=flightConstantx*tFlight) %>%
  dplyr::mutate(DEEkJ_forage=forageConstantx*tForage) %>%
  dplyr::mutate(DEEkJ_restland=ifelse(air_random <= LCT_air, (beta_land - air_random*TCx_air)*tLand, landConstantx*tLand)) %>%
  dplyr::mutate(DEEkJ_restland2=landConstantx*tLand) %>%
  dplyr::mutate(DEEkJ=DEEkJ_rest + DEEkJ_flight + DEEkJ_forage + DEEkJ_restland) %>%
  dplyr::mutate(DEEkJ_col=DEEkJ_rest_col + DEEkJ_flight + DEEkJ_forage + DEEkJ_restland) %>%
  dplyr::mutate(weight=weightG)

# 5: PREPARE OUTPUT
  
  return(energySub2)
  
}

calculateEnergetics_BLK_daily_map<-function(data, weightG, sstVals) {
  
# cf is caloric conversion factor of 20.1 J per mL O2 (Schmidt-Nielsen 1997)
  cf<-20.1
  
  # We determine RMR
  RMR<-data$RMR[1]
  
  # Rest coef is generated from Tremblay et al. sample size is 50
  restCoef<-data$c4[1]
  restCoef<-restCoef/24 # kJ.g.hr
  
  # forage coef is generated from  Tremblay et al 2024 & is a mix of flapping & swim
  forageCoef<-data$c2[1]/24 # kJ.g.hr
  
  # land coef taken from Tremblay et al. 
  landCoef<-data$c3[1]
  landCoef<-landCoef/24 # kJ.g.hr
  
  # Flap coef are from Tremblay et al. 
  flightCoef<-data$c1[1]
  flightCoef<-flightCoef/24 # kJ.g.hr
  
  # We will make a fake error distribution for beta & TC based on mean errors for other variables which is 29%
  betaCoef<-data$Beta_rest[1]
  TCCoef_water<-data$TC_water[1]
  TCCoef_air<-data$TC_air[1]
  
  # Account for change in constants
  flightConstantx<-((flightCoef*450)/450^0.717)*weightG^0.717
  restConstant2x<-((restCoef*450)/450^0.717)*weightG^0.717
  forageConstantx<-((forageCoef*450)/450^0.717)*weightG^0.717
  landConstantx<-((landCoef*450)/450^0.717)*weightG^0.717
  betax<-(((betaCoef*cf/1000)*365)/365^0.717)*weightG^0.717
  TCx_water<-(((TCCoef_water*cf/1000)*365)/365^0.717)*weightG^0.717
  TCx_air<-(((TCCoef_air*cf/1000)*365)/365^0.717)*weightG^0.717
  
  # Adjust beta so that beta-SST*TC is equal to rest constant 2 at LCT
  LCT_water<-data$LCT_water[1] # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1474-919X.2006.00618.x
  LCT_air<-data$LCT_air[1] # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1474-919X.2006.00618.x
  restConstant1x<-(LCT_water*TCx_water + restConstant2x)
  beta_land<-landConstantx + LCT_air*TCx_air
  
  # Turn SST raster into a data frame
  sst<-raster::subset(sstVals, 1)
  temp<-raster::subset(sstVals, 2)
  sstDf<-as.data.frame(sst, xy=TRUE)
  airDf<-as.data.frame(temp, xy=TRUE)
  colnames(sstDf)<-c("x", "y", "sst")
  colnames(airDf)<-c("x", "y", "temp")
  sstDf$temp<-airDf$temp
  
  # Calculate energy for every cell according to sst in that cell
  sstDf$DEEkJ_active=0
  sstDf$DEEkJ_rest=ifelse(sstDf$sst <=LCT_water, (restConstant1x - TCx_water*sstDf$sst)*data$tRestWater_month, restConstant2x*data$tRestWater_month)
  sstDf$DEEkJ_flight=flightConstantx*data$tFlight_month
  sstDf$DEEkJ_forage=forageConstantx*data$tForage_month
  sstDf$DEEkJ_restland=ifelse(sstDf$temp <= LCT_air, (beta_land - sstDf$temp*TCx_air)*data$tLand_month, landConstantx*data$tLand_month)
  sstDf$DEEkJ=sstDf$DEEkJ_rest + sstDf$DEEkJ_flight + sstDf$DEEkJ_forage + sstDf$DEEkJ_restland
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, sst, temp, DEEkJ) 
  
  return(sstDf_final)


}

##### Northern fulmar #####

calculateEnergetics_NF_daily<-function(data, weightG) {

# PURPOSE: Calculate daily energy expenditure for Northern fulmars from
# activity budgets, body mass, and environmental temperature.
  
# Broad workflow:
# 1. Extract activity-specific energetic and thermoregulatory parameters.
# 2. Convert and scale energetic parameters to the focal body mass.
# 3. Define temperature-dependent energetic costs below lower critical temperatures.
# 4. Calculate energetic costs separately for each activity.
# 5. Sum activity-specific costs to estimate total daily energy expenditure.
  
# INPUTS
# data - data frame containing daily activity budgets, environmental
# conditions, and energetic parameters
# weightG - body mass used to scale energetic costs (g)
  
# RETURNS
# A data frame containing the original daily activity data plus activity-specific
# and total daily energy expenditure estimates (kJ).
  
# 1: EXTRACT ENERGETIC PARAMETERS
  
# Caloric conversion factor used to convert oxygen consumption to energetic
# expenditure (20.1 J per mL O2; Schmidt-Nielsen 1997).
cf<-20.1
  
# Extract resting metabolic rate, sampled from values reported by
# Gabrielsen et al. (1988).
RMR<-data$RMR[1] # ml O2 g-1 hr-1
  
# Extract the energetic cost of resting on water from Bevan et al. (1997).
restCoef<-data$c4[1]
  
# Foraging cost is derived from Tremblay et al. (2024)
# Convert from daily to hourly cost.
forageCoef<-data$c2[1]/24 # kJ.g.hr
  
# Extract the energetic cost of being on land from Bevan et al. (1997).
landCoef<-data$c3[1]
  
# Extract the flight cost from Bevan et al. (1997).
flightCoef<-data$c1[1]
  
# Extract thermoregulatory parameters for resting in water and air.
TCCoef_water<-data$TC_water[1] # kJ g-1 hr-1 C-1
TCCoef_air<-data$TC_air[1] # kJ g-1 hr-1 C-1

# 2: SCALE ENERGETIC PARAMETERS TO BODY MASS

# Convert coefficients based on RMR and oxygen consumption to kJ and scale
# them from the 651-g reference body mass to the focal body mass using a
# mass-scaling exponent of 0.765.

flightConstantx<-(((flightCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
restConstant2x<-(((restCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
landConstantx<-(((landCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765

# Foraging costs originate from a different source and are scaled from a
# reference body mass of 450 g.

forageConstantx<-((forageCoef*450)/450^0.717)*weightG^0.765

# Convert and scale the resting intercept and thermal-conductance parameters
# using the 651-g reference body mass.

TCx_water<-(((TCCoef_water*cf/1000)*651)/651^0.765)*weightG^0.765
TCx_air<-(((TCCoef_air*cf/1000)*651)/651^0.765)*weightG^0.765
  
# 3: DEFINE TEMPERATURE-DEPENDENT RESTING COSTS

# Below the lower critical temperature (LCT), energetic expenditure increases
# as environmental temperature decreases. Separate LCTs are used for birds
# resting on water and on land.

LCT_water<-data$LCT_water # Gabrielsen et al. 1988
LCT_air<-data$LCT_air # Gabrielsen et al. 1988

# Set the intercepts of the temperature-dependent relationships so that the
# energetic cost at each LCT equals the corresponding thermoneutral resting
# or land cost.

restConstant1x<-(LCT_water*TCx_water + restConstant2x)
beta_land<-landConstantx + TCx_air*LCT_air
  
# 4: CALCULATE DAILY ENERGY EXPENDITURE

# Calculate energetic expenditure separately for each activity. Resting costs
# on water and land increase below their respective lower critical temperatures,
# whereas flight and foraging costs are independent of temperature here.
	
energySub2<-data %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(DEEkJ_active=0, DEEkJ_active_col=0) %>%
  dplyr::mutate(DEEkJ_rest=ifelse(sst_random <=LCT_water, (restConstant1x - TCx_water*sst_random)*tRestWater, restConstant2x*tRestWater)) %>%
  dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <=LCT_water, (restConstant1x - TCx_water*sst_random_colony)*tRestWater, restConstant2x*tRestWater)) %>%
  dplyr::mutate(DEEkJ_flight=flightConstantx*tFlight) %>%
  dplyr::mutate(DEEkJ_forage=forageConstantx*tForage) %>%
  dplyr::mutate(DEEkJ_restland=ifelse(air_random<=LCT_air, (beta_land - air_random*TCx_air)*tLand ,landConstantx*tLand)) %>%
  dplyr::mutate(DEEkJ_restland2=landConstantx*tLand) %>%
  dplyr::mutate(DEEkJ=DEEkJ_rest + DEEkJ_flight + DEEkJ_forage + DEEkJ_restland) %>%
  dplyr::mutate(DEEkJ_col=DEEkJ_rest_col + DEEkJ_flight + DEEkJ_forage + DEEkJ_restland) %>%
  dplyr::mutate(weight=weightG)
  
# 5: PREPARE OUTPUT

  return(energySub2)
  
}

calculateEnergetics_NF_daily_map<-function(data, weightG, sstVals) {
  
  # Calculate thermo-regulation in water based on same method used by Elliott & Gaston for guillemots
  # In Jodice et al. where SST = 14.1 in summer, tRest = 1.1*RMR
  # So we adujst the thermo equation from Geir so that at 14.1, beta-TC*SST is equal to 1.1*RMR
  
  # cf is caloric conversion factor of 20.1 J per mL O2 (Schmidt-Nielsen 1997)
  cf<-20.1
  
  # Set RMR : sample from a uniform distribution (Gabrielsen et al. 1988 - see table 1) 
  RMR<-data$RMR[1]
  
  # Rest coef is generated from Bevan et al. 1997
  restCoef<-data$c4[1]
  
  # forage coef is generated from Tremblay et al 2024. It's an average of flight & swim coefs
  #flappingCoef<-data$c1[1]
  #swimmingCoef<-data$c2[1]
  forageCoef<-data$c2[1]/24 # kJ.g.hr
  
  # land coef taken from bevan et al. 1997
  landCoef<-data$c3[1]
  
  # Flap & glide coefs are from Bevan et al. 1997
  flightCoef<-data$c1[1]
  
  # betaCoef - we assume error is equal to average error of others which in this case is 24%
  betaCoef<-data$Beta_rest[1]
  
  # TCCoef - we assume error is equal to average error of others which in this case is 24%
  TCCoef_water<-data$TC_water[1]
  TCCoef_air<-data$TC_air[1]
  
  # Account for change in constants
  #flightConstantx<-(7.3*RMR*cf/1000)*weightG
  flightConstantx<-(((flightCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  restConstant2x<-(((restCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  forageConstantx<-((forageCoef*450)/450^0.717)*weightG^0.765
  landConstantx<-(((landCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  betax<-(((betaCoef*cf/1000)*651)/651^0.765)*weightG^0.765
  TCx_water<-(((TCCoef_water*cf/1000)*651)/651^0.765)*weightG^0.765
  TCx_air<-(((TCCoef_air*cf/1000)*651)/651^0.765)*weightG^0.765
  
  # Adjust beta so that beta-SST*TC is equal to rest constant 2 at LCT
  LCT_water<-data$LCT_water # Gabrielsen et al. 1988
  LCT_air<-data$LCT_air # Gabrielsen et al. 1988
  restConstant1x<-(LCT_water*TCx_water + restConstant2x)
  beta_land<-landConstantx + TCx_air*LCT_air
  
  # Turn SST raster into a data frame
  sst<-subset(sstVals, 1)
  temp<-subset(sstVals, 2)
  sstDf<-as.data.frame(sst, xy=TRUE)
  airDf<-as.data.frame(temp, xy=TRUE)
  colnames(sstDf)<-c("x", "y", "sst")
  colnames(airDf)<-c("x", "y", "temp")
  sstDf$temp<-airDf$temp
  
  # Calculate energy for every cell according to sst in that cell
  sstDf$DEEkJ_active=0
  sstDf$DEEkJ_rest=ifelse(sstDf$sst <=LCT_water, (restConstant1x - TCx_water*sstDf$sst)*data$tRestWater_month, restConstant2x*data$tRestWater_month)
  sstDf$DEEkJ_flight=flightConstantx*data$tFlight_month
  sstDf$DEEkJ_forage=forageConstantx*data$tForage_month
  sstDf$DEEkJ_restland=ifelse(sstDf$temp <= LCT_air, (beta_land - sstDf$temp*TCx_air)*data$tLand_month, landConstantx*data$tLand_month)
  sstDf$DEEkJ=sstDf$DEEkJ_rest + sstDf$DEEkJ_flight + sstDf$DEEkJ_forage + sstDf$DEEkJ_restland
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, sst, temp, DEEkJ) 
  
  return(sstDf_final)
  
  
  
  
}

##### Common guillemot #####

calculateEnergetics_CoGu_daily<-function(data, CostDivider, weightG) {
  
# PURPOSE: Calculate daily energy expenditure for Common guillemots from
# activity budgets, body mass, and environmental temperature.
  
# Broad workflow:
# 1. Extract activity-specific energetic and thermoregulatory parameters.
# 2. Convert energetic coefficients to the units required by the model.
# 3. Scale energetic parameters from the reference body mass to the focal body mass.
# 4. Define temperature-dependent energetic costs below lower critical temperatures.
# 5. Calculate energetic costs separately for each activity.
# 6. Sum activity-specific costs to estimate total daily energy expenditure.
  
# INPUTS
# data - data frame containing daily activity budgets, environmental
# conditions, and energetic parameters
# CostDivider - reference body mass used for allometric scaling (g)
# weightG  - body mass used to scale energetic costs (g)
  
# RETURNS
# A data frame containing the original daily activity data plus activity-specific
# and total daily energy expenditure estimates (kJ).
  
# 1: EXTRACT & CONVERT ENERGETIC PARAMETERS
  
# Extract the energetic cost of flight and convert to kJ per hour
  
flightCoef<-data$c1[1]
flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
# Extract the intercepts used to describe temperature-dependent energetic costs
# for active and resting birds on water, together with thermal conductance in
# water and air.

activeCoef<-data$Beta_active[1] # kj.hr
restCoef1<-data$Beta_rest[1] # kj.hr
TC_water<-data$TC_water[1] # kJ.hr
TC_air<-data$TC_air[1] # ml O2 hr
  
# Extract the thermoneutral resting-water and land costs and convert them to
# kJ per hour.

restCoef2<-data$c4[1]
restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
landCoef<-data$c3[1]
landCoef_kj<-landCoef*3.6 # kj.hr 
  
# 2: SCALE ENERGETIC PARAMETERS TO BODY MASS

# Calculate the conversion factor used to scale energetic costs from the
# reference body mass to the focal body mass using an exponent of 0.689.
  
convf<-1/CostDivider^0.689 
  
# Apply allometric scaling to activity-specific costs and water
# thermoregulatory parameters.

flightConstant<-(flightCoef_kj*convf)*weightG^0.689
activeConstant<-(activeCoef*convf)*weightG^0.689
landConstant<-(landCoef_kj*convf)*weightG^0.689
restConstant1<-(restCoef1*convf)*weightG^0.689
restConstant2<-(restCoef2_kj*convf)*weightG^0.689
TC_water_Constant<-(TC_water*convf)*weightG^0.689

# Convert thermal conductance in air from oxygen consumption to kJ and scale
# from its reference body mass to the focal body mass.

TC_air_Constant<-TC_air*20.1/1000*819.3/819.3^0.689*weightG^0.689
  
# 3: DEFINE TEMPERATURE-DEPENDENT ACTIVITY COSTS

# Below the lower critical temperature (LCT), energetic costs on water increase
# as SST decreases. Above the LCT, costs are held at their value at the LCT.

LCT_water<-data$LCT_water[1]

# Calculate the active and resting costs at the water LCT. These provide the
# constant thermoneutral costs used when SST is above the LCT.

restConstant_adjust<-restConstant1 - LCT_water*TC_water_Constant # So that rest is equal to land when sst > LCT
activeConstant_adjust<-activeConstant - LCT_water*TC_water_Constant

# Define the corresponding temperature-dependent relationship for birds on
# land. The intercept is set so that energetic cost at the air LCT equals the
# thermoneutral land cost.

LCT_air<-data$LCT_air[1]
beta_land<-landConstant + TC_air_Constant*LCT_air
  
# 4: CALCULATE DAILY ENERGY EXPENDITURE

# Calculate energetic expenditure separately for each activity. Active and
# resting costs on water increase below the water LCT, while land costs increase
# below the air LCT.
  
energySub2<-data %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT_water, (activeConstant-sst_random*TC_water_Constant)*tActive, activeConstant_adjust*tActive)) %>%
  dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT_water, (activeConstant-sst_random_colony*TC_water_Constant)*tActive, activeConstant_adjust*tActive)) %>%
  dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT_water, (restConstant1 - TC_water_Constant*sst_random)*tRestWater , restConstant_adjust*tRestWater)) %>%
  dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT_water, (restConstant1 - TC_water_Constant*sst_random_colony)*tRestWater , restConstant_adjust*tRestWater)) %>%
  dplyr::mutate(DEEkJ_flight=flightConstant*tFlight) %>%
  dplyr::mutate(DEEkJ_forage=0) %>%
  dplyr::mutate(DEEkJ_restland=ifelse(air_random <= LCT_air, (beta_land - air_random*TC_air_Constant)*tLand, landConstant*tLand)) %>%
  dplyr::mutate(DEEkJ_restland2=landConstant*tLand) %>%
  dplyr::mutate(DEEkJ=DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland) %>%
  dplyr::mutate(DEEkJ_col=DEEkJ_active_col + DEEkJ_rest_col + DEEkJ_flight + DEEkJ_restland) %>%
  dplyr::mutate(weight=weightG)

# 5: PREPARE OUTPUT
  
return(energySub2)
  
}

calculateEnergetics_CoGu_daily_map<-function(data, CostDivider, weightG, sstVals) {
  
  # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC_water<-data$TC_water[1] # kJ.hr
  TC_air<-data$TC_air[1] # ml O2 hr
  
  # Rest & land coefs
  restCoef2<-data$c4[1]
  restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
  landCoef<-data$c3[1]
  landCoef_kj<-landCoef*3.6 # kj.hr -> should be the same as the rest coefficient
  
  # Active coef when thermoneutral (from kyle's paper)
  #activeCoef_neut<-data$c5[1]
  #activeCoef_neut_kj<-activeCoef_neut*3.6 # kJ.hr
  
  #if (activeCoef_neut_kj > activeCoef) stop (print("Error with active coefs")) # Because the coef when thermoneutral must be lower...
  
  # Here is a conversion factor to transform these to g.kj which incorporates allometric scaling
  convf<-1/CostDivider^0.689 # no wing loading
  
  # Account for change in constants
  flightConstant<-(flightCoef_kj*convf)*weightG^0.689
  activeConstant<-(activeCoef*convf)*weightG^0.689
 # activeConstant2<-(activeCoef_neut_kj*convf)*weightG^0.689
  restConstant1<-(restCoef1*convf)*weightG^0.689
  restConstant2<-(restCoef2_kj*convf)*weightG^0.689
  TC_water_Constant<-(TC_water*convf)*weightG^0.689
  TC_air_Constant<-TC_air*20.1/1000*819.3/819.3^0.689*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT_water<-data$LCT_water[1]
  restConstant_adjust<-restConstant1 - LCT_water*TC_water_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT_water*TC_water_Constant
  
  LCT_air<-data$LCT_air[1]
  beta_land<-landConstant + TC_air_Constant*LCT_air
  
  # Turn SST raster into a data frame
  sst<-subset(sstVals, 1)
  temp<-subset(sstVals, 2)
  sstDf<-as.data.frame(sst, xy=TRUE)
  airDf<-as.data.frame(temp, xy=TRUE)
  colnames(sstDf)<-c("x", "y", "sst")
  colnames(airDf)<-c("x", "y", "temp")
  sstDf$temp<-airDf$temp
  
  # Calculate energy for every cell according to sst in that cell
  sstDf$DEEkJ_active=ifelse(sstDf$sst <=LCT_water, (activeConstant-sstDf$sst*TC_water_Constant)*data$tActive_month, activeConstant_adjust*data$tActive_month)
  sstDf$DEEkJ_rest=ifelse(sstDf$sst <= LCT_water, (restConstant1 - TC_water_Constant*sstDf$sst)*data$tRestWater_month , restConstant_adjust*data$tRestWater_month)
  sstDf$DEEkJ_flight=flightConstant*data$tFlight_month
  sstDf$DEEkJ_forage=0
  sstDf$DEEkJ_restland=ifelse(sstDf$temp <= LCT_air, (beta_land - sstDf$temp*TC_air_Constant)*data$tLand_month, landConstant*data$tLand_month)
  sstDf$DEEkJ=sstDf$DEEkJ_rest + sstDf$DEEkJ_flight + sstDf$DEEkJ_active + sstDf$DEEkJ_restland
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, sst, temp, DEEkJ) 
  
  return(sstDf_final)
  
}

##### Brunnich's guillemot #####

calculateEnergetics_BrGu_daily<-function(data, CostDivider, weightG) {

# PURPOSE: Calculate daily energy expenditure for Common guillemots from
# activity budgets, body mass, and environmental temperature.
  
# Broad workflow:
# 1. Extract activity-specific energetic and thermoregulatory parameters.
# 2. Convert energetic coefficients to the units required by the model.
# 3. Scale energetic parameters from the reference body mass to the focal body mass.
# 4. Define temperature-dependent energetic costs below lower critical temperatures.
# 5. Calculate energetic costs separately for each activity.
# 6. Sum activity-specific costs to estimate total daily energy expenditure.
  
# INPUTS
# data - data frame containing daily activity budgets, environmental
# conditions, and energetic parameters
# CostDivider - reference body mass used for allometric scaling (g)
# weightG  - body mass used to scale energetic costs (g)
  
# RETURNS
# A data frame containing the original daily activity data plus activity-specific
# and total daily energy expenditure estimates (kJ).
  
# 1: EXTRACT & CONVERT ENERGETIC PARAMETERS
  
# Extract the energetic cost of flight and convert to kJ per hour

flightCoef<-data$c1[1]
flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
# Extract the intercepts used to describe temperature-dependent energetic costs
# for active and resting birds on water, together with thermal conductance in
# water and air.

activeCoef<-data$Beta_active[1] # kj.hr
restCoef1<-data$Beta_rest[1] # kj.hr
TC_water<-data$TC_water[1] # kJ.hr
TC_air<-data$TC_air[1] # ml O2 hr
  
# Extract the thermoneutral resting-water and land costs and convert them to
# kJ per hour.

restCoef2<-data$c4[1]
restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
landCoef<-data$c3[1]
landCoef_kj<-landCoef*3.6 # kj.hr 
  
# 2: SCALE ENERGETIC PARAMETERS TO BODY MASS

# Calculate the conversion factor used to scale energetic costs from the
# reference body mass to the focal body mass using an exponent of 0.689.

convf<-1/CostDivider^0.689 # no wing loading
  
# Apply allometric scaling to activity-specific costs and water
# thermoregulatory parameters.

flightConstant<-(flightCoef_kj*convf)*weightG^0.689
activeConstant<-(activeCoef*convf)*weightG^0.689
restConstant1<-(restCoef1*convf)*weightG^0.689
restConstant2<-(restCoef2_kj*convf)*weightG^0.689
landConstant<-(landCoef_kj*convf)*weightG^0.689
TC_water_Constant<-(TC_water*convf)*weightG^0.689

# Convert thermal conductance in air from oxygen consumption to kJ and scale
# from its reference body mass to the focal body mass.

TC_air_Constant<-TC_air*20.1/1000*819.3/819.3^0.689*weightG^0.689

# 3: DEFINE TEMPERATURE-DEPENDENT ACTIVITY COSTS

# Below the lower critical temperature (LCT), energetic costs on water increase
# as SST decreases. Above the LCT, costs are held at their value at the LCT.  

LCT_water<-data$LCT_water[1]

# Calculate the active and resting costs at the water LCT. These provide the
# constant thermoneutral costs used when SST is above the LCT.

restConstant_adjust<-restConstant1 - LCT_water*TC_water_Constant # So that rest is equal to land when sst > LCT
activeConstant_adjust<-activeConstant - LCT_water*TC_water_Constant

# Define the corresponding temperature-dependent relationship for birds on
# land. The intercept is set so that energetic cost at the air LCT equals the
# thermoneutral land cost.

LCT_air<-data$LCT_air[1]
beta_land<-landConstant + TC_air_Constant*LCT_air
  
# 4: CALCULATE DAILY ENERGY EXPENDITURE

# Calculate energetic expenditure separately for each activity. Active and
# resting costs on water increase below the water LCT, while land costs increase
# below the air LCT.
  
energySub2<-data %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT_water, (activeConstant-sst_random*TC_water_Constant)*tActive, activeConstant_adjust*tActive)) %>%
  dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT_water, (activeConstant-sst_random_colony*TC_water_Constant)*tActive, activeConstant_adjust*tActive)) %>%
  dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT_water, (restConstant1 - TC_water_Constant*sst_random)*tRestWater , restConstant_adjust*tRestWater)) %>%
  dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT_water, (restConstant1 - TC_water_Constant*sst_random_colony)*tRestWater , restConstant_adjust*tRestWater)) %>%
  dplyr::mutate(DEEkJ_flight=flightConstant*tFlight) %>%
  dplyr::mutate(DEEkJ_forage=0) %>%
  dplyr::mutate(DEEkJ_restland=ifelse(air_random <= LCT_air, (beta_land - air_random*TC_air_Constant)*tLand, landConstant*tLand)) %>%
  dplyr::mutate(DEEkJ_restland2=landConstant*tLand) %>%
  dplyr::mutate(DEEkJ=DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland) %>%
  dplyr::mutate(DEEkJ_col=DEEkJ_active_col + DEEkJ_rest_col + DEEkJ_flight + DEEkJ_restland) %>%
  dplyr::mutate(weight=weightG)
  
# 5: PREPARE OUTPUT

return(energySub2)
  
  
}

calculateEnergetics_BrGu_daily_map<-function(data, CostDivider, weightG, sstVals) {
  
  # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC_water<-data$TC_water[1] # kJ.hr
  TC_air<-data$TC_air[1] # ml O2 hr
  
  # Rest & land coefs
  restCoef2<-data$c4[1]
  restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
  landCoef<-data$c3[1]
  landCoef_kj<-landCoef*3.6 # kj.hr -> should be the same as the rest coefficient
  
  # Active coef when thermoneutral (from kyle's paper)
  #activeCoef_neut<-data$c5[1]
  #activeCoef_neut_kj<-activeCoef_neut*3.6 # kJ.hr
  
  #if (activeCoef_neut_kj > activeCoef) stop (print("Error with active coefs")) # Because the coef when thermoneutral must be lower...
  
  # Here is a conversion factor to transform these to g.kj which incorporates allometric scaling
  convf<-1/CostDivider^0.689 # no wing loading
  
  # Account for change in constants
  flightConstant<-(flightCoef_kj*convf)*weightG^0.689
  activeConstant<-(activeCoef*convf)*weightG^0.689
 # activeConstant2<-(activeCoef_neut_kj*convf)*weightG^0.689
  restConstant1<-(restCoef1*convf)*weightG^0.689
  restConstant2<-(restCoef2_kj*convf)*weightG^0.689
  TC_water_Constant<-(TC_water*convf)*weightG^0.689
  TC_air_Constant<-TC_air*20.1/1000*819.3/819.3^0.689*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT_water<-data$LCT_water[1]
  restConstant_adjust<-restConstant1 - LCT_water*TC_water_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT_water*TC_water_Constant
  
  LCT_air<-data$LCT_air[1]
  beta_land<-landConstant + TC_air_Constant*LCT_air
  
  # Turn SST raster into a data frame
  sst<-subset(sstVals, 1)
  temp<-subset(sstVals, 2)
  sstDf<-as.data.frame(sst, xy=TRUE)
  airDf<-as.data.frame(temp, xy=TRUE)
  colnames(sstDf)<-c("x", "y", "sst")
  colnames(airDf)<-c("x", "y", "temp")
  sstDf$temp<-airDf$temp
  
  # Calculate energy for every cell according to sst in that cell
  sstDf$DEEkJ_active=ifelse(sstDf$sst <=LCT_water, (activeConstant-sstDf$sst*TC_water_Constant)*data$tActive_month, activeConstant_adjust*data$tActive_month)
  sstDf$DEEkJ_rest=ifelse(sstDf$sst <= LCT_water, (restConstant1 - TC_water_Constant*sstDf$sst)*data$tRestWater_month , restConstant_adjust*data$tRestWater_month)
  sstDf$DEEkJ_flight=flightConstant*data$tFlight_month
  sstDf$DEEkJ_forage=0
  sstDf$DEEkJ_restland=ifelse(sstDf$temp <= LCT_air, (beta_land - sstDf$temp*TC_air_Constant)*data$tLand_month, landConstant*data$tLand_month)
  sstDf$DEEkJ=sstDf$DEEkJ_rest + sstDf$DEEkJ_flight + sstDf$DEEkJ_active + sstDf$DEEkJ_restland
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, sst, temp, DEEkJ) 
  
  return(sstDf_final)
  
  
}

##### Little auk #####

calculateEnergetics_LiA_daily<-function(data, CostDivider,  weightG) {

# PURPOSE: Calculate daily energy expenditure for Little auks from activity
# budgets, body mass, and environmental temperature.
  
# Broad workflow:
# 1. Extract activity-specific energetic and thermoregulatory parameters.
# 2. Convert energetic coefficients to the units required by the model.
# 3. Scale energetic parameters from reference body masses to the focal body mass.
# 4. Define temperature-dependent energetic costs below lower critical temperatures.
# 5. Calculate energetic costs separately for each activity.
# 6. Sum activity-specific costs to estimate total daily energy expenditure.
  
# INPUTS
# data- data frame containing daily activity budgets, environmental
#                 conditions, and energetic parameters
# CostDivider - reference body mass used for allometric scaling (g)
# weightG     - body mass used to scale energetic costs (g)
  
# RETURNS
# A data frame containing the original daily activity data plus activity-specific
# and total daily energy expenditure estimates (kJ).
  
# 1: EXTRACT & CONVERT ENERGETIC PARAMETERS
  
# Extract the mass-specific daily flight cost, convert it to an hourly cost,
# and multiply by the reference body mass (150.95 g) to obtain kJ per hour.
  
flightCoef<-data$c1[1]
flightCoef_kj.hr.g<-flightCoef/24 # kJ.hr.g
flightCoef_kj.hr<-flightCoef_kj.hr.g*150.95 # kJ.hr-1
  
# Extract the intercepts used to describe temperature-dependent energetic costs
# for active and resting birds on water, together with thermal conductance in
# water and air.

activeCoef<-data$Beta_active[1] # kj.hr
restCoef1<-data$Beta_rest[1] # kj.hr
TC_water<-data$TC_water[1] # kJ.hr
TC_air<-data$TC_air[1] # ml O2 hr
  
# Extract the thermoneutral resting-water and land costs and convert them to
# kJ per hour.

restCoef2<-data$c4[1]
restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
landCoef<-data$c3[1]
landCoef_kj<-landCoef*3.6 # kj.hr -> should be the same as the rest coefficient
  

# 2: SCALE ENERGETIC PARAMETERS TO BODY MASS

# Calculate the conversion factor used to scale energetic costs from the
# reference body mass to the focal body mass using an exponent of 0.689.

convf<-1/CostDivider^0.689 # no wing loading
  
# Scale flight cost from its 150.95-g reference body mass to the focal body mass.

flightConstant<-(flightCoef_kj.hr/150.95^0.689)*weightG^0.689

# Apply allometric scaling to active and resting energetic costs using
# CostDivider as the reference body mass.

activeConstant<-(activeCoef*convf)*weightG^0.689
restConstant1<-(restCoef1*convf)*weightG^0.689
restConstant2<-(restCoef2_kj*convf)*weightG^0.689
TC_water_Constant<-(TC_water*convf)*weightG^0.689
landConstant<-(landCoef_kj*convf)*weightG^0.689

# Convert thermal conductance in air from oxygen consumption to kJ and scale
# from its 163.7-g reference body mass to the focal body mass.

TC_air_Constant<-TC_air*20.1/1000*163.7/163.7^0.689*weightG^0.689
  
# 3: DEFINE TEMPERATURE-DEPENDENT ACTIVITY COSTS

# Below the lower critical temperature (LCT), energetic costs on water increase
# as SST decreases. Above the LCT, costs are held at their value at the LCT.

LCT_water<-data$LCT_water[1]

# Calculate the active and resting costs at the water LCT. These provide the
# constant thermoneutral costs used when SST is above the LCT.

restConstant_adjust<-restConstant1 - LCT_water*TC_water_Constant # So that rest is equal to land when sst > LCT
activeConstant_adjust<-activeConstant - LCT_water*TC_water_Constant

# Define the corresponding temperature-dependent relationship for birds on
# land. The intercept is set so that energetic cost at the air LCT equals the
# thermoneutral land cost.

LCT_air<-data$LCT_air[1]
beta_land<-landConstant + LCT_air*TC_air_Constant
  
# 4: CALCULATE DAILY ENERGY EXPENDITURE

# Calculate energetic expenditure separately for each activity. Active and
# resting costs on water increase below the water LCT, while land costs increase
# below the air LCT.
  
energySub2<-data %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT_water, (activeConstant-sst_random*TC_water_Constant)*tActive, activeConstant_adjust*tActive)) %>%
  dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT_water, (activeConstant-sst_random_colony*TC_water_Constant)*tActive, activeConstant_adjust*tActive)) %>%
  dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT_water, (restConstant1 - TC_water_Constant*sst_random)*tRestWater , restConstant_adjust*tRestWater)) %>%
  dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT_water, (restConstant1 - TC_water_Constant*sst_random_colony)*tRestWater , restConstant_adjust*tRestWater)) %>%
  dplyr::mutate(DEEkJ_flight=flightConstant*tFlight) %>%
  dplyr::mutate(DEEkJ_forage=0) %>%
  dplyr::mutate(DEEkJ_restland=ifelse(air_random<=LCT_air, (beta_land - air_random*TC_air_Constant)*tLand, landConstant*tLand)) %>%
  dplyr::mutate(DEEkJ_restland2=landConstant*tLand) %>%
  dplyr::mutate(DEEkJ=DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland) %>%
  dplyr::mutate(DEEkJ_col=DEEkJ_active_col + DEEkJ_rest_col + DEEkJ_flight + DEEkJ_restland) %>%
  dplyr::mutate(weight=weightG)

# 5: PREPARE OUTPUT

return(energySub2)
  
  
}

calculateEnergetics_LiA_daily_map<-function(data, CostDivider, weightG, sstVals) {
  
   # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj.hr.g<-flightCoef/24 # kJ.hr.g
  flightCoef_kj.hr<-flightCoef_kj.hr.g*150.95 # kJ.hr-1
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC_water<-data$TC_water[1] # kJ.hr
  TC_air<-data$TC_air[1] # ml O2 hr
  
  # Rest & land coefs
  restCoef2<-data$c4[1]
  restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
  landCoef<-data$c3[1]
  landCoef_kj<-landCoef*3.6 # kj.hr -> should be the same as the rest coefficient
  
  # Active coef when thermoneutral (from kyle's paper)
  #activeCoef_neut<-data$c5[1]
  #activeCoef_neut_kj<-activeCoef_neut*3.6 # kJ.hr
  
  #if (activeCoef_neut_kj > activeCoef) stop (print("Error with active coefs")) # Because the coef when thermoneutral must be lower...
  
  # Here is a conversion factor to transform these to g.kj which incorporates allometric scaling
  convf<-1/CostDivider^0.689 # no wing loading
  
  # Account for change in constants
  flightConstant<-(flightCoef_kj.hr/150.95^0.689)*weightG^0.689
  activeConstant<-(activeCoef*convf)*weightG^0.689
  #activeConstant2<-(activeCoef_neut_kj*convf)*weightG^0.689
  restConstant1<-(restCoef1*convf)*weightG^0.689
  restConstant2<-(restCoef2_kj*convf)*weightG^0.689
  #forageConstant<-(3.64*convother)*weightG^0.689
  TC_water_Constant<-(TC_water*convf)*weightG^0.689
  TC_air_Constant<-TC_air*20.1/1000*163.7/163.7^0.689*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT_water<-data$LCT_water[1]
  restConstant_adjust<-restConstant1 - LCT_water*TC_water_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT_water*TC_water_Constant
  
  LCT_air<-data$LCT_air[1]
  beta_land<-landConstant + LCT_air*TC_air_Constant
  
  # Turn SST raster into a data frame
  sst<-subset(sstVals, 1)
  temp<-subset(sstVals, 2)
  sstDf<-as.data.frame(sst, xy=TRUE)
  airDf<-as.data.frame(temp, xy=TRUE)
  colnames(sstDf)<-c("x", "y", "sst")
  colnames(airDf)<-c("x", "y", "temp")
  sstDf$temp<-airDf$temp
  
  # Calculate energy for every cell according to sst in that cell
  sstDf$DEEkJ_active=ifelse(sstDf$sst <=LCT_water, (activeConstant-sstDf$sst*TC_water_Constant)*data$tActive_month, activeConstant_adjust*data$tActive_month)
  sstDf$DEEkJ_rest=ifelse(sstDf$sst <= LCT_water, (restConstant1 - TC_water_Constant*sstDf$sst)*data$tRestWater_month , restConstant_adjust*data$tRestWater_month)
  sstDf$DEEkJ_flight=flightConstant*data$tFlight_month
  sstDf$DEEkJ_forage=0
  sstDf$DEEkJ_restland=ifelse(sstDf$temp <= LCT_air, (beta_land - sstDf$temp*TC_air_Constant)*data$tLand_month, landConstant*data$tLand_month)
  sstDf$DEEkJ=sstDf$DEEkJ_rest + sstDf$DEEkJ_flight + sstDf$DEEkJ_active + sstDf$DEEkJ_restland
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, sst, temp, DEEkJ) 
  
  return(sstDf_final)
  
}

##### Atlantic puffin #####

calculateEnergetics_AP_daily<-function(data, CostDivider, weightG) {

# PURPOSE: Calculate daily energy expenditure for Atlantic puffins from
# activity budgets, body mass, and environmental temperature.
  
# Broad workflow:
# 1. Extract activity-specific energetic and thermoregulatory parameters.
# 2. Convert energetic coefficients to the units required by the model.
# 3. Scale energetic parameters from reference body masses to the focal body mass.
# 4. Define temperature-dependent energetic costs below lower critical temperatures.
# 5. Calculate energetic costs separately for each activity.
# 6. Sum activity-specific costs to estimate total daily energy expenditure.
  
# INPUTS
# data  - data frame containing daily activity budgets, environmental
# conditions, and energetic parameters
# CostDivider - reference body mass used for allometric scaling (g)
# weightG - body mass used to scale energetic costs (g)
  
# RETURNS
# A data frame containing the original daily activity data plus activity-specific
# and total daily energy expenditure estimates (kJ).
  
# 1: EXTRACT & CONVERT ENERGETIC PARAMETERS
  
# Extract the energetic cost of flight and convert to kJ per hour.

flightCoef<-data$c1[1]
flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
# Extract the intercepts used to describe temperature-dependent energetic costs
# for active and resting birds on water, together with thermal conductance in water.

activeCoef<-data$Beta_active[1] # kj.hr
restCoef1<-data$Beta_rest[1] # kj.hr
TC_water<-data$TC_water[1] # kJ.hr
  
# Extract thermal conductance in air and the lower critical temperature for
# birds on land.

TC_air<-data$TC_air[1] # ml O2 hr
LCT_air<-data$LCT_air[1] # Degrees c
  
# Extract the thermoneutral resting-water and land costs and convert them to
# kJ per hour.

restCoef2<-data$c4[1]
restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
landCoef<-data$c3[1]
landCoef_kj<-landCoef*3.6 # kj.hr -> should be the same as the rest coefficient
  
# 2: SCALE ENERGETIC PARAMETERS TO BODY MASS

# Calculate the conversion factor used to scale energetic costs from the
# reference body mass to the focal body mass using an exponent of 0.689.

  convf<-1/CostDivider^0.689 # no wing loading
  
# Apply allometric scaling to activity-specific costs and water
# thermoregulatory parameters.
  
flightConstant<-(flightCoef_kj*convf)*weightG^0.689
activeConstant<-(activeCoef*convf)*weightG^0.689
restConstant1<-(restCoef1*convf)*weightG^0.689
restConstant2<-(restCoef2_kj*convf)*weightG^0.689
TC_water_Constant<-(TC_water*convf)*weightG^0.689
landConstant<-(landCoef_kj*convf)*weightG^0.689

# Convert thermal conductance in air from oxygen consumption to kJ and scale
# from its 819.3-g reference body mass to the focal body mass.

TC_air_Constant<-(TC_air*20.1/1000)*819.3/(819.3^0.689)*weightG^0.689
  
# 3: DEFINE TEMPERATURE-DEPENDENT ACTIVITY COSTS

# Below the lower critical temperature (LCT), energetic costs on water increase
# as SST decreases. Above the LCT, costs are held at their value at the LCT.

LCT_water<-data$LCT_water[1]

# Calculate the active and resting costs at the water LCT. These provide the
# constant thermoneutral costs used when SST is above the LCT.

restConstant_adjust<-restConstant1 - LCT_water*TC_water_Constant # So that rest is equal to land when sst > LCT
activeConstant_adjust<-activeConstant - LCT_water*TC_water_Constant
  
# Define the temperature-dependent relationship for birds on land. The
# intercept is set so that energetic cost at the air LCT equals the
# thermoneutral land cost.

beta_land<-landConstant + LCT_air*TC_air_Constant

# 4: CALCULATE DAILY ENERGY EXPENDITURE

# Calculate energetic expenditure separately for each activity. Active and
# resting costs on water increase below the water LCT, while land costs increase
# below the air LCT.
  
energySub2<-data %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT_water, (activeConstant-sst_random*TC_water_Constant)*tActive, activeConstant_adjust*tActive)) %>%
  dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT_water, (activeConstant-sst_random_colony*TC_water_Constant)*tActive, activeConstant_adjust*tActive)) %>%
  dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT_water, (restConstant1 - TC_water_Constant*sst_random)*tRestWater , restConstant_adjust*tRestWater)) %>%
  dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT_water, (restConstant1 - TC_water_Constant*sst_random_colony)*tRestWater , restConstant_adjust*tRestWater)) %>%
  dplyr::mutate(DEEkJ_flight=flightConstant*tFlight) %>%
  dplyr::mutate(DEEkJ_forage=0) %>%
  dplyr::mutate(DEEkJ_restland=ifelse(air_random <=LCT_air, (beta_land - TC_air_Constant*air_random)*tLand ,landConstant*tLand)) %>%
  dplyr::mutate(DEEkJ_restland2=landConstant*tLand) %>%
  dplyr::mutate(DEEkJ=DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland) %>%
  dplyr::mutate(DEEkJ_col=DEEkJ_active_col + DEEkJ_rest_col + DEEkJ_flight + DEEkJ_restland) %>%
  dplyr::mutate(weight=weightG)

# 5: PREPARE OUTPUT

  return(energySub2)
  
}

calculateEnergetics_AP_daily_map<-function(data,  CostDivider, weightG, sstVals) {
  
# First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC_water<-data$TC_water[1] # kJ.hr
  
  # Extract TC on Land
  TC_air<-data$TC_air[1] # ml O2 hr
  LCT_air<-data$LCT_air[1] # Degrees c
  
  # Rest & land coefs
  restCoef2<-data$c4[1]
  restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
  landCoef<-data$c3[1]
  landCoef_kj<-landCoef*3.6 # kj.hr -> should be the same as the rest coefficient
  
  # Active coef when thermoneutral (from kyle's paper)
  #activeCoef_neut<-data$c5[1]
  #activeCoef_neut_kj<-activeCoef_neut*3.6 # kJ.hr
  
  #if (activeCoef_neut_kj > activeCoef) stop (print("Error with active coefs")) # Because the coef when thermoneutral must be lower...
  
  # Here is a conversion factor to transform these to g.kj which incorporates allometric scaling
  convf<-1/CostDivider^0.689 # no wing loading
  
  # Account for change in constants
  flightConstant<-(flightCoef_kj*convf)*weightG^0.689
  activeConstant<-(activeCoef*convf)*weightG^0.689
  #activeConstant2<-(activeCoef_neut_kj*convf)*weightG^0.689
  restConstant1<-(restCoef1*convf)*weightG^0.689
  restConstant2<-(restCoef2_kj*convf)*weightG^0.689
  TC_water_Constant<-(TC_water*convf)*weightG^0.689
  TC_air_Constant<-(TC_air*20.1/1000)*819.3/(819.3^0.689)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT_water<-data$LCT_water[1]
  restConstant_adjust<-restConstant1 - LCT_water*TC_water_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT_water*TC_water_Constant
  
  # Calculate a beta_land
  beta_land<-landConstant + LCT_air*TC_air_Constant
  
  # Turn SST raster into a data frame
  sst<-subset(sstVals, 1)
  temp<-subset(sstVals, 2)
  sstDf<-as.data.frame(sst, xy=TRUE)
  airDf<-as.data.frame(temp, xy=TRUE)
  colnames(sstDf)<-c("x", "y", "sst")
  colnames(airDf)<-c("x", "y", "temp")
  sstDf$temp<-airDf$temp
  
  # Calculate energy for every cell according to sst in that cell
  sstDf$DEEkJ_active=ifelse(sstDf$sst <=LCT_water, (activeConstant-sstDf$sst*TC_water_Constant)*data$tActive_month, activeConstant_adjust*data$tActive_month)
  sstDf$DEEkJ_rest=ifelse(sstDf$sst <= LCT_water, (restConstant1 - TC_water_Constant*sstDf$sst)*data$tRestWater_month , restConstant_adjust*data$tRestWater_month)
  sstDf$DEEkJ_flight=flightConstant*data$tFlight_month
  sstDf$DEEkJ_forage=0
  sstDf$DEEkJ_restland=ifelse(sstDf$temp <= LCT_air, (beta_land - sstDf$temp*TC_air_Constant)*data$tLand_month, landConstant*data$tLand_month)
  sstDf$DEEkJ=sstDf$DEEkJ_rest + sstDf$DEEkJ_flight + sstDf$DEEkJ_active + sstDf$DEEkJ_restland
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, sst, temp, DEEkJ) 
  
  return(sstDf_final)
  
  
}
