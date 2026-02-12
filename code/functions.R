# Various functions.... # 

#### Splitting light data into years tracked ####

# years is number of unique year of light data in dataset
# lux.id is a session_id of light data
# lux.id is a session_id of light data

splitData_yearTracked<-function(lux.id) {
  
  lux.id.new<-lux.id %>%
    dplyr::mutate(month=as.numeric(substr(date_time, 6, 7))) %>%
    dplyr::mutate(day=as.numeric(substr(date_time, 9, 10))) %>%
    dplyr::mutate(track_year=ifelse(month<=5 & day<=31, paste0(year-1, "_", substr(year, 3, 4)), paste0(year, "_", substr(year+1, 3, 4))))
  
  return(lux.id.new)
  
}

#### Activity functions ####

# Function to calculate time spent in different activities
# Data is immersion data standardized to 0-1 (with 1= 100% wet and 0 = 100% dry)
# Function returns a list with two objects:
# Object 1 is the dataset annotated with different behaviours for exploratory use
# Object 2 is time spent in activity summarized per day
# Method can be equal to Lila or Caitlin for the auks... 

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

# Loop through these 

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

##### Kittiwake #####

# At the moment I am just applying Don-Jean's equation from 
# https://europepmc.org/article/ppr/ppr648557
# Max speed is 90 km.hr-1 according to the NINA geolocation processing algorithm
# Flight bout threshold is used to make a cut-off between what is likely flight & what is likely colony attendance
# Iterations is how many reps per individual we calculate

# The one below is in prep
calculateTimeInActivity_BLK<-function(data, irmaData){
  
  # Here we assign three behaviours: flight, forage & rest   
  dataCalc<-data %>%
    dplyr::filter(!is.na(col_lon)) %>%
    dplyr::ungroup() %>%
    rename(new_cond=new.cond) %>%
    dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
    dplyr::mutate(Activity=ifelse(new_cond<Th1, "Forage", "RestWater")) %>%
    dplyr::mutate(Activity=ifelse(new_cond<=Th2, "Dry", Activity)) %>%
    dplyr::mutate(Activity=ifelse(Th2==0 & new_cond==0, "Dry", Activity)) %>%
    dplyr::mutate(MaxDistColKm=max(distColonyKm)) 
  
  # Add value of next DistColKm - as we assume a bird could arrive at the colony if its next location is within 250 km of the colony 
  distances_next<-dataCalc %>%
    dplyr::group_by(distColonyKm) %>%
    dplyr::slice(1) %>%
    arrange(date_time) %>%
    dplyr::select(date_time, distColonyKm) %>%
    ungroup() %>%
    dplyr::mutate(distColonyKm_next=lead(distColonyKm)) %>%
    dplyr::mutate(distColonyKm_next=ifelse(is.na(distColonyKm_next), distColonyKm, distColonyKm_next)) %>%
    dplyr::select(-distColonyKm) %>%
	dplyr::mutate(date_time=as.character(date_time))
  
  # Now we assign bout numbers to dry bouts & determine whether it's in flight or on land 
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
	
 # Make sure no NAs in the date time
 nas_date<-subset(FlightBouts, is.na(date_time))
 
 if(nrow(nas_date)>0) {stop(print("Error: nas in date time")) }
  
  # Here we join the numbering & estimate the length of different flight bouts 
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
  
  # Scenario 1: Here we assume same rules for daylight and darkness
  
  activityAdjust1<-FlightBoutLengths %>%
    ungroup() %>%
    dplyr::full_join(distances_next, by=c("date_time")) %>%
    arrange(date_time) %>%
    fill(distColonyKm_next, .direction=c("down")) %>%
    dplyr::group_by(BoutNo) %>%
    # Determine distance to land & ice concentration
    dplyr::mutate(ice_random=ifelse(!is.na(BoutNo), ice_mean, NA)) %>%
    dplyr::mutate(ice_random=ifelse(!is.na(BoutNo) & ice_random<0, 0, ice_random)) %>%
    dplyr::mutate(ice_random=ifelse(!is.na(BoutNo) & ice_random>1, 1, ice_random)) %>%
    dplyr::mutate(distCoast_random=ifelse(!is.na(BoutNo), rnorm(mean=distCoastKm_mean, sd=distCoastKm_sd, n=n_distinct(BoutNo)),NA)) %>%
    dplyr::mutate(distCoast_random=ifelse(!is.na(BoutNo) & distCoast_random<0, 0, distCoast_random)) %>%
    dplyr::mutate(PossLand=ifelse( ice_random >0 | distColonyKm <= dist_colony | distCoast_random <= dist_colony | distColonyKm_next <= dist_colony, 1, 0)) %>%
    #dplyr::mutate(NewActivity=ifelse(Period=="Darkness" & Activity=="Dry" & PossLand==0, "Flight", NA)) %>%
    #dplyr::mutate(NewActivity=ifelse(Period=="Darkness" & Activity=="Dry" & PossLand==1, "Land", NewActivity)) %>%
	dplyr::mutate(NewActivity=ifelse(Activity=="Dry" & flightLengthMins > L1, "Land", NA)) %>%
	dplyr::mutate(NewActivity=ifelse(Activity=="Dry" & flightLengthMins <= L1, "Flight", NewActivity)) %>%
    dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity))
  
  # Check for remaining dry bouts & stop if there are some as error
  dryBouts<-subset(activityAdjust1, Activity=="Dry")
  
  if (nrow(dryBouts)>0) {
    stop(print("Error: remaining dry bouts"))
  }
  
  # Set sequence to resample from #
  uniqueValues<-unique(data$L1_colony_max[1], data$L1_colony_min[1])
  if (length(uniqueValues)>1){
 
  # Then I will re-allocate some land bouts to flight
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
    dplyr::mutate(LagMins=ifelse(Activity=="Land" & !PrevFlight %in% c(1) & firstLand==1, sample(c(seq(data$L1_colony_min[1], data$L1_colony_max[2], 10)), replace=TRUE), 0)) %>% # Here I determine a random number of 10-minute bouts to re-allocate from land to flight 
    dplyr::mutate(LagMinsAdj=ifelse(LagMins>=DurationLandMins, DurationLandMins-10, LagMins)) %>% # and here I make sure they are not longer than the actual land bout
    dplyr::group_by(LandBoutNo) %>%
    dplyr::mutate(LandBoutRow=ifelse(Activity=="Land", row_number()*10, 0)) %>% # Here i make a crazy system to re-allocate a certain number of rows...
    replace_na((list(LagMinsAdj=0))) %>%
    dplyr::mutate(NewActivity=ifelse(Activity=="Land" & first(LagMinsAdj)>0 & LandBoutRow<=first(LagMinsAdj) & !is.na(LandBoutRow) & first(PrevFlight) %in% c(0), "Flight", NA)) %>%
    dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity)) 
	
	} else {
	
	# Just for the sensitivity analysis #
	
	# Then I will re-allocate some land bouts to flight
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
	
	
  
  # Re-calculate final flight bout lengths 
  
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
  
  # Here we join the numbering & estimate the length of different flight bouts 
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
		dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tRestWater=sum(RestWater)*10/60, tLand=sum(Land)*10/60,
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
  
  # Calculate time spent per day doing different things
  dataCalcDay<-FlightBoutLengths_final %>%
    dplyr::ungroup() %>%
    dplyr::mutate(date=substr(date_time, 1, 10)) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
    dplyr::filter(Duration>=1430) %>% # Only keep full days %>%
    dplyr::group_by(species, colony, session_id, date) %>%
    dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                  Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0), sstRandom=sst_random_start) %>%
    dplyr::mutate(sstRandom=ifelse(sstRandom < -1.9, 1.9, sstRandom)) %>%
    dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tRestWater=sum(RestWater)*10/60, tLand=sum(Land)*10/60,
                     tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, MaxDistColKm=max(MaxDistColKm),
                     sst_random=mean(sstRandom), ice_random=mean(ice_random, na.rm=TRUE), distColonyKm_mean=mean(distColonyKm, na.rm=TRUE), mean.lon=mean(lon), mean.lat=mean(lat)) %>%
    ungroup() %>%
    dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
    dplyr::mutate(dayLengthHrs=tDaylight) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(DurationTot=sum(tForage, tRestWater, tLand, tFlight)) %>%
    dplyr::left_join(daily_max_flightBout, by=c("date")) %>%
	dplyr::left_join(dataCalcDay_period, by=c("date"))
  
  # Save results from Both
  actResults<-dataCalcDay
  boutResults<-dataCalcDay_period
  
  allResults<-list(actResults, boutResults)
  return(allResults)
  
}

##### Northern fulmar #####

calculateTimeInActivity_NF<-function(data, irmaData){
  
  # Set up parameter space to choose from  
  data$date_time<-as.character(data$date_time)
  
  # Here we assign three behaviours: flight, forage & rest   
  dataCalc<-data %>%
    dplyr::filter(!is.na(col_lon)) %>%
    dplyr::ungroup() %>%
    rename(new_cond=new.cond) %>%
    dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
    dplyr::mutate(Activity=ifelse(new_cond<Th1, "Forage", "RestWater")) %>%
    dplyr::mutate(Activity=ifelse(new_cond<Th2, "Dry", Activity)) %>%
    dplyr::mutate(Activity=ifelse(Th2==0 & new_cond==0, "Dry", Activity)) %>%
    dplyr::mutate(MaxDistColKm=max(distColonyKm)) 
  
  # Add value of next DistColKm - as we assume a bird could arrive at the colony if its next location is within 250 km of the colony 
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
  
  # Now we assign bout numbers to dry bouts & determine whether it's in flight or on land 
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
  
  # Here we join the numbering & estimate the length of different flight bouts 
  FlightBoutLengths<-dataCalc %>%
  #dplyr::mutate(date_time=as.character(date_time)) %>%
  #dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
    dplyr::left_join(FlightBouts, by=c("date_time")) %>%
    dplyr::group_by(Activity) %>%
    fill(BoutNo, .direction=c("down")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(BoutNo, distColonyKm) %>%
    dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
    dplyr::mutate(flightLengthMins=ifelse(Activity=="Dry", flightLengthMins, NA)) %>%
    ungroup()
  
  # Scenarios: No matter whether its day or dark, allocation is based on bout length & possibility of being on land
  
  activityAdjust1<-FlightBoutLengths %>%
    ungroup() %>%
	dplyr::mutate(date_characters=nchar(date_time)) %>%
	dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
	dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
    dplyr::full_join(distances_next, by=c("date_time")) %>%
    arrange(date_time) %>%
    fill(distColonyKm_next, .direction=c("down")) %>%
    dplyr::group_by(BoutNo, distColonyKm) %>%
    dplyr::mutate(ice_random=ifelse(!is.na(BoutNo), first(ice_random), NA)) %>%
    dplyr::mutate(ice_random=ifelse(ice_random<0, 0, ice_random)) %>%
    dplyr::mutate(ice_random=ifelse(ice_random>1, 1, ice_random)) %>%
    dplyr::mutate(distCoast_random=ifelse(!is.na(BoutNo), rnorm(mean=distCoastKm_mean, sd=distCoastKm_sd, n=n_distinct(BoutNo)),NA)) %>%
    dplyr::mutate(distCoast_random=ifelse(distCoast_random<0, 0, distCoast_random)) %>%
    dplyr::mutate(PossLand=ifelse(ice_random >0 | distColonyKm <= dist_colony | distColonyKm_next <=dist_colony, 1, 0)) %>%
    dplyr::mutate(p=runif(n=n_distinct(BoutNo))) %>%
    dplyr::mutate(p2=runif(n=n_distinct(BoutNo))) %>%
    dplyr::mutate(pLand=ifelse(p <= pLand_prob | p2 <= ice_random, 1, 0)) %>%
	replace_na(list(flightLengthMins=0)) %>%
    dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins >= L1 & PossLand==0, "Flight", NA)) %>% # If no land in sight, then it has to be flight
    dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins >= L1 & PossLand==1, "Land", NewActivity)) %>% # otherwise it's land
    dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==0, "Flight", NewActivity)) %>% # If it's short & no land then it's flight
    dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm <= dist_colony & pLand==1, "Land", NewActivity)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
    dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm <= dist_colony & pLand==0 , "Flight", NewActivity)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
    dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm > dist_colony & pLand==1, "Land", NewActivity)) %>% # In instances where bird could be on land or ice, then 50% prob of land overrules if ice concentration < 50%
    dplyr::mutate(NewActivity=ifelse( Activity=="Dry" & flightLengthMins < L1 & PossLand==1 & distColonyKm > dist_colony & pLand==0 , "Flight", NewActivity)) %>%
	# Create new column where I count how many rows this concerns
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
  
  # Determine how many unique values of L1_colony
  uniqueVals<-unique(data$L1_colony_min[1], data$L1_colony_max[1])
  if (length(uniqueVals)>1) {
	
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
    dplyr::mutate(LagMins=ifelse(Activity=="Land" & !PrevFlight %in% c(1) & firstLand==1, sample(seq(data$L1_colony_min[1], data$L1_colony_max[1], 10), replace=TRUE), 0)) %>% # Here I determine a random number of 10-minute bouts to re-allocate from land to flight 
    dplyr::mutate(LagMinsAdj=ifelse(LagMins>=DurationLandMins, DurationLandMins-10, LagMins)) %>% # and here I make sure they are not longer than the actual land bout
    dplyr::group_by(LandBoutNo) %>%
    dplyr::mutate(LandBoutRow=ifelse(Activity=="Land", row_number()*10, 0)) %>% # Here i make a crazy system to re-allocate a certain number of rows...
    replace_na((list(LagMinsAdj=0))) %>%
    dplyr::mutate(NewActivity=ifelse(Activity=="Land" & first(LagMinsAdj)>0 & LandBoutRow<=first(LagMinsAdj) & !is.na(LandBoutRow) & first(PrevFlight) %in% c(0), "Flight", NA)) %>%
    dplyr::mutate(Activity=ifelse(!is.na(NewActivity), NewActivity, Activity))
	
	
	} else {
	
	# Specific to sensitivity analysis
	
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
	
	# Make a new way of estimating bout no
	
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
	
	# We re-calculate flight length without varying distance to colony. Now if we have flight lengths which are crazy long, we change those which are within distance to land to land
	
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
	
  # Calculate time spent doing different activities per day so I can adjust surface-foraging
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
	  
		  # Attach max flight bout length
		  daily_max_flightBout<-FlightBoutLengths_final %>%
			dplyr::ungroup() %>%
			dplyr::mutate(date=substr(date_time, 1, 10)) %>%
			dplyr::filter(Activity=="Flight") %>%
			dplyr::group_by(date, BoutNo) %>%
			dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
			ungroup() %>%
			dplyr::group_by(date) %>%
			dplyr::summarise(maxFlightBoutsMins=max(flightLengthMins))
			
		  daily_max_flightBout_dark<- dataCalcDay_period%>%
			dplyr::ungroup() %>%
			#dplyr::filter(Activity=="Flight") %>%
			dplyr::group_by(date) %>%
			dplyr::summarise(maxFlightBoutsMins_dark=max(flightTimeMins_dark)) %>%
			ungroup() 
	  
	  # Calculate time spent per day doing different things
	  dataCalcDay<-FlightBoutLengths_final %>%
		dplyr::ungroup() %>%
		dplyr::mutate(date=substr(date_time, 1, 10)) %>%
		dplyr::group_by(date) %>%
		dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
		dplyr::filter(Duration>=1430) %>%  # Only keep full days %>%
		dplyr::group_by(species, colony, session_id, date) %>%
		dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
					  Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0), sstRandom=sst_random_start) %>%
		dplyr::mutate(sstRandom=ifelse(sstRandom < -1.9, 1.9, sstRandom)) %>%
		dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tRestWater=sum(RestWater)*10/60, tLand=sum(Land)*10/60,
						 tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, MaxDistColKm=max(MaxDistColKm),
						 sst_random=mean(sstRandom), ice_random=mean(ice_random, na.rm=TRUE), distColonyKm_mean=mean(distColonyKm, na.rm=TRUE), mean.lon=mean(lon),
						 mean.lat=mean(lat), boutsRandom=(sum(RandomAllocation)/144)) %>%
		ungroup() %>%
		dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
		#dplyr::filter(Breeding==0) %>% # We only include full day 
		dplyr::mutate(dayLengthHrs=tDaylight) %>%
		dplyr::group_by(date) %>%
		dplyr::mutate(DurationTot=sum(tForage, tRestWater, tLand, tFlight)) %>%
		dplyr::left_join(daily_max_flightBout, by=c("date")) %>%
		dplyr::left_join(dataCalcDay_period, by=c("date")) 
	
  
  # Create final dataset which has prop daylight/twilight etc doign various things
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
  
  # Now I will extract individual flight bouts for fulmars only to save for supplementary materials analysis (distribution of flight bout lengths during all times
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
  
  
  # Save results from Both
  actResults<-dataCalcDay
  boutResults<-dataCalcDay_period2
  
  allResults<-list(actResults, boutResults)
  return(allResults)
  
  
}

##### Methods for activity budgets - auks ####

methodCaitlin<-function(data, irmaData, speciesLatin) {
  
  # Here we assign day of year & filter to non-breeding period 
  dataCalc<-data %>%
    dplyr::ungroup() %>%
	#dplyr::group_by(session_id) %>%
    dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
    dplyr::mutate(MaxDistColKm=max(distColonyKm)) %>%
	arrange(date) %>%
	ungroup()  
  
  # Add value of next DistColKm - as we assume a bird could arrive at the colony if its next location is within 250 km of the colony 
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
  
  # Here we assign bouts to dry or other
  dataCalc$Activity<-ifelse(dataCalc$new.cond<=data$Th2[1], "Dry", "Other")
  
  # Now we allocate dry bouts to TFlight, TLand or TRest
  
  # First we assign bout numbers to dry bouts  
  FlightBouts<-dataCalc %>%
    dplyr::filter(Activity=="Dry") %>%
    ungroup() %>%
	dplyr::mutate(date_characters=nchar(date_time)) %>%
	dplyr::mutate(date_time=ifelse(date_characters<19, paste(date_time, "00:00:00", sep=" "), date_time)) %>%
	dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
    arrange(individ_id, session_id, date_time) %>%
	#dplyr::group_by(session_id) %>%
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
  
  # Here we join the numbering & estimate the length of different flight bouts 
  FlightBoutLengths<-dataCalc %>%
	dplyr::mutate(date_time=as.character(date_time)) %>%
  #dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
    dplyr::left_join(FlightBouts, by=c("date_time")) %>%
    dplyr::group_by(Activity) %>%
    fill(BoutNo, .direction=c("down")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(BoutNo) %>%
    dplyr::mutate(flightLengthMins=n_distinct(date_time)*10) %>%
    dplyr::mutate(flightLengthMins=ifelse(Activity=="Dry", flightLengthMins, 0)) %>%
    ungroup()
  
  # Scenario 1: Darkness -> the bird can either be on land if it's dry the whole night or resting on the water (leg-tucked) based on whether it's possible to be on land or not
  
  # Distance to sea-ice is only applicable to little auks & Brunnich's guillemots 
  FlightBoutLengths$use_seaice<-ifelse(dataCalc$species[1] %in% c("Little auk", "Brünnich's guillemot"), 1, 0)
  
  # Determine whether the bird is dry the whole night or not - first we number every night 
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
  
  # Next we determine the proportion of every night that is dry to identify periods which are entirely dry 
  darkness_propDry<-dataCalc %>%
  dplyr::mutate(date_time=as.character(date_time)) %>%
   #dplyr::mutate(date_time=as.POSIXct(date_time, format=c("%Y-%m-%d %H:%M:%S"), tz="UTC")) %>%
    dplyr::left_join(darkness_boutNo, by=c("date_time")) %>%
    dplyr::group_by(Period) %>%
    fill(BoutNo_dark, .direction=c("down")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Period, BoutNo_dark) %>%
    dplyr::mutate(totDry=sum(1-new.cond)*10, darknessLength=n_distinct(date_time)*10) %>%
    dplyr::mutate(propDarkDry=ifelse(Period=="Darkness", totDry/darknessLength, NA)) %>%
    ungroup() %>%
	dplyr::select(date_time, BoutNo_dark, propDarkDry)
  
  # Now we determine whether a bout during darkness if TLand or TRest based on these rules
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
    dplyr::mutate(distCoast_random=ifelse(!is.na(BoutNo), rnorm(mean=distCoastKm_mean, sd=distCoastKm_sd, n=n_distinct(BoutNo)),NA)) %>%
    dplyr::mutate(distCoast_random=ifelse(!is.na(BoutNo) & distCoast_random<0, 0, distCoast_random)) %>%
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
  
  # Scenario 2: It's daylight & bout length is long or short (these can be allocated to tFlight, TLand or TRest))
  
  # First we re-number the dry bouts that have yet to be allocated 
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
  
  # Then we re-calculate their duration
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
  
  activityAdjust2<-FlightBoutLengths2 %>%
    dplyr::group_by(BoutNo, distColonyKm) %>%
    # Determine distance to land & ice concentration
    dplyr::mutate(ice_random=ifelse(!is.na(BoutNo), ice_mean, NA)) %>%
    dplyr::mutate(ice_random=ifelse(ice_random<0, 0, ice_random)) %>%
    dplyr::mutate(ice_random=ifelse(ice_random>1, 1, ice_random)) %>%
	dplyr::mutate(ice_random=ifelse(use_seaice==1, ice_random, 0)) %>% # Change this to zero if species is incorrect
    dplyr::mutate(distCoast_random=ifelse(!is.na(BoutNo), rnorm(mean=distCoastKm_mean, sd=distCoastKm_sd, n=n_distinct(BoutNo)),NA)) %>%
    dplyr::mutate(distCoast_random=ifelse(distCoast_random<0, 0, distCoast_random)) %>%
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
  
  # Assign tActive & tInactive
  dry_final<-activityAdjust2
  dry_final$Activity<-ifelse(dry_final$new.cond<data$Th1[1] & dry_final$new.cond>0, "RestWater",dry_final$Activity)
  dry_final$Activity<-ifelse(dry_final$new.cond>=data$Th1[1], "Active", dry_final$Activity)
  
  # Add time in flight for time on land not preceded by tflight
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
  
  # Basically what happens below is that the start of a 'Land' section is re-adjusted to flight if
  # it is not preceded by flight. I'm allowing a change of up to flight bout threshold length at random
  # if the duration of the land attendance allows it 
  
  # Determine how many unique values of L1_colony
  uniqueValues<-unique(c(data$L1_colony_min, data$L1_colony_max))
  
  if (length(uniqueValues)>1) {
  
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
	
	# Sensitivity analysis Scenario
	
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
  
  # Now we check whether there are flight bouts which are too long as they will have to be re-classified to land or rest (with leg-tucked)
  
  ### Start by re-classfiying as land if possible ###
  
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
  
  # Here we join the numbering & estimate the length of different flight bouts 
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
	
  # If there are bouts which are too long and some of it could be on land -> then we re-assign to land
  FlightBoutLengths_lastcheck_reclassify<-FlightBoutLengths_lastcheck %>%
  ungroup() %>%
  dplyr::mutate(Activity=ifelse(Activity=="Flight" & flightLengthMins > L1 & PossLand==1 & Period=="Daylight", "Land", Activity)) 
  
  ### Now we re-classify as leg-tucked if possible ###
  
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
  
  # Here we join the numbering & estimate the length of different flight bouts 
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
	
  # If there are bouts which are too long and some of it could be on land -> then we re-assign to land
  FlightBoutLengths_lastcheck_reclassify2<-FlightBoutLengths_lastcheck2 %>%
  ungroup() %>%
  dplyr::group_by(BoutNo) %>%
  dplyr::mutate(index=row_number()) %>%
  dplyr::mutate(Activity_change=ifelse(Activity=="Flight" & flightLengthMins > L1, 1, 0)) %>% # Determine whether parts of a bout need to be re-assigned
  dplyr::mutate(Activity_change_duration=ifelse(Activity_change==1, (flightLengthMins-L1)/10, 0)) %>% # Calculate how many rows need to be changed
  dplyr::mutate(Activity=ifelse(Activity_change==1 & index <=ceiling(Activity_change_duration), "RestWater", Activity )) %>%
  dplyr::select(-c(Activity_change, Activity_change_duration)) 
  
  # No we do a final calculation of flight time & make an error if there are still bouts which are too long
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
  
  # Here we join the numbering & estimate the length of different flight bouts 
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
  
  # Attach max flight bout length
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
  
  ## Summarize per day time in activity & re-adjust time active using a coefficient varying between 1.5-3 (c)
  
  dataCalcDay<-FlightBoutLengths_lastcheck3 %>%
    dplyr::ungroup() %>%
    dplyr::mutate(date=substr(date_time, 1, 10)) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(Duration=n_distinct(date_time)*10) %>%
    dplyr::filter(Duration>=1430) %>% # This allows for full days including those which have a weird thing where they finish at midnight...
    dplyr::group_by(species, colony, session_id, date) %>%
    dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Active=ifelse(Activity=="Active", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                  Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0), sstRandom=sst_random_start) %>%
    dplyr::mutate(sstRandom=ifelse(sstRandom < -1.9, 1.9, sstRandom)) %>%
    dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater1=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tActive=sum(Active)*10/60, tLand=sum(Land)*10/60,
                     tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, MaxDistColKm=max(MaxDistColKm),
                     sst_random=mean(sstRandom), ice_random=mean(ice_random, na.rm=TRUE), distColonyKm_mean=mean(distColonyKm, na.rm=TRUE), mean.lon=mean(lon), mean.lat=mean(lat),
					 boutsRandom=(sum(RandomAllocation)/144)) %>%
    ungroup() %>%
    dplyr::mutate(doy=floor(as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01"))), unit=c("days"))) + 1) %>%
    dplyr::group_by(date) %>%
	dplyr::mutate(tRestWater2=tRestWater1*data$c[1], amountSub=tRestWater2-tRestWater1) %>%
	dplyr::mutate(tActive2=ifelse(amountSub > tActive, 0, tActive-amountSub), tRestWater2=ifelse(amountSub>tActive, tRestWater1 + tActive, tRestWater2)) %>%
	#dplyr::mutate(tActive2=ifelse(amountSub2 > tActive, 0, tActive-amountSub2), tForage2=ifelse(amountSub2 > tForage, 0, tForage-amountSub2), tRestWater2=ifelse(amountSub> (tActive + tForage), tRestWater1 + tActive + tForage, tRestWater2)) %>%
    #dplyr::mutate(tRestWater2=ifelse(tRestWater1*c>tActive, tRestWater1 + tActive, tRestWater1*c), tActive2=tActive-(tRestWater2 - tRestWater1)) %>%
    dplyr::select(-c(tActive, tRestWater1)) %>%
    dplyr::rename(tActive=tActive2, tRestWater=tRestWater2) %>%
    dplyr::mutate(tActive=ifelse(tActive<0, 0, tActive)) %>%
    dplyr::mutate(DurationTot=sum(tFlight, tForage, tActive, tRestWater, tLand)) %>%
    dplyr::mutate(dayLengthHrs=tDaylight) %>%
    dplyr::left_join(daily_max_flightBout, by=c("date")) 
  
  # Calculate what is going per daylight vs darkness just for verification purposes
  dataCalcDay_period2<-reAllocate %>%
    dplyr::ungroup() %>%
    dplyr::mutate(date=substr(date_time, 1, 10)) %>%
    dplyr::mutate(Period=ifelse(Period %in% c("Daylight", "Twilight"), "Daylight", "Darkness")) %>%
    dplyr::group_by(species, colony, session_id, date, Period) %>%
    dplyr::mutate(Forage=ifelse(Activity=="Forage", 1, 0), RestWater=ifelse(Activity=="RestWater", 1, 0), Active=ifelse(Activity=="Active", 1, 0), Flight=ifelse(Activity=="Flight", 1, 0), Land=ifelse(Activity=="Land", 1, 0), Daylight=ifelse(Period=="Daylight", 1, 0),
                  Darkness=ifelse(Period=="Darkness", 1, 0), Twilight=ifelse(Period=="Twilight", 1, 0), sstRandom=sst_random_start) %>%
    dplyr::mutate(sstRandom=ifelse(sstRandom < -1.9, 1.9, sstRandom)) %>%
    dplyr::summarise(tForage=sum(Forage)*10/60, tRestWater1=sum(RestWater)*10/60, tFlight=sum(Flight)*10/60, tActive=sum(Active)*10/60, tLand=sum(Land)*10/60,
                     tDaylight=sum(Daylight)*10/60, tDarkness=sum(Darkness)*10/60, tTwilight=sum(Twilight)*10/60, Duration=n_distinct(date_time)*10, MaxDistColKm=max(MaxDistColKm),
                     sst_random=mean(sstRandom), ice_random=mean(ice_random, na.rm=TRUE), distColonyKm_mean=mean(distColonyKm, na.rm=TRUE), mean.lon=mean(lon), mean.lat=mean(lat),
                     c=sample(seq(1.5, 3, 0.1), 1)) %>%
    rename(tRestWater=tRestWater1)
  
  # Save results from Both
  actResults<-dataCalcDay
  boutResults<-dataCalcDay_period2
  
  allResults<-list(actResults, boutResults)
  return(allResults)
  
}

##### Atlantic puffin #####

# At the moment I am just applying Fayet's approach : doi:10.1093/beheco/arw013 
# Max speed is 45 km.hr-1 according to the NINA geolocation processing algorithm

calculateTimeInActivity_AP<-function(data, irmaData){
    
    speciesLatin<-"Fratercula_arctica"    
    actResults<-methodCaitlin(data, irmaData, speciesLatin) 
  
  return(actResults)
  
  
}

##### Little auk #####

# At the moment I will use the puffin one???  
# Max speed is 45 km.hr-1 according to the NINA geolocation processing algorithm

calculateTimeInActivity_LiA<-function(data,  irmaData){
    
    speciesLatin<-"Alle_alle"    
    actResults<-methodCaitlin(data,  irmaData, speciesLatin) 
  
  return(actResults)
  
  
}

##### Common guillemot #####

# At the moment I use a very simplified version of Lila's equation 
# Max speed is 45 km.hr-1 according to the NINA geolocation processing algorithm

calculateTimeInActivity_CoGu<-function(data,  irmaData){
  
    speciesLatin<-"Uria_aalge"    
    actResults<-methodCaitlin(data,  irmaData, speciesLatin) 
  
  return(actResults)
  
  
}

##### Br guillemot #####

# At the moment I use a very simplified version of Lila's equation 
# Max speed is 45 km.hr-1 according to the NINA geolocation processing algorithm
# method refers to Lila's approach vs. (Caitlin)

calculateTimeInActivity_BrGu<-function(data, irmaData){
  
    speciesLatin<-"Uria_lomvia"    
    actResults<-methodCaitlin(data, irmaData, speciesLatin) 
  
  return(actResults)
  
}

##### Assign day.night #####

assignTimePeriod<-function(data) {
  
  data$date<-as.Date(data$date)
  timeDay<-getSunlightTimes(data=data, keep=c("nightEnd", "sunrise", "sunset", "night"))
  #timeDay<-timeDay %>%
    #dplyr::select(nightEnd, sunrise, sunset, night)
  #data2<-data %>%
    #dplyr::bind_cols(timeDay)
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
# Data is data frame containing activity budgets
# Function returns energy spent per day per species & colony...
# Type is daily or monthly resolution
# Type 2 is using sst from pop maps or individual locations ("pop", "ind")

calculateEnergetics<-function(species, data, colonySub, sstVals, type, type2) {
 
# Print randomized values
print(paste0("RMR = ", data$RMR[1]))
print(paste0("c1 = ", data$c1[1]))
print(paste0("c2 = ", data$c2[1]))
print(paste0("c3 = ", data$c3[1]))
print(paste0("c4 = ", data$c4[1]))
#print(paste0("c5 = ", data$c5[1]))
print(paste0("TC = ", data$TC[1]))
print(paste0("Beta_active = ", data$Beta_active[1]))
print(paste0("Beta_rest = ", data$Beta_rest[1]))

# Determine sessions to loop through
sessionNo<-unique(data$session_id)

# List to save results in
energyAll<-list() # Make a list to save results in

# Loop through the sessions
for (session in 1:length(sessionNo)) {

print(paste0("Calculating energy for session ", session))

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
  
  energySpent$weight<-weightG 
  energySpent$session_id<-sessionNo[session]
  energyAll<-rbind(energyAll, energySpent)
  
  }
  
  return(energyAll)    
  
}

##### Kittiwake #####

# sst is in degrees is fixed right now to 10 degrees but will eventually be pulled from a map
# RMR is resting metabolic rate & equal to 1.64 mL O2 g-1 h-1 Gabrielsen et al. 1988
# beta is intercept of RMR at sst = 0 degrees, 
# TC is thermal conductance in wate: 0.1000 mL O2 g-1 h-1 °C-1 (Gabrielsen et al. 1988)
# Weight is weight in g - fixed to 365 g at the moment (Gabrielsen et al., 1988)
# cf is caloric conversion factor of 20.1 J per mL O2 (Schmidt-Nielsen 1997)

calculateEnergetics_BLK_daily<-function(data, weightG) {
  
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
  TCCoef<-data$TC[1]
  
  # Account for change in constants
  flightConstantx<-((flightCoef*450)/450^0.717)*weightG^0.717
  restConstant2x<-((restCoef*450)/450^0.717)*weightG^0.717
  forageConstantx<-((forageCoef*450)/450^0.717)*weightG^0.717
  landConstantx<-((landCoef*450)/450^0.717)*weightG^0.717
  betax<-(((betaCoef*cf/1000)*365)/365^0.717)*weightG^0.717
  TCx<-(((TCCoef*cf/1000)*365)/365^0.717)*weightG^0.717
  
  # Adjust beta so that beta-SST*TC is equal to rest constant 2 at LCT
  LCT<-12.5 # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1474-919X.2006.00618.x
  restConstant1x<-(LCT*TCx + restConstant2x)
  
  # Calculate energetics
  energySub2<-data %>%
    dplyr::group_by(date) %>%
    #dplyr::mutate(DEEkJ=ifelse(sst_random <= LCT_newRest, flightConstant*tFlight + forageConstant*tForage + (restConstant1 - TC*sst_random)*tRestWater +  landConstant*tLand, flightConstant*tFlight + forageConstant*tForage + restConstant3*tRestWater + landConstant*tLand)) %>%
    dplyr::mutate(DEEkJ_active=0, DEEkJ_active_col=0) %>%
	dplyr::mutate(DEEkJ_rest=ifelse(sst_random <=LCT, (restConstant1x - TCx*sst_random)*tRestWater, restConstant2x*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <=LCT, (restConstant1x - TCx*sst_random_colony)*tRestWater, restConstant2x*tRestWater)) %>%
	dplyr::mutate(DEEkJ_flight=flightConstantx*tFlight) %>%
    dplyr::mutate(DEEkJ_forage=forageConstantx*tForage) %>%
	dplyr::mutate(DEEkJ_restland=landConstantx*tLand) %>%
	dplyr::mutate(DEEkJ=DEEkJ_rest + DEEkJ_flight + DEEkJ_forage + DEEkJ_restland) %>%
	dplyr::mutate(DEEkJ_col=DEEkJ_rest_col + DEEkJ_flight + DEEkJ_forage + DEEkJ_restland) %>%
    dplyr::mutate(weight=weightG)
  
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
  TCCoef<-data$TC[1]
  
  # Account for change in constants
  flightConstantx<-((flightCoef*450)/450^0.717)*weightG^0.717
  restConstant2x<-((restCoef*450)/450^0.717)*weightG^0.717
  forageConstantx<-((forageCoef*450)/450^0.717)*weightG^0.717
  landConstantx<-((landCoef*450)/450^0.717)*weightG^0.717
  betax<-(((betaCoef*cf/1000)*365)/365^0.717)*weightG^0.717
  TCx<-(((TCCoef*cf/1000)*365)/365^0.717)*weightG^0.717
  
  # Adjust beta so that beta-SST*TC is equal to rest constant 2 at LCT
  LCT<-12.5 # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1474-919X.2006.00618.x
  restConstant1x<-(LCT*TCx + restConstant2x)
  
  # Turn SST raster into a data frame
  sstDf<-as.data.frame(sstVals, xy=TRUE)
  
  # Calculate energy for every cell according to sst in that cell
  sstDf$DEEkJ<-ifelse(sstDf$layer <= LCT, flightConstantx*data$tFlight_month + forageConstantx*data$tForage_month +
  (restConstant1x - TCx*sstDf$layer)*data$tRestWater_month +  landConstantx*data$tLand_month, flightConstantx*data$tFlight_month + 
  forageConstantx*data$tForage_month + restConstant2x*data$tRestWater_month + landConstantx*data$tLand_month)
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, layer, DEEkJ) %>%
  rename(sst=layer) 
  
  return(sstDf_final)


}

##### Northern fulmar #####

# RMR is resting metabolic rate & equal to 1.00 mL O2 g-1 h-1 Gabrielsen et al. 1988
# beta is intercept of RMR at sst = 0 degrees, 
# TC is thermal conductance in wate: 0.04 mL O2 g-1 h-1 °C-1 (Gabrielsen et al. 1988)
# Weight is weight in g - fixed to 651 g at the moment (Gabrielsen et al., 1988)
# cf is caloric conversion factor of 20.1 J per mL O2 (Schmidt-Nielsen 1997)

calculateEnergetics_NF_daily<-function(data, weightG) {
  
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
  TCCoef<-data$TC[1]
  
  # Account for change in constants
  #flightConstantx<-(7.3*RMR*cf/1000)*weightG
  flightConstantx<-(((flightCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  restConstant2x<-(((restCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  forageConstantx<-((forageCoef*450)/450^0.717)*weightG^0.765
  landConstantx<-(((landCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  betax<-(((betaCoef*cf/1000)*651)/651^0.765)*weightG^0.765
  TCx<-(((TCCoef*cf/1000)*651)/651^0.765)*weightG^0.765
  
  # Adjust beta so that beta-SST*TC is equal to rest constant 2 at LCT
  LCT<-9 # Gabrielsen et al. 1988
  restConstant1x<-(LCT*TCx + restConstant2x)
  
  # Calculate energetics
  #energySub2<-data %>%
    #dplyr::group_by(date) %>%
    #dplyr::mutate(DEEkJ=ifelse(sst_random <= LCT_newRest, flightConstant*tFlight + forageConstant*tForage + (restConstant1 - TC*sst_random)*tRestWater +  landConstant*tLand, flightConstant*tFlight + forageConstant*tForage + restConstant3*tRestWater + landConstant*tLand)) %>%
   # dplyr::mutate(DEEkJ=ifelse(sst_random <= LCT, flightConstantx*tFlight + forageConstantx*tForage + (restConstant1x - TCx*sst_random)*tRestWater +  landConstantx*tLand, flightConstantx*tFlight + forageConstantx*tForage + restConstant2x*tRestWater + landConstantx*tLand)) %>%
   # dplyr::mutate(DEEkJ_thermo=ifelse(sst_random <=LCT, (restConstant1x - TCx*sst_random)*tRestWater, 0)) %>%
   # dplyr::mutate(DEEkJ_flight=flightConstantx*tFlight) %>%
   # dplyr::mutate(DEEkJ_forage=forageConstantx*tForage) %>%
   # dplyr::mutate(DEEkJ_restwater=ifelse(sst_random <=LCT, (restConstant1x - TCx*sst_random)*tRestWater, restConstant2x*tRestWater)) %>%
   # dplyr::mutate(DEEkJ_restland=landConstantx*tLand) %>%
    #dplyr::mutate(weight=weightG)
	
	energySub2<-data %>%
    dplyr::group_by(date) %>%
    #dplyr::mutate(DEEkJ=ifelse(sst_random <= LCT_newRest, flightConstant*tFlight + forageConstant*tForage + (restConstant1 - TC*sst_random)*tRestWater +  landConstant*tLand, flightConstant*tFlight + forageConstant*tForage + restConstant3*tRestWater + landConstant*tLand)) %>%
    dplyr::mutate(DEEkJ_active=0, DEEkJ_active_col=0) %>%
	dplyr::mutate(DEEkJ_rest=ifelse(sst_random <=LCT, (restConstant1x - TCx*sst_random)*tRestWater, restConstant2x*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <=LCT, (restConstant1x - TCx*sst_random_colony)*tRestWater, restConstant2x*tRestWater)) %>%
	dplyr::mutate(DEEkJ_flight=flightConstantx*tFlight) %>%
    dplyr::mutate(DEEkJ_forage=forageConstantx*tForage) %>%
	dplyr::mutate(DEEkJ_restland=landConstantx*tLand) %>%
	dplyr::mutate(DEEkJ=DEEkJ_rest + DEEkJ_flight + DEEkJ_forage + DEEkJ_restland) %>%
	dplyr::mutate(DEEkJ_col=DEEkJ_rest_col + DEEkJ_flight + DEEkJ_forage + DEEkJ_restland) %>%
    dplyr::mutate(weight=weightG)
  
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
  TCCoef<-data$TC[1]
  
  # Account for change in constants
  #flightConstantx<-(7.3*RMR*cf/1000)*weightG
  flightConstantx<-(((flightCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  restConstant2x<-(((restCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  forageConstantx<-((forageCoef*450)/450^0.717)*weightG^0.765
  landConstantx<-(((landCoef*RMR*cf/1000)*651)/651^0.765)*weightG^0.765
  betax<-(((betaCoef*cf/1000)*651)/651^0.765)*weightG^0.765
  TCx<-(((TCCoef*cf/1000)*651)/651^0.765)*weightG^0.765
  
  # Adjust beta so that beta-SST*TC is equal to rest constant 2 at LCT
  LCT<-9 # Gabrielsen et al. 1988
  restConstant1x<-(LCT*TCx + restConstant2x)
  
  # Turn SST raster into a data frame
  sstDf<-as.data.frame(sstVals, xy=TRUE)
  
  # Calculate energy for every cell according to sst in that cell
  sstDf$DEEkJ<-ifelse(sstDf$layer <= LCT, flightConstantx*data$tFlight_month + forageConstantx*data$tForage_month + 
  (restConstant1x - TCx*sstDf$layer)*data$tRestWater_month +  landConstantx*data$tLand_month, flightConstantx*data$tFlight_month + 
  forageConstantx*data$tForage_month + restConstant2x*data$tRestWater_month + landConstantx*data$tLand_month)
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, layer, DEEkJ) %>%
  rename(sst=layer) 
  
  return(sstDf_final)
  
  
  
  
}

##### Common guillemot #####

# CostDivider is mean weight in g of all available weight data for this species (from database). It will be used to
# calculate energy cost per g
# beta is intercept of RMR at sst = 0 degrees & here we have two for different activities -> see Kyle & Gaston 2014 J per hr
# TC is thermal conductance in water: 2.75 -> Kyle & Gaston 2014 J/hr
# Weight is weight in g - it will be a mean from colony-specific values
# LCT is lower critical tolerance -> 14.3 (following Buckingham et al. - in review)
# otherwise whole energetic equation is from Kyle & Gaston 2014/Buckingham et al. in review

calculateEnergetics_CoGu_daily<-function(data, CostDivider, weightG) {
  
  # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC<-data$TC[1] # kJ.hr
  
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
  TC_Constant<-(TC*convf)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT<-14.18
  restConstant_adjust<-restConstant1 - LCT*TC_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT*TC_Constant
  
  # And now we calculate new LCT?
  #LCT_newActive<-(activeConstant-restConstant2)/TC_Constant
  #LCT_newActive<-(activeConstant-activeConstant2)/TC_Constant
  #LCT_newRest<-(restConstant1-restConstant2)/TC_Constant
  
  # Determine whether any temps are above LCT active
  #temps<-subset(data, sst_random > LCT_newActive)
  
 # if (nrow(temps)>0) {
 #   stop(print("Active LCT has an issue"))
 # }
  
  energySub2<-data %>%
    dplyr::group_by(date) %>%
	#dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT_newActive, (activeConstant-sst_random*TC_Constant)*tActive, activeConstant2*tActive)) %>%
	dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT, (activeConstant-sst_random*TC_Constant)*tActive, activeConstant_adjust*tActive)) %>%
	#dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT_newActive, (activeConstant-sst_random_colony*TC_Constant)*tActive, activeConstant2*tActive)) %>%
	dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT, (activeConstant-sst_random_colony*TC_Constant)*tActive, activeConstant_adjust*tActive)) %>%
	#dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT_newRest, (restConstant1 - TC_Constant*sst_random)*tRestWater , restConstant2*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT, (restConstant1 - TC_Constant*sst_random)*tRestWater , restConstant_adjust*tRestWater)) %>%
	#dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT_newRest, (restConstant1 - TC_Constant*sst_random_colony)*tRestWater , restConstant2*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT, (restConstant1 - TC_Constant*sst_random_colony)*tRestWater , restConstant_adjust*tRestWater)) %>%
    dplyr::mutate(DEEkJ_flight=flightConstant*tFlight) %>%
    dplyr::mutate(DEEkJ_forage=0) %>%
    dplyr::mutate(DEEkJ_restland=landConstant*tLand) %>%
	dplyr::mutate(DEEkJ=DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland) %>%
	dplyr::mutate(DEEkJ_col=DEEkJ_active_col + DEEkJ_rest_col + DEEkJ_flight + DEEkJ_restland) %>%
    dplyr::mutate(weight=weightG)
	#dplyr::mutate(LCT_active=LCT_newActive) %>%
	#dplyr::mutate(LCT_rest=LCT_newRest)
  
  return(energySub2)
  
}

calculateEnergetics_CoGu_daily_map<-function(data, CostDivider, weightG, sstVals) {
  
  # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC<-data$TC[1] # kJ.hr
  
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
  TC_Constant<-(TC*convf)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT<-14.18
  restConstant_adjust<-landConstant + LCT*TC_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT*TC_Constant
  
  # Turn SST raster into a data frame
  sstDf<-as.data.frame(sstVals, xy=TRUE)
  
  # Compute cost of being active & rest seperately
  sstDf$DEEkJ_active<-ifelse(sstDf$layer <= LCT, (activeConstant-sstDf$layer*TC_Constant)*data$tActive_month, activeConstant_adjust*data$tActive_month)
  sstDf$DEEkJ_rest<-ifelse(sstDf$layer <= LCT, (restConstant1-sstDf$layer*TC_Constant)*data$tRestWater_month, restConstant2*data$tRestWater_month)
  sstDf$DEEkJ<-sstDf$DEEkJ_active + sstDf$DEEkJ_rest + flightConstant*data$tFlight_month + data$tLand_month*landConstant
  sstDf<-sstDf %>%
  dplyr::select(-c(DEEkJ_active, DEEkJ_rest))
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, layer, DEEkJ) %>%
  rename(sst=layer) 
  
  return(sstDf_final)
  
}

##### Brunnich's guillemot #####

# CostDivider is mean weight in g of all available weight data for this species (from database). It will be used to
# calculate energy cost per g
# beta is intercept of RMR at sst = 0 degrees & here we have two for different activities -> see Kyle & Gaston 2014 J per hr
# TC is thermal conductance in water: 2.75 -> Kyle & Gaston 2014 J/hr
# Weight is weight in g - it will be a mean from colony-specific values
# LCT is lower critical tolerance -> 14.3 (following Buckingham et al. - in review)
# otherwise whole energetic equation is from Kyle & Gaston 2014/Buckingham et al. in review

calculateEnergetics_BrGu_daily<-function(data, CostDivider, weightG) {
  
  # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC<-data$TC[1] # kJ.hr
  
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
  TC_Constant<-(TC*convf)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT<-14.18
  restConstant_adjust<-restConstant1 - LCT*TC_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT*TC_Constant
  
  # And now we calculate new LCT?
  #LCT_newActive<-(activeConstant-restConstant2)/TC_Constant
  #LCT_newActive<-(activeConstant-activeConstant2)/TC_Constant
  #LCT_newRest<-(restConstant1-restConstant2)/TC_Constant
  
  # Determine whether any temps are above LCT active
  #temps<-subset(data, sst_random > LCT_newActive)
  
 # if (nrow(temps)>0) {
 #   stop(print("Active LCT has an issue"))
 # }
  
  energySub2<-data %>%
    dplyr::group_by(date) %>%
	#dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT_newActive, (activeConstant-sst_random*TC_Constant)*tActive, activeConstant2*tActive)) %>%
	dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT, (activeConstant-sst_random*TC_Constant)*tActive, activeConstant_adjust*tActive)) %>%
	#dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT_newActive, (activeConstant-sst_random_colony*TC_Constant)*tActive, activeConstant2*tActive)) %>%
	dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT, (activeConstant-sst_random_colony*TC_Constant)*tActive, activeConstant_adjust*tActive)) %>%
	#dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT_newRest, (restConstant1 - TC_Constant*sst_random)*tRestWater , restConstant2*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT, (restConstant1 - TC_Constant*sst_random)*tRestWater , restConstant_adjust*tRestWater)) %>%
	#dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT_newRest, (restConstant1 - TC_Constant*sst_random_colony)*tRestWater , restConstant2*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT, (restConstant1 - TC_Constant*sst_random_colony)*tRestWater , restConstant_adjust*tRestWater)) %>%
    dplyr::mutate(DEEkJ_flight=flightConstant*tFlight) %>%
    dplyr::mutate(DEEkJ_forage=0) %>%
    dplyr::mutate(DEEkJ_restland=landConstant*tLand) %>%
	dplyr::mutate(DEEkJ=DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland) %>%
	dplyr::mutate(DEEkJ_col=DEEkJ_active_col + DEEkJ_rest_col + DEEkJ_flight + DEEkJ_restland) %>%
    dplyr::mutate(weight=weightG)
	#dplyr::mutate(LCT_active=LCT_newActive) %>%
	#dplyr::mutate(LCT_rest=LCT_newRest)
  
  return(energySub2)
  
  
}

calculateEnergetics_BrGu_daily_map<-function(data, CostDivider, weightG, sstVals) {
  
  # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC<-data$TC[1] # kJ.hr
  
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
  TC_Constant<-(TC*convf)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT<-14.18
  restConstant_adjust<-landConstant + LCT*TC_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT*TC_Constant
  
  # Turn SST raster into a data frame
  sstDf<-as.data.frame(sstVals, xy=TRUE)
  
  # Compute cost of being active & rest seperately
  sstDf$DEEkJ_active<-ifelse(sstDf$layer <= LCT, (activeConstant-sstDf$layer*TC_Constant)*data$tActive_month, activeConstant_adjust*data$tActive_month)
  sstDf$DEEkJ_rest<-ifelse(sstDf$layer <= LCT, (restConstant1-sstDf$layer*TC_Constant)*data$tRestWater_month, restConstant2*data$tRestWater_month)
  sstDf$DEEkJ<-sstDf$DEEkJ_active + sstDf$DEEkJ_rest + flightConstant*data$tFlight_month + data$tLand_month*landConstant
  sstDf<-sstDf %>%
  dplyr::select(-c(DEEkJ_active, DEEkJ_rest))
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, layer, DEEkJ) %>%
  rename(sst=layer) 
  
  return(sstDf_final)
  
  
}

##### Little auk #####

# CostDivider is mean weight in g of all available weight data for this species (from database). It will be used to
# calculate energy cost per g
# beta is intercept of RMR at sst = 0 degrees & here we have two for different activities -> see Kyle & Gaston 2014 J per hr
# TC is thermal conductance in water: 2.75 -> Kyle & Gaston 2014 J/hr
# Weight is weight in g - it will be a mean from colony-specific values
# LCT is lower critical tolerance -> 14.3 (following Buckingham et al. - in review)
# otherwise whole energetic equation is from Kyle & Gaston 2014/Buckingham et al. in review

calculateEnergetics_LiA_daily<-function(data, CostDivider,  weightG) {
  
  # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj.hr.g<-flightCoef/24 # kJ.hr.g
  flightCoef_kj.hr<-flightCoef_kj.hr.g*150.95 # kJ.hr-1
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC<-data$TC[1] # kJ.hr
  
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
  TC_Constant<-(TC*convf)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT<-14.18
  restConstant_adjust<-restConstant1 - LCT*TC_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT*TC_Constant
  
  # And now we calculate new LCT?
  #LCT_newActive<-(activeConstant-restConstant2)/TC_Constant
  #LCT_newActive<-(activeConstant-activeConstant2)/TC_Constant
  #LCT_newRest<-(restConstant1-restConstant2)/TC_Constant
  
  # Determine whether any temps are above LCT active
  #temps<-subset(data, sst_random > LCT_newActive)
  
 # if (nrow(temps)>0) {
 #   stop(print("Active LCT has an issue"))
 # }
  
  energySub2<-data %>%
    dplyr::group_by(date) %>%
	#dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT_newActive, (activeConstant-sst_random*TC_Constant)*tActive, activeConstant2*tActive)) %>%
	dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT, (activeConstant-sst_random*TC_Constant)*tActive, activeConstant_adjust*tActive)) %>%
	#dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT_newActive, (activeConstant-sst_random_colony*TC_Constant)*tActive, activeConstant2*tActive)) %>%
	dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT, (activeConstant-sst_random_colony*TC_Constant)*tActive, activeConstant_adjust*tActive)) %>%
	#dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT_newRest, (restConstant1 - TC_Constant*sst_random)*tRestWater , restConstant2*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT, (restConstant1 - TC_Constant*sst_random)*tRestWater , restConstant_adjust*tRestWater)) %>%
	#dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT_newRest, (restConstant1 - TC_Constant*sst_random_colony)*tRestWater , restConstant2*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT, (restConstant1 - TC_Constant*sst_random_colony)*tRestWater , restConstant_adjust*tRestWater)) %>%
    dplyr::mutate(DEEkJ_flight=flightConstant*tFlight) %>%
    dplyr::mutate(DEEkJ_forage=0) %>%
    dplyr::mutate(DEEkJ_restland=landConstant*tLand) %>%
	dplyr::mutate(DEEkJ=DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland) %>%
	dplyr::mutate(DEEkJ_col=DEEkJ_active_col + DEEkJ_rest_col + DEEkJ_flight + DEEkJ_restland) %>%
    dplyr::mutate(weight=weightG)
	#dplyr::mutate(LCT_active=LCT_newActive) %>%
	#dplyr::mutate(LCT_rest=LCT_newRest)
  
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
  TC<-data$TC[1]
  
  # Rest & land coefs
  restCoef2<-data$c4[1]
  restCoef2_kj<-restCoef2*3.6 #kJ.hr
  
  landCoef<-data$c3[1]
  landCoef_kj<-landCoef*3.6 # kj.hr -> should be the same as the rest coefficient
  
  # Here is a conversion factor to transform these to g.kj which incorporates allometric scaling
  convf<-1/CostDivider^0.689 # no wing loading
  
   # Active coef when thermoneutral (from kyle's paper)
  #activeCoef_neut<-data$c5[1]
  #activeCoef_neut_kj<-activeCoef_neut*3.6 # kJ.hr
  
  #if (activeCoef_neut_kj > activeCoef) stop (print("Error with active coefs")) # Because the coef when thermoneutral must be lower...
  
  # Account for change in constants
  flightConstant<-(flightCoef_kj.hr/150.95^0.689)*weightG^0.689
  activeConstant<-(activeCoef*convf)*weightG^0.689
  #activeConstant2<-(activeCoef_neut_kj*convf)*weightG^0.689
  restConstant1<-(restCoef1*convf)*weightG^0.689
  restConstant2<-(restCoef2_kj*convf)*weightG^0.689
  TC_Constant<-(TC*convf)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT<-14.18
  restConstant_adjust<-landConstant + LCT*TC_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT*TC_Constant
  
  # Turn SST raster into a data frame
  sstDf<-as.data.frame(sstVals, xy=TRUE)
  
  # Compute cost of being active & rest seperately
  sstDf$DEEkJ_active<-ifelse(sstDf$layer <= LCT, (activeConstant-sstDf$layer*TC_Constant)*data$tActive_month, activeConstant_adjust*data$tActive_month)
  sstDf$DEEkJ_rest<-ifelse(sstDf$layer <= LCT, (restConstant1-sstDf$layer*TC_Constant)*data$tRestWater_month, restConstant2*data$tRestWater_month)
  sstDf$DEEkJ<-sstDf$DEEkJ_active + sstDf$DEEkJ_rest + flightConstant*data$tFlight_month + data$tLand_month*landConstant
  sstDf<-sstDf %>%
  dplyr::select(-c(DEEkJ_active, DEEkJ_rest))
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, layer, DEEkJ) %>%
  rename(sst=layer) 
  
  return(sstDf_final)
  
}

##### Atlantic puffin #####

# CostDivider is mean weight in g of all available weight data for this species (from database). It will be used to
# calculate energy cost per g
# beta is intercept of RMR at sst = 0 degrees & here we have two for different activities -> see Kyle & Gaston 2014 J per hr
# TC is thermal conductance in water: 2.75 -> Kyle & Gaston 2014 J/hr
# Weight is weight in g - it will be a mean from colony-specific values
# LCT is lower critical tolerance -> 14.3 (following Buckingham et al. - in review)
# otherwise whole energetic equation is from Kyle & Gaston 2014/Buckingham et al. in review


calculateEnergetics_AP_daily<-function(data, CostDivider, weightG) {
  
  # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC<-data$TC[1] # kJ.hr
  
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
  TC_Constant<-(TC*convf)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT<-14.18
  restConstant_adjust<-restConstant1 - LCT*TC_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT*TC_Constant
  
  # And now we calculate new LCT?
  #LCT_newActive<-(activeConstant-restConstant2)/TC_Constant
  #LCT_newActive<-(activeConstant-activeConstant2)/TC_Constant
  #LCT_newRest<-(restConstant1-restConstant2)/TC_Constant
  
  # Determine whether any temps are above LCT active
  #temps<-subset(data, sst_random > LCT_newActive)
  
 # if (nrow(temps)>0) {
 #   stop(print("Active LCT has an issue"))
 # }
  
  energySub2<-data %>%
    dplyr::group_by(date) %>%
	#dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT_newActive, (activeConstant-sst_random*TC_Constant)*tActive, activeConstant2*tActive)) %>%
	dplyr::mutate(DEEkJ_active=ifelse(sst_random <=LCT, (activeConstant-sst_random*TC_Constant)*tActive, activeConstant_adjust*tActive)) %>%
	#dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT_newActive, (activeConstant-sst_random_colony*TC_Constant)*tActive, activeConstant2*tActive)) %>%
	dplyr::mutate(DEEkJ_active_col=ifelse(sst_random_colony <=LCT, (activeConstant-sst_random_colony*TC_Constant)*tActive, activeConstant_adjust*tActive)) %>%
	#dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT_newRest, (restConstant1 - TC_Constant*sst_random)*tRestWater , restConstant2*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest=ifelse(sst_random <= LCT, (restConstant1 - TC_Constant*sst_random)*tRestWater , restConstant_adjust*tRestWater)) %>%
	#dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT_newRest, (restConstant1 - TC_Constant*sst_random_colony)*tRestWater , restConstant2*tRestWater)) %>%
	dplyr::mutate(DEEkJ_rest_col=ifelse(sst_random_colony <= LCT, (restConstant1 - TC_Constant*sst_random_colony)*tRestWater , restConstant_adjust*tRestWater)) %>%
    dplyr::mutate(DEEkJ_flight=flightConstant*tFlight) %>%
    dplyr::mutate(DEEkJ_forage=0) %>%
    dplyr::mutate(DEEkJ_restland=landConstant*tLand) %>%
	dplyr::mutate(DEEkJ=DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland) %>%
	dplyr::mutate(DEEkJ_col=DEEkJ_active_col + DEEkJ_rest_col + DEEkJ_flight + DEEkJ_restland) %>%
    dplyr::mutate(weight=weightG)
	#dplyr::mutate(LCT_active=LCT_newActive) %>%
	#dplyr::mutate(LCT_rest=LCT_newRest)
  
  return(energySub2)
  
  
}

calculateEnergetics_AP_daily_map<-function(data,  CostDivider, weightG, sstVals) {
  
 # First I set up the activity cost multipliers which are from Elliott & Gaston et al. & Buckingham (in review)
  flightCoef<-data$c1[1]
  flightCoef_kj<-flightCoef*3.6 # kJ.hr
  
  # Add a fake error based on error distribution of other terms 
  activeCoef<-data$Beta_active[1] # kj.hr
  restCoef1<-data$Beta_rest[1] # kj.hr
  TC<-data$TC[1] # kJ.hr
  
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
  TC_Constant<-(TC*convf)*weightG^0.689
  landConstant<-(landCoef_kj*convf)*weightG^0.689
  
  # LCT is 14.18 (Buckingham et al. 2025)
  LCT<-14.18
  restConstant_adjust<-landConstant + LCT*TC_Constant # So that rest is equal to land when sst > LCT
  activeConstant_adjust<-activeConstant - LCT*TC_Constant
  
  # Turn SST raster into a data frame
  sstDf<-as.data.frame(sstVals, xy=TRUE)
  
  # Compute cost of being active & rest seperately
  sstDf$DEEkJ_active<-ifelse(sstDf$layer <= LCT, (activeConstant-sstDf$layer*TC_Constant)*data$tActive_month, activeConstant_adjust*data$tActive_month)
  sstDf$DEEkJ_rest<-ifelse(sstDf$layer <= LCT, (restConstant1-sstDf$layer*TC_Constant)*data$tRestWater_month, restConstant2*data$tRestWater_month)
  sstDf$DEEkJ<-sstDf$DEEkJ_active + sstDf$DEEkJ_rest + flightConstant*data$tFlight_month + data$tLand_month*landConstant
  sstDf<-sstDf %>%
  dplyr::select(-c(DEEkJ_active, DEEkJ_rest))
  
  # Add weight for converting later
  sstDf$weight<-weightG
  
  # Add other important information
  sstDf$individ_id<-data$individ_id[1]
  sstDf$species<-data$species[1]
  sstDf$colony<-data$colony[1]
  sstDf$rep<-data$rep[1]
  
  # Change order of columns
  sstDf_final<-sstDf %>%
  dplyr::select(rep, species, colony, individ_id, weight, x, y, layer, DEEkJ) %>%
  rename(sst=layer) 
  
  return(sstDf_final)
  
  
}

#### SST FUNCTIONS ####

#### Extract appropriate SST values (current) ####
# type is current vs. future

extractSSTmaps<-function(speciesSub, colonySub, type, metadataMaps, mapLox, sstLox) {
  
  # Selecting one of the six species
  species_name<-c("Fratercula_arctica","Uria_lomvia","Uria_aalge","Alle_alle",
                  "Fulmarus_glacialis","Rissa_tridactyla")
  
  # Selecting one of the six species
  common_name<-c("Atlantic puffin","Brünnich's guillemot","Common guillemot","Little auk",
                 "Northern fulmar","Black-legged kittiwake")  
  
  # set up a key to subset pop maps with
  key<-data.frame(latin=species_name, species=common_name)
  keysub<-subset(key, species==speciesSub)
  
  # List population maps
  pop.maps<-list.files(mapLox, full.names=TRUE)
  #pop.maps<-list.files("./data/intermediate/population_maps/", full.names=TRUE)
  
  # Subset to species of interest
  pop.maps.species<-pop.maps[grepl(keysub$latin, pop.maps)]
  
  # Figure out the colonies I am interest in based on model colony
  colonies<-subset(metadataMaps, modelcolonies==colonySub & species==speciesSub)
  
  # Make sure only one of each
  colonies<-colonies %>%
    dplyr::group_by(colonies) %>%
    dplyr::slice(1) 
  
  #random<-sample(1:nrow(colonies), 1) # We will just focus on one right now... to speed things up
  
  colonies<-colonies %>%
    ungroup() %>%
    dplyr::slice(1)
  
  # Create list for saving information in...
  sstAllModels<-list()
  
  # Open up sst data
  sst.data.location<-list.files(paste0(sstLox, type), full.names = TRUE)
  #sst.data.location<-list.files(paste0("./data/intermediate/sst/", type), full.names = TRUE)
  
  #for (i in seq_along(sst.data.location)) {
  
  for (i in 1) {
    
    #start.time <- Sys.time() 
    
    # open folder x
    sst.data.sub<-list.files(sst.data.location[i], full.names = TRUE)
    sstMean<-rast(sst.data.sub[grepl("mean", sst.data.sub)])
    sstSd<-rast(sst.data.sub[grepl("sd", sst.data.sub)])
    
    # Generate random sst value & then turn into data frame
    sstRandom<-sstMean
    #sstRandom<-calc(sstMean, function(x)rnorm(ncell(x), x, values(sstSd)))
    values(sstRandom)<-rnorm(ncell(sstRandom)*12, values(sstMean), values(sstSd))
    names(sstRandom)<-1:12
    sstRandomVals<-as.data.frame(sstRandom, xy=TRUE)
    sstRandomVals<-sstRandomVals %>% 
      gather(month, sstRandom, -c(x, y), na.rm = T) 
    colnames(sstRandomVals)<-c("x", "y", "month", "sstRandom")
    #sstRandomVals$month<-substr(sstRandomVals$month, 2, nchar(sstRandomVals$month))
    
    # Create list to save results in
    sstAllCols<-list()
    
    for (j in 1:nrow(colonies)) {
      
      start.time <- Sys.time()
      
      #for (i in 1:2) {
      
      print(paste0("Extracting sst for colony", j, "/", nrow(colonies), "..."))
      
      # subset sub-colony map  
      mapcolonySub<-colonies[j,]
      pop.maps.colony<-pop.maps.species[grepl(paste0("_", mapcolonySub$colonies[1], "_"), pop.maps.species)]
      
      # Make function fail is incorrect number of maps (should be 12 for all the months)
      if(!length(pop.maps.colony) %in% c(12)) {
        
        # Try & grep with shorter version of name
        pop.maps.colony<-pop.maps.species[grepl(substr(mapcolonySub$colonies[1], 1, 10), pop.maps.species)]
        
        if(!length(pop.maps.colony) %in% c(12)) {
          
          print("Number of pop maps weird")
          break
          
        }
        
      }
      
      #Stack the 12 months of data
      pop.maps.stack<-rast(pop.maps.colony)
      
      # pop data into data frame
      names(pop.maps.stack)<-1:12
      popDf<-as.data.frame(pop.maps.stack, xy=TRUE)
      popDf<-popDf %>% 
        gather(month, NoBirds, -c(x, y), na.rm = T) 
      colnames(popDf)<-c("x", "y", "month", "NoBirds")
      
      # Join values & change column names...
      joined<-popDf %>%
        #dplyr::filter(NoBirds>0) %>%
        arrange(month, x, y) %>% # otherwise file is enormous
        dplyr::left_join(sstRandomVals, by=c("x", "y", "month"))
      joined$species<-speciesSub
      joined$colony<-mapcolonySub$colonies[1]
      #joined$month<-substr(joined$month, 2, nchar(joined$month))
      joined$model<-sst.data.location[i]
      joined$sst_scenario<-type
      joined<-joined %>%
        #replace_na(list(sstRandom=999)) %>%
        dplyr::mutate(sstRandom=if_else(sstRandom< -2 & !is.na(sstRandom), -2, sstRandom))
      
      # Save all models
      sstAllCols<-rbind(sstAllCols, joined)
      
    }
    
    # Save models for all colonies
    
    sstAllCols$colonies<-mapcolonySub$colonies[1]
    sstAllCols$modelcolonies<-mapcolonySub$modelcolonies[1]
    sstAllModels<-rbind(sstAllModels, sstAllCols)
    
  }
  
  return(sstAllModels)
  
}


### Extract sst values for locations ####

extractSSTlox<-function(coords, type, sstLox) {
  
  # Establish standard projection
  projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # Flat projection
  
  # Open up sst data
  sst.data.location<-list.files(paste0(sstLox, type), full.names = TRUE)
  
  coordsModels<-list()
  
  #for (j in 1:length(sst.data.location)) {
  for (j in 1) {  
    
    # open folder x
    sst.data.sub<-list.files(sst.data.location[j], full.names = TRUE)
    sstMean<-raster::stack(sst.data.sub[grepl("mean", sst.data.sub)])
    sstSd<-raster::stack(sst.data.sub[grepl("sd", sst.data.sub)])
    
    # assign month to coordinate data
    coords$Month<-as.numeric(substr(coords$date, 6, 7))
    
    # Determine number of unique months
    months<-unique(coords$Month)
    
    # Make list to save monthly values
    coordsMonth<-list()
    
    for (k in 1:length(months)){
      
      # SST Sub
      sstMeanSub<-subset(sstMean, months[k])
      sstSdSub<-subset(sstSd, months[k])
      
      # Subset dataset
      coordsSub<-subset(coords, Month==months[k])
      
      # Extract sst at locations
      coordinates(coordsSub)<-~mean.lon + mean.lat
      proj4string(coordsSub)<-projection_84
      projMap<-proj4string(sstMeanSub)
      coordsProject<-spTransform(coordsSub, projMap)
      #loxSf<-st_as_sf(coordsProject)
      #bufferPoints <- st_buffer(loxSf, dist=200000)
      coordsProject$sst_mean<-raster::extract(sstMeanSub, coordsProject)
      coordsProject$sst_sd<-raster::extract(sstSdSub, coordsProject)
      
      # Turn back into data frame
      coordsDfSub<-data.frame(coordsProject) %>%
        dplyr::group_by(date) %>%
        dplyr::mutate(sst_random=ifelse(!is.na(sst_mean) & !is.na(sst_sd), rnorm(1, mean=sst_mean, sd=sst_sd), NA))
      
      # save monthly result
      coordsMonth<-rbind(coordsMonth, coordsDfSub)
      
    }  
    
    # Save for all models
    coordsMonth$sstModel<-sst.data.location[j]
    coordsMonth$time<-type
    coordsModels<-rbind(coordsModels, coordsMonth)
    
  }
  
  return(coordsModels)
  
  
}

############## Plotting activity & energy ######

scale_0_1 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Plot population-specific activity #

mapActivity_pop<-function(data, extentdf, behaviour, species) {

# Determine unique months we will need to plot #
months<-c(9, 10, 11, 12, 1, 2, 3, 4)

# prepare empty PDF plot
pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/popMaps/activity/", species, "_", data$colony_og[1], "_", behaviour, ".pdf"), width=10)

print("Mapping")

for (i in 1:length(months)) {

# Determine behaviour that will be plotted
behaviorplot<-paste0("time", behaviour)
behaviorplot2<-paste0("time", behaviour, "_sd")

# Determine possible maximum values
dataSub<-subset(data, month==months[i])
maxBeh1<-max(dataSub[[behaviorplot]])
maxBeh2<-max(dataSub[[behaviorplot2]])
maxBeh<-ifelse(maxBeh1>maxBeh2, maxBeh1, maxBeh2)

dataSub$MeanvsSD<-dataSub[[behaviorplot]]-dataSub[[behaviorplot2]]
minBeh1<-min(dataSub[[behaviorplot]])
minBeh2<-min(dataSub$MeanvsSD)
minBeh<-ifelse(minBeh1<minBeh2, minBeh1, minBeh2)

# plot Mean behaviour
plotMean<-ggplot() +
  geom_tile(data=subset(data, month==months[i]), aes_string(x="x", y="y", fill=paste0(behaviorplot))) +
  geom_text(data=subset(data, month==months[i] & NoBirds==1), aes_string(x="x", y="y", label="NoBirds"), cex=1.5) +
  #scale_alpha(range = c(0.5, 1)) +
  scale_fill_gradientn(paste0(behaviorplot), colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$coords.x1)-400000, max(extentdf$coords.x1) + 400000), ylim=c(min(extentdf$coords.x2)-400000, max(extentdf$coords.x2) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
  
# plot SD behaviour
plotSD<-ggplot() +
  geom_tile(data=subset(data, month==months[i]), aes_string(x="x", y="y", fill=paste0(behaviorplot2))) +
   geom_text(data=subset(data, month==months[i] & NoBirds==1), aes_string(x="x", y="y", label="NoBirds"), cex=1.5) +
  scale_fill_gradientn(paste0(behaviorplot2), colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$coords.x1)-400000, max(extentdf$coords.x1) + 400000), ylim=c(min(extentdf$coords.x2)-400000, max(extentdf$coords.x2) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
  
# Plot areas with high mean & low SD

data$meanVSsd<-data[[behaviorplot]] - data[[behaviorplot2]]
dataSub$mean_scaled<-scale_0_1(dataSub[[behaviorplot]])
dataSub$sd_scaled<-scale_0_1(dataSub[[behaviorplot2]])
dataSub$meanVSsd_scaled<-dataSub$mean_scaled- dataSub$sd_scaled

plotMeanSD<-ggplot() +
  geom_tile(data=subset(dataSub, month==months[i] & NoBirds >1), aes_string(x="x", y="y", fill="MeanvsSD")) +
  #geom_text(data=subset(dataSub, month==months[i] & NoBirds>1), aes_string(x="x", y="y", label="NoBirds"), cex=1.5) +
  scale_fill_gradientn('Mean - SD', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"),limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$coords.x1)-400000, max(extentdf$coords.x1) + 400000), ylim=c(min(extentdf$coords.x2)-400000, max(extentdf$coords.x2) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
  
grid.arrange(plotMean, plotSD, plotMeanSD, nrow=1)
  
  }
  
  print("Saving")
  
  dev.off()


}

# Plot population-specific activity #

mapActivity_species<-function(data_mean, data_sd, extentdf, behaviour, species) {

# Determine unique months we will need to plot #
months<-c(9, 10, 11, 12, 1, 2, 3, 4)

# Make some bins for plotting SD diff vs. colony number #
# Determine behaviour that will be plotted
behaviorplot<-paste0("time", behaviour, "Mean")
behaviorplot2<-paste0("time", behaviour, "_sd_diff")
behaviorplot3<-paste0("time", behaviour, "_sd")
data_sd$sd_diff<-data_sd[[behaviorplot2]] # give column standard name so easy to calculate quantiles
#data_sd_sub<-subset(data_sd, colonies>0) # remove 0 data
#data_sd_sub$sd_diff_bin<- cut(data_sd_sub$sd_diff, breaks = c(min(data_sd_sub$sd_diff),0, median(data_sd_sub$sd_diff), max(data_sd_sub$sd_diff)), include.lowest = TRUE)
#medCol<-median(data_sd_sub$colonies)
#medColFinal<-ifelse(medCol==2, 3, medCol)
#data_sd_sub$colony_bin<- cut(data_sd_sub$colonies, breaks = c(0, 1, medCol, max(data_sd_sub$colonies)), include.lowest = TRUE)
#data_sd_bi<-bi_class(data_sd_sub, x = sd_diff_bin, y = colony_bin, style = "quantile", dim = 3)

# prepare empty PDF plot
pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/activity/", species, "_", behaviour, ".pdf"), width=10, height=5)

print("Mapping")

for (i in 1:length(months)) {

# Scale differences for plotting
data_sd_sub<-subset(data_sd, month==months[i])
data_sd_sub$sd_diff2<-ifelse(data_sd_sub$sd_diff>=2, 2, data_sd_sub$sd_diff) # So that we can see the other differences better

# Make a mean vs. sd plot to see if any areas with high mean and low between-pop sd
data_sd_sub$mean_scaled<-scale_0_1(data_sd_sub[[behaviorplot]])
data_sd_sub$sd_scaled<-scale_0_1(data_sd_sub[[behaviorplot2]])
data_sd_sub$meanVSsd_scaled<-data_sd_sub[[behaviorplot]]- data_sd_sub[[behaviorplot3]]

# Determine possible maximum values
#dataSub<-subset(data, month==months[i])
data_mean_sub<-subset(data_mean, month==months[i])
maxBeh1<-max(data_mean_sub[[behaviorplot]])
maxBeh2<-max(data_sd_sub[[behaviorplot3]])
maxBeh<-ifelse(maxBeh1>maxBeh2, maxBeh1, maxBeh2)

minBeh1<-min(data_mean_sub[[behaviorplot]])
minBeh2<-min(data_sd_sub$meanVSsd_scaled)
minBeh<-ifelse(minBeh1<minBeh1, minBeh1, minBeh2)

# plot Mean behaviour
plotMean<-ggplot() +
  geom_tile(data=subset(data_mean, month==months[i]), aes_string(x="x", y="y", fill=paste0(behaviorplot))) +
  geom_text(data=subset(data_sd, month==months[i] & colonies<2 & colonies >0), aes_string(x="x", y="y", label="colonies"), cex=1.5) +
  #scale_alpha(range = c(0.5, 1)) +
  scale_fill_gradientn(paste0(behaviorplot), colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$x)-400000, max(extentdf$x) + 400000), ylim=c(min(extentdf$y)-400000, max(extentdf$y) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
  
# plot SD behaviour

# Calculate quantiles for SD
# Calculate quantiles for colonie #

plotSD<-ggplot() +
geom_tile(data=subset(data_sd_sub ,month==months[i] & colonies >1), aes_string(x="x", y="y", fill="sd_diff2"), color="black") +
#geom_tile(data=data_sd_sub, aes_string(x="x", y="y", fill="sd_diff"), color="black") +
#geom_text(data=filter(data_sd_sub, colonies <2), aes_string(x="x", y="y", label="colonies"), cex=1.5) +
  #geom_tile(data=subset(data_sd, month==months[i]), aes_string(x="x", y="y", fill=paste0(behaviorplot2))) +
  #geom_text(data=subset(data_sd_bi, month==months[i]), aes_string(x="x", y="y", label="colonies"), cex=1.5) +
  scale_fill_gradientn(
  colours = c(
    "#364B9A", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#A50026"
  ),
  limits = c(min(data_sd_sub$sd_diff), max(data_sd_sub$sd_diff2)),
  values = scales::rescale(c(min(data_sd_sub$sd_diff), 1, max(data_sd_sub$sd_diff2))),
  na.value = "darkblue"
) +
  #scale_fill_gradientn(paste0(behaviorplot2), colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #bi_scale_fill(pal = "DkCyan", dim = 3) + 
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$x)-400000, max(extentdf$x) + 400000), ylim=c(min(extentdf$y)-400000, max(extentdf$y) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  labs(fill="Between/within-pop SD") +
  theme(legend.position="bottom")

plotMeanSD<-ggplot() +
  geom_tile(data=subset(data_sd_sub, month==months[i] & NoBirds >1), aes_string(x="x", y="y", fill="meanVSsd_scaled")) +
  #geom_text(data=subset(dataSub, month==months[i] & NoBirds>1), aes_string(x="x", y="y", label="NoBirds"), cex=1.5) +
  scale_fill_gradientn('Mean - SD', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"),limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$x)-400000, max(extentdf$x) + 400000), ylim=c(min(extentdf$y)-400000, max(extentdf$y) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
 
grid.arrange(plotMean, plotSD, plotMeanSD, nrow=1)
 
# Arrange the first two plots side by side
#top_row <- arrangeGrob(plotMean, plotSD, nrow = 1)

# Combine the top row with the bottom plot
#final_layout <- arrangeGrob(top_row, rangeSD, nrow = 2, heights = c(2, 1))

#plot(final_layout)
  
  }
  
  print("Saving")
  
  dev.off()


}

##### Compute Moran's i - pop-level

computeMoran_pop<-function(df_mean, behaviour) {

# Determine unique months to loop through #
months<-unique(df_mean$month)

# Assign behaviour to a specific column #
behaviorplot<-paste0("time", behaviour)
df_mean$meanVal<-df_mean[[behaviorplot]] # give column standard name so easy to calculate quantiles

# Make a list to save these values in
moranRes<-list()

for (i in 1:length(months)) {

print(i)

# Subset data frame to month i
data_mean_i<-subset(df_mean, month==months[i])

data_mean_i_arrange<-data_mean_i %>%
ungroup() %>%
arrange(x, y) %>%
dplyr::group_by(index) %>%
dplyr::slice(1)

# Calculate range of values #
minVal<-min(data_mean_i$meanVal)
maxVal<-max(data_mean_i$meanVal)
if(abs(maxVal - minVal) < 1e-6) {
next}

coords <- cbind(data_mean_i$x, data_mean_i$y) # Make a coordinate grid to choose from
# 4 nearest neighbors
knn <- knearneigh(coords, k = 4)
nb <- knn2nb(knn)
# Row-standardized weights
lw <- nb2listw(nb, style = "W")

# moran test
moran_mean <- moran.test(data_mean_i$meanVal, lw) # Looks at spatial clustering at-species level
#moran_mean_pop<-moran.test(data_sd_i$sd_diff, lw)

# Make a df with results
results<-data.frame(month=months[i])
resultName<-paste0("Moran_", behaviour)
results[[resultName]]<-moran_mean$estimate[1]

# Save results
moranRes<-rbind(moranRes, results)

}

return(moranRes)

}

##### Compute Moran's i - species level #####

computeMoran<-function(df_mean, df_sd, behaviour) {

# Determine unique months to loop through #
months<-unique(df_mean$month)

# Assign behaviour to a specific column #
behaviorplot<-paste0("time", behaviour, "Mean")
behaviorplot2<-paste0("time", behaviour, "_sd_diff")
data_mean<-df_mean
data_sd<-df_sd
data_mean$meanVal<-data_mean[[behaviorplot]] # give column standard name so easy to calculate quantiles
data_sd$sd_diff<-data_sd[[behaviorplot2]] # give column standard name so easy to calculate quantiles

# Make a list to save these values in
moranRes<-list()

for (i in 1:length(months)) {

# Subset data frame to month i
data_mean_i<-subset(data_mean, month==months[i])
data_sd_i<-subset(data_sd, month==months[i])

# Group coordinates together
data_mean_i<-data_mean_i %>%
ungroup() %>%
dplyr::arrange(x, y)

data_sd_i<-data_sd_i %>%
ungroup() %>%
dplyr::arrange(x, y)

coords <- cbind(data_mean_i$x, data_mean_i$y)

# 4 nearest neighbors
knn <- knearneigh(coords, k = 4)
nb <- knn2nb(knn)

# Row-standardized weights
lw <- nb2listw(nb, style = "W")

# moran test
moran_mean <- moran.test(data_mean_i$meanVal, lw) # Looks at spatial clustering at-species level
moran_mean_pop<-moran.test(data_sd_i$sd_diff, lw)

# Make a df with results
results<-data.frame(month=months[i], moranI_speciesMean=moran_mean$estimate[1], moranI_popSD=moran_mean_pop$estimate[1], behaviour=behaviour)

# Save results
moranRes<-rbind(moranRes, results)

}

return(moranRes)

}

### Plot for calculating correlations by month between mean & SD activity ###

corrMonth<-function(data) {

# Determine months to loop through #
months<-unique(data$month)

# Make an empty list to save the results in
dataCorRes<-list()

for (i in 1:length(months)) {

# Subset to month i #
dataMonth_i<-subset(data, month==months[i])

# Determine correlations #
corrFlight<-cor(dataMonth_i$timeFlight, dataMonth_i$timeFlight_sd, method = "pearson")
corrRestWater<-cor(dataMonth_i$timeRestWater, dataMonth_i$timeRestWater_sd, , method = "pearson")
corrLand<-cor(dataMonth_i$timeLand, dataMonth_i$timeLand_sd,  method = "pearson")
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
corrForage<-cor(dataMonth_i$timeForage, dataMonth_i$timeForage_sd, method = "pearson")
} else {
corrActive<-cor(dataMonth_i$timeActive, dataMonth_i$timeActive_sd, method = "pearson")
}

# Add these values back to the main dataset
dataMonth_i$corFlight<-corrFlight
dataMonth_i$corRestWater<-corrRestWater
dataMonth_i$corLand<-corrLand
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
dataMonth_i$corForage<-corrForage
dataMonth_i$corActive<-NA
} else {
dataMonth_i$corForage<-NA
dataMonth_i$corActive<-corrActive
}

# Save monthly values
dataMonth_i_res<-dataMonth_i %>%
dplyr::select(month, corFlight, corRestWater, corForage, corLand, corActive) %>%
dplyr::group_by(month) %>%
dplyr::slice(1)

dataCorRes<-rbind(dataCorRes, dataMonth_i_res)

}

return(dataCorRes)

}

# At species level
corrMonth_species<-function(data) {

# Determine months to loop through #
months<-unique(data$month)

# Make an empty list to save the results in
dataCorRes<-list()

for (i in 1:length(months)) {

# Subset to month i #
dataMonth_i<-subset(data, month==months[i])

# Determine correlations #
corrFlight<-cor(dataMonth_i$timeFlightMean, dataMonth_i$timeFlight_sd, method = "pearson")
corrRestWater<-cor(dataMonth_i$timeRestWaterMean, dataMonth_i$timeRestWater_sd, , method = "pearson")
corrLand<-cor(dataMonth_i$timeLandMean, dataMonth_i$timeLand_sd,  method = "pearson")
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
corrForage<-cor(dataMonth_i$timeForageMean, dataMonth_i$timeForage_sd, method = "pearson")
} else {
corrActive<-cor(dataMonth_i$timeActiveMean, dataMonth_i$timeActive_sd, method = "pearson")
}

# Add these values back to the main dataset
dataMonth_i$corFlight<-corrFlight
dataMonth_i$corRestWater<-corrRestWater
dataMonth_i$corLand<-corrLand
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
dataMonth_i$corForage<-corrForage
dataMonth_i$corActive<-NA
} else {
dataMonth_i$corForage<-NA
dataMonth_i$corActive<-corrActive
}

# Save monthly values
dataMonth_i_res<-dataMonth_i %>%
dplyr::select(month, corFlight, corRestWater, corForage, corLand, corActive) %>%
dplyr::group_by(month) %>%
dplyr::slice(1)

dataCorRes<-rbind(dataCorRes, dataMonth_i_res)

}

return(dataCorRes)

}

### Calculate map extent for spatial comparisons ###

# This function will find the indices of the xy coordinates that cover the minimum square which comprises areas where there are birds accross all populations
# The reason for this is that I need a map of similar extent to conduct analyses like moran's I or correlations that are comparable accross pops

mapExtent<-function(actFiles) {

modelColonies<-list()

print("Starting for loop")

for (i in 1:length(actFiles)) {

print(paste0("Mapping file", i, "/", length(actFiles), "..."))

#### Open file i ###
actSub<-fread(actFiles[i])
actSub<-subset(actSub, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))

# Determine the model colony
metaSub<-subset(meta, colonies==actSub$colony[1])
modelColSub<-metaSub$modelcolonies[1]
print(modelColSub)

# Does the list already have this model colony?
if (modelColSub %in% modelColonies) {
print("Next")
next
}

# subset to data frame with only birds
actSubIndex<-actSub %>%
ungroup() %>%
dplyr::group_by(month) %>%
arrange(x, y) %>%
dplyr::mutate(index=row_number()) 
actSub_crop<-subset(actSubIndex, NoBirds >0)

if (i ==1) {

# Save results
minX<-min(actSub_crop$x)
maxX<-max(actSub_crop$x)
minY<-min(actSub_crop$y)
maxY<-max(actSub_crop$y)

# Figure out which indices these are #
index1<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  dplyr::filter(near(x, minX), near(y, minY)) %>%
  dplyr::mutate(type="corner1") %>%
  dplyr::select(index, type) 
index1<-index1$index
  
index2<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, minX), near(y, maxY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner2")
index2<-index2$index
  
index3<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, minY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner3")
index3<-index3$index
  
index4<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, maxY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner4")
index4<-index4$index

}

if (i>1) {

minX_2<-min(actSub_crop$x)
maxX_2<-max(actSub_crop$x)
minY_2<-min(actSub_crop$y)
maxY_2<-max(actSub_crop$y)

if(minX_2<minX) {

minX<-minX_2
index1<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  dplyr::filter(near(x, minX), near(y, minY)) %>%
  dplyr::mutate(type="corner1") %>%
  dplyr::select(index, type) 
index1<-index1$index
  
index2<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, minX), near(y, maxY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner2")
index2<-index2$index
}

if(maxX_2>maxX) {

maxX<-ifelse(maxX_2>maxX, maxX_2, maxX)
index3<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, minY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner3")
index3<-index3$index
  
index4<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, maxY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner4")
index4<-index4$index

}

if(minY_2<minY) {

minY<-ifelse(minY_2<minY, minY_2, minY)
index1<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  dplyr::filter(near(x, minX), near(y, minY)) %>%
  dplyr::mutate(type="corner1") %>%
  dplyr::select(index, type) 
index1<-index1$index
  
index3<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, minY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner3")
index3<-index3$index

}

if (maxY_2>maxY) {

maxY<-ifelse(maxY_2>maxY, maxY_2, maxY)

index2<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, minX), near(y, maxY_2)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner2")
index2<-index2$index
  
index4<-actSubIndex %>%
#dplyr::filter(month==1) %>%
  ungroup() %>%
  dplyr::filter(near(x, maxX), near(y, maxY_2)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner4")
index4<-index4$index

}

} 

extent<-data.frame(minX, maxX, minY, maxY, index1, index2, index3, index4)

# Otherwise add the modelcolony name to the list
modelColonies<-append(modelColonies, modelColSub)

}

# Join the indices together #

extent2<-extent %>%
dplyr::rename(index=index1) %>%
dplyr::mutate(type="corner1") %>%
dplyr::select(index, type)

extent3<-extent %>%
dplyr::rename(index=index2) %>%
dplyr::mutate(type="corner2") %>%
dplyr::select(index, type)

extent4<-extent %>%
dplyr::rename(index=index3) %>%
dplyr::mutate(type="corner3") %>%
dplyr::select(index, type)

extent5<-extent %>%
dplyr::rename(index=index4) %>%
dplyr::mutate(type="corner4") %>%
dplyr::select(index, type)

extentAll<-rbind(extent2, extent3, extent4, extent5)

return(extentAll)


}

#### Calculate pair-wise population correlations ####

popcorr<-function(data, dataFiles, meta, behaviour) {

# Create behaviour we will select
behaviourSelect<-paste0("time", behaviour)

# Create an empty list with colonies 
colonies<-list()
poptopopcor<-list()

# Determine current colony
metaSub<-subset(meta, colonies==data$colony[1])
modelColSub<-metaSub$modelcolonies[1]

# Add it to empty list
colonies<-rbind(colonies, modelColSub)

# Now loop through act files & conduct pair-wise analyses

for (i in 1:length(actFiles)) {

#### Open file i ###
actSub<-fread(actFiles[i])
actSub<-subset(actSub, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))

# Add an index for easier joining
actSub<-actSub %>%
  ungroup() %>%
  arrange(month, x, y) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(index=row_number()) %>%
  ungroup() 

# Determine which colony it is
metaSub<-subset(meta, colonies==actSub$colony[1])
modelColSub<-metaSub$modelcolonies[1]

# Check to see if it exists already #
if (modelColSub %in% colonies) {
print("Next")
next
}

# otherwise we proceed to month by month comparisons
corRes<-list()

# Months loop
monthsLoop<-c(9, 10, 11, 12, 1, 2, 3, 4)

print("Conducting pair-wise comparison")

for (j in 1:length(monthsLoop)) {

# Subset dataset to month j #
actMonth<-subset(actSub, month==monthsLoop[j])

# Subset to desired extent
indices<-data.frame(index=unique(data$index))

actMonthExt<-actMonth %>%
dplyr::inner_join(indices, by=c("index"))

# Subset dataset a 
dataMonth<-subset(data, month==monthsLoop[j])

# Now I find the intersection of both (i.e. I reduce to cells with birds in both)
actCompare_birds<-subset(actMonthExt, NoBirds>0)
actOriginal_birds<-subset(dataMonth, NoBirds>0)

actCompareSub<-actCompare_birds %>%
dplyr::filter(index %in% c(unique(actOriginal_birds$index)))

actOriginalSub<-actOriginal_birds %>%
dplyr::filter(index %in% c(unique(actCompare_birds$index)))

if (nrow(actCompareSub)==0 | nrow(actOriginalSub)==0) {
next
}

# Conduct pair-wise correlation for every behaviour #
corBeh<-cor(actCompareSub[[behaviourSelect]], actOriginalSub[[behaviourSelect]], method = "pearson")

# Save result for month & population #
results<-data.frame(behaviour, cor=corBeh, month=monthsLoop[j], colony=modelColSub)
corRes<-rbind(corRes, results)

}

# Save overall results
poptopopcor<-rbind(poptopopcor, corRes)

}

poptopopcor$colony_og<-colonies[1]

return(poptopopcor)

}

### Function to extract value from within mean +/- SE for plotting ###

energyEstimate<-function(data, iterations) {

print("Summing energy iteratively...")

allEnergy<-list()

for (i in 1:iterations) {

# Generate random energy expenditure number #
data<-data %>%
dplyr::filter(!is.na(energyPopkJ_mean)) %>%
ungroup()
data$energyRandom<-rnorm(n=nrow(data), mean=data$energyPopkJ_mean, sd=data$energyPopkJ_sd)

# Sum the values per month
dataMonth<-data %>%
dplyr::group_by(month) %>%
dplyr::summarise(energyTot=sum(energyRandom))

# Add iteration No
dataMonth$iteration<-i

# Save results #
allEnergy<-rbind(allEnergy, dataMonth)

}

return(allEnergy)


}

### Function to grid behaviour for seabirds ###

gridBeh<-function(data, map) {

print("Gridding behaviour..")

# Turn map into a grid #
grid<-as.data.frame(map, xy=TRUE)

# Add columns that I want 
grid$timeFlight<-0 # time spent in flight in hours at colony 1
grid$timeForage<-0 # time spent foraging in hours at colony 1
grid$timeActive<-0 # time spent foraging in hours at colony 1
grid$timeLand<-0 # time spent foraging in hours at colony 2
grid$timeRestWater<-0 # time spent foraging in hours at colony 2
grid$timeTotal<-0 # total time spent
grid<-grid %>%
arrange(x, y) 

# Determine number of months to loop through
Months<-unique(data$month)

# Make an an empty list to save the results in
BirdActivityFinal_All<-list()

print("Month...")

for (i in 1:length(Months)) {

print(paste0(i))

# Subset to month i
monthSub<-Months[i]

# Cut actBudgets
actMonth<-subset(data, month==monthSub)

# Calculate total time in different behaviours per month
totalTime<-actMonth %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id, month) %>%
dplyr::summarise(totFlight=sum(tFlight), totActive=sum(tActive), totRest=sum(tRestWater), totLand=sum(tLand), totForage=sum(tForage), totHrs=sum(totFlight, totActive, totRest, totLand, totForage), .groups = "drop")

# Rename our grid
gridMonth<-grid %>%
dplyr::filter(x>min(actMonth$mean.lon) - (res(map) + 1) & x<max(actMonth$mean.lon) + res(map) + 1 & y > min(actMonth$mean.lat) - (res(map) +1) & y < max(actMonth$mean.lat) + res(map) +1)

for (m in 1:nrow(grid)) {

# Subset grid x
gridSub<-gridMonth[m,]

resx<-res(map)[1]
resy<-res(map)[2]

# Subset coordinates which fit #
loxSub1<-subset(actMonth,   mean.lon > gridSub$x & mean.lon < gridSub$x + resx & mean.lat > gridSub$y & mean.lat < gridSub$y + resy)

if (nrow(loxSub1)>0) {

# Calculate time spent in flight
timeFlight<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tFlight), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent foraging
timeForage<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tForage), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeLand<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tLand), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeActive<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tActive), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeRest<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tRestWater), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate total time
timeTot<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tFlight, tRestWater, tForage, tActive, tLand), .groups = "drop") %>%
dplyr::full_join(totalTime, by=c("species", "colony", "individ_id")) %>%
replace_na(list(totTime=0)) %>%
dplyr::mutate(propTime=totTime/totHrs)

#nas<-subset(timeLand, is.na(propTime))
#if(nrow(nas)>0) {break}

# Save results in grid for plotting
gridMonth$timeFlight[m]<-mean(timeFlight$propTime)
gridMonth$timeForage[m]<-mean(timeForage$propTime)
gridMonth$timeActive[m]<-mean(timeActive$propTime)
gridMonth$timeLand[m]<-mean(timeLand$propTime)
gridMonth$timeRestWater[m]<-mean(timeRest$propTime)
gridMonth$timeTotal[m]<-mean(timeTot$propTime)

# Replace NAs with 0
gridMonth[is.na(gridMonth)] <- 0

}

}

# Add number of locations
gridMonth$month<-Months[i]
gridMonth$individ_id<-data$individ_id[1]
gridMonth$rep<-data$rep[1]

# Save results
BirdActivityFinal_All<-rbind(BirdActivityFinal_All, gridMonth)

}

return(BirdActivityFinal_All)

}

### Same function as above but much faster I think ###

gridBeh2<-function(data, map) {

print("Gridding behaviour...")
  
  # Get raster
  r <- map
  resx <- res(r)[1]
  resy <- res(r)[2]
  
  # Assign each point to a raster cell
  data <- data %>%
    mutate(cell = cellFromXY(r, cbind(mean.lon, mean.lat)))
  
  # Aggregate behavior times per cell, month, species, colony, individual
  agg <- data %>%
    group_by(month, species, colony, individ_id, cell) %>%
    summarise(
      timeFlight = sum(tFlight, na.rm = TRUE),
      timeForage = sum(tForage, na.rm = TRUE),
      timeActive = sum(tActive, na.rm = TRUE),
      timeLand   = sum(tLand, na.rm = TRUE),
      timeRestWater = sum(tRestWater, na.rm = TRUE),
      timeTotal = sum(tFlight + tForage + tActive + tLand + tRestWater, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add raster coordinates for plotting
  coords <- xyFromCell(r, agg$cell)
  agg <- agg %>%
    mutate(x = coords[,1],
           y = coords[,2])
  
  # Replace NA with 0 just in case
  agg <- agg %>%
    mutate(across(c(timeFlight, timeForage, timeActive, timeLand, timeRestWater, timeTotal), ~replace_na(., 0)))
	
	agg2 <- agg %>%
  mutate(
    propFlight = timeFlight / timeTotal,
    propForage = timeForage / timeTotal,
    propActive = timeActive / timeTotal,
    propLand   = timeLand / timeTotal,
    propRestWater = timeRestWater / timeTotal
  )
  
  return(agg2)
}

### Function to carry out one type of spatial interpolation ###

run_gam_grid <- function(ActMonthlyMap, rast.template, month, beh) {

  dat <- subset(ActMonthlyMap, month == month)
  if (nrow(dat) == 0) return(NULL)

  dat$individ_id <- factor(dat$individ_id)

  xyCoords <- as.data.frame(rast.template, xy = TRUE)
  xyCoords$individ_id <- dat$individ_id[1]

  # Prepare behavior to plot #
  behPredict<-paste0("prop", beh)

  # behaviour column used dynamically
  f <- as.formula(paste0(behPredict, " ~ s(x, y) + s(individ_id, bs='re')"))

  gam_mod <- mgcv::gam(
    formula = f,
    data = dat,
    family = betar,
    weights = timeTotal
  )

  preds <- predict(gam_mod, newdata = xyCoords, exclude = "s(individ_id)", type = "response")
  xyCoords$predictions <- preds

  coordinates(xyCoords) <- ~ x + y
  pred_raster <- rasterize(xyCoords, rast.template, field = "predictions", fun = mean)

  return(pred_raster)
}

run_idw <- function(gridBehBird, rast.template, month, beh) {

  grid_sub <- subset(gridBehBird, month == month)
  if (nrow(grid_sub) == 0) return(NULL)

  # behaviour field = beh
  names(grid_sub)[names(grid_sub) == beh] <- "beh_value"

  coordinates(grid_sub) <- ~x + y
  rastDF <- as.data.frame(rast.template, xy = TRUE)
  coordinates(rastDF) <- ~x + y
  gridded(rastDF) <- TRUE

  idw_result <- gstat::idw(
    formula = beh_value ~ 1,
    locations = grid_sub,
    newdata = rastDF,
    idp = 1
  )

  idw_raster <- raster(idw_result)
  return(idw_raster)
}

run_gam_points <- function(monthRandom, rast.template, beh) {

  dat <- monthRandom
  if (!beh %in% colnames(dat)) return(NULL)

  dat_sub <- dat[dat[[beh]] > 0 & dat[[beh]] < 1, ]   # beta family restrictions
  if (nrow(dat_sub) == 0) return(NULL)

  dat_sub$individ_id <- factor(dat_sub$individ_id)

  f <- as.formula(paste0(beh, " ~ s(mean.lon, mean.lat) + s(individ_id, bs='re')"))

  gam_mod <- mgcv::gam(
    formula = f,
    data = dat_sub,
    family = betar
  )

  xyCoords <- as.data.frame(rast.template, xy = TRUE)
  xyCoords$individ_id <- dat_sub$individ_id[1]
  names(xyCoords)[names(xyCoords) == "x"] <- "mean.lon"
  names(xyCoords)[names(xyCoords) == "y"] <- "mean.lat"

  preds <- predict(gam_mod, newdata = xyCoords, exclude = "s(individ_id)", type = "response")

  pred_df <- data.frame(predictions = preds, mean.lon = xyCoords$mean.lon, mean.lat = xyCoords$mean.lat)
  coordinates(pred_df) <- ~ mean.lon + mean.lat

  pred_raster <- rasterize(pred_df, rast.template, field = "predictions", fun = mean)

  return(pred_raster)
}

run_gam_krig <- function(monthRandom, rast.template, projection_NA, beh) {

  dat <- monthRandom
  if (!beh %in% colnames(dat)) return(NULL)

  dat_sub <- dat[dat[[beh]] > 0, ]  
  if (nrow(dat_sub) == 0) return(NULL)

  dat_sub$tBeh <- dat_sub[[beh]] * 24
  coordinates(dat_sub) <- ~mean.lon + mean.lat
  proj4string(dat_sub) <- "+proj=longlat +datum=WGS84"

  dat_trans <- data.frame(spTransform(dat_sub, projection_NA))

  f <- as.formula(paste0("tBeh ~ s(coords.x1, coords.x2) + s(individ_id, bs='re')"))

  gam_mod <- mgcv::gam(f, data = dat_trans)

  rast_trans <- raster::projectRaster(rast.template, crs = projection_NA)
  xy <- as.data.frame(rast_trans, xy = TRUE)
  xy$individ_id <- dat_trans$individ_id[1]
  names(xy)[1:2] <- c("coords.x1", "coords.x2")

  preds <- predict(gam_mod, newdata = xy, exclude = "s(individ_id)", type = "response")

  pred_df <- data.frame(predictions = preds, coords.x1 = xy$coords.x1, coords.x2 = xy$coords.x2)
  coordinates(pred_df) <- ~ coords.x1 + coords.x2
  pred_raster <- rasterize(pred_df, rast_trans, "predictions")

  # Kriging
  sp_pts <- SpatialPointsDataFrame(coords = dat_trans[, c("coords.x1", "coords.x2")],
                                   data = dat_trans, proj4string = CRS(projection_NA))

  sp_pts$residLMM <- resid(gam_mod, type = "response")

  variogram_model <- vgm(model = "Sph", nugget = 0.1)
  v <- variogram(residLMM ~ 1, sp_pts)
  v_fit <- fit.variogram(v, model = variogram_model)

  rstPix <- as(rast_trans, "SpatialPixelsDataFrame")
  crs(rstPix) <- projection_NA

  krig_map <- krige(residLMM ~ 1, sp_pts, rstPix, model = v_fit)
  krig_raster <- raster(krig_map, layer = "var1.pred")

  final <- pred_raster + krig_raster

  return(final)
}

### Function for plotting where top % of data is ###

topEnergy<-function(df, percent, columnName) {

df$energy<-df[[columnName]]

# Remove missing values if any
df <- df[!is.na(df$energy), ]

# Sort descending by energy
df_sorted <- df[order(df$energy, decreasing = TRUE), ]

# Cumulative proportion of total energy
df_sorted$cum_energy <- cumsum(df_sorted$energy) / sum(df_sorted$energy)

# Keep cells contributing to 90% of total energy
df_10 <- df_sorted[df_sorted$cum_energy <= percent, ]

return(df_10)

}

## Function for plotting top % of data as contour lines ##

topEnergy_contour<-function(df, percent, columnName) {

# Change df into a xyz data frame for transforming into a raster
plot1<-df[,c("x", "y")] 
plot1$z<-df[[columnName]]

# Turn into raster
rast1<-raster::rasterFromXYZ(plot1)
v <- values(rast1)

# Rank cells by energy
ord <- order(v, decreasing = TRUE)
cum <- numeric(length(v))
cum[ord] <- cumsum(v[ord]) / sum(v)

# Put cumulative surface back into raster
r_cum <- rast1
values(r_cum) <- NA
values(r_cum)[!is.na(values(rast1))] <- cum
r_terra <- rast(r_cum)
contours <- as.contour(
  r_terra,
  levels = c(percent)
)

# Turn back into data frame
r_df<-as.data.frame(r_terra, xy=TRUE)

return(r_df)

}

### Function for extracting activity budget for a given model colony ###
extractBudget<-function(modelColonySub, speciesSub, monthSub) {

# List table with Ids & colony names #
idCatalogue<-read.csv("./results/tables/main/table1_idcatalogue.csv")

# Translate the model colony name from the text in Per's rasters to the ones in mine #
modelcoloniesTranslate<-data.frame(modelColony_Per=c("Kap_Hoegh","Bjoernoeya","Hornsund","Isfjorden", "Franz_Josef_Land", "Witless_Bay","Isle_of_May","Faroe_Islands", "Vestmannaeyjar","Papey","Grimsey","Runde_and_Aalesund"
,"Sklinna", "Roest", "Anda", "Hjelmsoeya","Hornoeya", "Skellig_Michael", "Little_Saltee" ,"Eynhallow", "Jan_Mayen","Breidafjordur" ,"Skjalfandi","Jarsteinen", "Alkefjellet",
"Sermilinnguaq", "Kippaku", "Isle_of_Canna", "Nuuk", "Langanes", "Kongsfjorden", "Anda", "Cape_Krutik", "Kara_Gate", "Russkaya_Gavan",
"Latrabjarg", "Cape_Gorodetskiy", "Coats_Island",  "Oranskie_Islands"),
modelColony_saved=c("kaphoegh","bjornoya","hornsund","isfjorden", "franzjosefland", "witlessbay","isleofmay","faroeislands", "vestmannaeyjar","papey","grimsey","rundeandaalesund"
,"sklinna", "rost", "anda", "hjelmsoya","hornoya", "skelligmichael", "littlesaltee" ,"eynhallow", "janmayen","breidafjordur" ,"skjalfandi","jarsteinen", "alkefjellet",
"sermilinnguaq", "kippaku", "isleofcanna", "nuuk", "langanes", "kongsfjorden", "anda", "capekrutik", "karagate", "russkayagavan",
"latrabjarg", "capegorodetskiy", "coatsisland",  "oranskieislands"),  
modelColony_database=c("Kap Höegh", "Bjørnøya", "Hornsund", "Isfjorden", "Franz Josef Land", "Witless Bay", "Isle of May", "Faroe Island", "Vestmannaeyjar", "Papey", "Grimsey", "Runde and Ålesund", "Sklinna", "Røst", "Anda",
"Hjemlsøya", "Hornøya", "Skellig Michael", "Little Saltee", "Eynhallow", "Jan Mayen", "Breidafjordur", "Skjalfandi", "Jarsteinen", "Alkefjellet",
"Sermilinnguaq", "Kippaku", "Isle of Canna", "Nuuk", "Langanes", "Kongsfjorden", "Anda", "Cape Krutik", "Kara Gate", "Russkaya Gavan",
"Latrabjarg", "Cape Gorodetskiy", "Coats Island",  "Oranskie Islands"))

# Find correct model colony translation
modelColMatch<-subset(modelcoloniesTranslate, modelColony_Per==modelColonySub)$modelColony_database[1]
modelSave<-subset(modelcoloniesTranslate, modelColony_Per==modelColonySub)$modelColony_saved[1]

if(speciesSub=="Northern fulmar" & modelSave=="jarsteinen") {
modelSave<-"karmoy"
}

# Make sure we were able to find a colony match #
if (length(modelColMatch)<1) stop(print("Error: no match with model colony name"))

# Transform species sub into how it's written in the saved files
speciesMatch<-gsub(" ", "", speciesSub)
speciesMatch<-gsub("-", "", speciesMatch)
speciesMatch<-gsub("ü", "u", speciesMatch)
speciesMatch<-gsub("'", "", speciesMatch)

# But then also filter for files
allBudgets<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3/", full.names=TRUE)
allBudgets_colony<-allBudgets[grep(paste0("_", modelSave, "_"), allBudgets)] # Subset to colony of interest
allBudgets_monthly<-allBudgets_colony[grep("monthly", allBudgets_colony)] # Subset to monthly data
allBudgets_species<-allBudgets_monthly[grep(speciesMatch, allBudgets_monthly)]

# make sure we identified only one budget file #
if (!length(allBudgets_species) %in% c(1)) stop (print("Error: wrong number of budget files"))

# Now we open it #
budgets<-read.csv(allBudgets_species[1])

# Filter to month of interest #
budgets_sub<-subset(budgets, month==monthSub)

# Filter out birds with not enough days & sum per month #
birdMonth<-budgets_sub %>%
ungroup() %>%
dplyr::group_by(species, colony, month) %>%
dplyr::summarise(FlightHrs_mean=mean(tFlightMean), LandHrs_mean=mean(tLandMean), RestWaterHrs_mean=mean(tRestWaterMean), ActiveHrs_mean=mean(tActiveMean),
ForageHrs_mean=mean(tForageMean), .groups="drop")

return(birdMonth)

}

### Function for making randomized polygons. They are made in a new location and rotated. It's so I can compare values in these fake ones to real ones ###

# Function made by chat gpt so let's see how this goes

st_rotate <- function(x, angle_deg, center = st_centroid(st_union(x))) {
  angle <- angle_deg * pi/180
  M <- matrix(c(cos(angle), -sin(angle),
                sin(angle),  cos(angle)), nrow = 2, byrow = TRUE)

  ctr <- st_coordinates(center)[1, 1:2]

  g <- st_geometry(x)
  g_rot <- (g - ctr) * M + ctr

  st_set_geometry(x, g_rot)
}

make_random_poly <- function(patchNo, numberPolys, patchProj, availTrans) {

print("Generating random patches...")

# PatchProj is projected patch #

for (i in 1:numberPolys) {

repeat{

# Rotate polygon at random #
patchProj <- st_make_valid(patchProj)
#rotateRandom <- rotate.polygon(patchProj, angle=sample(0:360, 1))
rotateRandom <- st_rotate(patchProj, sample(0:360, 1))

# Find out extent of terra raster
e <- ext(availTrans)

xmin <- e$xmin
xmax <- e$xmax
ymin <- e$ymin
ymax <- e$ymax

# Determine maximum distances this patch can move by #
distancex<-xmax-xmin
distancey<-ymax-ymin

# Choose random amount to move patch by
dx <-sample(0:distancex, 1)
dy <-sample(0:distancey, 1)   # meters south

# Move polygon
poly_moved <- rotateRandom
st_geometry(poly_moved) <- st_geometry(rotateRandom) + c(dx, dy)
st_crs(poly_moved)<-projection_NA # Specify projection

# Extract values from our availability raster
vpoly <- vect(poly_moved)   # convert sf -> terra vector
vals <- terra::extract(availTrans, vpoly)
vals[is.na(vals)] <- 0 # to make my ifelse statement below easier

# all(vals == TRUE)

if (all(vals == 1)) {
break
} 

}

# Create an sf to return?

if (i==1) {

ids<-c(paste0(patchNo, "_observed"), paste0(patchNo, "_control_1"))
polys <- rbind(patchProj, poly_moved)

} else {

newId<-paste0(patchNo, "_control_", i)
ids<-append(ids, newId)
polys<-rbind(polys, poly_moved)

}

}

# Give the polygons new names #
polysFinal<-polys %>%
dplyr::mutate(patches=ids)

return(polysFinal)

}
