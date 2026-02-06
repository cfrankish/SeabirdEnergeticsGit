### Here I conduct a sensitivity analysis  ###
### Basically energy is caculated while I vary the values of different parameters many times ###
### This is done for every individual ###
### The input files are individual wet-dry csv files (birddata_ind/speciesname/id.csv) ###
### The output files are individual daily energy files ('tmp2/id_energyDay.csv') ###

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
library(tibble)

### Step 0: Open input file ####

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
#print(input_file)
#birds<-list.files("./data/birddata_ind/atlanticpuffin/", full.names=TRUE)
#birds<-list.files("./data/birddata_ind/northernfulmar/", full.names=TRUE)
#birds<-list.files("./data/birddata_ind/blackleggedkittiwake/", full.names=TRUE)
#birds<-list.files("./data/birddata_ind/commonguillemot/", full.names=TRUE)
#birds<-list.files("./data/birddata_ind/brunnichsguillemot/", full.names=TRUE)
#birds<-list.files("./data/birddata_ind/littleauk/", full.names=TRUE)
#dataSpeciesIdSub<-fread(birds[1])
dataSpeciesIdSub <- fread(input_file)	

#### Step 1: assign location of files & functions ####

# Source all necessary functions
source("./scripts/functions.R")

# Determine location of processed locations (IRMA data) 
irma.files<-list.files("./data/positionsIRMA/", full.names=TRUE)
metadata<-load(irma.files[8])
irma.files<-irma.files[grepl("IRMAlocs", irma.files)]
irma.files.df<-data.frame(irma.files)
colnames(irma.files.df)<-c("FileName")
irma.files.df$species<-c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot")

# Specify location of other metadata (information on body weights)
list.activity.meta<-list.files("./data/metadata/", full.names=TRUE)

# Determine model parameters to choose from #
speciesNo<-6
paramNo<-15
modelParams<-tibble(species=rep(c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot"), paramNo))
modelParams$parameter<-rep(c("L1", "Th1", "Th2", "L1_colony", "dist_colony", "pLand_prob", "c",
"RMR", "c1", "c2", "c3", "c4", "TC", "Beta_active", "Beta_rest"), each=6)
modelParams$values<-list(
c(240, 240-240*0.1, 240+240*0.1), c(810, 810-0.1*810, 810 + 0.1*810), c(90, 90-0.1*90, 90+0.1*90), c(134, 134-0.1*134, 134+0.1*134), c(88, 88-0.1*88, 88+0.1*88), c(88, 88-0.1*88, 88+0.1*88), # Species-specific flight bout duration (minutes)
c(0.965, 0.95, 0.98), c(0.965, 0.95, 0.98), c(0.885, 0.84, 0.93), c(0.885, 0.84, 0.93), c(0.885, 0.84, 0.93), c(0.885, 0.84, 0.93), # Th1 % wet threshold for differenciating between behaviors
c(0.05, 0, 0.10), c(0.05, 0, 0.10), c(0.05, 0, 0.10), c(0.05, 0, 0.10), c(0.05, 0, 0.10), c(0.05, 0, 0.10), # Th2 #% wet for determing dry vs. intermediate
c(0, 240), c(0, 810), c(0, 90), c(0, 134), c(0, 88), c(0, 88), # L1_colony: duration of dry bouts at start of tLand that can re-allocated to flight
c(500, 500-500*0.1, 500 + 500*0.1), c(500, 500-500*0.1, 500 + 500*0.1), c(500, 500-500*0.1, 500 + 500*0.1), c(500, 500-500*0.1, 500 + 500*0.1), c(500, 500-500*0.1, 500 + 500*0.1), c(500, 500-500*0.1, 500 + 500*0.1),  # Distance to colony (km) below which it is considered possible to be on land
c(0.5, 0.45, 0.55), c(0.5, 0.45, 0.55), c(0.5, 0.45, 0.55), c(0.5, 0.45, 0.55), c(0.5, 0.45, 0.55), c(0.5, 0.45, 0.55), # pland: probability of dry being land or something else
c(2.25, 1.5, 3), c(2.25, 1.5, 3), c(2.25, 1.5, 3), c(2.25, 1.5, 3), c(2.25, 1.5, 3), c(2.25, 1.5, 3), # coefficient for adjusting for leg-tucking
c(1.64, 1.60, 1.68), c(1.00, 0.96, 1.04), c(0), c(0), c(0), c(0), # RMR is just for a few of the species...BlKi / NoFu (Resting metabolic rate)
c(3.88, 2, 5.7), c(2.2, 2.2-2.2*0.1, 2.2+2.2*0.1), c(141, 106, 176), c(8.45, 5.4, 11.4), c(141, 106, 176), c(141, 106, 176), # c1 is the cost of Flight
list(c(2.06, 3.88) ,c(0.3, 2),c(3.8, 5.7)), list(c(2.06, 3.88) ,c(0.3, 2),c(3.8, 5.7)), c(0), c(0), c(0), c(0), # c2 is the cost of foraging...
c(0.56, 0.1, 1.1), c(0.8, 0.8-0.8*0.1, 0.8+0.8*0.1), c(9.2, 2.7, 15.3), c(9.2, 2.7, 15.3), c(9.2, 2.7, 15.3), c(9.2, 2.7, 15.3), # c3 is cost of being on land...
c(1.24, 0.1, 2.8), c(2, 2-2*0.1, 2+2*0.1), c(9.2, 2.7, 15.3), c(9.2, 2.7, 15.3), c(9.2, 2.7, 15.3), c(9.2, 2.7, 15.3), # c4 is the cost of resting on the water...
c(0.0997, 0.0997-0.1*0.0997, 0.0997 + 0.1*0.0997), c(0.07, 0.07-0.1*0.07, 0.07 + 0.1*0.07), # TC is thermal conductivity (code is TC - 1.96*TCError)
c(2.75, 2.75-0.1*2.75, 2.75 + 0.1*2.75), c(2.75, 2.75-0.1*2.75, 2.75 + 0.1*2.75),
c(2.75, 2.75-0.1*2.75, 2.75 + 0.1*2.75), c(2.75, 2.75-0.1*2.75, 2.75 + 0.1*2.75),
c(0), c(0), c(118, 118-118*0.1, 118 + 118*0.1),  # Intercepts of resting metabolic rate at 0°C during different activities (active)
c(118, 118-118*0.1, 118 + 118*0.1), c(118, 118-118*0.1, 118 + 118*0.1), c(118, 118-118*0.1, 118 + 118*0.1),
c(1.87, 1.87-1.87*0.1, 1.87 + 1.87*0.1), c(1.34, 1.34-1.34*0.1, 1.34+1.34*0.1),# # Intercepts of resting metabolic rate at 0°C during different activities (rest)
c(72.2, 72.2-72.2*0.1, 72.2 + 72.2*0.1), c(72.2, 72.2-72.2*0.1, 72.2 + 72.2*0.1),
c(72.2, 72.2-72.2*0.1, 72.2 + 72.2*0.1), c(72.2, 72.2-72.2*0.1, 72.2 + 72.2*0.1))

#### Step 2: estimate winter activity & energy budgets ####

# Set up lists to save results
energyDay<-list() # Daily energy estimates

# Determine Species
speciesSub<-dataSpeciesIdSub$species[1]

# Determine colony
colonySub<-dataSpeciesIdSub$colony[1]

# Open species-specific IRMA data (all individuals)
irma.file.sub<-subset(irma.files.df, species==speciesSub)
irmaFile<-readRDS(irma.file.sub$FileName)

# Open up species-specific metadata 
speciesName<-speciesSub
speciesName<-gsub(" ", "", speciesName)
speciesName<-gsub("ü", "u", speciesName)
metaSub<-readRDS(list.activity.meta[grepl(speciesName, list.activity.meta)])
    
# Determine bird id
idSub<-dataSpeciesIdSub$individ_id[1]  
      
# Here i make some modifications so that id in my dataset matches the IRMA dataset
idIRMA<-dataSpeciesIdSub$individ_id_IRMA[1]
        
# Subset corresponding IRMA file
irmaSub<-subset(irmaFile, individ_id==idIRMA)
        
	# If no match then stop as problem
	if (nrow(irmaSub)<1) {
			  
	stop(print("No matching IRMA file"))
			  
	}

# Determine sample size
sampleSize<-dataSpeciesIdSub %>%
dplyr::mutate(month=as.numeric(substr(date_time, 6, 7))) %>%
ungroup() %>%
dplyr::group_by(session_id) %>%
dplyr::summarise(months=n_distinct(month))

# Remove sessions with < 3 months (because probably means there is no full month available)
sampleSize_select<-subset(sampleSize, months>=3)

# Subset dataset to these sessions (just to reduce on storage space)
sessions<-subset(dataSpeciesIdSub, session_id %in% sampleSize_select$session_id)

### Make sure this bird has data within study period ###

# Determine study period
dates<-data.frame(dateKeep=seq(as.Date("2021-09-15"), as.Date("2022-04-15"), 1))
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
sessions$month<-as.numeric(substr(sessions$date, 6, 7))
sessions$year<-as.numeric(substr(sessions$date, 1, 4))
sessions$track_year<-ifelse(sessions$month<day1_month, paste0(sessions$year-1, "_", substr(sessions$year, 3, 4)), paste0(sessions$year, "_", substr(sessions$year + 1, 3, 4)))
sessions$session_year<-paste0(sessions$session_id, "_", sessions$track_year)
 
# Extract all unique days from our dataset
dates_sessions<-sessions %>%
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

if (nrow(dates_sessions) > 0) {

# Subset model parameters to our species of interest 
modelParamsSub<-subset(modelParams, species==dataSpeciesIdSub$species[1])

# Loop through parameter names
paramNames<-unique(modelParams$parameter)
# Add three extra parameters that are handled differently
paramNamesFinal<-c(paramNames, "year", "sst", "ice")

for (i in 1:length(paramNamesFinal)) {
      
print(paste("Parameter", i, "/", length(paramNamesFinal), "..."))  

if (paramNamesFinal[i] %in% c("year", "sst", "ice")) {

if (paramNamesFinal[i]=="year") {

# I figure out how many 'track years' I can cycle through... # 
ParamValues<-n_distinct(dates_sessions$session_year)

} else {


ParamValues<-3

}


} else {

if (paramNamesFinal[i]=="c5") {

# Subset specific parameter
ParamsSub<-subset(modelParamsSub, parameter==paramNamesFinal[i])

# Determine how many values to loop through
ParamValues<-n_distinct(ParamsSub$values[[1]][[1]])

} else {

# Subset specific parameter
ParamsSub<-subset(modelParamsSub, parameter==paramNamesFinal[i])

# Determine how many values to loop through
ParamValues<-n_distinct(ParamsSub$values[[1]])

}

}

# And then we loop through the number of values available for a parameter

for (j in 1:ParamValues) {

### Add model parameters to our dataset ###
print("Randomizing activity parameters...")
sessions$L1<-subset(modelParamsSub, parameter=="L1")$values[[1]][1] # Fixed value
sessions$Th1<-subset(modelParamsSub, parameter=="Th1")$values[[1]][1] # Sampled at random
sessions$Th2<-subset(modelParamsSub, parameter=="Th2")$values[[1]][1] # Fixed value
sessions$L1_colony_min<-subset(modelParamsSub, parameter=="L1_colony")$values[[1]][[1]] # Fixed value
sessions$L1_colony_max<-subset(modelParamsSub, parameter=="L1_colony")$values[[1]][[1]] # Fixed value
sessions$ice_random<-sessions$ice_mean # I will change this later to the first value of every dry bout (happens in the activity functions)
sessions$dist_colony<-subset(modelParamsSub, parameter=="dist_colony")$values[[1]][[1]] # Fixed value
sessions$pLand_prob<-subset(modelParamsSub, parameter=="pLand_prob")$values[[1]][[1]] # Fixed value
sessions$c<-subset(modelParamsSub, parameter=="c")$values[[1]][[1]] # Sampled at random
sessions$sst_random_start<-sessions$sst_lox_mean

# How to deal with most activity parameters #

if (paramNamesFinal[i] %in% c("L1", "Th1", "Th2", "dist_colony", "pLand_prob", "c")) {

# Cycle through the different values we are trying out!
sessions[[paramNamesFinal[i]]]<-subset(modelParamsSub, parameter==paramNamesFinal[i])$values[[1]][j]

}

if (paramNamesFinal[i] %in% c("L1_colony")) {

# Cycle through the different values we are trying out!
sessions$L1_colony_min<-subset(modelParamsSub, parameter==paramNamesFinal[i])$values[[1]][j]
sessions$L1_colony_max<-subset(modelParamsSub, parameter==paramNamesFinal[i])$values[[1]][j]

}

# Dealing with sst

if (paramNamesFinal[i] %in% c("sst")) {

if (j==1) {

# Cycle through the different values we are trying out!
sessions$sst_random_start<-sessions$sst_lox_mean

} 

if (j==2) {

sessions$sst_random_start<-sessions$sst_lox_mean

}

if (j==3) {

sessions$sst_random_start<-sessions$sst_lox_mean

}

}

# Dealing with ice

if (paramNamesFinal[i] %in% c("ice")) {

if (j==1) {

# Cycle through the different values we are trying out!
sessions$ice_random<-sessions$ice_mean

} 

if (j==2) {

sessions$ice_random<-sessions$ice_mean

}

if (j==3) {

sessions$ice_random<-sessions$ice_mean

}

}

# Dealing with year

if (paramNamesFinal[i] %in% c("year")) {

sessionsSelect<-unique(dates_sessions$session_year)
sessionSub<-sessionsSelect[j]
sessionRandom<-subset(sessions, session_year %in% sessionSub)
yearSave<-sessionsSelect[j]

}

# Otherwise we set the session to the first one in the list for all rounds

if (!paramNamesFinal[i] %in% c("year")) {

sessionsSelect<-unique(dates_sessions$session_year)
sessionSub<-sessionsSelect[1]
sessionRandom<-subset(sessions, session_year %in% sessionSub)
yearSave<-sessionsSelect[1]

}

### Calculate time in activity ####
print("Calculating Activity...")

act<-calculateTimeInActivity(species=sessions$species[1], data=sessionRandom, irmaData=irmaSub)

# Extract daily results, summarize by month & save 
if (is.data.frame(act)==TRUE) {
          
actRes<-act  
          
  } else {
          
actRes<-act[[1]]
actRes_prop<-act[[2]]

  }

nasAct<-subset(actRes, is.na(tLand))

if(nrow(nasAct)>0) {
  
  
  #stop(print("STOP"))
  next
}
        
# Remove incomplete days
maxDuration<-max(actRes$DurationTot)
actRes<-subset(actRes, !DurationTot<=23)
        
actRes$Month<-as.numeric(substr(actRes$date, 6, 7))
colonySub<-actRes$colony[1]
actRes$colony<-colonySub
actRes$individ_id<-idSub

# Add temperature at colony
tempColony<-dataSpeciesIdSub %>%
dplyr::select(date_time, sst_col_mean, sst_col_sd) %>%
dplyr::group_by(date_time) %>%
dplyr::mutate(sst_random_colony=sst_col_mean) %>%
ungroup() %>%
dplyr::mutate(date=substr(date_time, 1, 10)) %>%
dplyr::group_by(date) %>%
dplyr::summarise(sst_random_colony=mean(sst_random_colony)) %>%
dplyr::mutate(sst_random_colony=ifelse(sst_random_colony < (-1.9), -1.9, sst_random_colony))

actRes<-actRes %>%
dplyr::left_join(tempColony, by=c("date"))
        
print("Removing incomplete months...")
		
sstDay_current<-actRes
		
# Determine start month/year
start.month<-substr(min(sstDay_current$date), 6, 7)
start.year<-substr(min(sstDay_current$date), 1, 4)
			
# Determine end month/year
end.month<-substr(max(sstDay_current$date), 6, 7) 
end.year<-substr(max(sstDay_current$date), 1, 4)
end.date<-ifelse(end.month %in% c("01", "03", "05", "07", "08", "10", "12"), "31", "30")
end.date<-ifelse(end.month=="02", "28", end.date)
			
# Make a sequence of dates
seq.dates<-data.frame(date=seq(as.Date(paste0(start.year, "-", start.month, "-01")), as.Date(paste0(end.year, "-", end.month, "-", end.date)), 1))
sstDay<-rbind(sstDay_current)
sstDay$date<-as.Date(sstDay$date)
sstDay$sstModel<-"ECMWF_reanalysis"
sstDay$time<-"current"
sstJoin<-seq.dates %>%
	dplyr::full_join(sstDay, by=c("date")) %>%
	dplyr::mutate(nas=ifelse(is.na(species), 1, 0)) %>%
	dplyr::mutate(Month=as.numeric(substr(date, 6, 7))) %>%
	dplyr::mutate(Year=as.numeric(substr(date, 1, 4))) %>%
	dplyr::ungroup() %>%
	dplyr::group_by(Year, Month) %>%
	dplyr::mutate(nasSum=sum(nas)) %>%
	dplyr::filter(nasSum==0) 
        
if (nrow(sstJoin) <1 ) {
          
 stop(print("No data"))      
          
}
        
### Summarize activity & sst at monthly level ####
        
actResMonth<-sstJoin %>%
  dplyr::group_by(species, individ_id, session_id, colony, sstModel, time, Month, Year) %>%
  dplyr::mutate(tActive=ifelse(species %in% c("Black-legged kittiwake", "Northern fulmar"), 0, tActive)) %>%
  dplyr::mutate(EnergyDiveTot=0) %>%
  replace_na(list(maxFlightBoutsMins=0)) %>%
  dplyr::summarise(tForage_month=sum(tForage), tLand_month=sum(tLand), tFlight_month=sum(tFlight), tRestWater_month=sum(tRestWater), tActive_month=sum(tActive), sstMonth=mean(sst_random, na.rm=TRUE), dayLength=mean(dayLengthHrs), MaxDistColKm=max(MaxDistColKm),
                   EnergyDiveTotMonth=sum(EnergyDiveTot), totDays=n_distinct(date), PropForage=tForage_month/(totDays*24), PropFlight=tFlight_month/(totDays*24),
                   PropActive=tActive_month/(totDays*24), PropLand=tLand_month/(totDays*24), PropRest=tRestWater_month/(totDays*24), maxFlight=max(maxFlightBoutsMins, na.rm=TRUE),
                   iceMonth=mean(ice_random, na.rm=TRUE), distColKmMonth=mean(distColonyKm_mean)) %>%
  ungroup() %>%
  dplyr::mutate(propTot=PropLand + PropActive + PropForage + PropRest + PropFlight)
      
### Add model parameters to our dataset ###
print("Randomizing energy parameters...")
actRes$RMR<-subset(modelParamsSub, parameter=="RMR")$values[[1]][[1]] # choose at random from uniform distribution
actRes$c1<-subset(modelParamsSub, parameter=="c1")$values[[1]][[1]] # choose at random from uniform distribution
actRes$c2<-mean(subset(modelParamsSub, parameter=="c2")$values[[1]][[1]]) # choose at random from uniform distribution
actRes$c3<-subset(modelParamsSub, parameter=="c3")$values[[1]][[1]] # choose at random from uniform distribution
actRes$c4<-subset(modelParamsSub, parameter=="c4")$values[[1]][[1]] # Choose at random from uniform distribution
actRes$TC<-subset(modelParamsSub, parameter=="TC")$values[[1]][[1]] # Choose at random from uniform distribution
actRes$Beta_active<-subset(modelParamsSub, parameter=="Beta_active")$values[[1]][[1]] # Choose at random from uniform distribution
actRes$Beta_rest<-subset(modelParamsSub, parameter=="Beta_rest")$values[[1]][[1]] # Choose at random from uniform distribution

if (paramNamesFinal[i] %in% c("RMR", "c1", "c2", "c3", "c4", "TC", "Beta_active", "Beta_rest")) {

# Cycle through the different values we are trying out!
actRes[[paramNamesFinal[i]]]<-mean(subset(modelParamsSub, parameter==paramNamesFinal[i])$values[[1]][j][[1]])

}

### Calculate energetics (monthly) ####
print("Calculating energetics")

### Calculate energetics (daily) ####
energyDaily<-calculateEnergetics(species=speciesSub, data=actRes, colonySub=actRes$colony[1], sstVals=actRes, type="daily", type2="ind")
energyDaily$doy<-floor(as.numeric(difftime(energyDaily$date, paste0(substr(energyDaily$date, 1, 4), "-01-01"),  unit=c("days"))) + 1)#energyDay<-rbind(energyDay, energyDaily)

# Add all parameter values before Saving
energyDaily$L1<-sessions$L1[1]
energyDaily$Th1<-sessions$Th1[1]
energyDaily$Th2<-sessions$Th2[1]
energyDaily$L1_colony_min<-sessions$L1_colony_min[1]
energyDaily$L1_colony_max<-sessions$L1_colony_max[1]
energyDaily$dist_colony<-sessions$dist_colony[1]
energyDaily$pLand_prob<-sessions$pLand_prob[1]
energyDaily$c<-sessions$c[1]
energyDaily$RMR<-actRes$RMR[1]
energyDaily$c1<-actRes$c1[1]
energyDaily$c2<-actRes$c2[1]
energyDaily$c3<-actRes$c3[1]
energyDaily$c4<-actRes$c4[1] # This is equal to c3 for the auks... 
energyDaily$TC<-actRes$TC[1]
energyDaily$Beta_active<-actRes$Beta_active[1]
energyDaily$Beta_rest<-actRes$Beta_rest[1]
energyDaily$year<-yearSave

# Add modelling information
energyDaily$Parameter<-paramNamesFinal[i] # Parameter being tested
energyDaily$Run<-j # Iteration value which should cross-reference with my manuscript (run 1 is the control, and run 2-j are the test values)

energyDay<-rbind(energyDay, energyDaily)

} # End of j (parameter values)

} # end of i (parameter names)

# Save result per bird

#### Save outputs ####

print("Saving file...")
output_file1 <- args[2]

print("1")
write.csv(energyDay, file = output_file1, row.names = FALSE) # Daily energy data 
print("Finished")

} else {

print("data doesn't fall within study period...")

output_file1 <- args[2]
energyDay<-data.frame()

print("1")
write.csv(energyDay, file = output_file1, row.names = FALSE) # Daily energy data 
print("Finished")

}