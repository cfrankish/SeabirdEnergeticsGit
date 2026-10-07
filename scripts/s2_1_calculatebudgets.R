### Here I calculate time-activity budgets and energetics ###
### This is done many times for all session_ids at a daily and monthly scale ###
### Input file is individual csv wet-dry file stored in "./data/birddata_ind/speciesname/" ##
### Output file is an individual csv with activity and energy calculated for every day of tracking data ##
### This script also outputs some plots showing how this varies accross the tracking period (just for checking purposes)- stored in 'results/tempPlots'##

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

### Step s2_1_0: Open input file ####

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
dataSpeciesIdSub <- fread(input_file)	

# Set-up number of iterations...
overall.iterations<-100 # how many times this is calculated per individual

#### Step s2_1_1: assign location of files & functions ####

# Source all necessary functions for calculating activity budgets and energy expenditure
source("./scripts/functions.R")

# Determine location of processed locations (IRMA data) 
irma.files<-list.files("./data/positionsIRMA/", full.names=TRUE)
metadata<-load(irma.files[8])
irma.files<-irma.files[grepl("IRMAlocs", irma.files)]
irma.files.df<-data.frame(irma.files)
colnames(irma.files.df)<-c("FileName")
irma.files.df$species<-c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot")

# Determine model parameters to choose from #
# Note: for activity cost multipliers with a range of values, these correspond to lower and upper 95% confidence intervals (see Table S3 in publication)
speciesNo<-6
paramNo<-19
modelParams<-tibble(species=rep(c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot"), paramNo))
modelParams$parameter<-rep(c("L1", "Th1", "Th2", "L1_colony", "dist_colony", "pLand", "c",
"RMR", "c1", "c2", "c3", "c4", "c5", "TC_water", "Beta_active", "Beta_rest", "LCT_water", "TC_air", "LCT_air"), each=speciesNo)
modelParams$values<-list(
240, 810, 90, 134, 88, 88, # Species-specific flight bout duration (minutes)
0.95, 0.95, c(0.85, 0.9), c(0.85, 0.9), c(0.85, 0.9), c(0.85, 0.9), # Th1 % wet threshold for differenciating between behaviors
0, 0, 0, 0, 0, 0, # Th2 #% wet for determing dry vs. intermediate
c(0, 240), c(0, 810), c(0, 90), c(0, 134), c(0, 88), c(0, 88), # L1_colony: duration of dry bouts at start of tLand that can be re-allocated to flight in minutes
500, 500, 500, 500, 500, 500,  # Distance to colony (km) below which it is considered possible to be on land
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, # pland: probability of dry being land or something else
c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), c(seq(1.5, 3, 0.1)), # coefficient for adjusting for leg-tucking
c(0), c(seq(0.96, 1.04, 0.01)), c(0), c(0), c(0), c(0), # RMR is just for NoFu and stands for Resting metabolic rate in ml O2 g-1 hr-1
c(seq(2, 5.7, 0.1)), c(2.2), c(seq(106, 176, 1)), c(seq(5.4, 11.5, 0.1)), c(seq(106, 176, 1)), c(seq(106, 176, 1)), # c1 is the cost of Flight in source units (see Table S3 in publication)
c(seq(0.3, 3.8, 0.1)), c(seq(0.3, 3.8, 0.1)), c(0), c(0), c(0), c(0), # c2 is the cost of foraging in source units (see Table S3 in publication). Note, here is shown cost of swimming. c2 for Blki & Nofu is an average of these values and cost of flapping (c1). An averaging occurs later in the script.  
c(seq(0.05, 1.1, 0.1)), c(0.8), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), # c3 is cost of being on land in source units (see Table S3 in publication). 
c(seq(0.1, 2.8, 0.1)), c(2), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), c(seq(3.1, 15.3, 0.1)), # c4 is the cost of resting on the water in source units (see Table S3 in publication). Note the lower limit of Blki was adjusted to 0.1 to avoid being negative, and it represents cost of rest (see Tremblay et al. 2024). c4 is an average of swimming and rest and the averaging occurs later in the script. 
c(0), c(0), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), # c5 is the cost of being active on water when thermoneutral...(these values are no longer used, c5 is made in the code)
c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), c(list(seq(26 - 1.96*6, 29, 0.01), seq(27 - 1.96*12, 29, 0.01))), 
c(0.0997), c(0.07), # TC is thermal conductivity in water in kJ g-1 hr-1 C-1 for Blki and NoFu & kJ hr-1 C-1 for auks
c(2.75), c(2.75),
c(2.75), c(2.75),
c(0), c(0), c(118),  # Intercepts of resting metabolic rate at 0°C during different activities (active) in kJ hr-1
c(118), c(118), c(118),
c(0), c(0),# # Intercepts of resting metabolic rate at 0°C during different activities (rest) in kJ hr-1
c(72.2), c(72.2),
c(72.2), c(72.2),
12.5, 9, 14.18, 14.18, 14.18, 14.18, # LCT_water is lower critical temperature in water in degrees C
0.0466, 0.0336, 0.0282, 0.05, 0.0282, 0.0282, # TC is thermal conductivity in air in ml O2 g-1 hr-1 C-1
4.5, 9, 5.72, 4.5, 2, 2) # LCT_air is lower critical temperature in air in degrees C

#### Step s2_1_2: estimate winter activity & energy budgets ####

# Set up lists to save results
energyMonth<-list() # Monthly energy estimates
energyDay<-list() # Daily energy estimates
actRes_iterations<-list() # For plotting purposes
activityMonth<-list() # For plotting purposes - time in activity per month

# Determine species of interest
speciesSub<-dataSpeciesIdSub$species[1]

# Determine colony of interest
colonySub<-dataSpeciesIdSub$colony[1]

# Open species-specific IRMA data (all individuals)
irma.file.sub<-subset(irma.files.df, species==speciesSub)
irmaFile<-readRDS(irma.file.sub$FileName)

# Open up species-specific metadata 
speciesName<-speciesSub
speciesName<-gsub(" ", "", speciesName)
speciesName<-gsub("ü", "u", speciesName)
    
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

# Subset model parameters to our species of interest 
modelParamsSub<-subset(modelParams, species==dataSpeciesIdSub$species[1])

# Loop through iterations 

for (i in 1:overall.iterations) {
      
print(paste("iteration", i, "/", overall.iterations, "..."))  

### Add model parameters to our dataset ###
print("Randomizing activity parameters...")
sessions$L1<-subset(modelParamsSub, parameter=="L1")$values # Fixed value
sessions$Th1<-sample(subset(modelParamsSub, parameter=="Th1")$values[[1]], 1) # Sampled at random
sessions$Th2<-subset(modelParamsSub, parameter=="Th2")$values # Fixed value
sessions$L1_colony_min<-min(subset(modelParamsSub, parameter=="L1_colony")$values[[1]]) # Fixed value
sessions$L1_colony_max<-max(subset(modelParamsSub, parameter=="L1_colony")$values[[1]]) # Fixed value
sessions$ice_random<-sessions$ice_mean # I will change this later to the first value of every dry bout (happens in the activity functions)
sessions$dist_colony<-subset(modelParamsSub, parameter=="dist_colony")$values # Fixed value
sessions$pLand_prob<-subset(modelParamsSub, parameter=="pLand")$values # Fixed value
sessions$c<-sample(subset(modelParamsSub, parameter=="c")$values[[1]], 1) # Sampled at random
sessions$sst_random_start<-sessions$sst_lox_mean
sessions$sst_random_start<-ifelse(sessions$sst_random_start < (-1.9), -1.9, sessions$sst_random_start) # Make sure no un-realistic value

### Calculate time in activity ####
print("Calculating Activity...")

act<-calculateTimeInActivity(species=sessions$species[1], data=sessions, irmaData=irmaSub)

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

# Add temperature at colony to main dataset
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

actRes$rep<-i
actRes_iterations<-rbind(actRes_iterations, actRes)
        
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
        
### Summarize activity & sst at monthly level for plotting purposes ####
        
actResMonth<-sstJoin %>%
  dplyr::group_by(species, individ_id, session_id, colony, sstModel, time, Month, Year) %>%
  dplyr::mutate(tActive=ifelse(species %in% c("Black-legged kittiwake", "Northern fulmar"), 0, tActive)) %>%
  dplyr::mutate(EnergyDiveTot=0) %>%
  replace_na(list(maxFlightBoutsMins=0)) %>%
  dplyr::summarise(tForage_month=sum(tForage), tLand_month=sum(tLand), tFlight_month=sum(tFlight), tRestWater_month=sum(tRestWater), tActive_month=sum(tActive), sstMonth=mean(sst_random, na.rm=TRUE), dayLength=mean(dayLengthHrs), MaxDistColKm=max(MaxDistColKm),
                   EnergyDiveTotMonth=sum(EnergyDiveTot), totDays=n_distinct(date), PropForage=tForage_month/(totDays*24), PropFlight=tFlight_month/(totDays*24),
                   PropActive=tActive_month/(totDays*24), PropLand=tLand_month/(totDays*24), PropRest=tRestWater_month/(totDays*24), maxFlight=max(maxFlightBoutsMins, na.rm=TRUE),
                   iceMonth=mean(ice_random, na.rm=TRUE), airMonth=mean(air_random, na.rm=TRUE), distColKmMonth=mean(distColonyKm_mean)) %>%
  ungroup() %>%
  dplyr::mutate(propTot=PropLand + PropActive + PropForage + PropRest + PropFlight)

# Save result
actResMonth$rep<-i
activityMonth<-rbind(activityMonth, actResMonth)
      
### Add model parameters to our dataset ###
print("Randomizing energy parameters...")

if (actRes$species[1] %in% c("Northern fulmar", "Black-legged kittiwake")) {

actRes$c2<-mean(sample(subset(modelParamsSub, parameter=="c2")$values[[1]], 1), sample(c(seq(2, 5.7, 0.1)) ,1))
actRes$c5<-0 # No active behaviour for surface-foragers 

} else {

actRes$c2<-sample(subset(modelParamsSub, parameter=="c2")$values[[1]], 1) # choose at random from uniform distribution
actRes$c5<-mean(sample(subset(modelParamsSub, parameter=="c5")$values[[1]][[1]], 1), sample(subset(modelParamsSub, parameter=="c5")$values[[1]][[2]], 1)) # mean of two numbers

}

actRes$RMR<-sample(subset(modelParamsSub, parameter=="RMR")$values[[1]], 1) # choose at random from uniform distribution
actRes$c1<-sample(subset(modelParamsSub, parameter=="c1")$values[[1]], 1) # choose at random from uniform distribution
actRes$c3<-sample(subset(modelParamsSub, parameter=="c3")$values[[1]], 1) # choose at random from uniform distribution
actRes$c4<-sample(subset(modelParamsSub, parameter=="c4")$values[[1]], 1) # Choose at random from uniform distribution
actRes$TC_air<-subset(modelParamsSub, parameter=="TC_air")$values[[1]] # Choose one value
actRes$TC_water<-subset(modelParamsSub, parameter=="TC_water")$values[[1]] # Choose one value
actRes$Beta_active<-subset(modelParamsSub, parameter=="Beta_active")$values[[1]] # Choose one value
actRes$Beta_rest<-subset(modelParamsSub, parameter=="Beta_rest")$values[[1]] # Choose one value
actRes$LCT_air<-subset(modelParamsSub, parameter=="LCT_air")$values[[1]] # Choose one value
actRes$LCT_water<-subset(modelParamsSub, parameter=="LCT_water")$values[[1]] # Choose one value
	  
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
energyDaily$c5<-actRes$c5[1]
energyDaily$TC_water<-actRes$TC_water[1]
energyDaily$TC_air<-actRes$TC_air[1]
energyDaily$Beta_active<-actRes$Beta_active[1]
energyDaily$Beta_rest<-actRes$Beta_rest[1]

energyDay<-rbind(energyDay, energyDaily)

#### Calculate energetics monthly  - based on locations ####
print("Calculating energetics monthly")

datesJoin_month<-seq.dates%>%
dplyr::mutate(Month=as.numeric(substr(date, 6, 7))) %>%
dplyr::mutate(year=as.numeric(substr(date, 1, 4))) %>%
dplyr::group_by(year, Month) %>%
dplyr::summarise(days=n_distinct(date))

energyInd1<-energyDaily %>%
dplyr::mutate(year=as.numeric(substr(date, 1, 4))) %>%
dplyr::group_by(species, colony, individ_id, year, Month, weight) %>%
dplyr::summarise(days=n_distinct(date), monthly_DEEkJ=sum(DEEkJ), monthly_sst=mean(sst_random)) %>%
dplyr::left_join(datesJoin_month, by=c("year", "Month")) %>%
dplyr::filter(days.y==days.x) # Here I remove incomplete months of data...

energyIndAll<-energyInd1
energyIndAll$rep<-i
energyMonth<-rbind(energyMonth, energyIndAll)

}

# Save some plots as PDFs (temporary to see what it looks ike #

# Calculate average activity
activityPlot_prep<-activityMonth %>%
ungroup() %>%
 dplyr::group_by(Month) %>%
 dplyr::summarise(tForageMean=mean(PropForage), sdForage=sd(PropForage),
                     tLandMean=mean(PropLand), sdLand=sd(PropLand),
                     tRestMean=mean(PropRest), sdRest=sd(PropRest),
                     tFlightMean=mean(PropFlight), sdFlight=sd(PropFlight),
                     tActiveMean=mean(PropActive), sdActive=sd(PropActive), repNo=n_distinct(rep)) %>%
	dplyr::filter(Month %in% c(7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)) %>%
	dplyr::filter(repNo==overall.iterations)
	
order<-c(7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)
order_filtered <- order[order %in% activityPlot_prep$Month] # make sure it only contains values present in our dataset otherwise ggplot looks crazy
	
activityPlot<-activityPlot_prep %>%
	dplyr::mutate(Month=factor(Month, levels=order_filtered)) %>%
    ggplot(aes(x=Month, y=tForageMean)) +
    geom_point(aes(colour="Forage")) +
    geom_linerange(aes(ymin=tForageMean-sdForage, ymax=tForageMean + sdForage, colour="Forage")) +
    geom_line(aes(colour="Forage", y=tForageMean, x=as.numeric(Month))) +
    geom_ribbon(aes(ymin=tForageMean-sdForage, ymax=tForageMean + sdForage, fill="Forage", x=as.numeric(Month)), alpha=0.2) +
    geom_point(aes(x=Month, y=tLandMean, colour="Land")) +
    geom_linerange(aes(x=Month, y=tLandMean, ymin=tLandMean - sdLand, ymax=tLandMean + sdLand,colour="Land")) +
    geom_line(aes(y=tLandMean, colour="Land", x=as.numeric(Month))) +
    geom_ribbon(aes(x=as.numeric(Month), y=tLandMean, ymin=tLandMean - sdLand, ymax=tLandMean + sdLand,fill="Land"), alpha=0.2) +
    geom_point(aes(x=Month, y=tRestMean, colour="Rest water")) +
    geom_linerange(aes(x=Month, y=tRestMean, ymin=tRestMean - sdRest, ymax=tRestMean + sdRest,colour="Rest water")) +
    geom_line(aes(as.numeric(Month), y=tRestMean, colour="Rest water")) +
    geom_ribbon(aes(as.numeric(Month), y=tRestMean, ymin=tRestMean - sdRest, ymax=tRestMean + sdRest,fill="Rest water"), alpha=0.2) +
    geom_point(aes(x=Month, y=tFlightMean, colour="Flight")) +
    geom_linerange(aes(x=Month, y=tFlightMean, ymin=tFlightMean - sdFlight, ymax=tFlightMean + sdFlight,colour="Flight")) +
    geom_line(aes(as.numeric(Month), y=tFlightMean, colour="Flight")) +
    geom_ribbon(aes(as.numeric(Month), y=tFlightMean, ymin=tFlightMean - sdFlight, ymax=tFlightMean + sdFlight,fill="Flight"), alpha=0.1) +
    geom_point(aes(x=Month, y=tActiveMean, colour="Active")) +
    geom_linerange(aes(x=Month, y=tActiveMean, ymin=tActiveMean - sdActive, ymax=tActiveMean + sdActive,colour="Active")) +
    geom_line(aes(as.numeric(Month), y=tActiveMean, colour="Active")) +
    geom_ribbon(aes(as.numeric(Month), y=tActiveMean, ymin=tActiveMean - sdActive, ymax=tActiveMean + sdActive,fill="Active"), alpha=0.1) +
    theme_bw() +
    ylab("Time in activity (hours.month)")  +
    ggtitle("Monthly estimates") +
    labs(colour="Behaviour", fill="Behaviour") +
	ylim(0, 1)

# Same plot for energetics

energyPlot_prep<-energyMonth %>%
ungroup() %>%
 dplyr::group_by(Month) %>%
 dplyr::summarise(meanDEE=mean(monthly_DEEkJ), days.y=mean(days.y), sdDEE=sd(monthly_DEEkJ), meanSST=mean(monthly_sst), sdSST=sd(monthly_sst), repNo=n_distinct(rep)) %>%
	dplyr::filter(Month %in% c(7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)) %>%
	dplyr::filter(repNo==overall.iterations)
	
order<-c(7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)
order_filtered <- order[order %in% energyPlot_prep$Month] # make sure it only contains values present in our dataset otherwise ggplot looks crazy

coef<-(max(energyMonth$monthly_DEEkJ)/30)/max(energyMonth$monthly_sst)
	
energyPlot<-energyPlot_prep %>%
	dplyr::mutate(Month=factor(Month, levels=order_filtered)) %>%
    ggplot(aes(x=Month, y=meanDEE/days.y)) +
    geom_point(aes(colour="Energy expenditure")) +
    geom_linerange(aes(ymin=(meanDEE-sdDEE)/days.y, ymax=(meanDEE + sdDEE)/days.y, colour="Energy expenditure")) +
    geom_line(aes(colour="Energy expenditure", y=meanDEE/days.y, x=as.numeric(Month))) +
    geom_ribbon(aes(ymin=(meanDEE-sdDEE)/days.y, ymax=(meanDEE + sdDEE)/days.y, x=as.numeric(Month), fill="Energy expenditure"), alpha=0.2) +
    theme_bw() +
    ylab("Energy expenditure (kJ.month)")  +
	geom_point(aes(x=Month, y=meanSST*coef, colour="SST")) +
    geom_linerange(aes(ymin=(meanSST-sdSST)*coef, ymax=(meanSST + sdSST)*coef, colour="SST")) +
    geom_line(aes(colour="SST", y=meanSST*coef, x=as.numeric(Month))) +
    geom_ribbon(aes(ymin=(meanSST-sdSST)*coef, ymax=(meanSST + sdSST)*coef, fill="SST", x=as.numeric(Month)), alpha=0.2) +
    
    # Custom the Y scales:
    scale_y_continuous(
      
      # Features of the first axis
      name = "Energy expenditure (kJ.month)",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis( trans=~./coef, name="SST")) +
    labs(colour="Metric", fill="Metric") +
    ggtitle("Energy vs. SST") 
	

pdf(paste0("./results/tempPlots/", activityMonth$species[1], "_", activityMonth$colony[1], "_", activityMonth$individ_id[1], ".pdf"))
grid.arrange(activityPlot, energyPlot, nrow=2)
dev.off()

# Save result per bird

#### Save outputs ####

print("Saving file...")
output_file1 <- args[2]

print("1")
print(paste0("Sessions: ", n_distinct(energyDay$session_id)))
write.csv(energyDay, file = output_file1, row.names = FALSE) # Daily energy data 
print("Finished")
