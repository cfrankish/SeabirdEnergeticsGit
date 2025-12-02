### This script extracts immersion data for all species of interest from the SEATRACK database. It is parallelized by species ###
### The data is matched with IRMA data (by ID ) - only tracks with both immersion & positional data are retained ###
### Immersion data with unusual sampling regime are removed ###
### Juveniles are also removed ###
### Input will be a list of species to iterate through ###
### Main output file will be individual bird files (saved in folder named './data/wetdry_raw/' - as 'species_wetdry_date.csv') #### 
### There will also be a table showing number of birds per colony & species for every processing step ('species_sampleSizes.rds') ###
### There is also one table showing possible mismatches between IDs that need to be checked ('species_mismatchIds.rds') and one showing the sample size of logger types ('species_wetdry_type.rds') ###

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

### Step 1: Source wet-dry data ###

print("Step 1: extracting data from database...")

# Connect to database
print("Connecting...")
connectSeatrack(Username = "testreader", Password = "testreader",
                host = "seatrack.nina.no")
print("Connected")

# Download wet-dry data from species of interest				
wetdry_recording<-getRecordings(type="activity", species=speciesName)

# getLoggerInfo(asTibble = T)

# Find the metadata for this species (colony)
info <-getIndividInfo()

# Subset relevant columns
infoMeta<-info %>%
  dplyr::select(species, session_id, individ_id, colony, data_responsible) %>%
  group_by(session_id, individ_id, colony) %>%
  dplyr::slice(1)

# Join the metadata so we know which colony this corresponds to
wetdry_recording_meta<-wetdry_recording %>%
  dplyr::mutate(individ_id=ifelse(individ_id=="ISR_4110394", "ISR_4122906", individ_id)) %>% # There is a mismatch between IDs between posDat & wetdry entery for this bird (AtPu)
  dplyr::mutate(individ_id=ifelse(individ_id=="NOO_CJ_13025", "NOO_KA00991", individ_id)) %>% # There is a mismatch between IDs between posDat & wetdry entery for this bird (BlKi)
  dplyr::mutate(individ_id=ifelse(individ_id=="RUM_PS21120", "RUM_PS21220", individ_id)) %>% # There is a mismatch between IDs between posDat & wetdry entery for this bird (BlKi)
  dplyr::left_join(infoMeta, by=c("session_id", "individ_id")) 
  
# Make sure no NA Ids #
naIds<-subset(wetdry_recording_meta, is.na(species))

if (nrow(naIds) > 0) {

print("Na Ids")

break
}
  
# Record sample size (species, colony, individual ID) #
sampleSize1<-wetdry_recording_meta %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(birds_r1_rawact=n_distinct(individ_id)) 
				 
### Step 2: Match it with IRMA data (by session_id) ###

print("Step 2: Matching with IRMA data...")

# Load IRMA data (all)
IRMAdata<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")

# Load IRMA data by species
IRMAfiles<-list.files("./data/positionsIRMA/", full.names=TRUE) # Compile list of all IRMA data
IRMAfiles<-IRMAfiles[1:6] # Subset to species-specific files
IRMAfiles.df<-data.frame(files=IRMAfiles) # Turn into a data frame to make it easy to Subset
IRMAfiles.df$species<-c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot") # Assign species names

# Subset species-specific file
IRMAsub<-subset(IRMAfiles.df, species==speciesName)
IRMAsub<-readRDS(IRMAsub$files) # Open rds file

# Create a list of common and latin species names to help matching between datasets
speciesNames<-c("Northern fulmar", "Atlantic puffin", "Common guillemot", "Little auk", "Black-legged kittiwake", "Brünnich's guillemot")
latinNames<-c("Fulmarus_glacialis", "Fratercula_arctica", "Uria_aalge", "Alle_alle", "Rissa_tridactyla", "Uria_lomvia")
matchingNames<-data.frame(species=speciesNames, latin=latinNames) # make a dataframe out of this
matchingNamesSub<-subset(matchingNames, species==speciesName) # Subset to species of interest

# Create a list with the IDs in the IRMA dataset & make a key for matching between ids in seatrack database and in IRMA dataset
idsIRMA<-IRMAdata %>%
  ungroup() %>%
  dplyr::filter(species %in% matchingNamesSub$latin) %>%
  dplyr::select(colony, individ_id) %>%
  dplyr::rename(colony_IRMA=colony) %>%
  dplyr::rename(individ_id_IRMA=individ_id) %>%
  dplyr::mutate(individ_id=sub("-", "_", individ_id_IRMA)) %>%
  dplyr::mutate(individ_id=gsub("KAPH", "KapH", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("MK14", "mk14", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("MK13", "mk13", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("HORNØYA", "Hornøya", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("CAPEKRUTIK", "CapeKrutik", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("BTO_SKELLIG-1", "BTO_Skellig_1", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("BTO_SKELLIG-2", "BTO_Skellig_2", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("BTO_SKELLIG-9", "BTO_Skellig_9", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("GBT_PUFFIN-8793-IOM", "GBT_puffin_8793_IOM", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("CANADA_99691670", "CANADA_98691670", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("NOS_1116-LAT250A-2010", "NOS_1116_LAT250A_2010", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("NOS_1135-LAT250A-2010", "NOS_1135_LAT250A_2010", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("NOS_1136-LAT250A-2010", "NOS_1136_LAT250A_2010", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("NOS_1170-LAT250A-2010", "NOS_1170_LAT250A_2010", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("NOS_1174-LAT250A-2010", "NOS_1174_LAT250A_2010", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("NOS_1201-LAT250A-2010", "NOS_1201_LAT250A_2010", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-BK-NUUK-", "_BK_NUUK_", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-LA-KapH", "_LA_KapH", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-2010-mk13", "_2010_mk13", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("BK-CapeKrutik-", "BK_CapeKrutik_", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-BLK-Hornøya-", "_BLK_Hornøya_", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-Blue-RUM", "_Blue_RUM", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-BLK-", "_BLK_", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-2010-mk14", "_2010_mk14", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-2008-mk14", "_2008_mk14", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("BLUE", "Blue", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("DK-HORNOYA", "DK_HORNOYA", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-HORNOYA", "_HORNOYA", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-2009", "_2009", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-2021", "_2021", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-Blue-RUM", "_blue_RUM", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-LAT250A-", "_LAT250A_", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("-Blue", "_blue", individ_id)) %>%
  dplyr::mutate(individ_id=gsub("GBT_EA35285", "GBT_C2513_2018_IoM_BK", individ_id)) %>%
  dplyr::mutate(individ_id=ifelse(individ_id=="ISR_4110455", "ISR_4122905", individ_id)) %>%
  dplyr::filter(!individ_id %in% c("ISR_B09248")) %>% # Unable to re-construct track according to Vegard
  dplyr::filter(!individ_id %in% c("DKC_7159186")) %>% # Unable to re-construct track according to Vegard
  dplyr::filter(!individ_id %in% c("RUM_ES024446")) # Unable to re-construct track according to Vegard

# Inner join with wet-dry data
wetdry_subset<-wetdry_recording_meta %>%
  dplyr::inner_join(idsIRMA, by=c("individ_id")) %>%
  dplyr::mutate(colony=ifelse(is.na(colony), colony_IRMA, colony))
  
# Do some checks to see if anything is odd

# See which ids don't match (they have positional data but not immersion data)
idsDatabase<-wetdry_recording_meta %>%
  dplyr::select(individ_id, session_id) %>%
  dplyr::group_by(individ_id ,session_id) %>%
  dplyr::slice(1)

idsNoMatch1<-idsIRMA %>%
  dplyr::anti_join(idsDatabase, by=c("individ_id")) %>% # There are 2 (now 4) for fulmars, 36 for puffins, 71 for CoGu, 43 for LiAu, 83 for Blk, 165 for BrGu
  dplyr::mutate(type="Lox_noAct") %>%
  dplyr::select(individ_id, type)
 
  
# Look for non-matching ids in positional data
connectSeatrack(Username = "testreader", Password = "testreader",
                host = "seatrack.nina.no")
posDat<-getPosdata(individId = idsNoMatch1$individ_id)

# Do we find the same number of files?
posDataIds<-unique(posDat$individ_id)

if (!n_distinct(posDat$individ_id) == nrow(idsNoMatch1)) {
  
  print("Issue with non-matching IDs")
  break
}

# Look for logger status so I can share with Vegard
connectSeatrack(Username = "testreader", Password = "testreader",
                host = "seatrack.nina.no") 
loggerInfo<-getLoggerInfo(asTibble = T)
sessions<-posDat$session_id
loggerSub<-subset(loggerInfo, session_id %in% sessions)
loggerSubInfo<-loggerSub %>%
dplyr::select(session_id, download_type)
posDatInfo<-posDat %>%
dplyr::group_by(individ_id, session_id) %>%
dplyr::slice(1) %>%
dplyr::select(individ_id, session_id) %>%
dplyr::inner_join(loggerSubInfo, by=c("session_id"))
idsNoMatch1_final<-idsNoMatch1 %>%
dplyr::inner_join(posDatInfo, by=c("individ_id"))

# Finally, which wet-dry data have no matching positional data?
idsNoMatch2<-idsDatabase %>%
  dplyr::anti_join(idsIRMA, by=c("individ_id")) %>% # There are 2 (now 4) for fulmars, 36 for puffins, 71 for CoGu, 43 for LiAu, 83 for Blk, 165 for BrGu
  dplyr::mutate(type="Act_noLox") 
loggerSub2<-subset(loggerInfo, session_id %in% idsNoMatch2$session_id)
loggerSubInfo2<-loggerSub2 %>%
dplyr::select(session_id, download_type)
idsNoMatch2_final<-idsNoMatch2 %>%
dplyr::inner_join(loggerSubInfo2, by=c("session_id"))
   
# Save these outputs to analyze later
idsNoMatch<-rbind(idsNoMatch1_final, idsNoMatch2_final)
idsNoMatch$species<-speciesName
write.table(idsNoMatch, file=paste0("./data/wetdry_raw/", speciesNameSave, "_mismatchIds.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Record sample size (species, colony, individual ID) #
sampleSize_intermediate<-wetdry_subset %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(birds_r2_matchIRMAids=n_distinct(individ_id)) 

# Attach to main sample size dataset to keep tracking of changing number of individuals
sampleSize2<-sampleSize1 %>%
ungroup() %>%
dplyr::full_join(sampleSize_intermediate, by=c("species", "colony")) 

### Step 3: Assign sampling interval & remove those that don't fit ###

print("Step 3: Remove tracks with rare sampling interval...")

# Summarise max conductivity reading per session id & mean sampling interval
act.characteristics<-wetdry_subset %>%
    dplyr::group_by(individ_id, session_id) %>%
    arrange(individ_id, date_time) %>%
    dplyr::mutate(time.interval=as.numeric(difftime(date_time, lag(date_time), unit=c("mins")))) %>%
    dplyr::mutate(max.cond=max(conductivity)) %>%
    ungroup() %>%
    dplyr::group_by(individ_id, session_id, max.cond) %>%
    dplyr::count(time.interval) %>%
    dplyr::filter(!is.na(time.interval)) %>%
    dplyr::filter(n==max(n))
	
# Join with main dataset
 act.logger.type<-wetdry_subset %>%
    dplyr::left_join(act.characteristics, by=c("individ_id", "session_id"))
 
# Save complete sample size by logger type
wetdryType<-act.logger.type %>%
ungroup() %>%
dplyr::group_by(species, max.cond) %>%
dplyr::summarise(birds=n_distinct(individ_id))

saveRDS(wetdryType, file=paste0("./data/wetdry_raw/", speciesNameSave, "_wetdry_type.rds"))
 
# Remove logger types which are rare & not really comparable (max conduct = 480 & 50) # following Amelineau et al. 2021
# https://doi.org/10.3354/meps13872 
act.remove<-act.logger.type %>%
 dplyr::filter(max.cond %in% c(20, 200)) # These are the ones I am keeping
  
# Standardize to between 0 & 1 by dividing by max conductivity reading

print("Step 3: Standardizing immersion data to 0-1...")

act.standard<-act.remove %>%
 dplyr::mutate(new.cond=conductivity/max.cond)
  
# Check this worked
min.cond<-min(act.standard$new.cond)
max.cond<-max(act.standard$new.cond)
  
if(!min.cond %in% c(0) | !max.cond %in% c(1)) {
    
    print ("Range is outside of 0-1")
    break
  }

# Record sample size (species, colony, individual ID) #
sampleSize_intermediate2<-act.standard %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(birds_r3_removeRareLoggers=n_distinct(individ_id)) 

# Attach to main sample size dataset to keep tracking of changing number of individuals
sampleSize3<-sampleSize2 %>%
ungroup() %>%
dplyr::full_join(sampleSize_intermediate2, by=c("species", "colony")) 

### Step 4: Remove juveniles ###

print("Step 4: Removing juveniles...")

#find age at deployment
Individs<-info[info$eventType%in%"Deployment",]

#limit to individuals deployed on as chicks/pullus
sessions_juveniles<-Individs[(Individs$status_age%in%c("pullus","chick")),]

# Remove these sessions from the main dataset
act.adults<-subset(act.standard, !session_id %in% sessions_juveniles$session_id)

# Record sample size (species, colony, individual ID) #
sampleSize_intermediate3<-act.adults %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(birds_r4_removeJuveniles=n_distinct(individ_id)) 

# Attach to main sample size dataset to keep tracking of changing number of individuals
sampleSize4<-sampleSize3 %>%
ungroup() %>%
dplyr::full_join(sampleSize_intermediate3, by=c("species", "colony")) %>%
dplyr::mutate(date=Sys.Date()) # add date to know when this happened

# Save sample size Output
saveRDS(sampleSize4, file=paste0("./data/wetdry_raw/", speciesNameSave, "_sampleSizes.rds"))

### Step 5: Save outputs ###

print("Step 5: Saving outputs...") 

# Intermediate 'species' files
act.final<-act.adults %>%
dplyr::ungroup() %>%
dplyr::select(species, colony, colony_IRMA, individ_id, individ_id_IRMA, session_id, date_time, max.cond, new.cond, data_responsible) %>% # select most important columns to clean up
dplyr::mutate(extractionDate=Sys.Date())

# Number 1
output_file1 <- args[3]
write.csv(act.final, file = output_file1, row.names = FALSE) # Wet-dry data by species

#write.csv(act.final, file=paste0("./data/wetdry_raw/", speciesNameSave, "_wetdry_2025-08-08.csv")) # Save output file

print("Done")
