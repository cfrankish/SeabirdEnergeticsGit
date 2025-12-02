# This script starts by creating an ID catalogue of all birds for which we have energy files to make the following steps easier #
# Otherwise, what it is really doing is carrying out ana analysis to determine how many iterations we need before estimates of energy expenditure stabilize #
# This number is carried through the rest of the analysis #
# Input files are start & end dates of the study period #
# Output file is a table with total non-breeding energy expenditure vs. number of iterarions as well as a figure showing this (Figure S11): "./results/tables/supplementary/table1_miniteration.csv" #

# load functions
library(ggplot2)
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
#library(adehabitatHR)

# Readin arguments
args <- commandArgs(trailingOnly = TRUE)
startDate<-args[1] # Read-in start of study period
endDate<-args[2] # Read-in end date of study period

### Step 1: Load all individual Ids from result files ###

# This is so we can loop through them later #

allResults<-list.files("./tmp", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# energyRes_day<-sample(energyRes_day, 100, replace=FALSE)

energyAll<-list() 

print("Step 1: making catalogue of energy files...")

energyAll<-list() 

print("Step 1: making catalogue of energy files...")

for (i in 1:length(energyRes_day)) {
  
  print(i)
  
  # Subset individual file
  energySub<-fread(energyRes_day[i], nrow=2)
  
  energySub2<-energySub %>%
  ungroup() %>%
  dplyr::group_by(species, colony, individ_id) %>%
  dplyr::slice(1)
  
  # save results
  energyAll<-rbind(energyAll, energySub2)
 
  
}

# Save output files #
output_file1 <- args[3]
print("Saving output file 1")
write.csv(energyAll, file = output_file1, row.names = FALSE) # directory of all ids that will be re-used after

### Step 2: Calculate minimum iteration number ###

print("Step 2: determine minimum sample size per pop")

# First we make a list of populations

colony.summary<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")

speciesMatch<-data.frame(speciesLatin=c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle"), 
                         species=c("Brünnich's guillemot", "Black-legged kittiwake", "Common guillemot", "Atlantic puffin", "Northern fulmar", "Litte auk"))

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
                          "RU", "NO", "GR", "RU", "RU", "NO", "NO", "NO", "CA", "IC", "IC", "IC", "IC", "IC", "IC", 
                          "GR", "NO", "IC", "IC", "GR", "IC", "IC", "NO", "IC", "CA", "CA", "NO", "FA", "GR", "UK", "NO", "UK", "UK", "DK",
                          "UK", "UK", "UK", "IR", "IR", "CA", "IR", "IR", "IR", "UK", "CA", "CA")) %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(colonyNo=row_number()) %>%
  dplyr::mutate(colonyName=paste0(country, colonyNo)) # Here i make a special naming system with country first and colony second. It is ordered from North to South
  
 colonies_lox<-colonyMatch %>%
 dplyr::left_join(colonyNames, by=c("colony")) %>%
 ungroup() %>%
 arrange(colonyName) %>%
 dplyr::group_by(colonyName) %>%
 dplyr::slice(1)

# First we subset the dataset to large populations (this might become smaller based on whether birds have the entire study period or not) 

largePops<-energyAll %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::mutate(birds=n_distinct(individ_id)) %>%
dplyr::filter(birds>=5) %>%
dplyr::left_join(colonyNames, by=c("colony"))

# Determine location of input files
allResults<-list.files("./tmp/", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

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
  
# We will use this to turn DEE into kj.g
species<-data.frame(species=c("Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot", "Little auk", "Atlantic puffin"))
species$allometryCoef<-c(0.765, 0.717, 0.689, 0.689, 0.689, 0.689)

# Determine number of species to loop through
speciesList<-unique(largePops$species)

# Make a list to save the results in
minInterationNoRes<-list()

for (i in 1:length(speciesList)) {

# Subset to species i
speciesSub<-speciesList[i]
largePops_speciesSub<-subset(largePops, species==speciesSub)

# Determine list of colonies to loop through (one population per country group at random)

randomSample<-largePops_speciesSub %>%
ungroup() %>%
dplyr::group_by(country) %>%
dplyr::count(colony) %>%
dplyr::slice_sample(n=1)

Pops<-unique(randomSample$colony)

# make a list to save results in
minIterationNo_colRes<-list()

for (j in 1:length(Pops)) {

# Subset to population j
colonySub<-Pops[j]
largePops_colSub<-subset(largePops_speciesSub, colony==colonySub)

# make a list of ids to cycle through
ids<-unique(largePops_colSub$individ_id)
ids<-sample(ids, length(ids), replace=FALSE) # re-sample at random so birds are really chosen at random

# make a list to save results in
minSampleSize_idRes<-list()

print(paste0("Assembling birds, Species ", i, " Pop", j))

for (k in 1:length(ids)) {

print(k)

# Subset to bird k
birdSub<-ids[k]

# Change hyphen to underscore if it's there
birdSub<-gsub("-", "_", birdSub)
birdSub<-gsub("ø", "o", birdSub)

# Open id k
birdSub_csv<-read.csv(energyRes_day[grepl(birdSub, energyRes_day)])

# Split data into session/track year combinations
birdSub_csv$month<-as.numeric(substr(birdSub_csv$date, 6, 7))
birdSub_csv$year<-as.numeric(substr(birdSub_csv$date, 1, 4))
birdSub_csv$track_year<-ifelse(birdSub_csv$month<day1_month, paste0(birdSub_csv$year-1, "_", substr(birdSub_csv$year, 3, 4)), paste0(birdSub_csv$year, "_", substr(birdSub_csv$year + 1, 3, 4)))
birdSub_csv$session_year<-paste0(birdSub_csv$session_id, "_", birdSub_csv$track_year)
      
# Figure out whether data covers study period
dates_sessions<-birdSub_csv %>%
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

if (nrow(dates_sessions)<1) {
next}

# Otherwise we choose a random session_id_year every rep 
maxiteration<-max(birdSub_csv$rep) # Determine how many times
randomSelect<-sample(dates_sessions$session_year, maxiteration, replace=TRUE) # Select at random
selectDf<-data.frame(rep=1:maxiteration, session_year=randomSelect)

# Now we sub-select our full dataset % match with correct dates
birdSub_csv$day<-as.numeric(substr(birdSub_csv$date, 9, 10))
birdSub_csv$month<-as.numeric(substr(birdSub_csv$date, 6, 7))
birdSub_rep1<-birdSub_csv %>%
dplyr::inner_join(selectDf, by=c("rep", "session_year")) %>%
ungroup() %>%
dplyr::inner_join(dates_weekly2, by=c("month", "day")) %>%
dplyr::group_by(rep) %>%
dplyr::mutate(totalDays=n_distinct(date)) # make sure correct number of days has been selected

# Make sure the number of days is correct
print("Subsetting to study period")
studyPeriod<-unique(birdSub_rep1$totalDays)
if (length(studyPeriod)>1 || ! studyPeriod %in% c(nrow(dates_weekly2))) {
 (print("Error: number of days off"))
 next
}
      
# Save results
minSampleSize_idRes<-rbind(minSampleSize_idRes, birdSub_rep1)

}

# Determine list of ids
idsPop<-unique(minSampleSize_idRes$individ_id)

if (length(idsPop) <1) {
next

}

# Choose bird at random
idChoose<-sample(idsPop, 1)

# Here we look at mean DEE for an increasing number of iterarions
print("Calculating number of iterations...")

maxReps<-max(minSampleSize_idRes$rep)

# make a list to save the results in
increasingIterationNumber<-list()

for (m in 1:max(minSampleSize_idRes$rep)){

print(paste0("Calculating energy for species ", i, " Pop ", j, "rep", m, "/", max(minSampleSize_idRes$rep)))

# Subset iterations 1:m
repsSub<-1:m

# Subset dataset to these birds
birdsAnalysisSub<-subset(minSampleSize_idRes, rep %in% repsSub & individ_id %in% c(idChoose))

# Calculate mean energy expenditure per bird
energyCost<-birdsAnalysisSub %>%
dplyr::left_join(species, by=c("species")) %>%
dplyr::mutate(DEEg=DEEkJ/weight^allometryCoef) %>%
dplyr::ungroup() %>%
dplyr::group_by(rep, species, colony, individ_id) %>%
dplyr::summarise(totenergy=sum(DEEg)) %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totenergy_ind_iterations=mean(totenergy), sdenergy=sd(totenergy)) %>%
dplyr::mutate(reps=m)

# save results
increasingIterationNumber<-rbind(increasingIterationNumber, energyCost)

}

# Save colony results
minIterationNo_colRes<-rbind(minIterationNo_colRes, increasingIterationNumber)
#minSampleSize_colRes<-rbind(minSampleSize_colRes, increasingSampleSize)

}

# Save species results
minInterationNoRes<-rbind(minInterationNoRes, minIterationNo_colRes)
#minSampleSizeRes<-rbind(minSampleSizeRes, minSampleSize_colRes)

}

# Make a plot
FigureS11<-minInterationNoRes %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
												 ggplot(aes(x=reps, y=totenergy_ind_iterations)) +
												 geom_line(aes(group=individ_id, color=individ_id)) +
												 geom_ribbon(aes(ymin=totenergy_ind_iterations - sdenergy, ymax=totenergy_ind_iterations + sdenergy, group=individ_id, fill=individ_id), alpha=0.2) +
												 facet_wrap(~species) +
												 ylab("NB energy expenditure (mean +/- SD)") +
												 xlab("Number of iterations") +
												 guides(color=FALSE, fill=FALSE) +
												 theme_bw()
												 
pdf("./results/figures/supplementary/FigureS11.pdf")
plot(FigureS11)
dev.off()

# Save output file 2
output_file2 <- args[4]
print("Saving output file 2")
write.csv(minInterationNoRes, file = output_file2, row.names = FALSE) # Table for calculation min number of iterations
 

