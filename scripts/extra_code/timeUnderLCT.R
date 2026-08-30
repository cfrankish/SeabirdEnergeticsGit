library(dplyr)
library(data.table)
library(sp)
library(raster)
library(sf)
library(terra)
library(ncdf4)
library(ggplot2)

# Source all necessary functions
source("functions.R")

#### How much time below LCT? or above? #### + degrees below

lcts<-data.frame(species=c("Black-legged kittiwake", "Northern fulmar", "Little auk", "Common guillemot", "Brünnich's guillemot", "Atlantic puffin"),
LCT_water=c(12.5, 9, 14.18, 14.18, 14.18, 14.18), LCT_air=c(4.5, 9, 4.5, 2, 2, 5.72))

# Open each individual file cry
idCatalogue<-read.csv("../results/tables/main/table1_idcatalogue.csv")

# Determine where daily files are
allResults<-list.files("../tmp/", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# Determine study period
startDate<-"2021-09-15" # Read-in start of study period
endDate<-"2022-04-15" # Read-in end date of study period

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

# Make a species list to iterate through #
speciesUnique<-unique(idCatalogue$species)

# Make a list to save results
species_lct_stats<-list()

for (i in 2:length(speciesUnique)) {

#for (i in 1:1) {

# subset to species i
speciesSub<-speciesUnique[i]

# make a list of ids to iterate through
idSub<-subset(idCatalogue, species==speciesSub)
idList<-unique(idSub$individ_id)

# list to save results
bird_lct_stats<-list()

for (j in 1:length(idList)) {

# Subset bird j
    birdSub<-idList[j]
	birdSub<-gsub("-", "_", birdSub)
	birdSub<-gsub("ø", "o", birdSub)
	birdSub<-gsub("6T_", "6T-", birdSub)
	birdSub<-gsub("4F_", "4F-", birdSub)
    
    # Find csv file
    print(paste0("Opening file..Species ", i, " Bird ", j, "/", length(idList)))
    birdSub<-fread(energyRes_day[grep(birdSub, energyRes_day)])
	
	# Determine month and year
	birdSub$day<-as.numeric(substr(birdSub$date, 9, 10))
	birdSub$month<-as.numeric(substr(birdSub$date, 6, 7))
	birdSub$year<-as.numeric(substr(birdSub$date, 1, 4))
	birdSub$track_year<-ifelse(birdSub$month<day1_month, paste0(birdSub$year-1, "_", substr(birdSub$year, 3, 4)), paste0(birdSub$year, "_", substr(birdSub$year + 1, 3, 4)))
    birdSub$session_year<-paste0(birdSub$session_id, "_", birdSub$track_year)

    # Determine session_year combos with enough data	
    dates_sessions<-birdSub%>%
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
    
	if (nrow(dates_sessions) < 1) {
      next}
	
   # Subset to one location per date to reduce file size
    birdcsv_reducted<-birdSub %>%
      ungroup() %>%
	  dplyr::filter(session_year %in% c(dates_sessions$session_year)) %>%
      dplyr::inner_join(dates_weekly2, by=c("day", "month")) %>%
	  dplyr::filter(weekNo>0) %>%
	  dplyr::select(rep, species, colony, individ_id, track_year, date, session_year, mean.lon, mean.lat, sst_random, tFlight, tLand) 

   # Extract sst...
   birdcsv_airTemp<-extract_air_temp(birdcsv_reducted)
	
   # Attach LCT
   bird_lct<-birdcsv_airTemp %>%
   dplyr::left_join(lcts, by=c("species")) %>%
   dplyr::mutate(below_lct_water=ifelse(sst_random<LCT_water, 1, 0), below_lct_air=ifelse(airTemp<LCT_air, 1, 0)) %>%
   dplyr::group_by(rep, species, colony, individ_id, track_year) %>%
   dplyr::summarise(propBelow_lct_water=sum(below_lct_water)/n_distinct(date), propBelow_lct_air=sum(below_lct_air)/n_distinct(date)) %>%
   ungroup() %>%
   dplyr::group_by(species, colony, individ_id) %>%
   dplyr::summarise(meanpropBelow_lct_water=mean(propBelow_lct_water), meanpropBelow_lct_air=mean(propBelow_lct_air))
   
   # Compute degrees below when it is below (water)
   bird_degrees_below_water<-birdcsv_airTemp %>%
   dplyr::left_join(lcts, by=c("species")) %>%
   dplyr::mutate(below_lct_water=ifelse(sst_random<LCT_water, 1, 0)) %>%
   dplyr::filter(below_lct_water==1) %>%
   dplyr::mutate(degrees_below_lct_water=LCT_water - sst_random) %>%
   dplyr::group_by(rep, species, colony, individ_id, track_year) %>%
   dplyr::summarise(degrees_below_lct_water_mean=mean(degrees_below_lct_water)) %>%
   ungroup() %>%
   dplyr::group_by(species, colony, individ_id) %>%
   dplyr::summarise(degrees_below_lct_water=mean(degrees_below_lct_water_mean), degrees_below_lct_water_max=max(degrees_below_lct_water_mean), degrees_below_lct_water_min=min(degrees_below_lct_water_mean))
   
   bird_degrees_below_air<-birdcsv_airTemp %>%
   dplyr::left_join(lcts, by=c("species")) %>%
   dplyr::mutate(below_lct_air=ifelse(airTemp<LCT_air, 1, 0)) %>%
   dplyr::filter(below_lct_air==1) %>%
   dplyr::mutate(degrees_below_lct_air=LCT_air - airTemp) %>%
   dplyr::group_by(rep, species, colony, individ_id, track_year) %>%
   dplyr::summarise(degrees_below_lct_air_mean=mean(degrees_below_lct_air), totFlight=sum(tFlight), totLand=sum(tLand), totDry=sum(totFlight, totLand), totTime=n_distinct(date)*24, propDry=totDry/totTime) %>%
   ungroup() %>%
   dplyr::group_by(species, colony, individ_id) %>%
   dplyr::summarise(degrees_below_lct_air=mean(degrees_below_lct_air_mean), degrees_below_lct_air_max=max(degrees_below_lct_air_mean), degrees_below_lct_air_min=min(degrees_below_lct_air_mean), meanPropLand_below_air=mean(totLand/totTime))
   
   # Compute degrees below when it is below
   bird_degrees_above_water<-birdcsv_airTemp %>%
   dplyr::left_join(lcts, by=c("species")) %>%
   dplyr::mutate(below_lct_water=ifelse(sst_random<LCT_water, 1, 0)) %>%
   dplyr::filter(below_lct_water==0)

   bird_degrees_above_air<-birdcsv_airTemp %>%
   dplyr::left_join(lcts, by=c("species")) %>%
   dplyr::mutate(below_lct_air=ifelse(airTemp<LCT_air, 1, 0)) %>%
   dplyr::filter(below_lct_air==0)   
   
   if(nrow(bird_degrees_above_water)>0) {
   
   bird_degrees_above_final<-bird_degrees_above_water %>%
   dplyr::mutate(degrees_above_lct_water=sst_random - LCT_water) %>%
   ungroup() %>%
   dplyr::group_by(rep, species, colony, individ_id, track_year) %>%
   dplyr::summarise(degrees_above_lct_water_mean=mean(degrees_above_lct_water), totFlight=sum(tFlight), totLand=sum(tLand), totDry=sum(totFlight, totLand), totTime=n_distinct(date)*24, propDry=totDry/totTime) %>%
   ungroup() %>%
   dplyr::group_by(species, colony, individ_id) %>%
   dplyr::summarise(degrees_above_lct_water=mean(degrees_above_lct_water_mean), degrees_above_lct_water_max=max(degrees_above_lct_water_mean), degrees_above_lct_water_min=min(degrees_above_lct_water_mean), meanPropDry=mean(propDry))
   
   # Join these datasets together
   bird_lct_sum<-bird_lct %>%
   dplyr::left_join(bird_degrees_below_water, by=c("species", "colony", "individ_id")) %>%
   dplyr::left_join(bird_degrees_below_air, by=c("species", "colony", "individ_id")) %>%
   dplyr::left_join(bird_degrees_above_final,by=c("species", "colony", "individ_id"))
   
   } else {
   
   # Join these datasets together
   bird_lct_sum<-bird_lct %>%
   dplyr::left_join(bird_degrees_below_water, by=c("species", "colony", "individ_id")) %>%
   dplyr::left_join(bird_degrees_below_air, by=c("species", "colony", "individ_id"))

}

if (nrow(bird_degrees_above_air)>0) {

bird_degrees_above_final2<-bird_degrees_above_air %>%
   dplyr::mutate(degrees_above_lct_air=airTemp - LCT_air) %>%
   ungroup() %>%
   dplyr::group_by(rep, species, colony, individ_id, track_year) %>%
   dplyr::summarise(degrees_above_lct_air_mean=mean(degrees_above_lct_air), totFlight=sum(tFlight), totLand=sum(tLand), totDry=sum(totFlight, totLand), totTime=n_distinct(date)*24, propDry=totDry/totTime) %>%
   ungroup() %>%
   dplyr::group_by(species, colony, individ_id) %>%
   dplyr::summarise(degrees_above_lct_air=mean(degrees_above_lct_air_mean), degrees_above_lct_air_max=max(degrees_above_lct_air_mean), degrees_above_lct_air_min=min(degrees_above_lct_air_mean), meanPropLand_above_air=mean(totLand/totTime))

bird_lct_sum<-bird_lct_sum %>%
dplyr::left_join(bird_degrees_above_final2, by=c("species", "colony", "individ_id"))

}

 bird_lct_stats<-rbind(bird_lct_stats, bird_lct_sum)
 
 }
 
 print("Saving intermediate file...")
 write.csv(bird_lct_stats, file=paste0("../results/LCT_stats_", speciesSub, ".csv"))
 
 species_lct_stats<-rbind(species_lct_stats, bird_lct_stats)
 
 }
 
## Now we analyze the results ##

results_check_lct<-list.files("../results/", full.names=TRUE) # find folder with results
lct_files<-results_check_lct[grep("LCT", results_check_lct)] # Subset to LCT files

# make some lists to save the results as we loop through the species #
species_lct_water_all<-list()
species_lct_air_all<-list()
colony_lct_water_all<-list()
colony_lct_air_all<-list()

for (i in 1:length(lct_files)) {

print(i)

# open results file i 
lct_files_i<-read.csv(lct_files[i])

# LCT water first #

# Prop time below LCT_water - species
species_lct_water1<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(meanPropBelowLCT=mean(meanpropBelow_lct_water, na.rm=TRUE), birds=n_distinct(individ_id)) %>%
dplyr::filter(birds>=5) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::summarise(metric="prop_time_below_LCT_water",  metric_mean=mean(meanPropBelowLCT), metric_sd=sd(meanPropBelowLCT), metric_se=metric_sd/sqrt(n_distinct(colony)), type="species")

# Prop time below LCT_water - colony
colony_lct_water1<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(metric="prop_time_below_LCT_water", metric_mean=mean(meanpropBelow_lct_water), metric_sd=sd(meanpropBelow_lct_water), metric_se=metric_sd/sqrt(n_distinct(individ_id)), birds=n_distinct(individ_id), type="colony") %>%
dplyr::filter(birds>=5) 

# degrees below LCT_water - species
species_lct_water2<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(meanPropBelowLCT=mean(degrees_below_lct_water, na.rm=TRUE), birds=n_distinct(individ_id)) %>%
dplyr::filter(birds>=5) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::summarise(metric="degrees_below_LCT_water",  metric_mean=mean(meanPropBelowLCT), metric_sd=sd(meanPropBelowLCT), metric_se=metric_sd/sqrt(n_distinct(colony)), type="species")

# degrees below LCT_water - colony
colony_lct_water2<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(metric="degrees_below_LCT_water", metric_mean=mean(degrees_below_lct_water), metric_sd=sd(degrees_below_lct_water), metric_se=metric_sd/sqrt(n_distinct(individ_id)), birds=n_distinct(individ_id), type="colony") %>%
dplyr::filter(birds>=5) 

# Bind together #
species_lct_water<-rbind(species_lct_water1, species_lct_water2) 
colony_lct_water<-rbind(colony_lct_water1, colony_lct_water2)

# LCT air now #

species_lct_air1<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(meanPropBelowLCT=mean(meanpropBelow_lct_air), birds=n_distinct(individ_id)) %>%
dplyr::filter(birds>=5) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::summarise(metric="prop_time_below_LCT_air",  metric_mean=mean(meanPropBelowLCT), metric_sd=sd(meanPropBelowLCT), metric_se=metric_sd/sqrt(n_distinct(colony)), type="species")

# Prop time below LCT_air - colony
colony_lct_air1<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(metric="prop_time_below_LCT_air", metric_mean=mean(meanpropBelow_lct_air), metric_sd=sd(meanpropBelow_lct_air), metric_se=metric_sd/sqrt(n_distinct(individ_id)), birds=n_distinct(individ_id), type="colony") %>%
dplyr::filter(birds>=5) 

# degrees below LCT_air - species
species_lct_air2<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(meanPropBelowLCT=mean(degrees_below_lct_air, na.rm=TRUE), birds=n_distinct(individ_id)) %>%
dplyr::filter(birds>=5) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::summarise(metric="degrees_below_LCT_air",  metric_mean=mean(meanPropBelowLCT), metric_sd=sd(meanPropBelowLCT), metric_se=metric_sd/sqrt(n_distinct(colony)), type="species")

# degrees below LCT_air - colony
colony_lct_air2<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(metric="degrees_below_LCT_air", metric_mean=mean(degrees_below_lct_air, na.rm=TRUE), metric_sd=sd(degrees_below_lct_air, na.rm=TRUE), metric_se=metric_sd/sqrt(n_distinct(individ_id)), birds=n_distinct(individ_id), type="colony") %>%
dplyr::filter(birds>=5) 

# prop land - species
species_lct_air3<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(meanPropBelowLCT=mean(meanPropLand_below_air, na.rm=TRUE), birds=n_distinct(individ_id)) %>%
dplyr::filter(birds>=5) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::summarise(metric="prop_land",  metric_mean=mean(meanPropBelowLCT), metric_sd=sd(meanPropBelowLCT), metric_se=metric_sd/sqrt(n_distinct(colony)), type="species")

# prop land- colony
colony_lct_air3<-lct_files_i %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(metric="prop_land", metric_mean=mean(meanPropLand_below_air, na.rm=TRUE), metric_sd=sd(meanPropLand_below_air, na.rm=TRUE), metric_se=metric_sd/sqrt(n_distinct(individ_id)), birds=n_distinct(individ_id), type="colony") %>%
dplyr::filter(birds>=5) 

species_lct_air<-rbind(species_lct_air1, species_lct_air2, species_lct_air3) 
colony_lct_air<-rbind(colony_lct_air1, colony_lct_air2, colony_lct_air3)

# Save all results #
species_lct_water_all<-rbind(species_lct_water_all, species_lct_water)
species_lct_air_all<-rbind(species_lct_air_all, species_lct_air)
colony_lct_water_all<-rbind(colony_lct_water_all, colony_lct_water)
colony_lct_air_all<-rbind(colony_lct_air_all, colony_lct_air)
}

# now we plot the results 

species_lct_water_all<-species_lct_water_all %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(metric=factor(metric, levels=c("prop_time_below_LCT_water", "degrees_below_LCT_water")))

colony_lct_water_all<-colony_lct_water_all %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(metric=factor(metric, levels=c("prop_time_below_LCT_water", "degrees_below_LCT_water")))

plot1<-ggplot() +
geom_pointrange(data=colony_lct_water_all, aes(x=species, y=metric_mean, ymin=metric_mean-1.96*metric_se, ymax=metric_mean + 1.96*metric_se, color=species), position=position_dodge2(width=0.2), alpha=0.1) +
geom_pointrange(data=species_lct_water_all, aes(x=species, y=metric_mean, ymin=metric_mean-1.96*metric_se, ymax=metric_mean + 1.96*metric_se, color=species)) +
facet_wrap(~metric, scales="free_y") +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab("Metric value") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") 

species_lct_air_all2<-species_lct_air_all %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(metric=factor(metric, levels=c("prop_time_below_LCT_air", "degrees_below_LCT_air", "prop_land")))

colony_lct_air_all2<-colony_lct_air_all %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(metric=factor(metric, levels=c("prop_time_below_LCT_air", "degrees_below_LCT_air", "prop_land")))

plot2<-ggplot() +
geom_pointrange(data=colony_lct_air_all2, aes(x=species, y=metric_mean, ymin=metric_mean-1.96*metric_se, ymax=metric_mean + 1.96*metric_se, color=species), position=position_dodge2(width=0.2), alpha=0.1) +
geom_pointrange(data=species_lct_air_all2, aes(x=species, y=metric_mean, ymin=metric_mean-1.96*metric_se, ymax=metric_mean + 1.96*metric_se, color=species)) +
facet_wrap(~metric, scales="free_y") +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab("Metric value") +
xlab("") +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") 
 
