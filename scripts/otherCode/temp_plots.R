### These are extra plots I have made t understand what is happening in the results ###

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

### Step 1: Plot migratory distance per population ###

# Read in main dataset #
totalCosts<-read.csv("./results/tables/main/table9_totalNBCosts.csv")

# Make a population summary of migratory distances # 

minSampleSize<-5
reps<-50

colonyRes<-totalCosts %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::summarise(birds=n_distinct(individ_id), meanDEE=mean(totalDEE), minDEE=min(totalDEE), maxDEE=max(totalDEE), sdDEE=sd(totalDEE), seDEE=sdDEE/sqrt(birds),
  meanDev=mean(totDeviance), minDev=min(totDeviance), maxDev=max(totDeviance), sdDev=sd(totDeviance), seDev=sdDev/sqrt(birds),
  meanCov=mean(totCov), minCov=min(totCov), maxCov=max(totCov), sdCov=sd(totCov), seCov=sdCov/sqrt(birds),
  meanDEE2=mean(DEE_scale), minDEE2=min(DEE_scale), maxDEE2=max(DEE_scale), sdDEE2=sd(DEE_scale), seDEE2=sdDEE2/sqrt(birds),
  meanDev2=mean(deviance_scale), minDev2=min(deviance_scale), maxDev2=max(deviance_scale), sdDev2=sd(deviance_scale), seDev2=sdDev2/sqrt(birds),
  meanCov2=mean(cov_scale), minCov2=min(cov_scale), maxCov2=max(cov_scale), sdCov2=sd(cov_scale), seCov2=sdCov2/sqrt(birds), meanMigratory=mean(MigratoryDistKm), sdMigratory=sd(MigratoryDistKm), 
  seMigratory=sdMigratory/sqrt(birds))  %>% 
  dplyr::filter(birds>=minSampleSize) %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>%
  dplyr::summarise(meanDEEg=mean(meanDEE), sd=mean(sdDEE), seDEE=mean(seDEE), minDEE=mean(minDEE), maxDEE=mean(maxDEE),
  meanDev=mean(meanDev), sd2=mean(sdDev), seDev=mean(seDev), minDev=mean(minDev), maxDev=mean(maxDev),
  meanCov=mean(meanCov), sd3=mean(sdCov), seCov=mean(seCov), minCov=mean(minCov), maxCov=mean(maxCov),
  meanDEE2=mean(meanDEE2), seDEE2=mean(seDEE2),
  meanDev2=mean(meanDev2), seDev2=mean(seDev2), 
  meanCov2=mean(meanCov2), seCov2=mean(seCov2), meanBirds=mean(birds), meanMigr=mean(meanMigratory), sdMigr=mean(sdMigratory), seMigr=mean(seMigratory)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  ungroup() %>%
  dplyr::group_by(species) %>%
  dplyr::arrange(meanCov)%>%
  dplyr::mutate(colonySp=paste0(colony, "_", species))
  
# Order colonyRes populations according to value
colonyOrder<-unique(colonyRes$colonySp)
colonyRes$colonySp<-factor(colonyRes$colonySp, levels=colonyOrder)	

# open colony information & join so we have longitude & latitude & project into correct format
colony.summary<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")

speciesMatch<-data.frame(speciesLatin=c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle"), 
                         species=c("Brünnich's guillemot", "Black-legged kittiwake", "Common guillemot", "Atlantic puffin", "Northern fulmar", "Little auk"))

colonyMatch<-colony.summary %>%
 dplyr::filter(species %in% c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle")) %>%
  dplyr::group_by(colony) %>%
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
                          "RU", "NO", "GR", "RU", "RU", "NO", "RU", "NO", "CA", "IC", "IC", "IC", "IC", "IC", "IC", 
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
 dplyr::slice(1) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
												 dplyr::select(species, colony, col_lon, col_lat)
												 
colonies_lox2<-colonies_lox %>%
dplyr::select(-species)

colonyRes_lox<-colonyRes %>%
dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ungroup() %>%
dplyr::arrange(col_lat) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

colonyOrder<-unique(colonyRes_lox$colony)
colonyOrder2<-unique(colonyRes_lox$colonyName)

colonyRes_lox$colony<-factor(colonyRes_lox$colony, levels=colonyOrder)
colonyRes_lox$colonyName<-factor(colonyRes_lox$colonyName, levels=colonyOrder2)
												 

Figure_migr<-ggplot() +
 geom_pointrange(data=colonyRes_lox, aes(y=meanMigr, x=colonyName, ymin=meanMigr - 1.96*seMigr, ymax=meanMigr + 1.96*seMigr, color=species)) +
 facet_wrap(~species, scales="free_y", nrow=3) +
 scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+
  theme_bw() +
  ylab("Migratory distance (km)") +
  xlab("") +
  theme_bw() +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("./results/figures/supplementary/figure_migr.pdf", width=15, height=9)
plot(Figure_migr)
dev.off()

#### Step 2: PLOT WHERE POPS ARE WHEN WEE is HIGH #####

# First we open the appropriate files
allResults_deviation<-list.files("./results/tables/main/", full.names=TRUE)
deviation_weekly<-allResults_deviation[grepl("/weeklydeviance_", allResults_deviation)]

# make an empty list to save data in
devianceMean_weekly<-list()

for (m in 1:length(deviation_weekly)) {
  
  deviance_sub<-readRDS(deviation_weekly[m])
  
  devianceMean_weekly<-rbind(devianceMean_weekly, deviance_sub)
  
}

# Make a species-mean (weighted for varying sample size):
print("Calculating species weighted averages by week...")

devianceSpeciesRes<-list()

for (p in 1:reps) {

print(paste0("Rep ", p))

devianceSpeciesRep<-devianceMean_weekly %>%
  ungroup() %>%
  dplyr::filter(rep==p) %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::mutate(birds=n_distinct(individ_id)) %>%
  dplyr::filter(birds >=minSampleSize) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, individ_id) %>%
  dplyr::mutate(meanPropFlight=mean(propFlight), meanPropActive=mean(propActive), meanPropRest=mean(propRest), meanPropLand=mean(propLand), meanPropForage=mean(propForage), annualSST=mean(meanSST)) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, individ_id, weekNo) %>%
  dplyr::mutate(devianceFlight=(propFlight - meanPropFlight)/meanPropFlight, devianceActive=(propActive - meanPropActive)/meanPropActive, devianceRest=(propRest - meanPropRest)/meanPropRest, devianceLand=(propLand - meanPropLand)/meanPropLand, devianceForage=(propForage - meanPropForage)/meanPropForage, devianceSST=(meanSST - annualSST)) %>%
  dplyr::group_by(species, colony, individ_id) %>%
  dplyr::mutate(repsTot=n_distinct(rep)) %>%
  ungroup() %>%
  replace_na(list(devianceFlight=0, devianceActive=0, devianceLand=0)) %>%
  dplyr::group_by(rep, species, colony, weekNo) %>%
  dplyr::reframe(averageDeviance=mean(devianceDEE), averageDevianceFlight=mean(devianceFlight), averageDevianceActive=mean(devianceActive, na.rm=TRUE), 
  averageDevianceRest=mean(devianceRest), averageDevianceLand=mean(devianceLand), averageDevianceForage=mean(devianceForage), averageDevianceSST=mean(devianceSST), sampleSize=n_distinct(individ_id), weight=1/sampleSize, birds=n_distinct(individ_id)) %>%
  ungroup() %>%
  #dplyr::filter(birds>=minSampleSize) %>%
  dplyr::group_by(rep, species, weekNo) %>%
  dplyr::reframe(weightedDeviance=sum(averageDeviance*weight)/sum(weight), weightedDevianceFlight=sum(averageDevianceFlight*weight)/sum(weight), weightedDevianceActive=sum(averageDevianceActive*weight)/sum(weight),
  weightedDevianceRest=sum(averageDevianceRest*weight)/sum(weight), weightedDevianceLand=sum(averageDevianceLand*weight)/sum(weight), weightedDevianceForage=sum(averageDevianceForage*weight)/sum(weight),
  weightedDevianceSST=sum(averageDevianceSST*weight)/sum(weight))%>%
  ungroup()
  
  devianceSpeciesRes<-rbind(devianceSpeciesRes, devianceSpeciesRep)
  
  }

devianceSpecies<-devianceSpeciesRes %>%
  dplyr::group_by(species, weekNo) %>%
  dplyr::reframe(mean_weightedDeviance=mean(weightedDeviance), sd_weightedDeviance=sd(weightedDeviance), se=sd_weightedDeviance/sqrt(reps),
  mean_weightedDevianceFlight=mean(weightedDevianceFlight), sd_weightedDevianceFlight=sd(weightedDevianceFlight), seFlight=sd_weightedDevianceFlight/sqrt(reps),
  mean_weightedDevianceActive=mean(weightedDevianceActive), sd_weightedDevianceActive=sd(weightedDevianceActive), seActive=sd_weightedDevianceActive/sqrt(reps),
  mean_weightedDevianceRest=mean(weightedDevianceRest), sd_weightedDevianceRest=sd(weightedDevianceRest), seRest=sd_weightedDevianceRest/sqrt(reps),
  mean_weightedDevianceLand=mean(weightedDevianceLand), sd_weightedDevianceLand=sd(weightedDevianceLand), seLand=sd_weightedDevianceLand/sqrt(reps),
  mean_weightedDevianceForage=mean(weightedDevianceForage), sd_weightedDevianceForage=sd(weightedDevianceForage), seForage=sd_weightedDevianceForage/sqrt(reps),
  mean_weightedDevianceSST=mean(weightedDevianceSST), sd_weightedDevianceSST=sd(weightedDevianceSST), seSST=sd_weightedDevianceSST/sqrt(reps))%>%
  dplyr::mutate(scale="weekly") %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

print("Calculating colony weighted averages...")

devianceColonyRes<-list()

for (q in 1:reps) {

print(paste0("Rep ", q))

devianceColonyRep<-devianceMean_weekly %>%
  ungroup() %>%
  dplyr::filter(rep==q) %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::mutate(birds=n_distinct(individ_id)) %>%
  dplyr::filter(birds >=minSampleSize) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, individ_id) %>%
  dplyr::mutate(meanPropFlight=mean(propFlight), meanPropActive=mean(propActive), meanPropRest=mean(propRest), meanPropLand=mean(propLand), meanPropForage=mean(propForage), annualSST=mean(meanSST)) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, individ_id, weekNo) %>%
  dplyr::mutate(devianceFlight=(propFlight - meanPropFlight)/meanPropFlight, devianceActive=(propActive - meanPropActive)/meanPropActive, devianceRest=(propRest - meanPropRest)/meanPropRest, devianceLand=(propLand - meanPropLand)/meanPropLand, devianceForage=(propForage - meanPropForage)/meanPropForage, devianceSST=(meanSST - annualSST)) %>%
  dplyr::group_by(species, colony, individ_id) %>%
  dplyr::mutate(repsTot=n_distinct(rep)) %>%
  ungroup() %>%
  replace_na(list(devianceFlight=0, devianceActive=0, devianceLand=0)) %>%
  ungroup() %>%
  dplyr::group_by(rep, species, colony, weekNo) %>%
  dplyr::reframe(averageDeviance=mean(devianceDEE), averageDevianceFlight=mean(devianceFlight), averageDevianceActive=mean(devianceActive), 
  averageDevianceRest=mean(devianceRest), averageDevianceLand=mean(devianceLand), averageDevianceForage=mean(devianceForage), averageDevianceSST=mean(devianceSST), sampleSize=n_distinct(individ_id), weight=1/sampleSize, birds=n_distinct(individ_id)) %>%
  ungroup()  
  
  devianceColonyRes<-rbind(devianceColonyRes, devianceColonyRep)
  
  }


deviancePopulation<-devianceColonyRes %>%
  ungroup() %>%
  dplyr::group_by(species, colony, weekNo) %>%
  dplyr::summarise(mean_weightedDeviance=mean(averageDeviance), sd_weightedDeviance=sd(averageDeviance), se=sd_weightedDeviance/sqrt(reps),
  mean_weightedDevianceFlight=mean(averageDevianceFlight), sd_weightedDevianceFlight=sd(averageDevianceFlight), seFlight=sd_weightedDevianceFlight/sqrt(reps),
  mean_weightedDevianceActive=mean(averageDevianceActive), sd_weightedDevianceActive=sd(averageDevianceActive), seActive=sd_weightedDevianceActive/sqrt(reps),
  mean_weightedDevianceRest=mean(averageDevianceRest), sd_weightedDevianceRest=sd(averageDevianceRest), seRest=sd_weightedDevianceRest/sqrt(reps),
  mean_weightedDevianceLand=mean(averageDevianceLand), sd_weightedDevianceLand=sd(averageDevianceLand), seLand=sd_weightedDevianceLand/sqrt(reps),
  mean_weightedDevianceForage=mean(averageDevianceForage), sd_weightedDevianceForage=sd(averageDevianceForage), seForage=sd_weightedDevianceForage/sqrt(reps),
  mean_weightedDevianceSST=mean(averageDevianceSST), sd_weightedDevianceSST=sd(averageDevianceSST), seSST=sd_weightedDevianceSST/sqrt(reps)) %>%
  dplyr::mutate(scale="weekly") %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

#### Step 3: Manually determine key periods for plotting ####

# Determine species, weeks & pops of interest #

sub1<-subset(devianceSpecies, species=="Little auk") 

sub2<-subset(sub1, mean_weightedDeviance > 0.02) # Subset to highest deviances

sub3<-sub2 %>%
ungroup() %>%
dplyr::group_by(weekNo) %>%
dplyr::summarise( maxDEE=max(mean_weightedDeviance), minDEE=min(mean_weightedDeviance), minSST=min(mean_weightedDevianceSST), maxActive=max(mean_weightedDevianceActive), maxFlight=max(mean_weightedDevianceFlight))

# Common guillemot #

CoGu1_species<-subset(devianceSpecies, species=="Common guillemot" & weekNo %in% c(16, 17, 18, 19, 20, 21, 22)) 
#maxDeviance1<-max(CoGu1_species$mean_weightedDeviance) # Basically here I try and see which populations are driving the species-level trend
CoGu1<-subset(deviancePopulation, species=="Common guillemot" & weekNo %in% c(16, 17, 18, 19, 20, 21, 22)) 
CoGu1<-CoGu1%>%
dplyr::left_join(sub3, by=c("weekNo")) %>%
dplyr::left_join(colonyNames, by=c("colony")) %>%
dplyr::filter(mean_weightedDeviance > maxDEE)
CoGu1$period<-"Winter"

CoGu2_species<-subset(devianceSpecies, species=="Common guillemot" & weekNo %in% c(24, 25, 26, 27, 28, 29, 30)) 
#maxDeviance2<-max(CoGu2_species$mean_weightedDeviance) # Basically here I try and see which populations are driving the species-level trend
CoGu2<-subset(deviancePopulation, species=="Common guillemot" & weekNo %in% c(24, 25, 26, 27, 28, 29, 30)) 
CoGu2<-CoGu2%>%
dplyr::left_join(sub3, by=c("weekNo")) %>%
dplyr::left_join(colonyNames, by=c("colony")) %>%
dplyr::filter(mean_weightedDeviance > maxDEE)
CoGu2$period<-"Spring"

allCogu<-rbind(CoGu1, CoGu2)

# Brunnich's guillemot #

BrGu1_species<-subset(devianceSpecies, species=="Brünnich's guillemot" & weekNo %in% c(13, 14, 15, 16, 17, 18, 19, 20, 21)) 
maxDeviance1<-max(BrGu1_species$mean_weightedDeviance) # Basically here I try and see which populations are driving the species-level trend
BrGu1<-subset(deviancePopulation, species=="Common guillemot" & weekNo %in% c(13, 14, 15, 16, 17, 18, 19, 20, 21) & mean_weightedDeviance > maxDeviance) 
BrGu1$period<-"Winter"

BrGu2_species<-subset(devianceSpecies, species=="Common guillemot" & weekNo %in% c(23, 24, 25, 26, 27, 28, 29, 30)) 
maxDeviance2<-max(BrGu2_species$mean_weightedDeviance) # Basically here I try and see which populations are driving the species-level trend
BrGu2<-subset(deviancePopulation, species=="Common guillemot" & weekNo %in% c(23, 24, 25, 26, 27, 28, 29, 30) & mean_weightedDeviance > maxDeviance2) 
BrGu2$period<-"Spring"

allBrGu<-rbind(BrGu1, BrGu2)

### Step 4: Plot! ###

# Determine where locations are stored
allResults<-list.files("./tmp/", full.names=TRUE)
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

# Determine specific period-colony combinations
periodsPops<-allCogu %>%
ungroup() %>%
dplyr::group_by(period) %>%
dplyr::count(colony)

# Subset main dataset to get information on which IDs we need to extract
birdids<-subset(devianceMean_weekly, species=="Common guillemot" & colony %in% unique(allCogu$colony) & weekNo %in% unique(allCogu$weekNo))
colonies<-unique(birdids$colony)

# Loop through colonies & save in empty list below
allLoxPeriod<-list()

for (i in 1:nrow(periodsPops)) {

print(paste0("periodPop", i, "/", nrow(periodsPops)))

Subset<-periodsPops[i,]
colSub<-Subset$colony
periodSub<-Subset$period
periodSubWeeks<-subset(allCogu, period==periodSub)
colSubBirds<-subset(birdids, colony==colSub & weekNo %in% unique(periodSubWeeks$weekNo))
ids<-unique(colSubBirds$individ_id)

# Make a list to save locations in
allLox<-list()

print("Finding locations...")

for (j in 1:length(ids)) {

#for (j in 1:5) {

print(paste0("Bird", j, "/", length(ids)))

energyOpen<-fread(energyRes_day[grepl(ids[j], energyRes_day)])

# Subset to relevant weeks #
energyOpen$month<-as.numeric(substr(energyOpen$date, 6, 7))
energyOpen$day<-as.numeric(substr(energyOpen$date, 9, 10))
energyMean<-energyOpen %>%
dplyr::group_by(species, colony, date, day, month) %>%
dplyr::summarise(meanLon=mean(mean.lon), meanLat=mean(mean.lat)) %>%
dplyr::inner_join(dates_weekly, by=c("month", "day")) %>%
dplyr::filter(weekNo %in% unique(colSubBirds$weekNo)) 
energyMean$period<-periodSub

# Save results
allLox<-rbind(allLox, energyMean)

}

# Map where these are

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

# project locations
coordinates(allLox)<-~ meanLon + meanLat
proj4string(allLox)<-projection_84
allLoxTrans<-data.frame(spTransform(allLox, projection_NA))

# Make a polygon around these points
allLocations_polygon<-allLoxTrans %>%
  st_as_sf(coords = c("coords.x1", "coords.x2")) %>%
  group_by(species, colony, period) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_convex_hull() %>%
  st_set_crs(projection_NA)
  
# Convert to a line object so I can change thickness of outline...
polygon_lines <- st_cast(allLocations_polygon, "MULTILINESTRING")

# Save result
allLoxPeriod[[i]]<-polygon_lines

}

# Map all of these #

# Open country layers
print("Making map S1...")
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# open colony information & join so we have longitude & latitude & project into correct format
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
  
  colonies_lox<-colonyMatch %>%
  dplyr::left_join(colonyNames, by=c("colony")) %>%
  ungroup() %>%
  arrange(colonyName) %>%
  dplyr::group_by(colonyName) %>%
  dplyr::slice(1)
  
coordinates(colonies_lox)<-~col_lon + col_lat
proj4string(colonies_lox)<-projection_84
colonies_lox_trans<-data.frame(spTransform(colonies_lox, projection_NA))

# Somehow join the polygons
joined_lines <- do.call(rbind, allLoxPeriod)

# Figure out extent so I can crop my map
bbox <- st_bbox(joined_lines)
xmin <- bbox["xmin"]
xmax <- bbox["xmax"]
ymin <- bbox["ymin"]
ymax <- bbox["ymax"]

# Make on other which is subset to correct colonies
colonies_lox_trans_sub<-subset(colonies_lox_trans, colony %in% unique(joined_lines$colony))

# Make a plot
ggplot() +
  #geom_raster(data=subset(SSTProject_df, !is.na(layer)), aes(x=x, y=y, fill=layer), alpha=0.5, color=NA) +
  #scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(-2.5, 33)) +
  geom_sf(data=joined_lines, aes(color=species), fill=NA, size=4) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=colonies_lox_trans, aes(x=coords.x1, y=coords.x2),  cex=2, fill="yellow", shape=21) + 
  geom_point(data=colonies_lox_trans_sub, aes(x=coords.x1, y=coords.x2), fill="purple", cex=2, shape=21) + 
  scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=projection_NA, xlim=c(xmin-100000, xmax + 100000), ylim=c(ymin-100000, ymax + 100000)) +
  xlab("") +
  ylab("") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~period)
  
  ggplot() +
  #geom_raster(data=subset(SSTProject_df, !is.na(layer)), aes(x=x, y=y, fill=layer), alpha=0.5, color=NA) +
  #scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(-2.5, 33)) +
  geom_sf(data=joined_lines, aes(color=colony), fill=NA, size=4) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=colonies_lox_trans, aes(x=coords.x1, y=coords.x2),  cex=2, fill="yellow", shape=21) + 
  geom_point(data=colonies_lox_trans_sub, aes(x=coords.x1, y=coords.x2, fill=colony), cex=2, shape=21) + 
  #scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  #scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=projection_NA, xlim=c(xmin-100000, xmax + 100000), ylim=c(ymin-100000, ymax + 100000)) +
  xlab("") +
  ylab("") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~period)

