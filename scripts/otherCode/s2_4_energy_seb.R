# Script where, for every individual seabird from Svalbard I calculate mean +/ SD, SE, upper/lower energy expenditure & time spent in different activities
# This is done for every complete day of tracking data (output = csv file) # 
# I could for every colony-species combo, provide a plot showing change in energy, and time spent in different behaviors #
# and change in SST? #

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

### Step 1: Figure out which birds we are working with ####

args <- commandArgs(trailingOnly = TRUE)
print("Step 1: Load id catalogue")
input_file <- args[1]
energyAll <- read.csv(input_file)

# Subset to birds from Svalbard # (not Bjornoya?)
energySvalbard<-subset(energyAll, colony %in% c("Isfjorden", "Kongsfjorden", "Alkefjellet", "Jan Mayen", "Hornsund"))

### Step 2: Now we loop through every individual in this dataset ###

# We will use this to turn DEE into kj.g
species<-data.frame(species=c("Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot", "Little auk", "Atlantic puffin"))
species$allometryCoef<-c(0.765, 0.717, 0.689, 0.689, 0.689, 0.689)

print("Step 2: looping through energy files...")
energyFiles<-list.files("./tmp/", full.names=TRUE)
speciesList<-unique(energySvalbard$species)

# Make a list to save all results in
energySval<-list()

for (i in 1:length(speciesList)) {

# Subset species i
speciesSub<-speciesList[i]

# Subset dataset to different colonies #
speciesSubdf<-subset(energySvalbard, species==speciesSub)
colonies<-unique(speciesSubdf$colony)

# Now we loop through every colony & save results in a list
energyCol<-list()

for (j in 1:length(colonies)) {

# Subset species i
colonySub<-colonies[j]

# Subset dataset to different colonies #
colonySubdf<-subset(speciesSubdf, colony==colonySub)
birds<-unique(colonySubdf$individ_id)

# Now we loop through every bird & save results in a list
energyBird<-list()

for (k in 1:length(birds)) {

#for (k in 1:5) {

### Now I need to grep what's in my results & summarize by date ### 
print(paste0("Step2: calculating for species ", i, " colony ", j, ": Bird ", k))  
    
# Subset bird j
birdSub<-birds[k]
birdSub<-gsub("-", "_", birdSub)
birdSub<-gsub("ø", "o", birdSub)
    
# Find csv file
print("Opening file...")
birdSub<-fread(energyFiles[grep(birdSub, energyFiles)])

if (speciesSub %in% c("Black-legged kittiwake", "Northern fulmar")) {
birdSub$tActive<-0
}

# Summarize by day of the year what is going on
birdSum<-birdSub %>%
dplyr::left_join(species, by=c("species")) %>%
dplyr::mutate(DEEg=DEEkJ/(weight^allometryCoef)) %>%
dplyr::group_by(species, colony, individ_id, date) %>%
dplyr::summarise(DEEg_mean=mean(DEEg), DEEg_sd=sd(DEEg), DEEg_se=DEEg_sd/sqrt(100),
FlightHrs_mean=mean(tFlight), FlightHrs_sd=sd(tFlight), FlightHrs_se=FlightHrs_sd/sqrt(100),
LandHrs_mean=mean(tLand), LandHrs_sd=sd(tLand), LandHrs_se=LandHrs_sd/sqrt(100),
RestWaterHrs_mean=mean(tRestWater), RestWaterHrs_sd=sd(tRestWater), RestWaterHrs_se=RestWaterHrs_sd/sqrt(100),
ActiveHrs_mean=mean(tActive), ActiveHrs_sd=sd(tActive), ActiveHrs_se=ActiveHrs_sd/sqrt(100),
ForageHrs_mean=mean(tForage), ForageHrs_sd=sd(tForage), ForageHrs_se=ForageHrs_sd/sqrt(100),
SST_mean=mean(sst_random), SST_sd=sd(sst_random), SST_se=SST_sd/sqrt(100)
)

# Save result #
energyBird<-rbind(energyBird, birdSum)

}

# Save result
energyCol<-rbind(energyCol, energyBird)

}

# Save result
energySval<-rbind(energySval, energyCol)

}

# make some plots #
energySvalSum<-energySval %>%
ungroup() %>%
dplyr::group_by(individ_id) %>%
dplyr::mutate(doy=as.numeric(difftime(date, as.Date(paste0(substr(date, 1, 4), "-01-01")), unit=c("days")))+1) %>%
ungroup() %>%
dplyr::group_by(species, colony, doy) %>%
dplyr::summarise(meanDEE=mean(DEEg_mean), sdDEE=sd(DEEg_mean), birds=n_distinct(individ_id), seDEE=sdDEE/sqrt(birds),
meanForage=mean(ForageHrs_mean), sdForage=sd(ForageHrs_mean), birds=n_distinct(individ_id), seForage=sdForage/sqrt(birds),
meanFlight=mean(FlightHrs_mean), sdFlight=sd(FlightHrs_mean), birds=n_distinct(individ_id), seFlight=sdFlight/sqrt(birds),
meanActive=mean(ActiveHrs_mean), sdActive=sd(ActiveHrs_mean), birds=n_distinct(individ_id), seActive=sdActive/sqrt(birds),
meanRestWater=mean(RestWaterHrs_mean), sdRestWater=sd(RestWaterHrs_mean), birds=n_distinct(individ_id), seRestWater=sdRestWater/sqrt(birds),
meanLand=mean(LandHrs_mean), sdLand=sd(LandHrs_mean), birds=n_distinct(individ_id), seLand=sdLand/sqrt(birds),
meanSST=mean(SST_mean), sdSST=sd(SST_mean), birds=n_distinct(individ_id), seSST=sdSST/sqrt(birds)
)

# Create list of weeks to roll through
dates<-data.frame(dateKeep=seq(as.Date("2021-07-01"), as.Date("2022-06-30"), 1))
dates$doy2<-1:nrow(dates)
dates$month<-as.numeric(substr(dates$date, 6, 7))
dates$day<-as.numeric(substr(dates$date, 9, 10))
dates$doy<-as.numeric(difftime(dates$dateKeep, as.Date(paste0(substr(dates$dateKeep, 1, 4), "-01-01")), unit=c("days"))) + 1

energySvalSum2<-energySvalSum %>%
dplyr::left_join(dates, by=c("doy"))

datesSub<-dates %>%
dplyr::filter(day==1)

EnergyPlot<-energySvalSum2 %>%
ggplot(aes(x=doy2, y=meanDEE)) +
geom_smooth(aes(group=colony, color=colony, fill=colony)) +
facet_wrap(~species) +
theme_bw() +
xlab("") +
ylab("Daily energy expenditure (kJ.g)") +
scale_x_continuous(breaks=datesSub$doy2, labels=c("July", "August", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "June")) 

pdf("./results/figures/supplementary/seb_energy.pdf")
plot(EnergyPlot)
dev.off()

ActivityPlot<-energySvalSum2 %>%
ggplot(aes(x=doy2, y=meanForage)) +
geom_smooth(aes(group=colony, linetype=colony, color="Forage", fill="Forage")) +
geom_smooth(aes(x=doy2, y=meanFlight, group=colony, linetype=colony, color="Flight", fill="Flight")) +
geom_smooth(aes(x=doy2, y=meanActive, group=colony, linetype=colony, color="Active", fill="Active")) +
geom_smooth(aes(x=doy2, y=meanRestWater, group=colony, linetype=colony, color="RestWater", fill="RestWater")) +
geom_smooth(aes(x=doy2, y=meanLand, group=colony, linetype=colony, color="Land", fill="Land")) +
facet_wrap(~species) +
theme_bw() +
xlab("Day of the year") +
ylab("Time spent in activity (hours)") +
labs(color="Behavior", fill="Behavior") +
scale_x_continuous(breaks=datesSub$doy2, labels=c("July", "August", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "June")) 

pdf("./results/figures/supplementary/seb_activity.pdf")
plot(ActivityPlot)
dev.off()

SSTPlot<-energySvalSum2 %>%
ggplot(aes(x=doy2, y=meanSST)) +
geom_smooth(aes(group=colony, color=colony, fill=colony)) +
facet_wrap(~species) +
theme_bw() +
xlab("") +
ylab("SST") +
scale_x_continuous(breaks=datesSub$doy2, labels=c("July", "August", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "June")) 

pdf("./results/figures/supplementary/seb_sst.pdf")
plot(SSTPlot)
dev.off()

# Save output files
print("Saving output files...")

# Number 3
output_file1 <- args[2]
print("Saving output file 1")
write.csv(energySval, file = output_file1, row.names = FALSE) # Species-mean deviance