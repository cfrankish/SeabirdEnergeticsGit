## This script prepares population-specific activity distributions ##
## To do this, by species, the script will determine which model populations it needs to map ##
## Then it will search for individual birds ##
## Then it will estimate a distribution of behaviour by pop ##
## The output will be saved as 'species_population.csv ##

library(dplyr)
library(ncdf4)
library(raster)
library(sp)
library(sf)
library(lubridate)
library(data.table)
library(tidyr)
library(gridExtra)
library(ggplot2)

#### Step 0: setting up basic conditions ####

# Set-up number of iterations...
overall.iterations<-50 # how many times this is calculated per individual
print(paste0("Determining activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1] # This will read in a species-specific density map - input_file1<-"data/popdata_raw/SEATRACK_Abundance_Model_Rissa_tridactyla_Ver_3_1.nc"
input_file2<- args[2] # This is list of all bird Ids (id catalogue) - input_file2<-"results/tables/main/table1_idcatalogue.csv"

# open two input files
nc<-nc_open(input_file1) # nc<-nc_open("data/popdata_raw/SEATRACK_Abundance_Model_Rissa_tridactyla_Ver_3_1.nc")
idCatalogue<-read.csv(input_file2) # idCatalogue<-read.csv("results/tables/main/table1_idcatalogue.csv")

#### Step 1: Determine number of model colonies & other pops ####

print("Step 1: determine number of model colonies to loop through for...")

# Extract file name to determine which species this is
mapsFileName<-nc$filename
print(mapsFileName)
speciesLatinSub<-substr(mapsFileName, 91, nchar(mapsFileName)-11)
#speciesLatinSub<-substr(mapsFileName, 43, nchar(mapsFileName)-11)
speciesMatch<-data.frame(speciesLatin=c("Alle_alle", "Fratercula_arctica", "Fulmarus_glacialis", "Rissa_tridactyla", "Uria_lomvia", "Uria_aalge"), species=c("Little auk", "Atlantic puffin", "Northern fulmar",
"Black-legged kittiwake", "Brünnich's guillemot", "Common guillemot"))
speciesSub<-subset(speciesMatch, speciesLatin==speciesLatinSub)$species
print(paste0(speciesSub, "..."))

# Subset to relevant activity budget data files
idCatalogueSub<-subset(idCatalogue, species==speciesSub)

# Figure out number of model colonies & colonies for this species
modelcolonies <- (ncvar_get(nc,"SmcolName")) # Find list of populations
colonies <- (ncvar_get(nc,"colonyName")) # Find list of ALL colonies
colcode <- (ncvar_get(nc,"SmcolCode")) # Colony code
meta <- data.frame(modelcolonies, colonies, colcode) # Make some metadata so easier to understand raster structure
meta$species<-speciesSub # Add species names

#### Step 2: Summarize time in activity per population ####

# Determine unique number of model colonies
modelColonies<-unique(meta$modelcolonies)

# List to save all info in
meanActSpecies<-list()
meanActSpecies_day<-list()

# Loop through These
for (i in 1:length(modelColonies)) {

# Subset to model colony i
popSub<-modelColonies[i]
popSub<-gsub("_", " ", popSub)
popSub<-gsub("oeg", "öeg", popSub)
popSub<-gsub("oe", "ø", popSub)
popSub<-gsub("Aa", "Å", popSub)
#popSub<-gsub("ae", "æ", popSub)
#popSub<-gsub("oe", "ö", popSub)

popSub<-ifelse(popSub =="Farø Islands", "Faroe Islands", popSub)
popSub<-ifelse(popSub =="Jarsteinen", "Karmøy", popSub)

print(popSub)

# Determine individuals we will loop through #
idColonySub<-subset(idCatalogueSub, colony==popSub)
if (nrow(idColonySub)<1) stop(print("Error: no match between colonies"))
birds<-unique(idColonySub$individ_id)
#print(birds)

# Now we determine where the daily files are located 
actfiles<-list.files("./tmp/", full.names=TRUE)
dailyfiles<-actfiles[grepl("Day.csv", actfiles)]

# And loop through every bird #
meanActMonthly<-list()
meanActDaily<-list()

for (j in 1:length(birds)) {

# Print status update #
print(paste0("Summing activity for modelCol ", i, "/", length(modelColonies), " : Summarizig bird ", j, "/", length(birds)))

# Subset to bird j

ID<-birds[j]
ID<-gsub("-", "_", ID) 
ID<-gsub("Hornøya", "Hornoya", ID)

birdSub<-dailyfiles[grepl(ID, dailyfiles)]

if (length(birdSub)<1) {
next}

# open csv
birdCsv<-fread(birdSub)

# Subset to max number of repetitions to save memory
activitySubSub<-subset(birdCsv, rep<=overall.iterations)

# Figure out for which months we will be estimating time spent in activity #

# Determine number of complete months of data #
startMonth<-as.Date(paste0(substr(min(activitySubSub$date), 1, 7), "-01"))
endMonth<-as.Date(paste0(substr(max(activitySubSub$date), 1, 7), "-01"))
endMonthNew<-endMonth %m+% months(1) # Add another month so I get a full month of days
endMonthFinal<-endMonthNew-1 # and remove one day!
fullMonths<-data.frame(date=seq(from=startMonth, to=endMonthFinal, by="day")) # Choose a year with 28 days in Feb
fullMonths$day<-as.numeric(substr(fullMonths$date, 9, 10))
fullMonths$month<-as.numeric(substr(fullMonths$date, 6, 7))

# Attach activity sequence and determine which months are missing data
activityDates<-data.frame(date=as.Date(unique(activitySubSub$date)), type="myData")

fullMonthsLength<-fullMonths %>%
dplyr::mutate(year=as.numeric(substr(date, 1, 4))) %>%
dplyr::group_by(year, month) %>%
dplyr::mutate(days=n_distinct(day)) %>%
ungroup()

fullMonthsAttach<-fullMonthsLength %>%
dplyr::inner_join(activityDates, by=c("date")) %>%
ungroup() %>%
dplyr::group_by(year, month) %>%
dplyr::mutate(daysMerge=n_distinct(day)) %>%
dplyr::filter(daysMerge==days)

# Determine unique year-month combos
fullMonthsYear<-fullMonthsAttach %>%
ungroup() %>%
dplyr::group_by(year) %>%
dplyr::count(month)

# Select a month at random per year 
fullMonthsRandom<-fullMonthsYear %>%
ungroup() %>%
dplyr::group_by(month) %>%
dplyr::slice_sample(n=1) %>%
dplyr::mutate(month_year=paste0(month, "_", year))

# add t_active if needed
if (!any(colnames(activitySubSub) %in% c("tActive"))) {
activitySubSub$tActive<-0
}

# Now we subset activity to these months
activitySubSub$month<-as.numeric(substr(activitySubSub$date, 6, 7))
activitySubSub$year<-as.numeric(substr(activitySubSub$date, 1, 4))
activitySubSub$month_year<-paste0(activitySubSub$month, "_", activitySubSub$year)
meanBehMonths<-subset(activitySubSub, month_year %in% fullMonthsRandom$month_year) # Subset to random month-year combos

# Add time foraging / active if it doesn't exist
# add t_active if needed
if (!any(colnames(activitySubSub) %in% c("tActive"))) {
activitySubSub$tActive<-0
}

# Now we sum time spent in different behaviors per rep-month
sumBehMonth<-meanBehMonths %>%
ungroup() %>%
dplyr::group_by(rep, month) %>%
dplyr::summarise(tForage=sum(tForage), tActive=sum(tActive), tRestWater=sum(tRestWater), tLand=sum(tLand), tFlight=sum(tFlight), 
days=n_distinct(date), totHrs=days*24, totHrsAct=sum(tForage, tActive, tRestWater, tLand, tFlight), meanLon=mean(mean.lon), meanLat=mean(mean.lat)) 

# make sure these behaviors sum to correct number (otherwise there are incomplete months #
sumBehNoMatch<-subset(sumBehMonth, !all.equal(totHrs, totHrsAct, tolerance=0.2))
if (nrow(sumBehNoMatch) >0) stop(print("Error: activity hours weird"))

# Now we estimate mean time in behaviour with confidence intervals #
meanBehMonth<-sumBehMonth %>%
dplyr::ungroup() %>%
dplyr::group_by(month) %>%
dplyr::summarise(tForageMean=mean(tForage), tForageSD=sd(tForage), tForageSE=tForageSD/sqrt(overall.iterations), tForage_lower=tForageMean-1.96*tForageSE, tForage_upper=tForageMean+1.96*tForageSE,
tFlightMean=mean(tFlight), tFlightSD=sd(tFlight), tFlightSE=tFlightSD/sqrt(overall.iterations), tFlight_lower=tFlightMean-1.96*tFlightSE, tFlight_upper=tFlightMean+1.96*tFlightSE,
tActiveMean=mean(tActive), tActiveSD=sd(tActive), tActiveSE=tActiveSD/sqrt(overall.iterations), tActive_lower=tActiveMean-1.96*tActiveSE, tActive_upper=tActiveMean+1.96*tActiveSE,
tLandMean=mean(tLand), tLandSD=sd(tLand), tLandSE=tLandSD/sqrt(overall.iterations), tLand_lower=tLandMean-1.96*tLandSE, tLand_upper=tLandMean+1.96*tLandSE,
tRestWaterMean=mean(tRestWater), tRestWaterSD=sd(tRestWater), tRestWaterSE=tRestWaterSD/sqrt(overall.iterations), tRestWater_lower=tRestWaterMean-1.96*tRestWaterSE, tRestWater_upper=tRestWaterMean+1.96*tRestWaterSE,
totHrs=mean(totHrs), meanLon=mean(meanLon), meanLat=mean(meanLat))

# We also estimate time spent in behaviours by day #

# But first we create a key of day of the year for summarizing
yeardates<-data.frame(date=seq(as.Date("2021-01-01"), as.Date("2021-12-31"), 1))
yeardates$month<-as.numeric(substr(yeardates$date, 6, 7))
yeardates$day<-as.numeric(substr(yeardates$date,9, 10))

meanBehDay<-activitySubSub %>%
dplyr::ungroup() %>%
dplyr::mutate(day=as.numeric(substr(date, 9, 10))) %>%
dplyr::inner_join(yeardates, by=c("month", "day")) %>%
dplyr::group_by(month, day) %>%
dplyr::summarise(tForageMean=mean(tForage), tForageSD=sd(tForage), tForageSE=tForageSD/sqrt(overall.iterations), tForage_lower=tForageMean-1.96*tForageSE, tForage_upper=tForageMean+1.96*tForageSE,
tFlightMean=mean(tFlight), tFlightSD=sd(tFlight), tFlightSE=tFlightSD/sqrt(overall.iterations), tFlight_lower=tFlightMean-1.96*tFlightSE, tFlight_upper=tFlightMean+1.96*tFlightSE,
tActiveMean=mean(tActive), tActiveSD=sd(tActive), tActiveSE=tActiveSD/sqrt(overall.iterations), tActive_lower=tActiveMean-1.96*tActiveSE, tActive_upper=tActiveMean+1.96*tActiveSE,
tLandMean=mean(tLand), tLandSD=sd(tLand), tLandSE=tLandSD/sqrt(overall.iterations), tLand_lower=tLandMean-1.96*tLandSE, tLand_upper=tLandMean+1.96*tLandSE,
tRestWaterMean=mean(tRestWater), tRestWaterSD=sd(tRestWater), tRestWaterSE=tRestWaterSD/sqrt(overall.iterations), tRestWater_lower=tRestWaterMean-1.96*tRestWaterSE, 
tRestWater_upper=tRestWaterMean+1.96*tRestWaterSE, meanLon=mean(mean.lon), mean.lat=mean(mean.lat)) %>%
dplyr::mutate(og_colony=meanBehMonths$colony[1])

# Save results
meanBehMonth$species<-speciesSub
meanBehMonth$colony<-popSub
meanBehMonth$individ_id<-birds[j]

meanBehDay$species<-speciesSub
meanBehDay$colony<-popSub
meanBehDay$individ_id<-birds[j]

meanActMonthly<-rbind(meanActMonthly, meanBehMonth)
meanActDaily<-rbind(meanActDaily, meanBehDay)

}

meanActSpecies<-rbind(meanActSpecies, meanActMonthly)
meanActSpecies_day<-rbind(meanActSpecies_day, meanActDaily)

# Save output
speciesSave<-gsub(" ", "", speciesSub)
speciesSave<-gsub("-", "", speciesSave)
speciesSave<-gsub("ü", "u", speciesSave)
speciesSave<-gsub("'", "", speciesSave)
speciesSave<-gsub("`", "", speciesSave)

# Make a special case for a population name which is cray cray
startName<-meanActMonthly$colony[1]

if (startName =="Rodolph") {
popSave<-"rodolphisland"
} else {
popSave<-gsub("_", "", meanActMonthly$colony[1])
}
popSave<-gsub("ö", "o", popSave)
popSave<-gsub("ø", "o", popSave)
popSave<-gsub("Å", "Aa", popSave)
popSave<-gsub("æ", "ae", popSave)
popSave<-gsub(" ", "", popSave)
popSave<-gsub("-", "", popSave)
popSave<-gsub(".,", "", popSave)
popSave<-gsub(",", "", popSave)
popSave<-gsub("/", "", popSave)
popSave<-gsub("\\)", "", popSave)
popSave<-gsub("\\(", "", popSave)
popSave<-gsub("–", "", popSave)
popSave<-gsub("ý", "", popSave)
popSave<-gsub("ð", "", popSave)
popSave<-gsub("á", "", popSave)
popSave<-gsub("ó", "", popSave)
popSave<-gsub(":", "", popSave)
popSave<-gsub("í", "", popSave)
popSave<-gsub("`", "", popSave)
popSave<-gsub("ú", "", popSave)
popSave<-gsub("é", "", popSave)
popSave<-gsub("'", "", popSave)
popSave<-tolower(popSave)
print(popSave)

write.csv(meanActMonthly, file=paste0("tmp3/", speciesSave, "_", popSave, "_actbudgets_monthly.csv")) 
write.csv(meanActDaily, file=paste0("tmp3/", speciesSave, "_", popSave, "_actbudgets_daily.csv")) 

}

# Make a plot to summarize sample size per month by population #

pdf(paste0("results/figures/speciesMaps/energy/datadistribution_", speciesSave, ".pdf"))

monthlySampleSize<-meanActSpecies %>%
dplyr::group_by(species, colony) %>%
dplyr::count(month) %>%
dplyr::mutate(sampleSize=1) %>%
ggplot() +
geom_col(aes(x=factor(month), y=sampleSize, fill=colony)) +
theme_bw() +
ylab("Has monthly data") +
xlab("Month")

grid.arrange(monthlySampleSize, nrow=2)

# Make a plot to summarize which dates are available for incomplete months (possible to move forwards with daily approach... #

monthsAll<-data.frame(month=rep(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), each=31))
monthsAll$day<-rep(seq(1, 31, 1), 12)

incompleteColonies<-meanActSpecies %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(monthsUnique=n_distinct(month)) %>%
dplyr::filter(monthsUnique<12)

# Figure out which months are missing

for (i in 1:nrow(incompleteColonies)) {

# Subset to colony i
incompleteColonySub<-incompleteColonies[i,]

# subset main dataset to this
incompleteData<-meanActSpecies %>%
dplyr::filter(colony==incompleteColonies$colony[1]) 

# Match with month data
monthsIncomplete<-monthsAll %>%
dplyr::anti_join(incompleteData, by=c("month"))

# Subset daily dataset
missingDays<-meanActSpecies_day %>%
dplyr::filter(colony==incompleteColonies$colony[1]) %>%
dplyr::filter(month %in% c(monthsIncomplete$month)) %>%
dplyr::group_by(species, colony, month) %>%
dplyr::count(day) %>%
dplyr::full_join(monthsAll ,by=c("month", "day")) %>%
ungroup() %>%
fill(species) %>%
fill(colony)

missingDaysPlot<-missingDays %>%
dplyr::filter(month %in% c(monthsIncomplete$month)) %>%
ggplot(aes(x=day, y=n)) +
geom_col() +
facet_wrap(~species + colony + month) +
ylab("Sample size (track-years)") +
xlab("") +
scale_x_continuous(breaks=seq(1, 31, 1)) +
theme_bw()

grid.arrange(missingDaysPlot, nrow=2)

}

dev.off()

# Save metadata as main output file
output_file1 <- args[3]
write.csv(meanActSpecies, file=output_file1, row.names=FALSE)



