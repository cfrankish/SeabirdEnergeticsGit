library(dplyr)
library(ncdf4)
library(raster)
library(sp)
library(sf)
library(lubridate)
library(data.table)
library(tidyr)
library(ggplot2)
library(rnaturalearthdata)
library(rnaturalearth)
library(lme4)
library(mgcv)

### Step 1: Show an example fulmar population map ###

# Find matching species raster #
rasterList<-data.frame(files=list.files("data/popdata_raw/", full.names=TRUE), species=c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot"))
rasterSub<-subset(rasterList, species=="Northern fulmar") # Subset to correct species
nc<-nc_open(rasterSub$files[1]) # open the file

# Create metadata of this file so I can easily pull the correct raster from it
modelcolonies <- (ncvar_get(nc,"SmcolName")) # Find list of populations
colonies <- (ncvar_get(nc,"colonyName")) # Find list of ALL colonies
colcode <- (ncvar_get(nc,"SmcolCode")) # Colony code
meta <- data.frame(modelcolonies, colonies, colcode) # Make some metadata so easier to understand raster structure

# Pull an example one - Faroe Islands in December #
iMth<-12
icols<-which(ncvar_get(nc, "colonyName") %in% "Fareoes") # Determine where colony info is located
if (length(icols)>1) {icols<-icols[1]}

# Add number of birds & error

rast.mean <- raster(nc$filename, varname="PredictedAbundanceMean", band = icols, level = iMth) # Mean density
rast.ci.high<- raster(nc$filename, varname="PredictedAbundance95%CIHigh", band = icols, level = iMth) # 95% CI
rast.ci.low<- raster(nc$filename, varname="PredictedAbundance95%CILow", band = icols, level = iMth) # 95% CI

# Make a plot #
rast.mean.df<-as.data.frame(rast.mean, xy=TRUE)

popMap<-ggplot() +
geom_tile(data=rast.mean.df, aes(x=x, y=y, fill=Mean.predicted.abundance)) +
scale_fill_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  )) +
  theme_classic() +
  theme(legend.position="bottom") +
  ggtitle("Faroe Islands - December") 
  
### Step 2: Find all bird locations ###

idCatalogue<-read.csv("results/tables/main/table1_idcatalogue.csv")

# Now we loop through all populations for a given species #
budgetLox<-list.files("tmp3/", full.names=TRUE)
budgetLox_species<-budgetLox[grep("Northernfulmar", budgetLox)]
budgetLox_species_daily<-budgetLox_species[grep("daily", budgetLox_species)]

# Now we subset idCatalogue to the correct species & colony
actBudgets<-fread(budgetLox_species_daily[5])
idSelect<-subset(idCatalogue, species==actBudgets$species[1] & colony==actBudgets$colony[1])
birds<-unique(idSelect$individ_id)

# Now we determine where the daily files are located 
actfiles<-list.files("./tmp/", full.names=TRUE)
dailyfiles<-actfiles[grepl("Day.csv", actfiles)]

# Now we loop through all individuals, open their files & summarize monthly locations & behavioral budgets

ActTemporal<-list() # List to save temporal variation in behaviour for a given month
meanActMonthly<-list() # List to save daily behaviour in

# Assign some variables
monthChoice<-12
Behaviour<-"Forage"
input_file1<-"Northernfulmar"
i<-1
pop<-5
overall.iterations<-5
iterationsRandom<-sample(1:100, overall.iterations, replace=FALSE) # Determine randomly which iterations I will draw from 

for (j in 1:length(birds)) {

# Print status update #
print(paste0("Iteration ", i, ", pop ", pop, "/", length(budgetLox_species_daily), ", bird ", j, "/", length(birds)))

# Subset to bird j
ID<-birds[j]
ID<-gsub("-", "_", ID) 
ID<-gsub("Hornøya", "Hornoya", ID)

# Subset correct file and open
birdSub<-dailyfiles[grepl(ID, dailyfiles)]

# Temporary fix for the few troublesome birds #
if(length(birdSub)<1) {
next}

birdCsv<-fread(birdSub)

# Subset to max number of repetitions to save memory
activitySubSub<-subset(birdCsv, rep==iterationsRandom[i])

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

# Subset dataset to the month of choice
meanBehMonths_sub<-subset(meanBehMonths, month==monthChoice)
meanBehMonths_sub$variable<-meanBehMonths_sub[[paste0("t", Behaviour)]]/24
meanBehMonths_save<-meanBehMonths_sub %>% # Save just the info that matters so I can make a temporal plot
 dplyr::select(individ_id, date, variable)
ActTemporal<-rbind(ActTemporal, meanBehMonths_save)

# Turn behaviors into prop of day doing something
meanBeh<-meanBehMonths_sub %>%
ungroup() %>%
dplyr::mutate(propForage=tForage/24, propActive=tActive/24, propRestWater=tRestWater/24, propLand=tLand/24, propFlight=tFlight/24) %>%
dplyr::select(species, colony, rep, month, date, individ_id, mean.lon, mean.lat, propForage, propActive, propRestWater, propLand, propFlight)

# Save results #
meanActMonthly<-rbind(meanActMonthly, meanBeh)

} # End of pop loop

coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

rast.mean.df$newVal<-ifelse(is.na(rast.mean.df$Mean.predicted.abundance), NA, 1)

extentPop<-filter(rast.mean.df, Mean.predicted.abundance>0)

popLox<-ggplot() +
geom_tile(data=rast.mean.df, aes(x=x, y=y, fill=factor(newVal)), alpha=0.2) +
scale_fill_manual(values = c("1" = "lightgrey"), na.value = "white") +
geom_tile(data=filter(rast.mean.df, Mean.predicted.abundance>0), aes(x=x, y=y), fill="grey") +
geom_point(data=meanActMonthly, aes(x=mean.lon, y=mean.lat, color=propForage)) +
scale_color_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  )) +
  theme_classic() +
  theme(legend.position="bottom") +
  ggtitle("Faroe Islands - December") +
  labs(fill="", x="", y="") +
  guides(fill = "none") 
 # coord_sf(xlim=c(min(extentPop$x)-2, max(extentPop$x) + 2), ylim=c(min(extentPop$y)-5, max(extentPop$y) + 2))
  
 
### Step 3a: Show some interpolations using gam ###

# Choose iterations at random 
iterationsRandom<-sample(1:100, overall.iterations, replace=FALSE) # Determine randomly which iterations I will draw from 

# Make a list to save all populations in
ActMonthlyAll<-list() # All locations & budgets
ActTemporalAll<-list() # Summary of budgets by day

for (pop in 1:length(budgetLox_species_daily)) {

# Now we subset idCatalogue to the correct species & colony
actBudgets<-fread(budgetLox_species_daily[pop])
idSelect<-subset(idCatalogue, species==actBudgets$species[1] & colony==actBudgets$colony[1])
birds<-unique(idSelect$individ_id)

# Now we determine where the daily files are located 
actfiles<-list.files("./tmp/", full.names=TRUE)
dailyfiles<-actfiles[grepl("Day.csv", actfiles)]

# Now we loop through all individuals, open their files & summarize monthly locations & behavioral budgets

ActTemporal<-list() # List to save temporal variation in behaviour for a given month
meanActMonthly<-list() # List to save daily behaviour in

for (j in 1:length(birds)) {

# Print status update #
print(paste0("Iteration ", i, ", pop ", pop, "/", length(budgetLox_species_daily), ", bird ", j, "/", length(birds)))

# Subset to bird j
ID<-birds[j]
ID<-gsub("-", "_", ID) 
ID<-gsub("Hornøya", "Hornoya", ID)

# Subset correct file and open
birdSub<-dailyfiles[grepl(ID, dailyfiles)]

# Temporary fix for the few troublesome birds #
if(length(birdSub)<1) {
next}

birdCsv<-fread(birdSub)

# Subset to max number of repetitions to save memory
activitySubSub<-subset(birdCsv, rep==iterationsRandom[i])

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

# Subset dataset to the month of choice
meanBehMonths_sub<-subset(meanBehMonths, month==monthChoice)
meanBehMonths_sub$variable<-meanBehMonths_sub[[paste0("t", Behaviour)]]/24
meanBehMonths_save<-meanBehMonths_sub %>% # Save just the info that matters so I can make a temporal plot
 dplyr::select(individ_id, date, variable)
ActTemporal<-rbind(ActTemporal, meanBehMonths_save)

# Turn behaviors into prop of day doing something
meanBeh<-meanBehMonths_sub %>%
ungroup() %>%
dplyr::mutate(propForage=tForage/24, propActive=tActive/24, propRestWater=tRestWater/24, propLand=tLand/24, propFlight=tFlight/24) %>%
dplyr::select(species, colony, rep, month, date, individ_id, mean.lon, mean.lat, propForage, propActive, propRestWater, propLand, propFlight)

# Save results #
meanActMonthly<-rbind(meanActMonthly, meanBeh)

} # End of pop loop

### Save oher data frames ###
ActMonthlyAll<-rbind(ActMonthlyAll, meanActMonthly) # This list contains all individuals, all pops

} # end of all pops loop 

###### Build a GAM to predict behavior accross whole North Atlantic ######

# Now we prepare a map for griddig behaviors (use one of Per's at random) #
nc<-nc_open("data/popdata_raw/SEATRACK_Abundance_Model_Rissa_tridactyla_Ver_3_1.nc")
rast.gr <- raster(nc$filename, varname="PredictedAbundanceMean", band = 1, level = 1)
crs(rast.gr)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Non-projected
xyCoords<-as.data.frame(rast.gr, xy=TRUE)
xyCoords<-xyCoords %>%
rename(coords.x1=x, coords.x2=y) %>%
dplyr::mutate(individ_id=ActMonthlyAll$individ_id[1]) %>%
dplyr::mutate(colony=ActMonthlyAll$colony[1])

# Define response variable
responseVar<-paste0("prop", Behaviour) 
responseVar2<-paste0(Behaviour) 

gam_formula2 <- as.formula(
  paste0(responseVar2, " ~ s(coords.x1, coords.x2) + s(individ_id, bs='re')")
)

gam_formula3 <- as.formula(
  paste0(responseVar2, " ~  s(coords.x1, coords.x2, by=colony) + colony +  s(individ_id, bs='re')")
)

# Add a column to the dataset with weights of different populations (so that pops with less individuals are weighted more?)
ActMonthlyAllTrans<-ActMonthlyAll
ActMonthlyFit<-ActMonthlyAllTrans %>%
ungroup() %>%
dplyr::group_by(colony) %>%
dplyr::mutate(birds=n_distinct(individ_id), popWeights=1/birds) %>%
ungroup() %>%
dplyr::mutate(popWeightsScaled=popWeights/mean(popWeights)) %>%
dplyr::rename(coords.x1=mean.lon, coords.x2=mean.lat)

# Change values so they cannot be close to 0 or 1 (necessary for the family 'betar'
ActMonthlyFit[[responseVar]]<-ifelse(ActMonthlyFit[[responseVar]]==1, 0.99, ActMonthlyFit[[responseVar]]) # can't be 1 if beta reg family
ActMonthlyFitSub <- ActMonthlyFit[ ActMonthlyFit[[responseVar]] > 0, ] # can't be 0 either

# Factorize individ id and colony
ActMonthlyFitSub $individ_id<-factor(ActMonthlyFitSub $individ_id)
ActMonthlyFitSub $colony<-factor(ActMonthlyFitSub $colony)
ActMonthlyFitSub[[Behaviour]]<-ActMonthlyFitSub[[responseVar]]*24 # So that I can just use a gaussian gam #

# Fit two gams
gam2 <- gam(gam_formula2, family = gaussian(), data=ActMonthlyFitSub)
gam3 <- gam(gam_formula3, family = gaussian(), data=ActMonthlyFitSub)

# Make prediction outputs 

#### Model 2 ####

xyCoords<-as.data.frame(rast.gr, xy=TRUE)
xyCoords<-xyCoords %>%
rename(coords.x1=x, coords.x2=y) %>%
dplyr::mutate(individ_id=ActMonthlyAll$individ_id[1]) %>%
dplyr::mutate(colony=ActMonthlyAll$colony[1])

Predictions_Outputs <- predict(gam2, newdata=xyCoords, exclude = c("s(individ_id)"), type="response")
xyCoords$predictions<-Predictions_Outputs
coordinates(xyCoords) <- ~ coords.x1 + coords.x2
statPoints<-xyCoords
predictions_space <- rasterize(statPoints, rast.gr, field="predictions", fun=mean) # Species-level December map # 

#### Model 3 ####

# Select random 3 colonies for plotting #
coloniesAll<-unique(ActMonthlyFitSub$colony)
coloniesRandom<-"Faroe Islands"

xyCoords<-as.data.frame(rast.gr, xy=TRUE)
xyCoords<-xyCoords %>%
rename(coords.x1=x, coords.x2=y) %>%
dplyr::mutate(individ_id=ActMonthlyAll$individ_id[1]) %>%
dplyr::mutate(colony=coloniesRandom)

Predictions_Outputs <- predict(gam3, newdata=xyCoords, exclude = c("s(individ_id)"), type="response")
xyCoords$predictions<-Predictions_Outputs
coordinates(xyCoords) <- ~ coords.x1 + coords.x2
statPoints<-xyCoords
predictions_colony <- rasterize(statPoints, rast.gr, field="predictions", fun=mean)

# Then I will crop by non-zero values in the density map #
mask_r <- rast.mean
mask_r[mask_r == 0] <- NA

r1_masked <- mask(predictions_colony, mask_r)
r1_masked.df<-as.data.frame(r1_masked, xy=TRUE)

popLox_interpol<-ggplot() +
  geom_tile(
    data = filter(rast.mean.df, newVal == 1),
    aes(x = x, y = y),
    fill = "lightgrey",
    alpha = 0.2
  ) +
  geom_tile(
    data = filter(rast.mean.df, Mean.predicted.abundance > 0),
    aes(x = x, y = y),
    fill = "grey"
  ) +
  geom_tile(
    data = filter(r1_masked.df, layer > 0),
    aes(x = x, y = y, fill = layer/24)
  ) +
  scale_fill_gradientn(
    colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Faroe Islands - December") +
  labs(fill="Prop forage")
 
### Step 3b: Show some gridding ###

# grid size 
gridSize<-5 # degrees

# Now we prepare a map for griddig behaviors (use one of Per's at random) #
nc<-nc_open("data/popdata_raw/SEATRACK_Abundance_Model_Rissa_tridactyla_Ver_3_1.nc")
rast.gr <- raster(nc$filename, varname="PredictedAbundanceMean", band = 1, level = 1)
crs(rast.gr)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Now we loop through all individuals, open their files & summarize monthly locations & behavioral budgets
griddedBirds<-list()

for (j in 1:length(birds)) {

# Print status update #
print(paste0("Iteration ", i, ", pop ", pop, "/", length(budgetLox_species_daily), ", bird ", j, "/", length(birds)))

# Subset to bird j
ID<-birds[j]
ID<-gsub("-", "_", ID) 
ID<-gsub("Hornøya", "Hornoya", ID)

# Subset correct file and open
birdSub<-dailyfiles[grepl(ID, dailyfiles)]

# Temporary fix for the few troublesome birds #
if(length(birdSub)<1) {
next}

birdCsv<-fread(birdSub)

# Subset to max number of repetitions to save memory
activitySubSub<-subset(birdCsv, rep==iterationsRandom[i])

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

# Subset dataset to the month of choice
meanBehMonths_sub<-subset(meanBehMonths, month==monthChoice)
meanBehMonths_sub$variable<-meanBehMonths_sub[[paste0("t", Behaviour)]]/24
meanBehMonths_save<-meanBehMonths_sub %>% # Save just the info that matters so I can make a temporal plot
 dplyr::select(individ_id, date, variable)

# Turn behaviors into prop of day doing something
meanBeh<-meanBehMonths_sub %>%
ungroup() %>%
dplyr::mutate(propForage=tForage/24, propActive=tActive/24, propRestWater=tRestWater/24, propLand=tLand/24, propFlight=tFlight/24) %>%
dplyr::select(species, colony, rep, month, date, individ_id, mean.lon, mean.lat, propForage, propActive, propRestWater, propLand, propFlight) %>%
dplyr::mutate(coords.x1=mean.lon, coords.x2=mean.lat)

if (nrow(meanBeh)<1) {
next}

# Change projections of coordinates
#monthlyidw<-meanBeh
#coordinates(monthlyidw)<-~mean.lon + mean.lat
#proj4string(monthlyidw)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#ActMonthAllTrans1<-as.data.frame(spTransform(monthlyidw, projection_NA))
responseVar2<-paste0("prop", Behaviour)
ActMonthAllTrans1<-meanBeh
ActMonthAllTrans1$Beh<-ActMonthAllTrans1[[responseVar2]]*24

# Conduct IDW on behavior of choice & save result in list 
rast.coarse<-aggregate(rast.gr, fact=20)
rastDF<-as.data.frame(rast.coarse, xy=TRUE)
rastDF$propBeh<-0
rastDF$birds<-0
rastDF$index<-1:nrow(rastDF)

# Now we subset this raster, otherwise it will take too long to fill it up
rastDF_sub<-subset(rastDF, x> min(ActMonthAllTrans1$coords.x1) - gridSize*2 & x < max(ActMonthAllTrans1$coords.x1) + gridSize*2 & y> min(ActMonthAllTrans1$coords.x2) - gridSize*2 & y < max(ActMonthAllTrans1$coords.x2) + gridSize*2)

# Grid behaviour first #
for(cell in 1:nrow(rastDF_sub)) {

# Below we sum the total time in that behaviour by the number of days
rows<-subset(ActMonthAllTrans1, coords.x1>=rastDF_sub$x[cell] & coords.x1 < rastDF_sub$x[cell] + gridSize & coords.x2 >= rastDF_sub$y[cell] & coords.x2 < rastDF_sub$y[cell] + gridSize)

# If there are positions in that cell that calculate proportion of total time in that cell doing a specific behaviour
if (nrow(rows)>0) {
rastDF_sub$propBeh[cell]<-sum(rows$Beh)/(nrow(rows)*24) # Divided by time spent in that cell
rastDF_sub$birds[cell]<-1

}

}

# Now we add other information & subset to rows where a bird was present #
rastDF_sub$species<-input_file1
rastDF_sub$colony<-meanBeh$colony[1]
rastDF_sub$month<-monthChoice
rastDF_sub$behaviour<-Behaviour
rastDF_sub$individ_id<-meanBeh$individ_id[1]
rastDF_bird<-subset(rastDF_sub, birds>0)

griddedBirds<-rbind(griddedBirds, rastDF_bird)

} # End of individual loop

griddedPop<-griddedBirds %>%
ungroup() %>%
dplyr::group_by(x, y) %>%
dplyr::summarise(meanPropForage=mean(propBeh), totBirds=n_distinct(individ_id))

popLox<-ggplot() +
  geom_tile(data = filter(rast.mean.df, is.na(newVal)), aes(x = x, y = y), fill = "white", alpha = 0.2) +
  geom_tile(data = dplyr::filter(rast.mean.df, newVal== 1),
            aes(x = x, y = y), fill = "grey", alpha=0.2) +
  geom_tile(data=filter(rast.mean.df, Mean.predicted.abundance>1), aes(x=x, y=y), fill="grey") +
  geom_tile(data = griddedPop, aes(x = x, y = y, fill = meanPropForage), alpha = 0.5) +
  scale_fill_gradientn(colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Faroe Islands - December") +
  labs(fill = "", x = "", y = "")
  

### Step 4: Show population smooth ###

allPops<-readRDS("tmp4/Northernfulmar_Forage_m12_rasters_gam3_v2.rds")

# Stack together all the iterations for a given colony #
restacked <- lapply(seq_len(nlayers(allPops[[1]])), function(i) {
  stack(lapply(allPops, raster::subset, i))
})

faroes<-overlay(restacked[[5]], fun="mean")

# Then I will crop by non-zero values in the density map #
mask_r <- rast.mean
mask_r[mask_r == 0] <- NA

faroes_dis<-resample(faroes, mask_r) # make same extent etc.
r1_masked <- mask(faroes_dis, mask_r)
r1_masked.df<-as.data.frame(r1_masked, xy=TRUE)

coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

popLox_interpol2<-ggplot() +
  geom_tile(
    data = filter(rast.mean.df, newVal == 1),
    aes(x = x, y = y),
    fill = "white", alpha=0.2
  ) +
  geom_tile(
    data = filter(rast.mean.df, Mean.predicted.abundance > 0),
    aes(x = x, y = y),
    fill = "grey"
  ) +
  geom_tile(
    data = filter(r1_masked.df, layer > 0),
    aes(x = x, y = y, fill = layer)
  ) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
  scale_fill_gradientn(
    colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c") 
  ) +
  coord_sf(xlim=c(min(rast.mean.df$x), max(rast.mean.df$x)), ylim=c(min(rast.mean.df$y), max(rast.mean.df$y))) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Faroe Islands - December") +
  labs(fill="Prop forage")
  
### Step 5: Show multiple types of maps ###

rescale01 <- function(x, na.rm = TRUE) {
  rng <- range(x, na.rm = na.rm)
  (x - rng[1]) / (rng[2] - rng[1])
}

# Foraging layer - re-scale between 0 and 1
r1_masked.df_sub<-subset(r1_masked.df, layer>0)
r1_masked.df_sub$scaled<-rescale01(r1_masked.df_sub$layer)
r1_masked.df_sub$metric<-"PropForaging"
r1_masked.df_sub<-r1_masked.df_sub[c("x", "y", "scaled", "metric")]

# Population layer -re.scale 
mask_r.df<-as.data.frame(mask_r, xy=TRUE)
mask_r_sub<-subset(mask_r.df, Mean.predicted.abundance>0)
mask_r_sub$scaled<-rescale01(mask_r_sub$Mean.predicted.abundance)
mask_r_sub$metric<-"NoBirds"
mask_r_sub<-mask_r_sub[c("x", "y", "scaled", "metric")]

# Multiply the two? #
mask_r_sub_join<-mask_r_sub %>%
dplyr::left_join(r1_masked.df_sub, by=c("x", "y")) %>%
dplyr::mutate(scaled=scaled.x*scaled.y) %>%
dplyr::select(x, y, scaled) %>%
dplyr::mutate(metric="Forage x birds")

# Energy map #
energymap<-read.csv("tmp5/Northernfulmar_fareoes_energyMap_v2.csv")
energySub<-energymap %>%
dplyr::filter(month==12) %>%
dplyr::filter(NoBirds_mean>0) %>%
dplyr::mutate(scaled=rescale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric="EnergyDemand") %>%
dplyr::select(x, y, scaled, metric)

# Join all layers #
allLayers<-rbind(r1_masked.df_sub, mask_r_sub, mask_r_sub_join)
allLayers<-allLayers %>%
dplyr::mutate(metric=factor(metric, levels=c("NoBirds", "PropForaging", "Forage x birds")))

tiles<-ggplot() +
  geom_tile(
    data = allLayers,
    aes(x = x, y = y, fill = scaled)
  ) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
  scale_fill_gradientn(
    colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c") 
  ) +
  coord_sf(xlim=c(min(rast.mean.df$x), max(rast.mean.df$x)), ylim=c(min(rast.mean.df$y), max(rast.mean.df$y))) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Faroe Islands - December") +
  labs(fill="Scaled 0-1") +
  facet_wrap(~metric, nrow=1)
  
  pdf("tiles.pdf", width=10, height=7)
  plot(tiles)
  dev.off()
  
# Highlight top 5% of data #

thresh <- quantile(mask_r_sub_join$scaled, 0.99, na.rm = TRUE)
mask_r_sub_join$top5 <- mask_r_sub_join$scaled >= thresh

thresh2 <- quantile(mask_r_sub$scaled, 0.99, na.rm = TRUE)
mask_r_sub$top5 <- mask_r_sub$scaled >= thresh

allLayers2<-rbind(mask_r_sub, mask_r_sub_join)
allLayers2<-allLayers2 %>%
dplyr::mutate(metric=factor(metric, levels=c("NoBirds", "Forage x birds")))

tiles2<-ggplot() +
  geom_tile(
    data = allLayers2,
    aes(x = x, y = y, fill = top5)
  ) +
  scale_fill_manual(values=c("darkblue", "yellow")) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
  coord_sf(xlim=c(min(rast.mean.df$x), max(rast.mean.df$x)), ylim=c(min(rast.mean.df$y), max(rast.mean.df$y))) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Faroe Islands - December") +
  labs(fill="top 1% data?") +
  facet_wrap(~metric, nrow=1)
  
  pdf("tiles2.pdf", width=10, height=7)
  plot(tiles2)
  dev.off()
  
### Step 6: Show multiple types of maps - species-level ###

# Population layer -re.scale 

# Extract all population maps for December, stack and sum 
speciesDensity<-list()
for (i in 1:nrow(meta)) {
rast.mean <- raster(nc$filename, varname="PredictedAbundanceMean", band = i, level = iMth) # # Add number of birds & error
speciesDensity[[i]]<-rast.mean # save
}
speciesDensityAll<-stack(speciesDensity)
speciesDensitySum<-overlay(speciesDensityAll, fun="sum")
speciesDensityDf<-as.data.frame(speciesDensitySum, xy=TRUE)
colnames(speciesDensityDf)<-c("x", "y", "NoBirds")
speciesDensityDf_sub<-subset(speciesDensityDf, NoBirds>0)
speciesDensityDf_sub$scaled<-rescale01(speciesDensityDf_sub$NoBirds)
speciesDensityDf_sub$metric<-"NoBirds"
speciesDensityDf_sub<-speciesDensityDf_sub[c("x", "y", "scaled", "metric")]

# Foraging layer - re-scale between 0 and 1
foraging_species<-predictions_space
mask_r <- speciesDensitySum
mask_r[mask_r == 0] <- NA
foraging_species_resample<-resample(faroes, foraging_species) # make same extent etc.
r1_masked <- mask(foraging_species_resample, mask_r)
r1_masked.df<-as.data.frame(r1_masked, xy=TRUE)
r1_masked.df_sub<-subset(r1_masked.df, layer>0)
r1_masked.df_sub$scaled<-rescale01(r1_masked.df_sub$layer)
r1_masked.df_sub$metric<-"PropForaging"
r1_masked.df_sub<-r1_masked.df_sub[c("x", "y", "scaled", "metric")]

# Extract energy density maps #
allmaps<-list.files("tmp5/", full.names=TRUE)
allmaps_fulmars<-allmaps[grep("Northernfulmar", allmaps)] 
allmaps_fulmars_spatial<-allmaps_fulmars[grep("energyMap_v2", allmaps_fulmars)] 
allenergy<-list()
for(i in 1:length(allmaps_fulmars_spatial)) {
print(paste0("Opening map", i, length(allmaps_fulmars_spatial)))
mapi<-fread(allmaps_fulmars_spatial[i])
mapi_dec<-subset(mapi, month==12)
allenergy<-rbind(allenergy, mapi_dec)
}
allenergySum<-allenergy %>%
dplyr::group_by(x, y) %>%
dplyr::summarise(totenergy=sum(energyPopkJ_mean), totBirds=sum(NoBirds_mean)) %>%
dplyr::filter(totBirds>0)
allenergySum$scaled<-rescale01(allenergySum$totenergy)
allenergySum$metric<-"energyDemand"
allenergySum<-allenergySum[c("x", "y", "scaled", "metric")]

# Multiply the two? #
mask_r_sub_join<-speciesDensityDf_sub %>%
dplyr::left_join(r1_masked.df_sub, by=c("x", "y")) %>%
dplyr::mutate(scaled=scaled.x*scaled.y) %>%
dplyr::select(x, y, scaled) %>%
dplyr::mutate(metric="Forage x birds")

# Join all layers #
allLayers<-rbind(r1_masked.df_sub, speciesDensityDf_sub, mask_r_sub_join)
allLayers<-allLayers %>%
dplyr::mutate(metric=factor(metric, levels=c("NoBirds", "PropForaging", "Forage x birds")))

tiles3<-ggplot() +
  geom_tile(
    data = allLayers,
    aes(x = x, y = y, fill = scaled)
  ) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
  scale_fill_gradientn(
    colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c") 
  ) +
  coord_sf(xlim=c(min(rast.mean.df$x), max(rast.mean.df$x)), ylim=c(min(rast.mean.df$y), max(rast.mean.df$y))) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Faroe Islands - December") +
  labs(fill="Scaled 0-1") +
  facet_wrap(~metric, nrow=1)
  
  pdf("tiles3.pdf", width=10, height=7)
  plot(tiles3)
  dev.off()
  
# Highlight top 5% of data #

thresh <- quantile(mask_r_sub_join$scaled, 0.99, na.rm = TRUE)
mask_r_sub_join$top5 <- mask_r_sub_join$scaled >= thresh
mask_r_sub_join$top5<-ifelse(mask_r_sub_join$top5==TRUE, "Forage hotspot", mask_r_sub_join$top5)

thresh2 <- quantile(speciesDensityDf_sub$scaled, 0.99, na.rm = TRUE)
speciesDensityDf_sub$top5 <- speciesDensityDf_sub$scaled >= thresh
speciesDensityDf_sub$top5<-ifelse(speciesDensityDf_sub$top5==TRUE, "Bird hotspot", speciesDensityDf_sub$top5)

thresh3 <- quantile(allenergySum$scaled, 0.99, na.rm = TRUE)
allenergySum$top5 <- allenergySum$scaled >= thresh
allenergySum$top5<-ifelse(allenergySum$top5==TRUE, "Energy demand hotspot", allenergySum$top5)

allLayers2<-rbind(speciesDensityDf_sub, mask_r_sub_join, allenergySum)
allLayers2<-allLayers2 %>%
dplyr::mutate(metric=factor(metric, levels=c("NoBirds", "Forage x birds")))

tiles4<-ggplot() +
  geom_tile(
    data = filter(allLayers2, top5==FALSE), 
    aes(x = x, y = y, fill = "Species Distribution")
  ) +
  geom_tile(
    data = filter(allLayers2, top5=="Bird hotspot"), 
    aes(x = x, y = y, fill = "Bird hotspot")
  ) +
  geom_tile(
    data = filter(allLayers2, top5=="Forage hotspot"), 
    aes(x = x, y = y, fill = "Forage hotspot")
  ) +
  geom_tile(
    data = filter(allLayers2, top5=="Energy demand hotspot"), 
    aes(x = x, y = y, fill = "Energy hotspot"), alpha=0.9
  ) +
  scale_fill_manual(values=c("#ffffbf", "#fdae61", "#d7191c", "#2c7bb6", "#abd9e9")) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
  coord_sf(xlim=c(min(rast.mean.df$x), max(rast.mean.df$x)), ylim=c(min(rast.mean.df$y), max(rast.mean.df$y))) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("Species-level :December") +
  labs(fill="") 
  
  pdf("tiles4.pdf", width=10, height=7)
  plot(tiles4)
  dev.off()
  
  "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"


