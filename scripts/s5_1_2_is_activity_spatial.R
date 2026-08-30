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
library(lmerTest)
library(mgcv)
library(gstat)
library(gridExtra)

#### Step 0: setting up basic conditions ####

# Projections
# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Set-up number of iterations...
overall.iterations<-3 # how many times this is calculated per individual
print(paste0("Determining activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1] # This will read in a species-name 
input_file2<- args[2] # This is list of all bird Ids (id catalogue) - input_file2<-"results/tables/main/table1_idcatalogue.csv"  
Behaviour<-args[3] # This will be a behaviour like 'forage'
print(paste0("Behavior=", Behaviour))
monthChoice<-args[4] # This will be a month like 3
print(paste0("Month:", monthChoice))

# open input file 2
idCatalogue<-read.csv(input_file2) # idCatalogue<-read.csv("results/tables/main/table1_idcatalogue.csv")

# Now we prepare a map for griddig behaviors (use one of Per's at random) #
nc<-nc_open("data/popdata_raw/SEATRACK_Abundance_Model_Rissa_tridactyla_Ver_3_1.nc")
rast.gr <- raster(nc$filename, varname="PredictedAbundanceMean", band = 1, level = 1)
crs(rast.gr)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
rast.trans<-raster::projectRaster(rast.gr, crs=projection_NA, res=200000, method="bilinear")
rast.trans[is.na(rast.trans)] <- 0

# Now we loop through all populations for a given species #
budgetLox<-list.files("tmp3/", full.names=TRUE)
budgetLox_species<-budgetLox[grep(input_file1, budgetLox)]
budgetLox_species_daily<-budgetLox_species[grep("daily", budgetLox_species)]

# Make lists to save interpolations in #
predictionsAllGam_model1<-list()
predictionsAllGam_model2<-list()
predictionsAllGam_model3<-list()
predictionSummaries<-list()

for (i in 1:overall.iterations) {

startIteration <- Sys.time()

# Choose iterations at random 
iterationsRandom<-sample(1:100, overall.iterations, replace=FALSE) # Determine randomly which iterations I will draw from 

# Make a list to save all populations in
ActMonthlyAll<-list() # All locations & budgets
ActTemporalAll<-list() # Summary of budgets by day

for (pop in 1:length(budgetLox_species_daily)) {

startPop <- Sys.time()

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

EndPop <- Sys.time()

print("Pop loop took:")
print(EndPop-startPop)

} # end of all pops loop 

###### Build a GAM to predict behavior accross whole North Atlantic ######

# Make a dataset to predicton (Per's pop maps)# projected
#xyCoords<-as.data.frame(rast.trans, xy=TRUE)
#xyCoords<-xyCoords %>%
#rename(coords.x1=x, coords.x2=y) %>%
#dplyr::mutate(individ_id=ActMonthlyAll$individ_id[1]) %>%
#dplyr::mutate(colony=ActMonthlyAll$colony[1])

# Non-projected
xyCoords<-as.data.frame(rast.gr, xy=TRUE)
xyCoords<-xyCoords %>%
rename(coords.x1=x, coords.x2=y) %>%
dplyr::mutate(individ_id=ActMonthlyAll$individ_id[1]) %>%
dplyr::mutate(colony=ActMonthlyAll$colony[1])

# Define response variable
responseVar<-paste0("prop", Behaviour) 
responseVar2<-paste0(Behaviour) 

# Define formulas to test

gam_formula1<-as.formula(
  paste0(responseVar2, " ~  s(individ_id, bs='re')")
)

gam_formula2 <- as.formula(
  paste0(responseVar2, " ~ s(coords.x1, coords.x2, bs='gp', k=30) + s(individ_id, bs='re')")
)

gam_formula3 <- as.formula(
  paste0(responseVar2, " ~  s(coords.x1, coords.x2, by=colony, bs='gp', k=30) + colony +  s(individ_id, bs='re')")
)

# Transform the coordinates of the dataset #
#ActMonthlyAll1<-ActMonthlyAll
#coordinates(ActMonthlyAll1)<-~mean.lon + mean.lat
#proj4string(ActMonthlyAll1)<-projection_84
#ActMonthlyAllTrans<-as.data.frame(spTransform(ActMonthlyAll1, projection_NA))

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

# Fit gam
startGam<-Sys.time()

print("Fitting gams...")

gam1 <- bam(gam_formula1, family = gaussian(), data=ActMonthlyFitSub, method = "fREML", discrete = TRUE)
gam2 <- bam(gam_formula2, family = gaussian(), data=ActMonthlyFitSub, method = "fREML", discrete = TRUE)
gam3 <- bam(gam_formula3, family = gaussian(), data=ActMonthlyFitSub, method = "fREML", discrete = TRUE)

#loovc_1<-conduct_loovc(model_formula=gam_formula1, data=ActMonthlyFitSub, responseVar, testType="Time") # Conduct loovc

# Store model results
AICModels<-AIC(gam1, gam2, gam3) # Conduct AIC

# Add results to the dataset (so I also have locations)
ActMonthlyFitSub$AIC1<-AICModels$AIC[1]
ActMonthlyFitSub$AIC2<-AICModels$AIC[2]
ActMonthlyFitSub$AIC3<-AICModels$AIC[3]
ActMonthlyFitSub$dev1<-summary(gam1)$dev.expl
ActMonthlyFitSub$dev2<-summary(gam2)$dev.expl
ActMonthlyFitSub$dev3<-summary(gam3)$dev.expl

# Save these results
predictionSummaries<-rbind(predictionSummaries, ActMonthlyFitSub)

# Make prediction outputs 

#### Model 1 ####

Predictions_Outputs <- predict(gam1, newdata=xyCoords, exclude = c("s(individ_id)"), type="response")
xyCoords$predictions<-Predictions_Outputs
coordinates(xyCoords) <- ~ coords.x1 + coords.x2
statPoints<-xyCoords
predictions_temporal <- rasterize(statPoints, rast.gr, field="predictions", fun=mean)

# save prediction raster
predictionsAllGam_model1[[i]]<-predictions_temporal

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
predictions_space <- rasterize(statPoints, rast.gr, field="predictions", fun=mean)

# save prediction raster
predictionsAllGam_model2[[i]]<-predictions_space

#### Model 3 ####

# Save in a small stack
colonyMaps<-list()

# Select random 3 colonies for plotting #
coloniesAll<-unique(ActMonthlyFitSub$colony)
coloniesRandom<-sample(coloniesAll, 3, replace=FALSE)

for (col in 1:length(coloniesRandom)) {

xyCoords<-as.data.frame(rast.gr, xy=TRUE)
xyCoords<-xyCoords %>%
rename(coords.x1=x, coords.x2=y) %>%
dplyr::mutate(individ_id=ActMonthlyAll$individ_id[1]) %>%
dplyr::mutate(colony=coloniesRandom[col])

Predictions_Outputs <- predict(gam3, newdata=xyCoords, exclude = c("s(individ_id)"), type="response")
xyCoords$predictions<-Predictions_Outputs
coordinates(xyCoords) <- ~ coords.x1 + coords.x2
statPoints<-xyCoords
predictions_colony <- rasterize(statPoints, rast.gr, field="predictions", fun=mean)

colonyMaps<-append(colonyMaps, predictions_colony)

}

# Stack these
colonyMapsStack<-stack(colonyMaps)

# give names so I know which colonies
names(colonyMapsStack)<-coloniesRandom

# save prediction raster
predictionsAllGam_model3[[i]]<-colonyMapsStack

endIteration <- Sys.time()

print("Iteration took:")
print(endIteration-startIteration)


} # end of iteration looop

#### Make some summary plots of stats ####

statsSum<-predictionSummaries %>%
dplyr::group_by(species, month) %>%
dplyr::summarise(AIC1_mean=mean(AIC1), AIC1_sd=sd(AIC1), AIC2_mean=mean(AIC2), AIC2_sd=sd(AIC2), AIC3_mean=mean(AIC3), AIC3_sd=sd(AIC3),
dev1_mean=mean(dev1), dev1_sd=sd(dev1), dev2_mean=mean(dev2), dev2_sd=sd(dev2), dev3_mean=mean(dev3), dev3_sd=sd(dev3))

stats1<-statsSum %>%
dplyr::select(species, month, AIC1_mean, AIC1_sd, dev1_mean, dev1_sd) %>%
rename(AIC_mean=AIC1_mean, AIC_sd=AIC1_sd, dev_mean=dev1_mean, dev_sd=dev1_sd) %>%
dplyr::mutate(modelType="Null")

stats2<-statsSum %>%
dplyr::select(species, month, AIC2_mean, AIC2_sd, dev2_mean, dev2_sd) %>%
rename(AIC_mean=AIC2_mean, AIC_sd=AIC2_sd, dev_mean=dev2_mean, dev_sd=dev2_sd) %>%
dplyr::mutate(modelType="Spatial")

stats3<-statsSum %>%
dplyr::select(species, month, AIC3_mean, AIC3_sd, dev3_mean, dev3_sd) %>%
rename(AIC_mean=AIC3_mean, AIC_sd=AIC3_sd, dev_mean=dev3_mean, dev_sd=dev3_sd) %>%
dplyr::mutate(modelType="Spatial_pop")

stats_all<-rbind(stats1, stats2, stats3)

stats1<-ggplot() +
geom_pointrange(data=stats_all, aes(x=modelType, y=AIC_mean, ymin=AIC_mean-AIC_sd, ymax=AIC_mean + AIC_sd)) +
ylab("AIC (mean +/- SD") +
xlab("Model type") +
theme_bw()+ 
ggtitle(paste(stats_all$species[1], monthChoice, responseVar2))

stats2<-ggplot() +
geom_pointrange(data=stats_all, aes(x=modelType, y=dev_mean, ymin=dev_mean-dev_sd, ymax=dev_mean + dev_sd)) +
ylab("dev expl (mean +/- SD") +
xlab("Model type") +
theme_bw()+ 
ggtitle(paste(stats_all$species[1], monthChoice, responseVar2))

#### Make some summary plots of rasters ####

coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# First we plot mean species-level raster & plot locations on top with propForage #
predictionsAllGam_model2_stack<-stack(predictionsAllGam_model2) # Stack all rasters
predictionsAllGam_model2_mean<-raster::overlay(predictionsAllGam_model2_stack, fun="mean") # Calculate mean accross iterations
predictionsAllGam_model2_df<-as.data.frame(predictionsAllGam_model2_mean, xy=TRUE) # Turn into data frame
predictionSummaries$Beh<-predictionSummaries[[responseVar2]] # Add a column with the behaviour we are modelling
locations<-predictionSummaries %>%
ungroup() %>%
dplyr::group_by(species, colony, month, date) %>%
dplyr::summarise(mean.lon=mean(coords.x1), mean.lat=mean(coords.x2), meanBeh=mean(Beh)) # Summarize locations & activity budgets for plotting

# example: points as sf
pts_sf <- st_as_sf(locations, coords = c("mean.lon", "mean.lat"), crs = projection_84)

# convex hull
hull <- st_convex_hull(st_union(pts_sf))
hull_84 <- st_transform(hull, crs = projection_84)

maxRaw<-max(locations$meanBeh)
maxPred<-max(predictionsAllGam_model2_df$layer)
maxPlot<-ifelse(maxRaw>maxPred, maxRaw, maxPred)
maxPlot<-ifelse(maxPlot>24, 24, maxPlot)

minRaw<-min(locations$meanBeh)
minPred<-min(predictionsAllGam_model2_df$layer)
minPlot<-ifelse(minRaw>minPred, minPred, minRaw)
minPlot<-ifelse(minPlot<0, 0, minPlot)

plot1<-ggplot() +
geom_tile(data=predictionsAllGam_model2_df, aes(x=x, y=y, fill=layer)) +
geom_sf(data = hull, fill = NA, color = "red", linewidth = 1) +
geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
coord_sf(crs=projection_84, xlim=c(min(locations$mean.lon), max(locations$mean.lon)), ylim=c(min(locations$mean.lat), max(locations$mean.lat))) +
scale_fill_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
scale_color_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
  labs(color=paste0(Behaviour), fill=paste0(Behaviour)) +
  ggtitle(paste0(input_file1, monthChoice, ": species-level pred")) +
  xlab("") +
  ylab("") 
  
plot2<-ggplot() +
#geom_tile(data=predictionsAllGam_model2_df, aes(x=x, y=y, fill=layer)) +
geom_point(data=locations, aes(x=mean.lon, y=mean.lat, color=meanBeh)) +
geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
scale_fill_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
  scale_color_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
  coord_sf(crs=projection_84, xlim=c(min(locations$mean.lon), max(locations$mean.lon)), ylim=c(min(locations$mean.lat), max(locations$mean.lat))) +
  labs(color=paste0(Behaviour), fill=paste0(Behaviour)) +
  ggtitle(paste0(input_file1, monthChoice, ": raw data")) +
  xlab("") +
  ylab("")
  
pdf(paste0("results/", input_file1, "_", Behaviour, "_", monthChoice, ".pdf"))
grid.arrange(stats1, stats2, plot1, plot2, nrow=2)
#dev.off()
  
### Now we summarize & compare some of the population-level predictions ###

predictionsAllGam_model3_stack<-predictionsAllGam_model3[[1]] # Stack all rasters

for (i in 1:nlayers(predictionsAllGam_model3_stack)) {

colonySub<-predictionsAllGam_model3_stack[[i]]
predictionsAllGam_model3_df<-as.data.frame(colonySub, xy=TRUE) # Turn into data frame
predictionsAllGam_model3_df$predictions<-predictionsAllGam_model3_df[,3]
predictionSummaries$Beh<-predictionSummaries[[responseVar2]] # Add a column with the behaviour we are modelling
colonyNameMatch<-gsub("\\.", " ", names(colonySub)) # Extract colony name from the raster so I can subset the locations
locations<-predictionSummaries %>%
dplyr::filter(colony==colonyNameMatch) %>%
ungroup() %>%
dplyr::group_by(species, colony, month, date) %>%
dplyr::summarise(mean.lon=mean(coords.x1), mean.lat=mean(coords.x2), meanBeh=mean(Beh)) # Summarize locations & activity budgets for plotting

# example: points as sf
pts_sf <- st_as_sf(locations, coords = c("mean.lon", "mean.lat"), crs = projection_84)

# convex hull
hull <- st_convex_hull(st_union(pts_sf))
hull_84 <- st_transform(hull, crs = projection_84)

# Determine range to plot 
maxRaw<-max(locations$meanBeh)
maxPred<-max(predictionsAllGam_model3_df$predictions)
maxPlot<-ifelse(maxRaw>maxPred, maxRaw, maxPred)
maxPlot<-ifelse(maxPlot>24, 24, maxPlot)

minRaw<-min(locations$meanBeh)
minPred<-min(predictionsAllGam_model3_df$predictions)
minPlot<-ifelse(minRaw>minPred, minPred, minRaw)
minPlot<-ifelse(minPlot<0, 0, minPlot)

speciesPlot<-ggplot() +
geom_tile(data=predictionsAllGam_model2_df, aes(x=x, y=y, fill=layer)) +
geom_sf(data = hull, fill = NA, color = "red", linewidth = 1) +
geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
coord_sf(crs=projection_84, xlim=c(min(locations$mean.lon), max(locations$mean.lon)), ylim=c(min(locations$mean.lat), max(locations$mean.lat))) +
scale_fill_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
scale_color_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
  labs(color=paste0(Behaviour), fill=paste0(Behaviour)) +
  ggtitle(paste(input_file1, monthChoice, colonyNameMatch, ": species-level pred")) +
  xlab("") +
  ylab("") 

popPlot<-ggplot() +
geom_tile(data=predictionsAllGam_model3_df, aes(x=x, y=y, fill=predictions)) +
geom_sf(data = hull, fill = NA, color = "red", linewidth = 1) +
geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
coord_sf(crs=projection_84, xlim=c(min(locations$mean.lon), max(locations$mean.lon)), ylim=c(min(locations$mean.lat), max(locations$mean.lat))) +
scale_fill_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
scale_color_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
  labs(color=paste0(Behaviour), fill=paste0(Behaviour)) +
  ggtitle(paste(input_file1, monthChoice, colonyNameMatch, ": pop-level pred")) +
  xlab("") +
  ylab("") 
  
rawData<-ggplot() +
#geom_tile(data=predictionsAllGam_model2_df, aes(x=x, y=y, fill=layer)) +
geom_point(data=locations, aes(x=mean.lon, y=mean.lat, color=meanBeh)) +
geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
scale_fill_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
  scale_color_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  ),
  limits = c(minPlot, maxPlot)) +
  coord_sf(crs=projection_84, xlim=c(min(locations$mean.lon), max(locations$mean.lon)), ylim=c(min(locations$mean.lat), max(locations$mean.lat))) +
  labs(color=paste0(Behaviour), fill=paste0(Behaviour)) +
  ggtitle(paste(input_file1, monthChoice, colonyNameMatch, ": raw data")) +
  xlab("") +
  ylab("")
  
  grid.arrange(rawData, speciesPlot, popPlot, nrow=3)
  
  }
  
dev.off()


### Save one of the examples as the main output so it recognizes the script has finished ###

# Save output files
print("Saving output files")

# Number 2
output_file1 <- args[5]
print("Saving output file 1")
write.csv(predictionSummaries, file = output_file1, row.names = FALSE) # AIC, dev explained, and raw data points for plotting

# Number 2
output_file2 <- args[6]
print("Saving output file 2")
saveRDS(predictionsAllGam_model1, file = output_file2) # Prediction rasters from model 1

# Number 3
output_file3 <- args[7]
print("Saving output file 3")
saveRDS(predictionsAllGam_model2, file = output_file3) # Prediction rasters from model 2

# Number 4
output_file4<- args[8]
print("Saving output file 4")
saveRDS(predictionsAllGam_model3, file = output_file4) # Prediction rasters from model 3