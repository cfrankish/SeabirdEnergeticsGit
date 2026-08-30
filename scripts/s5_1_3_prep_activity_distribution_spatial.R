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

#### Step 0: setting up basic conditions ####

# Set-up number of iterations...
overall.iterations<-10 # how many times this is calculated per individual
print(paste0("Determining activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1] # This will read in a population-specific activity budget
input_file2<- args[2] # This is list of all bird Ids (id catalogue) - input_file2<-"results/tables/main/table1_idcatalogue.csv"  

# open input file 1
actBudgets<-read.csv(input_file1) # actBudgets<-read.csv("tmp3/Atlanticpuffin_isleofmay_actbudgets_daily.csv")
idCatalogue<-read.csv(input_file2) # idCatalogue<-read.csv("results/tables/main/table1_idcatalogue.csv")

# Now we subset idCatalogue to the correct species & colony
idSelect<-subset(idCatalogue, species==actBudgets$species[1] & colony==actBudgets$colony[1])
birds<-unique(idSelect$individ_id)

# Now we determine where the daily files are located 
actfiles<-list.files("./tmp/", full.names=TRUE)
dailyfiles<-actfiles[grepl("Day.csv", actfiles)]

### Step 1: Map behaviours somehow ###

for (i in 1:overall.iterations) {

# Choose iterations at random 
iterationsRandom<-sample(1:100, overall.iterations, replace=FALSE) # Determine randomly which iterations I will draw from 

# And loop through every bird #
meanActMonthly<-list()

for (j in 1:length(birds)) {

# Print status update #
print(paste0("Summing activity for bird ", j, "/", length(birds)))

# Subset to bird j
ID<-birds[j]
ID<-gsub("-", "_", ID) 
ID<-gsub("Hornøya", "Hornoya", ID)

birdSub<-dailyfiles[grepl(ID, dailyfiles)]

# open csv
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

# Turn behaviors into prop of day doing something
meanBeh<-meanBehMonths %>%
ungroup() %>%
dplyr::mutate(propForage=tForage/24, propActive=tActive/24, propRestWater=tRestWater/24, propLand=tLand/24, propFlight=tFlight/24) %>%
dplyr::select(species, colony, rep, month, date, individ_id, mean.lon, mean.lat, propForage, propActive, propRestWater, propLand, propFlight)

# Save results #
meanActMonthly<-rbind(meanActMonthly, meanBeh)

}

# Subset to month x to begin with #
monthRandom<-subset(meanActMonthly, month==12)
monthRandom<-monthRandom %>%
droplevels()

# Do this month by month? #

monthRandom$individ_id<-factor(monthRandom$individ_id)
lmm.out <- lmerTest::lmer(propActive ~(mean.lon)+ (mean.lat)+ 
                                            (1|individ_id),data=monthRandom, REML=T) 
											
# Check outputs
anova(lmm.out)
plot(lmm.out)
qqnorm(resid(lmm.out))
qqline(resid(lmm.out))
plot(lmm.out, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
#qqmath(lmm.out, id = 0.05)
plot(fixef(lmm.out))#plot fixe effects

monthRandom$predictions<-predict(lmm.out, newdata=monthRandom)

## Make a grid to plot predictions on #

proj4Str<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

statPoints<-SpatialPointsDataFrame(coords = monthRandom[, c("mean.lon", "mean.lat")], 
                                     data = monthRandom, 
                                     proj4string = CRS(proj4Str))

# Open one of Per's maps at random to get a master raster layer #
nc<-nc_open("data/popdata_raw/SEATRACK_Abundance_Model_Rissa_tridactyla_Ver_3_1.nc")
rast.gr <- raster(nc$filename, varname="PredictedAbundanceMean", band = 1, level = 1)

# Now I rasterize my predictions
coordinates(monthRandom) <- ~ mean.lon + mean.lat
gridded(monthRandom) <- FALSE
r_val <- rasterize(monthRandom, rast.gr, field="propActive", fun=mean)

# Now I run the krigging

# Define the kriging parameters and fit the variogram using OLS
formMod <- residLMM~ 1
mod<- vgm(model  = "Sph" )
  variog[[i]] <- variogram(formMod[[i]], statPointsTMP[[i]])
  variogFitOLS[[i]] <- fit.variogram(variog[[i]], model = mod[[i]])
}

# Plot the results
plot(variog[[10]], variogFitOLS[[10]], main="Semi-variogram of LMM residuals")
plot(variog[[1]], variogFitOLS[[1]], main="Semi-variogram of LMM residuals")

# store the results of the krigind
residKrigMap<-list()
for (i in 1:10) { 
residKrigMap[[i]] <- krige(formula = formMod[[i]] ,
                      locations = statPointsTMP[[i]], 
                      model = variogFitOLS[[i]],
                      newdata = rstPixDF[[i]])
}


# residKrigMap is a "SpatialPixelDataFrame". 
# the other data are Raster Layer
# So I need to change residKrigMap into a raster layer

residKrigMap_pred<-list() # Store the predictions of the kriging
residKrigMap_var<-list()  # store the variance of the krigins

for (i in 1:10) {
  residKrigMap_pred[[i]] <- rasterize(x=statPoints[[i]][, c("lon_medoids", "lat_medoids")], # lon-lat data
                                 rast, # raster object
                                 field=residKrigMap[[i]]$var1.pred) # vals to fill raster with
  
  residKrigMap_var[[i]] <- rasterize(x=statPoints[[i]][, c("lon_medoids", "lat_medoids")], # lon-lat data
                                     rast, # raster object
                                     field=residKrigMap[[i]]$var1.var) # vals to fill raster with
  
                            
}

par(mfrow=c(1,1))

nb.cols<-18
pal2 <- colorRampPalette(brewer.pal(8,"PuOr"))(nb.cols)
plot(residKrigMap_pred[[10]],main="[Hg]\n(LMM regression-kriging predictions)",col=pal2,
     xlab="Longitude", ylab="Latitude", cex.main=0.8, cex.axis=0.7, cex=0.8)
plot(residKrigMap_var[[1]],main="[Hg]\n(GAM regression-kriging variance)",col=rev(pal2),
     xlab="Longitude", ylab="Latitude", cex.main=0.8, cex.axis=0.7, cex=0.8)

# create a new list of raster stack to add the new data from the kriging to the original dataset
B<-stack()
B<-list()
for (i in 1:10) {
  B[[i]]<-stack(X[[i]], YEAR[[i]], Species[[i]], Sampling_site[[i]], Ring_number[[i]], 
                Hg_HF_BF_Winter_Spatial[[i]], lon_medoids[[i]],lat_medoids[[i]], predictions[[i]], 
                residKrigMap_pred[[i]], residKrigMap_var[[i]])
  names(B[[i]]) <- c("X", "YEAR", "Species", "Sampling_site", "Ring_number", 
                     "Hg_HF_BF_Winter_Spatial", "lon_medoids","lat_medoids", "predictions", "pred_krig", "se_krig")
  #change name layers
}

plot(B[[1]])



# change to raster layer to have the same format
# residKrigMap_Pred <- as(residKrigMap$var1.pred, "RasterLayer")
# We can now add the predictions of the LMM to the predictions of the kriging
for (i in 1:10) {
  residKrigMap[[i]]<-B[[i]]$pred_krig+B[[i]]$predictions
  stack(residKrigMap[[i]])
}
plot(residKrigMap[[10]])

#residKrigMap_se<-B$se_krig+B$s
#stack(residKrigMap_se)
# create a new list of raster stack to add the new data  "residKrigMap" to the previous list "B"
C<-stack()
C<-list()
for (i in 1:10) {
  C[[i]]<-stack(X[[i]], YEAR[[i]], Species[[i]], Sampling_site[[i]], Ring_number[[i]], 
                Hg_HF_BF_Winter_Spatial[[i]], lon_medoids[[i]],lat_medoids[[i]], predictions[[i]], 
                residKrigMap_pred[[i]], residKrigMap_var[[i]],residKrigMap[[i]])
  names(C[[i]]) <- c("X", "YEAR", "Species", "Sampling_site", "Ring_number", 
                     "Hg_HF_BF_Winter_Spatial", "lon_medoids","lat_medoids", "predictions", "pred_krig", "se_krig", "residKrigMap")
  #change name layers
}

plot(C[[10]]$residKrigMap)
plot(C[[1]]$residKrigMap)
plot(C[[5]]$residKrigMap)

#On peut repasser en spatialPixelDataFrame pour pouvoir utiliser spplot. 
#residKrigMap <- as(residKrigMap, "SpatialPixelsDataFrame")
#rstPredGAM_2 <- as(rstPredGAM, "SpatialPixelsDataFrame")
#gamKrigMap <- rstPredGAM_2+ residKrigMap
