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
#rast.trans2 <- aggregate(rast.trans, fact = 10) 
#rast.template <- rast.gr

# Now we loop through all populations for a given species #
budgetLox<-list.files("tmp3/", full.names=TRUE)
budgetLox_species<-budgetLox[grep(input_file1, budgetLox)]
budgetLox_species_daily<-budgetLox_species[grep("daily", budgetLox_species)]

# Make lists to save interpolations in #
predictionsAllGam_species<-list()
predictionsAllGam_pop<-list()

for (i in 1:overall.iterations) {

startIteration <- Sys.time()

# Choose iterations at random 
iterationsRandom<-sample(1:100, overall.iterations, replace=FALSE) # Determine randomly which iterations I will draw from 

# Make a list to save all populations in
ActMonthlyAll<-list() # All locations & budgets
ActTemporalAll<-list() # Summary of budgets by day
predictionsAllIDW<-list() # list for saving output from IDW

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

###### Build an IDW ####

# Change projections of coordinates
monthlyidw<-meanActMonthly
coordinates(monthlyidw)<-~mean.lon + mean.lat
proj4string(monthlyidw)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
ActMonthAllTrans1<-spTransform(monthlyidw, projection_NA)

# Conduct IDW on behavior of choice & save result in list 
rastDF<-as.data.frame(rast.trans, xy=TRUE)
coordinates(rastDF) <- ~x+y
gridded(rastDF) <- TRUE
crs(rastDF)<-projection_NA

# First we interpolate the behavior #
responseVar<-paste0("prop", Behaviour)
formula_idw <- as.formula(paste0(responseVar, " ~ 1"))
idw_result <- idw(formula = formula_idw, locations = ActMonthAllTrans1, newdata = rastDF, idp = 4) # idp = the bigger it is the smoother the interpolation 
idw_raster <- raster(idw_result)
idw_raster_weighted<-idw_raster*n_distinct(ActMonthAllTrans1$individ_id)
print(n_distinct(ActMonthAllTrans1$individ_id))

# save prediction raster
predictionsAllIDW[[pop]]<-idw_raster_weighted

### Save oher data frames ###
ActMonthlyAll<-rbind(ActMonthlyAll, meanActMonthly) # This list contains all individuals, all pops
ActTemporalAll<-rbind(ActTemporalAll, ActTemporal) # This list contains time spent doing different things by day (maybe summarize by pop?)

EndPop <- Sys.time()

print("Pop loop took:")
print(EndPop-startPop)

} # end of all pops loop 

###### stack IDW & divide by total numbers of birds ####

stackV1 <- raster::stack(predictionsAllIDW)
SumV1 <- raster::calc(stackV1, fun = sum, na.rm = TRUE)
predictionsWeighted<-SumV1/n_distinct(ActMonthlyAll$individ_id)

# save prediction raster
predictionsAllIDW_allpops[[i]]<-predictionsWeighted

###### Build a GAM to predict behavior accross whole North Atlantic ######

# Make a dataset to predicton (Per's pop maps)#
xyCoords<-as.data.frame(rast.trans, xy=TRUE)
xyCoords<-xyCoords %>%
rename(coords.x1=x, coords.x2=y) %>%
dplyr::mutate(individ_id=ActMonthlyAll$individ_id[1]) %>%
dplyr::mutate(colony=ActMonthlyAll$colony[1])

# Define response variable
responseVar<-paste0("prop", Behaviour) 

# Define formula
gam_formula <- as.formula(
  paste0(responseVar, " ~ s(coords.x1, coords.x2) + s(individ_id, bs='re') + s(colony, bs='re')")
)

# Transform the coordinates of the dataset #
ActMonthlyAll1<-ActMonthlyAll
coordinates(ActMonthlyAll1)<-~mean.lon + mean.lat
proj4string(ActMonthlyAll1)<-projection_84
ActMonthlyAllTrans<-as.data.frame(spTransform(ActMonthlyAll1, projection_NA))

# Add a column to the dataset with weights of different populations (so that pops with less individuals are weighted more?)
ActMonthlyFit<-ActMonthlyAllTrans %>%
ungroup() %>%
dplyr::group_by(colony) %>%
dplyr::mutate(birds=n_distinct(individ_id), popWeights=1/birds) %>%
ungroup() %>%
dplyr::mutate(popWeightsScaled=popWeights/mean(popWeights))

# Change values so they cannot be close to 0 or 1 (necessary for the family 'betar'
ActMonthlyFit[[responseVar]]<-ifelse(ActMonthlyFit[[responseVar]]==1, 0.99, ActMonthlyFit[[responseVar]]) # can't be 1 if beta reg family
ActMonthlyFitSub <- ActMonthlyFit[ ActMonthlyFit[[responseVar]] > 0, ] # can't be 0 either

# Factorize individ id and colony
ActMonthlyFitSub $individ_id<-factor(ActMonthlyFitSub $individ_id)
ActMonthlyFitSub $colony<-factor(ActMonthlyFitSub $colony)

# Fit gam
startGam<-Sys.time()
gam1 <- gam(gam_formula, family = betar, data=ActMonthlyFitSub, weights=popWeightsScaled)
endGam<-Sys.time()
print("Gam took")
print(endGam-startGam)
Predictions_Outputs <- predict(gam1, newdata=xyCoords, exclude = c("s(population)", "s(individ_id)"), type="response")
xyCoords$predictions<-Predictions_Outputs
coordinates(xyCoords) <- ~ coords.x1 + coords.x2
statPoints<-xyCoords
predictions <- rasterize(statPoints, rast.trans, field="predictions", fun=mean)

# save prediction raster
predictionsAllGam[[i]]<-predictions

# Save gam summary for later
s <- summary(gam1)
smooth_df <- as.data.frame(s$s.table)
smooth_df$term <- rownames(smooth_df)
smooth_df$p_value <- smooth_df$`p-value`
smooth_df <- smooth_df[, c("term", "edf", "Ref.df", "Chi.sq", "p_value")]
smooth_df$rep<-i # Add iteration number
smooth_df$species<-input_file1
smooth_df$month<-monthChoice
smooth_df$Behv<-Behaviour

# Save these results
predictionSummaries<-rbind(predictionSummaries, smooth_df)

###### Add krigging ###### (has to be continuous for this to work)

# Define response variable (has to be continous for the predictions to work #
#ActMonthlyFitSub$variable<-ActMonthlyFitSub[[responseVar]]*24

# Define formula
#gam_formula2 <- as.formula(
 # paste0("variable", " ~ s(coords.x1, coords.x2) + s(individ_id, bs='re') + s(colony, bs='re')")
#)

# Re-fit gam & make new predictions #

#xyCoords<-as.data.frame(rast.trans, xy=TRUE)
#xyCoords<-xyCoords %>%
#rename(coords.x1=x, coords.x2=y) %>%
#dplyr::mutate(individ_id=ActMonthlyAll$individ_id[1]) %>%
#dplyr::mutate(colony=ActMonthlyAll$colony[1])

#gam2 <- gam(gam_formula2, data=ActMonthlyFitSub, weights=popWeightsScaled)
#Predictions_Outputs2 <- predict(gam2, newdata=xyCoords, exclude = c("s(colony)", "s(individ_id)"), type="response")
#coordinates(xyCoords) <- ~ coords.x1 + coords.x2
#xyCoords$predictions<-Predictions_Outputs2
#statPoints<-xyCoords
#predictions2 <- rasterize(statPoints, rast.trans, field="predictions", fun=mean)

# Set-up variogram for krigging #
#statPoints<-SpatialPointsDataFrame(coords = ActMonthlyFitSub[, c("coords.x1", "coords.x2")], 
 #                                    data = ActMonthlyFitSub, 
  #                                   proj4string = CRS(projection_NA))
#statPointsTMP <- statPoints
#statPointsTMP@data<- cbind(statPointsTMP@data, residLMM = resid(gam2, type="response"))
#formMod <- residLMM~ 1
#mod<- vgm(model="Sph", nugget=0.1)
#variog <- variogram(formMod, statPointsTMP)
#variogFitOLS <- fit.variogram(variog, model = mod)
# Plot the results
#plot(variog, variogFitOLS, main="Semi-variogram of LMM residuals")
#rast.template[is.na(rast.template)] <- 0
#rstPixSPDF <- as(rast.trans, "SpatialPixelsDataFrame")
#crs(rstPixSPDF)<-projection_NA

#krigTime<-Sys.time()
#residKrigMap <- krige(formula = formMod ,
#                      locations = statPointsTMP, 
#                      model = variogFitOLS,
#                      newdata = rstPixSPDF)
#krigEnd<-Sys.time()
#print("Krigging took")
#print(krigEnd-krigTime)
									 
#residKrigRaster <- raster(residKrigMap, layer="var1.pred")
#finalPred <- predictions2 + residKrigRaster

#predictionsAllGamKrig[[i]]<-finalPred

endIteration <- Sys.time()

print("Iteration took:")
print(endIteration-startIteration)


} # end of iteration looop

### Now we make some plots for saving ####

#### Make summaries from all the iterations & save these files somewhere ###

# Prepare predictions from method #1
stackV1 <- raster::stack(predictionsAllGam)
meanV1 <- raster::calc(stackV1, fun = mean, na.rm = TRUE)
sdV1 <- raster::calc(stackV1, fun = sd, na.rm = TRUE)

# Prepare predictions from method #2
#stackV2 <- raster::stack(predictionsAllGamKrig)
#meanV2 <- raster::calc(stackV2, fun = mean, na.rm = TRUE)
#sdV2 <- raster::calc(stackV2, fun = sd, na.rm = TRUE)

# Prepare predictions from method #3
stackV3 <- raster::stack(predictionsAllIDW_allpops)
meanV3 <- raster::calc(stackV3, fun = mean, na.rm = TRUE)
sdV3 <- raster::calc(stackV3, fun = sd, na.rm = TRUE)

#### make some ggplots to compare these behaviours ####

coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

method1<-as.data.frame(meanV1, xy=TRUE)
method1_sd<-as.data.frame(sdV1, xy=TRUE)
method1$method<-"Gam"
method1$sd<-method1_sd$layer

#method2<-as.data.frame(meanV2, xy=TRUE)
#method2$layer<-method2$layer/24
#method2_sd<-as.data.frame(sdV2, xy=TRUE)
#method2$method<-"Gam + krig"
#method2$sd<-method2_sd$layer

method3<-as.data.frame(meanV3, xy=TRUE)
method3_sd<-as.data.frame(sdV3, xy=TRUE)
method3$method<-"IDW"
method3$sd<-method3_sd$layer

#methodsAll<-rbind(method1, method2, method3)
methodsAll<-rbind(method1, method3)

plot1<-ggplot() +
geom_tile(data=methodsAll, aes(x=x, y=y, fill=layer)) +
scale_fill_gradientn("Prop time", colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, 1)) +
geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast) +
coord_sf(xlim=c(min(method1$x), max(method1$x)), ylim=c(min(method1$y), max(method1$y)), crs=projection_NA) +
xlab("") +
ylab("") +
ggtitle(paste0(Behaviour, " month ", monthChoice)) +
facet_wrap(~method, nrow=2)

# now we add plots of raw points & gridded data for show #

ActMonthlyAllTrans$fill_value <- ActMonthlyAllTrans[[responseVar]]

plot3<-ggplot() +
stat_summary_hex(
    data = ActMonthlyAllTrans,
    aes(x = coords.x1, y = coords.x2, z = fill_value),
    bins = 20, fun = mean
  ) +
scale_fill_gradientn("Prop time", colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, 1)) +
geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
geom_sf(data=coast)+
coord_sf(xlim=c(min(method1$x), max(method1$x)), ylim=c(min(method1$y), max(method1$y)), crs=projection_NA) +
xlab("") +
ylab("") +
ggtitle(paste0(Behaviour, " month ", monthChoice, ": Daily locations"))  

# Final temporal plot #

ActMonthlyAll$variable<-ActMonthlyAll[[responseVar]]

ActTemporalAllSum<-ActMonthlyAll%>%
ungroup() %>%
dplyr::mutate(day=as.numeric(substr(date, 9, 10))) %>%
dplyr::group_by(colony, month) %>%
dplyr::summarise(meanVar=mean(variable), sdVar=sd(variable)) 

plot5<-ggplot() +
geom_pointrange(data=ActTemporalAllSum, aes(x=colony, y=meanVar, ymin=meanVar-sdVar, ymax=meanVar +sdVar)) +
#geom_smooth(data=ActTemporalAllSum, aes(x=day, y=meanVar)) +
ylab("Prop time Beh") +
xlab("") +
theme_bw()

### Save all these in a crazy plot ###

saveString<-paste0("./results/figures/speciesMaps/activityTrial/", actBudgets$species[1], "_", Behaviour, "_month", monthChoice, ".pdf")
saveString<-gsub(" ", "", saveString)
pdf(saveString)
grid.arrange(plot3, plot5, nrow=2)
plot1
dev.off()

### Save one of the examples as the main output so it recognizes the script has finished ###

# Save output files
print("Saving output files")

# Number 2
output_file1 <- args[5]
print("Saving output file 1")
write.csv(methodsAll, file = output_file1, row.names = FALSE) # species-specific budgets

# Number 2
output_file2 <- args[6]
print("Saving output file 1")
write.csv(predictionSummaries, file = output_file2, row.names = FALSE) # species-specific budgets




