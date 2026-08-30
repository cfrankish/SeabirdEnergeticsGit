# Estimate SST within polygons #
valsSST <- raster::extract(sstRast, poly_sp)

sstpolygons<-list()
for (j in 1:length(valsSST)){
sstPolygons<-data.frame(SST=valsSST[[j]])
sstPolygons$ID<-j
sstpolygons<-rbind(sstpolygons, sstPolygons)
}

# Estimate SST in all other tiles
sstRast_terra<-rast(sstRast)
crs(sstRast_terra)<-projection_84
out <- ifel(!is.na(remainingLandscape) & remainingLandscape, sstRast_terra, NA)
sstAll<-data.frame(SST=values(out))
colnames(sstAll)<-c("SST")
sstAll$ID<-"Whole landscape"

# Join values #
sstBoth<-rbind(sstpolygons, sstAll)

ggplot() +
geom_boxplot(data=sstBoth, aes(x=ID, y=SST)) +
theme_bw() +
ylab("SST")

#### Step 5: Extract characteristics of polygons ####

print("Step 5: Extracting characteristics of polygons...")

# First we start with... species/population composition #

# I will make a key identifying the months that need to be extracted for specific seasons #
seasonsMonth<-data.frame(month=c(9, 10, 11, 12, 1, 2, 3, 4), season=c("Autumn", "Autumn", "Winter", "Winter", "Winter", "Winter", "Spring", "Spring"))

# I will also make a key identifying SST from which to extract information #
sstMonthly<-list.files("./data/sstPopMaps", full.names=TRUE)
sstMaps<-sstMonthly[grepl("tif", sstMonthly)] # Choose the tif file
sstMap<-raster::stack(sstMaps[1]) # Stack all layers (lol)
sstKey<-data.frame(layerNo=1:nlayers(sstMap)) # Make a key with dates so sstMap is easy to access
# Set start and end dates
start_date <- as.Date("2012-01-01")
end_date   <- as.Date("2023-01-01")
date_seq <- seq(from = start_date, to = end_date, by = "month") # Generate sequence for the 1st of each month
# Annotate key
sstKey$dates<-date_seq
sstKey$Year<-substr(sstKey$dates, 1, 4)
sstKey$Month<-as.numeric(substr(sstKey$dates, 6, 7))

# Loop through each season
birdsPerPolygon_seasons<-list()

for (i in 1:length(seasons)) {

#for (i in 1:2) {

# Subset to season i 
seasonSub<-seasons[i]

# Determine where pop rasters are located #
rasterList<-data.frame(files=list.files("data/popdata_raw/", full.names=TRUE), species=c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot"))

# Loop through each raster
birdsPerPolygon_species<-list()
sstPerPolygon_months<-list()

for (j in 1:nrow(rasterList)) {

#for (j in 1:3) {

# Subset to raster j #
nc<-nc_open(rasterList[j,]$files[1]) # open the file

# Create metadata of this file so I can easily pull the correct raster from it
modelcolonies <- (ncvar_get(nc,"SmcolName")) # Find list of populations
colonies <- (ncvar_get(nc,"colonyName")) # Find list of ALL colonies
colcode <- (ncvar_get(nc,"SmcolCode")) # Colony code
colLat<-(ncvar_get(nc,"colonyLatitude"))
colLon<-(ncvar_get(nc,"colonyLongitude"))
meta <- data.frame(modelcolonies, colonies, colcode, colLon, colLat) # Make some metadata so easier to understand raster structure
meta$species<-rasterList[j,]$species[1] # Add species names

# Loop through every population, choose correct month(s) & extract pop numbers within polygons #
birdsPerPolygon_populations<-list()

#for (k in 1:nrow(meta)) {

for (k in 1:nrow(meta)) {

#for (k in 1:4) {

# Subset to population k 
metaSub<-meta[k,]

# Figure out which months I need #
monthsNeeded<-subset(seasonsMonth, season==seasonSub)

# Loop through these monts #
birdsPerPolygon_months<-list()

for (l in 1:nrow(monthsNeeded)) {

print(paste0("Extracting info; season ", i, ", species ", j, ", population ", k, ", month ", l))

# Subset to month l
monthSub<-monthsNeeded[l,]$month

### Extract number of birds in polygons ###
print("Extracting No of birds...")

# Figure out which layer this is in the raster layer
iMth<-monthSub
icols<-which(ncvar_get(nc, "colonyName") %in% metaSub$colonies[1]) # Determine where colony info is located
if (length(icols)>1) {icols<-icols[1]}

# Open up raster
rast.mean <- raster(nc$filename, varname="PredictedAbundanceMean", band = icols, level = iMth) # Mean density
rast.proj<-raster::projectRaster(rast.mean, crs=projection_NA) # Project to metric coordinate system

# Extract number of birds using me polygons #
seasonalHotspots_sub<-seasonalHotspots[[i]] # subset to correct polygons
poly_sp <- as(seasonalHotspots_sub	, "Spatial")
v <- raster::extract(rast.proj, poly_sp)  # list (one element per polygon)
birdsPerPolygon<- bind_rows(lapply(seq_along(v), function(i) {
  x <- v[[i]]
  if (is.null(x)) return(NULL)  # in case a polygon hits no cells

  # x is a vector (single-layer) or a matrix (multi-layer)
  out <- as.data.frame(x)
  out$poly_id <- i
  out
}))

# put ID first
birdsPerPolygon <- birdsPerPolygon %>% relocate(poly_id)
#birdsPerPolygon$poly_id<-seasonalHotspots_sub$patches

# Measure area of the polygons so we get number of birds per m2 or something
areaPatches_km2<-st_area(seasonalHotspots_sub)/ 1e6

# Add other identifiers
birdsPerPolygonSum<-birdsPerPolygon %>%
dplyr::mutate(species=rasterList[j,]$species[1], season=seasonSub, colony=metaSub$colonies[1], colonyLat=metaSub$colLat, colonyLon=metaSub$colLon, month=monthSub) %>%
dplyr::group_by(species, season ,colony, colonyLat, colonyLon, poly_id, month) %>%
dplyr::summarise(totalBirds=sum(x, na.rm=TRUE)) %>%
ungroup() %>%
dplyr::mutate(areaPatcheskm2=areaPatches_km2, birdsPerKm2=totalBirds/areaPatcheskm2)

birdsPerPolygonSum$poly_id<-seasonalHotspots_sub$patches

### Extract SST in polygons ###
print("Extracting SST...")

if (j==1 & k==1) {
# Prepare corresponding SST surface
sstKeySub<-subset(sstKey, Month == monthSub)
sstSub <- sstMap[[c(unique(sstKeySub$layerNo))]]
sstMean<-raster::calc(sstSub, fun = mean, na.rm = TRUE)
sstMean_proj<-raster::projectRaster(sstMean, crs=projection_NA)

# Extract polygon information #
sstVals <- raster::extract(sstMean, poly_sp)  # list (one element per polygon)
sstPolygon<- bind_rows(lapply(seq_along(sstVals), function(i) {
  x <- sstVals[[i]]
  if (is.null(x)) return(NULL)  # in case a polygon hits no cells

  # x is a vector (single-layer) or a matrix (multi-layer)
  out <- as.data.frame(x)
  out$poly_id <- i
  out
}))

# Summarize information
sstPolygon_mean<-sstPolygon %>%
dplyr::group_by(poly_id) %>%
dplyr::summarise(meanSST=mean(x, na.rm=TRUE)) 

sstPolygon_mean$poly_id<-seasonalHotspots_sub$patches
birdsPerPolygonSum$poly_id<-seasonalHotspots_sub$patches

# Add month information #
sstPolygon_mean$month<-monthSub

# save info
sstPerPolygon_months<-rbind(sstPerPolygon_months, sstPolygon_mean)

} # end of SST calc

### Extract behavioral budgets ###
print("Extracting behaviour..")

# First we identify the corresponding model colony #
modelColonySub<-metaSub$modelcolonies[1]
speciesSub<-rasterList[j,]$species[1]

# Now we extract mean activity budgets from raw data #
popBudget<-extractBudget(modelColonySub, speciesSub, monthSub) 
popBudget<-popBudget %>%
rename(modelColony=colony)

# join results to main data frame
birdsPerPolygonSum2<-birdsPerPolygonSum %>%
dplyr::left_join(popBudget, by=c("species", "month"))

# Save results
birdsPerPolygon_months<-rbind(birdsPerPolygon_months, birdsPerPolygonSum2)


} # End of month loop
birdsPerPolygon_populations<-rbind(birdsPerPolygon_populations, birdsPerPolygon_months)

} # end of population loop

birdsPerPolygon_species<-rbind(birdsPerPolygon_species, birdsPerPolygon_populations)

} # End of species loop 

# Attach SST information

birdsPerPolygon_species<-birdsPerPolygon_species %>%
dplyr::left_join(sstPerPolygon_months, by=c("month", "poly_id"))

birdsPerPolygon_seasons<-rbind(birdsPerPolygon_seasons, birdsPerPolygon_species)

} # End of seasonal loop

# Intermediate save #
write.csv(birdsPerPolygon_seasons, file="./results/tables/energyhotspots.csv")

#### Step 6: Plot observed vs. random polygon characteristics ####

print("Step 6: intepreting hotspots...")

# First we add some columns so we can differentiate between observed/control patches #
birdsPerPolygon_id<-birdsPerPolygon_seasons %>%
ungroup() %>%
dplyr::mutate(type=sub("^[^_]*_([^_]*)_.*$", "\\1", poly_id)) %>%
dplyr::mutate(type=sub("^[^_]*_", "", type)) %>%
dplyr::mutate(rep=sub("^[^_]*_[^_]*_", "", poly_id)) %>%
dplyr::mutate(rep=ifelse(type=="observed", 1, rep)) %>%
dplyr::mutate(patchNo=sub("_.*$", "", poly_id)) 

# Now we sum number of birds per patch? #
birdsPerPolygon_totBirds<-birdsPerPolygon_id %>%
dplyr::group_by(season, patchNo, type, rep) %>%
dplyr::summarise(totBirds=sum(totalBirds), area=mean(areaPatcheskm2), birdsPerKm2=totBirds/area) %>%
dplyr::ungroup() %>%
dplyr::group_by(season, patchNo, type) %>%
dplyr::summarise(totBirds_mean=mean(totBirds), totBirds_sd=sd(totBirds), birdsPerKm2_mean=mean(birdsPerKm2), birdsPerKm2_sd=sd(birdsPerKm2)) %>%
dplyr::mutate(type=factor(type, levels=c("observed", "control"))) %>%
dplyr::mutate(season=factor(season, levels=c("Autumn", "Winter", "Spring"))) %>%
replace_na(list(birdsPerKm2_sd=0))

charac1<-birdsPerPolygon_totBirds %>%
filter(season=="Winter") %>%
ggplot(aes(x=patchNo, y=birdsPerKm2_mean)) +
geom_pointrange(aes(ymin=birdsPerKm2_mean-birdsPerKm2_sd, ymax=birdsPerKm2_mean + birdsPerKm2_sd, color=type), position=position_dodge2(width=0.4)) +
facet_wrap(~season, scales="free", nrow=1) +
theme_bw() +
ylab("Birds per km2") +
xlab("Patch #")

# Now we compute SST #
birdsPerPolygon_sst<-birdsPerPolygon_id %>%
dplyr::group_by(season, patchNo, type, rep) %>%
dplyr::summarise(SST_mean=mean(meanSST), SST_sd=sd(meanSST)) %>%
dplyr::ungroup() %>%
dplyr::group_by(season, patchNo, type) %>%
dplyr::summarise(SST_mean2=mean(SST_mean), SST_sd=sd(SST_mean)) %>%
dplyr::mutate(type=factor(type, levels=c("observed", "control")))

charac2<-birdsPerPolygon_sst %>%
filter(season=="Winter") %>%
ggplot(aes(x=patchNo, y=SST_mean2)) +
geom_pointrange(aes(ymin=SST_mean2-SST_sd, ymax=SST_mean2 + SST_sd, color=type), position=position_dodge2(width=0.4)) +
facet_wrap(~season, scales="free", nrow=1) +
theme_bw() +
ylab("SST") +
xlab("Patch #")

# Now we compute proportion of species per patch #
all_species <- sort(unique(birdsPerPolygon_id$species))

birdsPerPolygon_prop <- birdsPerPolygon_id %>%
  group_by(season, patchNo, type, rep, species) %>%
  summarise(totBirds = sum(totalBirds), .groups = "drop") %>%
  group_by(season, patchNo, type, rep) %>%
  complete(species = all_species, fill = list(totBirds = 0)) %>%
  mutate(
    totBirdsAll = sum(totBirds),
    propBirds   = if_else(totBirdsAll > 0, totBirds / totBirdsAll, NA_real_)
  ) %>%
  ungroup() %>%
  group_by(season, patchNo, type, species) %>%
  summarise(propBirds_mean = mean(propBirds, na.rm = TRUE), .groups = "drop") %>%
  mutate(type = factor(type, levels = c("observed", "control")))
  
charac3<-birdsPerPolygon_prop %>%
filter(season=="Winter") %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
												 dplyr::filter(type=="observed") %>%
  ggplot(aes(x = type, y = propBirds_mean, fill = species)) +
  geom_col(aes(group = type), position = "stack") +
  scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
  facet_grid(season ~ patchNo, scales = "free_x") +  # clean separation by type
  ylab("Hotspot composition") +
  theme_bw()

# Where are the birds from? -> First we make a metadata file with all information #

rasterList<-data.frame(files=list.files("data/popdata_raw/", full.names=TRUE), species=c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot"))

metaAll<-list()

for (i in 1:nrow(rasterList)) {

nc<-nc_open(rasterList[i,]$files[1]) # open the file

# Create metadata of this file so I can easily pull the correct raster from it
modelColony <- (ncvar_get(nc,"SmcolName")) # Find list of populations
colony <- (ncvar_get(nc,"colonyName")) # Find list of ALL colonies
colcode <- (ncvar_get(nc,"SmcolCode")) # Colony code
colLat<-(ncvar_get(nc,"colonyLatitude"))
colLon<-(ncvar_get(nc,"colonyLongitude"))
country<-(ncvar_get(nc,"colonyCountry"))
meta <- data.frame(modelColony, colony, country) # Make some metadata so easier to understand raster structure
meta$species<-rasterList[i,]$species[1] # Add species names
meta <- distinct(meta)


# Save results #
metaAll<-rbind(metaAll, meta)

}

metaAll_unique<-distinct(metaAll)

naModelColonies<-data.frame(modelColony=c("Faroe Islands", "Hjelmsøya", "Hornøya", "Isle of May", "Runde and Ålesund", "Røst",
"Skellig Michael", "Witless Bay", "Bjørnøya", "Cape Krutik", "Franz Josef Land", "Isle of Canna", "Kara Gate", "Russkaya Gavan", "Cape Gorodetskiy", "Coats Island", "Jan Mayen", "Oranskie Islands", "Kap Höegh", "Karmøy",
"Little Saltee"), 
country=c("Faroes", "Norway", "Norway", "Scotland", "Norway", "Norway", "Republic of Ireland", "Canada", "Norway", "Russia", "Russia", "Scotland", "Russia", "Russia", "Russia", "Canada", "Norway", "Russia", 
"Greenland", "Norway", "Republic of Ireland"))

PropCountry<-birdsPerPolygon_id %>%
dplyr::filter(type=="observed" & season=="Winter") %>%
dplyr::group_by(season, patchNo, type, species, modelColony, colony) %>%
dplyr::summarise(totBirds=sum(totalBirds), area=mean(areaPatcheskm2), birdsPerKm2=totBirds/area) %>%
dplyr::left_join(metaAll, by=c("species", "modelColony", "colony")) %>%
dplyr::left_join(naModelColonies, by=c("modelColony")) %>%
dplyr::mutate(country=ifelse(is.na(country.x), country.y, country.x)) %>%
ungroup() %>%
dplyr::group_by(season, patchNo, species, country) %>%
dplyr::summarise(birdsPerKm2_tot=sum(birdsPerKm2)) %>%
ungroup() %>%
dplyr::group_by(season, patchNo, species, country) %>%
dplyr::summarise(birdsPerKm2_tot2=sum(birdsPerKm2_tot)) %>%
ungroup() %>%
dplyr::group_by(season, patchNo, species) %>%
dplyr::mutate(totBirds=sum(birdsPerKm2_tot2), propBirds=birdsPerKm2_tot2/totBirds) %>%
dplyr::mutate(country=ifelse(propBirds<0.1, "Other", country)) %>%
dplyr::ungroup() %>%
dplyr::group_by(season, patchNo, species, country) %>%
dplyr::summarise(propBirds=sum(propBirds))

PropCountry %>%
filter(season=="Winter") %>%
  ggplot(aes(x = patchNo, y = propBirds, fill = country)) +
  geom_col(position = "stack") +
  #scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
  ylab("Hotspot composition") +
  facet_wrap(~species)
  theme_bw() 

# Now we compute time spent in different behaviours # -> weighted by number of birds contributing to activity
birdsPerPolygon_beh<-birdsPerPolygon_id %>%
dplyr::group_by(season, patchNo, type, rep, species, colony) %>%
dplyr::summarise(birds=sum(totalBirds), Flight_tot=sum(FlightHrs_mean), Forage_tot=sum(ForageHrs_mean), Active_tot=sum(ActiveHrs_mean), Rest_tot=sum(RestWaterHrs_mean))  %>%
#dplyr::filter(season=="Spring" & type=="control" & rep==8 & patchNo==2) %>%
ungroup() %>%
dplyr::mutate(Flight_weighted=Flight_tot*birds, Forage_weighted=Forage_tot*birds, Active_weighted=Active_tot*birds, Rest_weighted=Rest_tot*birds) %>%
dplyr::group_by(season, patchNo, type, rep) %>%
dplyr::summarise(totBirds=sum(birds), Flight_mean_weighted=sum(Flight_weighted/totBirds), Forage_mean_weighted=sum(Forage_weighted/totBirds),
Active_mean_weighted=sum(Active_weighted/totBirds), Rest_mean_weighted=sum(Rest_weighted/totBirds)) %>%
ungroup() %>%
dplyr::group_by(season, patchNo, type) %>%
dplyr::summarise(Flight_mean=mean(Flight_mean_weighted, na.rm=TRUE), Flight_sd=sd(Flight_mean_weighted, na.rm=TRUE),
Forage_mean=mean(Forage_mean_weighted, na.rm=TRUE), Forage_sd=sd(Forage_mean_weighted, na.rm=TRUE),
Active_mean=mean(Active_mean_weighted, na.rm=TRUE), Active_sd=sd(Active_mean_weighted, na.rm=TRUE),
Rest_mean=mean(Rest_mean_weighted, na.rm=TRUE), Rest_sd=sd(Rest_mean_weighted, na.rm=TRUE))

charac4<-birdsPerPolygon_beh %>%
filter(season=="Winter") %>%
dplyr::mutate(type=factor(type, levels=c("observed", "control"))) %>%
dplyr::mutate(season=factor(season, levels=c("Autumn", "Winter", "Spring"))) %>%
ggplot(aes(x=patchNo, y=Flight_mean/2880)) +
geom_pointrange(aes(ymin=(Flight_mean-Flight_sd)/2880, ymax=(Flight_mean + Flight_sd)/2880, shape=type, color="Flight"), position=position_dodge2(width=0.4)) + 
scale_shape_manual(values=c(19, 21)) +
facet_wrap(~season, scales="free", nrow=1) +
theme_bw() +
ylab("Time spent in Flight (Hrs)") +
xlab("Patch #")

charac5<-birdsPerPolygon_beh %>%
filter(season=="Winter") %>%
dplyr::mutate(type=factor(type, levels=c("observed", "control"))) %>%
dplyr::mutate(season=factor(season, levels=c("Autumn", "Winter", "Spring"))) %>%
ggplot(aes(x=patchNo, y=Forage_mean/2880)) +
geom_pointrange(aes(ymin=(Forage_mean-Forage_sd)/2880, ymax=(Forage_mean + Forage_sd)/2880, shape=type, color="Forage"), position=position_dodge2(width=0.4)) + 
scale_shape_manual(values=c(19, 21)) +
facet_wrap(~season, scales="free", nrow=1) +
theme_bw() +
ylab("Time spent in Forage (Hrs)") +
xlab("Patch #")

charac6<-birdsPerPolygon_beh %>%
filter(season=="Winter") %>%
dplyr::mutate(type=factor(type, levels=c("observed", "control"))) %>%
dplyr::mutate(season=factor(season, levels=c("Autumn", "Winter", "Spring"))) %>%
ggplot(aes(x=patchNo, y=Rest_mean/2880)) +
geom_pointrange(aes(ymin=(Rest_mean-Rest_sd)/2880, ymax=(Rest_mean + Rest_sd)/2880, shape=type, color="Rest"), position=position_dodge2(width=0.4)) + 
scale_shape_manual(values=c(19, 21)) +
facet_wrap(~season, scales="free", nrow=1) +
theme_bw() +
ylab("Time resting on water (hrs)") +
labs(color="Behaviour") +
#ylab("Time spent in Rest (Hrs)") +
xlab("Patch #")

charac7<-birdsPerPolygon_beh %>%
filter(season=="Winter") %>%
dplyr::mutate(type=factor(type, levels=c("observed", "control"))) %>%
dplyr::mutate(season=factor(season, levels=c("Autumn", "Winter", "Spring"))) %>%
ggplot(aes(x=patchNo, y=Active_mean/2880)) +
geom_pointrange(aes(ymin=(Active_mean-Active_sd)/2880, ymax=(Active_mean + Active_sd)/2880, shape=type, color="Active"), position=position_dodge2(width=0.4)) + 
scale_shape_manual(values=c(19, 21)) +
facet_wrap(~season, scales="free", nrow=1) +
theme_bw() +
#ylab("Time resting on water (hrs)") +
labs(color="Behaviour")+
ylab("Time spent in Active (Hrs)") +
xlab("Patch #")
