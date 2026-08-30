#### SST FUNCTIONS ####

#### Extract appropriate SST values (current) ####
# type is current vs. future

extractSSTmaps<-function(speciesSub, colonySub, type, metadataMaps, mapLox, sstLox) {
  
  # Selecting one of the six species
  species_name<-c("Fratercula_arctica","Uria_lomvia","Uria_aalge","Alle_alle",
                  "Fulmarus_glacialis","Rissa_tridactyla")
  
  # Selecting one of the six species
  common_name<-c("Atlantic puffin","Brünnich's guillemot","Common guillemot","Little auk",
                 "Northern fulmar","Black-legged kittiwake")  
  
  # set up a key to subset pop maps with
  key<-data.frame(latin=species_name, species=common_name)
  keysub<-subset(key, species==speciesSub)
  
  # List population maps
  pop.maps<-list.files(mapLox, full.names=TRUE)
  #pop.maps<-list.files("./data/intermediate/population_maps/", full.names=TRUE)
  
  # Subset to species of interest
  pop.maps.species<-pop.maps[grepl(keysub$latin, pop.maps)]
  
  # Figure out the colonies I am interest in based on model colony
  colonies<-subset(metadataMaps, modelcolonies==colonySub & species==speciesSub)
  
  # Make sure only one of each
  colonies<-colonies %>%
    dplyr::group_by(colonies) %>%
    dplyr::slice(1) 
  
  #random<-sample(1:nrow(colonies), 1) # We will just focus on one right now... to speed things up
  
  colonies<-colonies %>%
    ungroup() %>%
    dplyr::slice(1)
  
  # Create list for saving information in...
  sstAllModels<-list()
  
  # Open up sst data
  sst.data.location<-list.files(paste0(sstLox, type), full.names = TRUE)
  #sst.data.location<-list.files(paste0("./data/intermediate/sst/", type), full.names = TRUE)
  
  #for (i in seq_along(sst.data.location)) {
  
  for (i in 1) {
    
    #start.time <- Sys.time() 
    
    # open folder x
    sst.data.sub<-list.files(sst.data.location[i], full.names = TRUE)
    sstMean<-rast(sst.data.sub[grepl("mean", sst.data.sub)])
    sstSd<-rast(sst.data.sub[grepl("sd", sst.data.sub)])
    
    # Generate random sst value & then turn into data frame
    sstRandom<-sstMean
    #sstRandom<-calc(sstMean, function(x)rnorm(ncell(x), x, values(sstSd)))
    values(sstRandom)<-rnorm(ncell(sstRandom)*12, values(sstMean), values(sstSd))
    names(sstRandom)<-1:12
    sstRandomVals<-as.data.frame(sstRandom, xy=TRUE)
    sstRandomVals<-sstRandomVals %>% 
      gather(month, sstRandom, -c(x, y), na.rm = T) 
    colnames(sstRandomVals)<-c("x", "y", "month", "sstRandom")
    #sstRandomVals$month<-substr(sstRandomVals$month, 2, nchar(sstRandomVals$month))
    
    # Create list to save results in
    sstAllCols<-list()
    
    for (j in 1:nrow(colonies)) {
      
      start.time <- Sys.time()
      
      #for (i in 1:2) {
      
      print(paste0("Extracting sst for colony", j, "/", nrow(colonies), "..."))
      
      # subset sub-colony map  
      mapcolonySub<-colonies[j,]
      pop.maps.colony<-pop.maps.species[grepl(paste0("_", mapcolonySub$colonies[1], "_"), pop.maps.species)]
      
      # Make function fail is incorrect number of maps (should be 12 for all the months)
      if(!length(pop.maps.colony) %in% c(12)) {
        
        # Try & grep with shorter version of name
        pop.maps.colony<-pop.maps.species[grepl(substr(mapcolonySub$colonies[1], 1, 10), pop.maps.species)]
        
        if(!length(pop.maps.colony) %in% c(12)) {
          
          print("Number of pop maps weird")
          break
          
        }
        
      }
      
      #Stack the 12 months of data
      pop.maps.stack<-rast(pop.maps.colony)
      
      # pop data into data frame
      names(pop.maps.stack)<-1:12
      popDf<-as.data.frame(pop.maps.stack, xy=TRUE)
      popDf<-popDf %>% 
        gather(month, NoBirds, -c(x, y), na.rm = T) 
      colnames(popDf)<-c("x", "y", "month", "NoBirds")
      
      # Join values & change column names...
      joined<-popDf %>%
        #dplyr::filter(NoBirds>0) %>%
        arrange(month, x, y) %>% # otherwise file is enormous
        dplyr::left_join(sstRandomVals, by=c("x", "y", "month"))
      joined$species<-speciesSub
      joined$colony<-mapcolonySub$colonies[1]
      #joined$month<-substr(joined$month, 2, nchar(joined$month))
      joined$model<-sst.data.location[i]
      joined$sst_scenario<-type
      joined<-joined %>%
        #replace_na(list(sstRandom=999)) %>%
        dplyr::mutate(sstRandom=if_else(sstRandom< -2 & !is.na(sstRandom), -2, sstRandom))
      
      # Save all models
      sstAllCols<-rbind(sstAllCols, joined)
      
    }
    
    # Save models for all colonies
    
    sstAllCols$colonies<-mapcolonySub$colonies[1]
    sstAllCols$modelcolonies<-mapcolonySub$modelcolonies[1]
    sstAllModels<-rbind(sstAllModels, sstAllCols)
    
  }
  
  return(sstAllModels)
  
}


### Extract sst values for locations ####

extractSSTlox<-function(coords, type, sstLox) {
  
  # Establish standard projection
  projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # Flat projection
  
  # Open up sst data
  sst.data.location<-list.files(paste0(sstLox, type), full.names = TRUE)
  
  coordsModels<-list()
  
  #for (j in 1:length(sst.data.location)) {
  for (j in 1) {  
    
    # open folder x
    sst.data.sub<-list.files(sst.data.location[j], full.names = TRUE)
    sstMean<-raster::stack(sst.data.sub[grepl("mean", sst.data.sub)])
    sstSd<-raster::stack(sst.data.sub[grepl("sd", sst.data.sub)])
    
    # assign month to coordinate data
    coords$Month<-as.numeric(substr(coords$date, 6, 7))
    
    # Determine number of unique months
    months<-unique(coords$Month)
    
    # Make list to save monthly values
    coordsMonth<-list()
    
    for (k in 1:length(months)){
      
      # SST Sub
      sstMeanSub<-subset(sstMean, months[k])
      sstSdSub<-subset(sstSd, months[k])
      
      # Subset dataset
      coordsSub<-subset(coords, Month==months[k])
      
      # Extract sst at locations
      coordinates(coordsSub)<-~mean.lon + mean.lat
      proj4string(coordsSub)<-projection_84
      projMap<-proj4string(sstMeanSub)
      coordsProject<-spTransform(coordsSub, projMap)
      #loxSf<-st_as_sf(coordsProject)
      #bufferPoints <- st_buffer(loxSf, dist=200000)
      coordsProject$sst_mean<-raster::extract(sstMeanSub, coordsProject)
      coordsProject$sst_sd<-raster::extract(sstSdSub, coordsProject)
      
      # Turn back into data frame
      coordsDfSub<-data.frame(coordsProject) %>%
        dplyr::group_by(date) %>%
        dplyr::mutate(sst_random=ifelse(!is.na(sst_mean) & !is.na(sst_sd), rnorm(1, mean=sst_mean, sd=sst_sd), NA))
      
      # save monthly result
      coordsMonth<-rbind(coordsMonth, coordsDfSub)
      
    }  
    
    # Save for all models
    coordsMonth$sstModel<-sst.data.location[j]
    coordsMonth$time<-type
    coordsModels<-rbind(coordsModels, coordsMonth)
    
  }
  
  return(coordsModels)
  
  
}

############## Plotting activity & energy ######

scale_0_1 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Plot population-specific activity #

mapActivity_pop<-function(data, extentdf, behaviour, species) {

# Determine unique months we will need to plot #
months<-c(9, 10, 11, 12, 1, 2, 3, 4)

# prepare empty PDF plot
pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/popMaps/activity/", species, "_", data$colony_og[1], "_", behaviour, ".pdf"), width=10)

print("Mapping")

for (i in 1:length(months)) {

# Determine behaviour that will be plotted
behaviorplot<-paste0("time", behaviour)
behaviorplot2<-paste0("time", behaviour, "_sd")

# Determine possible maximum values
dataSub<-subset(data, month==months[i])
maxBeh1<-max(dataSub[[behaviorplot]])
maxBeh2<-max(dataSub[[behaviorplot2]])
maxBeh<-ifelse(maxBeh1>maxBeh2, maxBeh1, maxBeh2)

dataSub$MeanvsSD<-dataSub[[behaviorplot]]-dataSub[[behaviorplot2]]
minBeh1<-min(dataSub[[behaviorplot]])
minBeh2<-min(dataSub$MeanvsSD)
minBeh<-ifelse(minBeh1<minBeh2, minBeh1, minBeh2)

# plot Mean behaviour
plotMean<-ggplot() +
  geom_tile(data=subset(data, month==months[i]), aes_string(x="x", y="y", fill=paste0(behaviorplot))) +
  geom_text(data=subset(data, month==months[i] & NoBirds==1), aes_string(x="x", y="y", label="NoBirds"), cex=1.5) +
  #scale_alpha(range = c(0.5, 1)) +
  scale_fill_gradientn(paste0(behaviorplot), colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$coords.x1)-400000, max(extentdf$coords.x1) + 400000), ylim=c(min(extentdf$coords.x2)-400000, max(extentdf$coords.x2) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
  
# plot SD behaviour
plotSD<-ggplot() +
  geom_tile(data=subset(data, month==months[i]), aes_string(x="x", y="y", fill=paste0(behaviorplot2))) +
   geom_text(data=subset(data, month==months[i] & NoBirds==1), aes_string(x="x", y="y", label="NoBirds"), cex=1.5) +
  scale_fill_gradientn(paste0(behaviorplot2), colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$coords.x1)-400000, max(extentdf$coords.x1) + 400000), ylim=c(min(extentdf$coords.x2)-400000, max(extentdf$coords.x2) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
  
# Plot areas with high mean & low SD

data$meanVSsd<-data[[behaviorplot]] - data[[behaviorplot2]]
dataSub$mean_scaled<-scale_0_1(dataSub[[behaviorplot]])
dataSub$sd_scaled<-scale_0_1(dataSub[[behaviorplot2]])
dataSub$meanVSsd_scaled<-dataSub$mean_scaled- dataSub$sd_scaled

plotMeanSD<-ggplot() +
  geom_tile(data=subset(dataSub, month==months[i] & NoBirds >1), aes_string(x="x", y="y", fill="MeanvsSD")) +
  #geom_text(data=subset(dataSub, month==months[i] & NoBirds>1), aes_string(x="x", y="y", label="NoBirds"), cex=1.5) +
  scale_fill_gradientn('Mean - SD', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"),limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$coords.x1)-400000, max(extentdf$coords.x1) + 400000), ylim=c(min(extentdf$coords.x2)-400000, max(extentdf$coords.x2) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
  
grid.arrange(plotMean, plotSD, plotMeanSD, nrow=1)
  
  }
  
  print("Saving")
  
  dev.off()


}

# Plot population-specific activity #

mapActivity_species<-function(data_mean, data_sd, extentdf, behaviour, species) {

# Determine unique months we will need to plot #
months<-c(9, 10, 11, 12, 1, 2, 3, 4)

# Make some bins for plotting SD diff vs. colony number #
# Determine behaviour that will be plotted
behaviorplot<-paste0("time", behaviour, "Mean")
behaviorplot2<-paste0("time", behaviour, "_sd_diff")
behaviorplot3<-paste0("time", behaviour, "_sd")
data_sd$sd_diff<-data_sd[[behaviorplot2]] # give column standard name so easy to calculate quantiles
#data_sd_sub<-subset(data_sd, colonies>0) # remove 0 data
#data_sd_sub$sd_diff_bin<- cut(data_sd_sub$sd_diff, breaks = c(min(data_sd_sub$sd_diff),0, median(data_sd_sub$sd_diff), max(data_sd_sub$sd_diff)), include.lowest = TRUE)
#medCol<-median(data_sd_sub$colonies)
#medColFinal<-ifelse(medCol==2, 3, medCol)
#data_sd_sub$colony_bin<- cut(data_sd_sub$colonies, breaks = c(0, 1, medCol, max(data_sd_sub$colonies)), include.lowest = TRUE)
#data_sd_bi<-bi_class(data_sd_sub, x = sd_diff_bin, y = colony_bin, style = "quantile", dim = 3)

# prepare empty PDF plot
pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/activity/", species, "_", behaviour, ".pdf"), width=10, height=5)

print("Mapping")

for (i in 1:length(months)) {

# Scale differences for plotting
data_sd_sub<-subset(data_sd, month==months[i])
data_sd_sub$sd_diff2<-ifelse(data_sd_sub$sd_diff>=2, 2, data_sd_sub$sd_diff) # So that we can see the other differences better

# Make a mean vs. sd plot to see if any areas with high mean and low between-pop sd
data_sd_sub$mean_scaled<-scale_0_1(data_sd_sub[[behaviorplot]])
data_sd_sub$sd_scaled<-scale_0_1(data_sd_sub[[behaviorplot2]])
data_sd_sub$meanVSsd_scaled<-data_sd_sub[[behaviorplot]]- data_sd_sub[[behaviorplot3]]

# Determine possible maximum values
#dataSub<-subset(data, month==months[i])
data_mean_sub<-subset(data_mean, month==months[i])
maxBeh1<-max(data_mean_sub[[behaviorplot]])
maxBeh2<-max(data_sd_sub[[behaviorplot3]])
maxBeh<-ifelse(maxBeh1>maxBeh2, maxBeh1, maxBeh2)

minBeh1<-min(data_mean_sub[[behaviorplot]])
minBeh2<-min(data_sd_sub$meanVSsd_scaled)
minBeh<-ifelse(minBeh1<minBeh1, minBeh1, minBeh2)

# plot Mean behaviour
plotMean<-ggplot() +
  geom_tile(data=subset(data_mean, month==months[i]), aes_string(x="x", y="y", fill=paste0(behaviorplot))) +
  geom_text(data=subset(data_sd, month==months[i] & colonies<2 & colonies >0), aes_string(x="x", y="y", label="colonies"), cex=1.5) +
  #scale_alpha(range = c(0.5, 1)) +
  scale_fill_gradientn(paste0(behaviorplot), colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$x)-400000, max(extentdf$x) + 400000), ylim=c(min(extentdf$y)-400000, max(extentdf$y) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
  
# plot SD behaviour

# Calculate quantiles for SD
# Calculate quantiles for colonie #

plotSD<-ggplot() +
geom_tile(data=subset(data_sd_sub ,month==months[i] & colonies >1), aes_string(x="x", y="y", fill="sd_diff2"), color="black") +
#geom_tile(data=data_sd_sub, aes_string(x="x", y="y", fill="sd_diff"), color="black") +
#geom_text(data=filter(data_sd_sub, colonies <2), aes_string(x="x", y="y", label="colonies"), cex=1.5) +
  #geom_tile(data=subset(data_sd, month==months[i]), aes_string(x="x", y="y", fill=paste0(behaviorplot2))) +
  #geom_text(data=subset(data_sd_bi, month==months[i]), aes_string(x="x", y="y", label="colonies"), cex=1.5) +
  scale_fill_gradientn(
  colours = c(
    "#364B9A", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#A50026"
  ),
  limits = c(min(data_sd_sub$sd_diff), max(data_sd_sub$sd_diff2)),
  values = scales::rescale(c(min(data_sd_sub$sd_diff), 1, max(data_sd_sub$sd_diff2))),
  na.value = "darkblue"
) +
  #scale_fill_gradientn(paste0(behaviorplot2), colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026")) +
  #bi_scale_fill(pal = "DkCyan", dim = 3) + 
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$x)-400000, max(extentdf$x) + 400000), ylim=c(min(extentdf$y)-400000, max(extentdf$y) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  labs(fill="Between/within-pop SD") +
  theme(legend.position="bottom")

plotMeanSD<-ggplot() +
  geom_tile(data=subset(data_sd_sub, month==months[i] & NoBirds >1), aes_string(x="x", y="y", fill="meanVSsd_scaled")) +
  #geom_text(data=subset(dataSub, month==months[i] & NoBirds>1), aes_string(x="x", y="y", label="NoBirds"), cex=1.5) +
  scale_fill_gradientn('Mean - SD', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"),limits=c(minBeh, maxBeh)) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(extentdf$x)-400000, max(extentdf$x) + 400000), ylim=c(min(extentdf$y)-400000, max(extentdf$y) + 400000)) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  ggtitle(paste0(species, ": month", months[i])) +
  theme(legend.position="bottom")
 
grid.arrange(plotMean, plotSD, plotMeanSD, nrow=1)
 
# Arrange the first two plots side by side
#top_row <- arrangeGrob(plotMean, plotSD, nrow = 1)

# Combine the top row with the bottom plot
#final_layout <- arrangeGrob(top_row, rangeSD, nrow = 2, heights = c(2, 1))

#plot(final_layout)
  
  }
  
  print("Saving")
  
  dev.off()


}

##### Compute Moran's i - pop-level

computeMoran_pop<-function(df_mean, behaviour) {

# Determine unique months to loop through #
months<-unique(df_mean$month)

# Assign behaviour to a specific column #
behaviorplot<-paste0("time", behaviour)
df_mean$meanVal<-df_mean[[behaviorplot]] # give column standard name so easy to calculate quantiles

# Make a list to save these values in
moranRes<-list()

for (i in 1:length(months)) {

print(i)

# Subset data frame to month i
data_mean_i<-subset(df_mean, month==months[i])

data_mean_i_arrange<-data_mean_i %>%
ungroup() %>%
arrange(x, y) %>%
dplyr::group_by(index) %>%
dplyr::slice(1)

# Calculate range of values #
minVal<-min(data_mean_i$meanVal)
maxVal<-max(data_mean_i$meanVal)
if(abs(maxVal - minVal) < 1e-6) {
next}

coords <- cbind(data_mean_i$x, data_mean_i$y) # Make a coordinate grid to choose from
# 4 nearest neighbors
knn <- knearneigh(coords, k = 4)
nb <- knn2nb(knn)
# Row-standardized weights
lw <- nb2listw(nb, style = "W")

# moran test
moran_mean <- moran.test(data_mean_i$meanVal, lw) # Looks at spatial clustering at-species level
#moran_mean_pop<-moran.test(data_sd_i$sd_diff, lw)

# Make a df with results
results<-data.frame(month=months[i])
resultName<-paste0("Moran_", behaviour)
results[[resultName]]<-moran_mean$estimate[1]

# Save results
moranRes<-rbind(moranRes, results)

}

return(moranRes)

}

##### Compute Moran's i - species level #####

computeMoran<-function(df_mean, df_sd, behaviour) {

# Determine unique months to loop through #
months<-unique(df_mean$month)

# Assign behaviour to a specific column #
behaviorplot<-paste0("time", behaviour, "Mean")
behaviorplot2<-paste0("time", behaviour, "_sd_diff")
data_mean<-df_mean
data_sd<-df_sd
data_mean$meanVal<-data_mean[[behaviorplot]] # give column standard name so easy to calculate quantiles
data_sd$sd_diff<-data_sd[[behaviorplot2]] # give column standard name so easy to calculate quantiles

# Make a list to save these values in
moranRes<-list()

for (i in 1:length(months)) {

# Subset data frame to month i
data_mean_i<-subset(data_mean, month==months[i])
data_sd_i<-subset(data_sd, month==months[i])

# Group coordinates together
data_mean_i<-data_mean_i %>%
ungroup() %>%
dplyr::arrange(x, y)

data_sd_i<-data_sd_i %>%
ungroup() %>%
dplyr::arrange(x, y)

coords <- cbind(data_mean_i$x, data_mean_i$y)

# 4 nearest neighbors
knn <- knearneigh(coords, k = 4)
nb <- knn2nb(knn)

# Row-standardized weights
lw <- nb2listw(nb, style = "W")

# moran test
moran_mean <- moran.test(data_mean_i$meanVal, lw) # Looks at spatial clustering at-species level
moran_mean_pop<-moran.test(data_sd_i$sd_diff, lw)

# Make a df with results
results<-data.frame(month=months[i], moranI_speciesMean=moran_mean$estimate[1], moranI_popSD=moran_mean_pop$estimate[1], behaviour=behaviour)

# Save results
moranRes<-rbind(moranRes, results)

}

return(moranRes)

}

### Plot for calculating correlations by month between mean & SD activity ###

corrMonth<-function(data) {

# Determine months to loop through #
months<-unique(data$month)

# Make an empty list to save the results in
dataCorRes<-list()

for (i in 1:length(months)) {

# Subset to month i #
dataMonth_i<-subset(data, month==months[i])

# Determine correlations #
corrFlight<-cor(dataMonth_i$timeFlight, dataMonth_i$timeFlight_sd, method = "pearson")
corrRestWater<-cor(dataMonth_i$timeRestWater, dataMonth_i$timeRestWater_sd, , method = "pearson")
corrLand<-cor(dataMonth_i$timeLand, dataMonth_i$timeLand_sd,  method = "pearson")
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
corrForage<-cor(dataMonth_i$timeForage, dataMonth_i$timeForage_sd, method = "pearson")
} else {
corrActive<-cor(dataMonth_i$timeActive, dataMonth_i$timeActive_sd, method = "pearson")
}

# Add these values back to the main dataset
dataMonth_i$corFlight<-corrFlight
dataMonth_i$corRestWater<-corrRestWater
dataMonth_i$corLand<-corrLand
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
dataMonth_i$corForage<-corrForage
dataMonth_i$corActive<-NA
} else {
dataMonth_i$corForage<-NA
dataMonth_i$corActive<-corrActive
}

# Save monthly values
dataMonth_i_res<-dataMonth_i %>%
dplyr::select(month, corFlight, corRestWater, corForage, corLand, corActive) %>%
dplyr::group_by(month) %>%
dplyr::slice(1)

dataCorRes<-rbind(dataCorRes, dataMonth_i_res)

}

return(dataCorRes)

}

# At species level
corrMonth_species<-function(data) {

# Determine months to loop through #
months<-unique(data$month)

# Make an empty list to save the results in
dataCorRes<-list()

for (i in 1:length(months)) {

# Subset to month i #
dataMonth_i<-subset(data, month==months[i])

# Determine correlations #
corrFlight<-cor(dataMonth_i$timeFlightMean, dataMonth_i$timeFlight_sd, method = "pearson")
corrRestWater<-cor(dataMonth_i$timeRestWaterMean, dataMonth_i$timeRestWater_sd, , method = "pearson")
corrLand<-cor(dataMonth_i$timeLandMean, dataMonth_i$timeLand_sd,  method = "pearson")
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
corrForage<-cor(dataMonth_i$timeForageMean, dataMonth_i$timeForage_sd, method = "pearson")
} else {
corrActive<-cor(dataMonth_i$timeActiveMean, dataMonth_i$timeActive_sd, method = "pearson")
}

# Add these values back to the main dataset
dataMonth_i$corFlight<-corrFlight
dataMonth_i$corRestWater<-corrRestWater
dataMonth_i$corLand<-corrLand
if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
dataMonth_i$corForage<-corrForage
dataMonth_i$corActive<-NA
} else {
dataMonth_i$corForage<-NA
dataMonth_i$corActive<-corrActive
}

# Save monthly values
dataMonth_i_res<-dataMonth_i %>%
dplyr::select(month, corFlight, corRestWater, corForage, corLand, corActive) %>%
dplyr::group_by(month) %>%
dplyr::slice(1)

dataCorRes<-rbind(dataCorRes, dataMonth_i_res)

}

return(dataCorRes)

}

### Calculate map extent for spatial comparisons ###

# This function will find the indices of the xy coordinates that cover the minimum square which comprises areas where there are birds accross all populations
# The reason for this is that I need a map of similar extent to conduct analyses like moran's I or correlations that are comparable accross pops

mapExtent<-function(actFiles) {

modelColonies<-list()

print("Starting for loop")

for (i in 1:length(actFiles)) {

print(paste0("Mapping file", i, "/", length(actFiles), "..."))

#### Open file i ###
actSub<-fread(actFiles[i])
actSub<-subset(actSub, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))

# Determine the model colony
metaSub<-subset(meta, colonies==actSub$colony[1])
modelColSub<-metaSub$modelcolonies[1]
print(modelColSub)

# Does the list already have this model colony?
if (modelColSub %in% modelColonies) {
print("Next")
next
}

# subset to data frame with only birds
actSubIndex<-actSub %>%
ungroup() %>%
dplyr::group_by(month) %>%
arrange(x, y) %>%
dplyr::mutate(index=row_number()) 
actSub_crop<-subset(actSubIndex, NoBirds >0)

if (i ==1) {

# Save results
minX<-min(actSub_crop$x)
maxX<-max(actSub_crop$x)
minY<-min(actSub_crop$y)
maxY<-max(actSub_crop$y)

# Figure out which indices these are #
index1<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  dplyr::filter(near(x, minX), near(y, minY)) %>%
  dplyr::mutate(type="corner1") %>%
  dplyr::select(index, type) 
index1<-index1$index
  
index2<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, minX), near(y, maxY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner2")
index2<-index2$index
  
index3<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, minY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner3")
index3<-index3$index
  
index4<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, maxY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner4")
index4<-index4$index

}

if (i>1) {

minX_2<-min(actSub_crop$x)
maxX_2<-max(actSub_crop$x)
minY_2<-min(actSub_crop$y)
maxY_2<-max(actSub_crop$y)

if(minX_2<minX) {

minX<-minX_2
index1<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  dplyr::filter(near(x, minX), near(y, minY)) %>%
  dplyr::mutate(type="corner1") %>%
  dplyr::select(index, type) 
index1<-index1$index
  
index2<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, minX), near(y, maxY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner2")
index2<-index2$index
}

if(maxX_2>maxX) {

maxX<-ifelse(maxX_2>maxX, maxX_2, maxX)
index3<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, minY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner3")
index3<-index3$index
  
index4<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, maxY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner4")
index4<-index4$index

}

if(minY_2<minY) {

minY<-ifelse(minY_2<minY, minY_2, minY)
index1<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  dplyr::filter(near(x, minX), near(y, minY)) %>%
  dplyr::mutate(type="corner1") %>%
  dplyr::select(index, type) 
index1<-index1$index
  
index3<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, maxX), near(y, minY)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner3")
index3<-index3$index

}

if (maxY_2>maxY) {

maxY<-ifelse(maxY_2>maxY, maxY_2, maxY)

index2<-actSubIndex %>%
dplyr::filter(month==1) %>%
  ungroup() %>%
  filter(near(x, minX), near(y, maxY_2)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner2")
index2<-index2$index
  
index4<-actSubIndex %>%
#dplyr::filter(month==1) %>%
  ungroup() %>%
  dplyr::filter(near(x, maxX), near(y, maxY_2)) %>%
  dplyr::select(index) %>%
  dplyr::mutate(type="corner4")
index4<-index4$index

}

} 

extent<-data.frame(minX, maxX, minY, maxY, index1, index2, index3, index4)

# Otherwise add the modelcolony name to the list
modelColonies<-append(modelColonies, modelColSub)

}

# Join the indices together #

extent2<-extent %>%
dplyr::rename(index=index1) %>%
dplyr::mutate(type="corner1") %>%
dplyr::select(index, type)

extent3<-extent %>%
dplyr::rename(index=index2) %>%
dplyr::mutate(type="corner2") %>%
dplyr::select(index, type)

extent4<-extent %>%
dplyr::rename(index=index3) %>%
dplyr::mutate(type="corner3") %>%
dplyr::select(index, type)

extent5<-extent %>%
dplyr::rename(index=index4) %>%
dplyr::mutate(type="corner4") %>%
dplyr::select(index, type)

extentAll<-rbind(extent2, extent3, extent4, extent5)

return(extentAll)


}

#### Calculate pair-wise population correlations ####

popcorr<-function(data, dataFiles, meta, behaviour) {

# Create behaviour we will select
behaviourSelect<-paste0("time", behaviour)

# Create an empty list with colonies 
colonies<-list()
poptopopcor<-list()

# Determine current colony
metaSub<-subset(meta, colonies==data$colony[1])
modelColSub<-metaSub$modelcolonies[1]

# Add it to empty list
colonies<-rbind(colonies, modelColSub)

# Now loop through act files & conduct pair-wise analyses

for (i in 1:length(actFiles)) {

#### Open file i ###
actSub<-fread(actFiles[i])
actSub<-subset(actSub, month %in% c(1, 2, 3, 4, 9, 10, 11, 12))

# Add an index for easier joining
actSub<-actSub %>%
  ungroup() %>%
  arrange(month, x, y) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(index=row_number()) %>%
  ungroup() 

# Determine which colony it is
metaSub<-subset(meta, colonies==actSub$colony[1])
modelColSub<-metaSub$modelcolonies[1]

# Check to see if it exists already #
if (modelColSub %in% colonies) {
print("Next")
next
}

# otherwise we proceed to month by month comparisons
corRes<-list()

# Months loop
monthsLoop<-c(9, 10, 11, 12, 1, 2, 3, 4)

print("Conducting pair-wise comparison")

for (j in 1:length(monthsLoop)) {

# Subset dataset to month j #
actMonth<-subset(actSub, month==monthsLoop[j])

# Subset to desired extent
indices<-data.frame(index=unique(data$index))

actMonthExt<-actMonth %>%
dplyr::inner_join(indices, by=c("index"))

# Subset dataset a 
dataMonth<-subset(data, month==monthsLoop[j])

# Now I find the intersection of both (i.e. I reduce to cells with birds in both)
actCompare_birds<-subset(actMonthExt, NoBirds>0)
actOriginal_birds<-subset(dataMonth, NoBirds>0)

actCompareSub<-actCompare_birds %>%
dplyr::filter(index %in% c(unique(actOriginal_birds$index)))

actOriginalSub<-actOriginal_birds %>%
dplyr::filter(index %in% c(unique(actCompare_birds$index)))

if (nrow(actCompareSub)==0 | nrow(actOriginalSub)==0) {
next
}

# Conduct pair-wise correlation for every behaviour #
corBeh<-cor(actCompareSub[[behaviourSelect]], actOriginalSub[[behaviourSelect]], method = "pearson")

# Save result for month & population #
results<-data.frame(behaviour, cor=corBeh, month=monthsLoop[j], colony=modelColSub)
corRes<-rbind(corRes, results)

}

# Save overall results
poptopopcor<-rbind(poptopopcor, corRes)

}

poptopopcor$colony_og<-colonies[1]

return(poptopopcor)

}

### Function to extract value from within mean +/- SE for plotting ###

energyEstimate<-function(data, iterations) {

print("Summing energy iteratively...")

allEnergy<-list()

for (i in 1:iterations) {

# Generate random energy expenditure number #
data<-data %>%
dplyr::filter(!is.na(energyPopkJ_mean)) %>%
ungroup()
data$energyRandom<-rnorm(n=nrow(data), mean=data$energyPopkJ_mean, sd=data$energyPopkJ_sd)

# Sum the values per month
dataMonth<-data %>%
dplyr::group_by(month) %>%
dplyr::summarise(energyTot=sum(energyRandom))

# Add iteration No
dataMonth$iteration<-i

# Save results #
allEnergy<-rbind(allEnergy, dataMonth)

}

return(allEnergy)


}

### Function to grid behaviour for seabirds ###

gridBeh<-function(data, map) {

print("Gridding behaviour..")

# Turn map into a grid #
grid<-as.data.frame(map, xy=TRUE)

# Add columns that I want 
grid$timeFlight<-0 # time spent in flight in hours at colony 1
grid$timeForage<-0 # time spent foraging in hours at colony 1
grid$timeActive<-0 # time spent foraging in hours at colony 1
grid$timeLand<-0 # time spent foraging in hours at colony 2
grid$timeRestWater<-0 # time spent foraging in hours at colony 2
grid$timeTotal<-0 # total time spent
grid<-grid %>%
arrange(x, y) 

# Determine number of months to loop through
Months<-unique(data$month)

# Make an an empty list to save the results in
BirdActivityFinal_All<-list()

print("Month...")

for (i in 1:length(Months)) {

print(paste0(i))

# Subset to month i
monthSub<-Months[i]

# Cut actBudgets
actMonth<-subset(data, month==monthSub)

# Calculate total time in different behaviours per month
totalTime<-actMonth %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id, month) %>%
dplyr::summarise(totFlight=sum(tFlight), totActive=sum(tActive), totRest=sum(tRestWater), totLand=sum(tLand), totForage=sum(tForage), totHrs=sum(totFlight, totActive, totRest, totLand, totForage), .groups = "drop")

# Rename our grid
gridMonth<-grid %>%
dplyr::filter(x>min(actMonth$mean.lon) - (res(map) + 1) & x<max(actMonth$mean.lon) + res(map) + 1 & y > min(actMonth$mean.lat) - (res(map) +1) & y < max(actMonth$mean.lat) + res(map) +1)

for (m in 1:nrow(grid)) {

# Subset grid x
gridSub<-gridMonth[m,]

resx<-res(map)[1]
resy<-res(map)[2]

# Subset coordinates which fit #
loxSub1<-subset(actMonth,   mean.lon > gridSub$x & mean.lon < gridSub$x + resx & mean.lat > gridSub$y & mean.lat < gridSub$y + resy)

if (nrow(loxSub1)>0) {

# Calculate time spent in flight
timeFlight<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tFlight), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent foraging
timeForage<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tForage), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeLand<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tLand), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeActive<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tActive), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate time spent on land
timeRest<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tRestWater), totTimeAll=sum(tFlight, tActive, tRestWater, tLand, tForage), .groups = "drop") %>%
dplyr::mutate(propTime=totTime/totTimeAll)

# Calculate total time
timeTot<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(tFlight, tRestWater, tForage, tActive, tLand), .groups = "drop") %>%
dplyr::full_join(totalTime, by=c("species", "colony", "individ_id")) %>%
replace_na(list(totTime=0)) %>%
dplyr::mutate(propTime=totTime/totHrs)

#nas<-subset(timeLand, is.na(propTime))
#if(nrow(nas)>0) {break}

# Save results in grid for plotting
gridMonth$timeFlight[m]<-mean(timeFlight$propTime)
gridMonth$timeForage[m]<-mean(timeForage$propTime)
gridMonth$timeActive[m]<-mean(timeActive$propTime)
gridMonth$timeLand[m]<-mean(timeLand$propTime)
gridMonth$timeRestWater[m]<-mean(timeRest$propTime)
gridMonth$timeTotal[m]<-mean(timeTot$propTime)

# Replace NAs with 0
gridMonth[is.na(gridMonth)] <- 0

}

}

# Add number of locations
gridMonth$month<-Months[i]
gridMonth$individ_id<-data$individ_id[1]
gridMonth$rep<-data$rep[1]

# Save results
BirdActivityFinal_All<-rbind(BirdActivityFinal_All, gridMonth)

}

return(BirdActivityFinal_All)

}

### Same function as above but much faster I think ###

gridBeh2<-function(data, map) {

print("Gridding behaviour...")
  
  # Get raster
  r <- map
  resx <- res(r)[1]
  resy <- res(r)[2]
  
  # Assign each point to a raster cell
  data <- data %>%
    mutate(cell = cellFromXY(r, cbind(mean.lon, mean.lat)))
  
  # Aggregate behavior times per cell, month, species, colony, individual
  agg <- data %>%
    group_by(month, species, colony, individ_id, cell) %>%
    summarise(
      timeFlight = sum(tFlight, na.rm = TRUE),
      timeForage = sum(tForage, na.rm = TRUE),
      timeActive = sum(tActive, na.rm = TRUE),
      timeLand   = sum(tLand, na.rm = TRUE),
      timeRestWater = sum(tRestWater, na.rm = TRUE),
      timeTotal = sum(tFlight + tForage + tActive + tLand + tRestWater, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add raster coordinates for plotting
  coords <- xyFromCell(r, agg$cell)
  agg <- agg %>%
    mutate(x = coords[,1],
           y = coords[,2])
  
  # Replace NA with 0 just in case
  agg <- agg %>%
    mutate(across(c(timeFlight, timeForage, timeActive, timeLand, timeRestWater, timeTotal), ~replace_na(., 0)))
	
	agg2 <- agg %>%
  mutate(
    propFlight = timeFlight / timeTotal,
    propForage = timeForage / timeTotal,
    propActive = timeActive / timeTotal,
    propLand   = timeLand / timeTotal,
    propRestWater = timeRestWater / timeTotal
  )
  
  return(agg2)
}

### Function to carry out one type of spatial interpolation ###

run_gam_grid <- function(ActMonthlyMap, rast.template, month, beh) {

  dat <- subset(ActMonthlyMap, month == month)
  if (nrow(dat) == 0) return(NULL)

  dat$individ_id <- factor(dat$individ_id)

  xyCoords <- as.data.frame(rast.template, xy = TRUE)
  xyCoords$individ_id <- dat$individ_id[1]

  # Prepare behavior to plot #
  behPredict<-paste0("prop", beh)

  # behaviour column used dynamically
  f <- as.formula(paste0(behPredict, " ~ s(x, y) + s(individ_id, bs='re')"))

  gam_mod <- mgcv::gam(
    formula = f,
    data = dat,
    family = betar,
    weights = timeTotal
  )

  preds <- predict(gam_mod, newdata = xyCoords, exclude = "s(individ_id)", type = "response")
  xyCoords$predictions <- preds

  coordinates(xyCoords) <- ~ x + y
  pred_raster <- rasterize(xyCoords, rast.template, field = "predictions", fun = mean)

  return(pred_raster)
}

run_idw <- function(gridBehBird, rast.template, month, beh) {

  grid_sub <- subset(gridBehBird, month == month)
  if (nrow(grid_sub) == 0) return(NULL)

  # behaviour field = beh
  names(grid_sub)[names(grid_sub) == beh] <- "beh_value"

  coordinates(grid_sub) <- ~x + y
  rastDF <- as.data.frame(rast.template, xy = TRUE)
  coordinates(rastDF) <- ~x + y
  gridded(rastDF) <- TRUE

  idw_result <- gstat::idw(
    formula = beh_value ~ 1,
    locations = grid_sub,
    newdata = rastDF,
    idp = 1
  )

  idw_raster <- raster(idw_result)
  return(idw_raster)
}

run_gam_points <- function(monthRandom, rast.template, beh) {

  dat <- monthRandom
  if (!beh %in% colnames(dat)) return(NULL)

  dat_sub <- dat[dat[[beh]] > 0 & dat[[beh]] < 1, ]   # beta family restrictions
  if (nrow(dat_sub) == 0) return(NULL)

  dat_sub$individ_id <- factor(dat_sub$individ_id)

  f <- as.formula(paste0(beh, " ~ s(mean.lon, mean.lat) + s(individ_id, bs='re')"))

  gam_mod <- mgcv::gam(
    formula = f,
    data = dat_sub,
    family = betar
  )

  xyCoords <- as.data.frame(rast.template, xy = TRUE)
  xyCoords$individ_id <- dat_sub$individ_id[1]
  names(xyCoords)[names(xyCoords) == "x"] <- "mean.lon"
  names(xyCoords)[names(xyCoords) == "y"] <- "mean.lat"

  preds <- predict(gam_mod, newdata = xyCoords, exclude = "s(individ_id)", type = "response")

  pred_df <- data.frame(predictions = preds, mean.lon = xyCoords$mean.lon, mean.lat = xyCoords$mean.lat)
  coordinates(pred_df) <- ~ mean.lon + mean.lat

  pred_raster <- rasterize(pred_df, rast.template, field = "predictions", fun = mean)

  return(pred_raster)
}

run_gam_krig <- function(monthRandom, rast.template, projection_NA, beh) {

  dat <- monthRandom
  if (!beh %in% colnames(dat)) return(NULL)

  dat_sub <- dat[dat[[beh]] > 0, ]  
  if (nrow(dat_sub) == 0) return(NULL)

  dat_sub$tBeh <- dat_sub[[beh]] * 24
  coordinates(dat_sub) <- ~mean.lon + mean.lat
  proj4string(dat_sub) <- "+proj=longlat +datum=WGS84"

  dat_trans <- data.frame(spTransform(dat_sub, projection_NA))

  f <- as.formula(paste0("tBeh ~ s(coords.x1, coords.x2) + s(individ_id, bs='re')"))

  gam_mod <- mgcv::gam(f, data = dat_trans)

  rast_trans <- raster::projectRaster(rast.template, crs = projection_NA)
  xy <- as.data.frame(rast_trans, xy = TRUE)
  xy$individ_id <- dat_trans$individ_id[1]
  names(xy)[1:2] <- c("coords.x1", "coords.x2")

  preds <- predict(gam_mod, newdata = xy, exclude = "s(individ_id)", type = "response")

  pred_df <- data.frame(predictions = preds, coords.x1 = xy$coords.x1, coords.x2 = xy$coords.x2)
  coordinates(pred_df) <- ~ coords.x1 + coords.x2
  pred_raster <- rasterize(pred_df, rast_trans, "predictions")

  # Kriging
  sp_pts <- SpatialPointsDataFrame(coords = dat_trans[, c("coords.x1", "coords.x2")],
                                   data = dat_trans, proj4string = CRS(projection_NA))

  sp_pts$residLMM <- resid(gam_mod, type = "response")

  variogram_model <- vgm(model = "Sph", nugget = 0.1)
  v <- variogram(residLMM ~ 1, sp_pts)
  v_fit <- fit.variogram(v, model = variogram_model)

  rstPix <- as(rast_trans, "SpatialPixelsDataFrame")
  crs(rstPix) <- projection_NA

  krig_map <- krige(residLMM ~ 1, sp_pts, rstPix, model = v_fit)
  krig_raster <- raster(krig_map, layer = "var1.pred")

  final <- pred_raster + krig_raster

  return(final)
}

### Function for plotting where top % of data is ###

topEnergy<-function(df, percent, columnName) {

df$energy<-df[[columnName]]

# Remove missing values if any
df <- df[!is.na(df$energy), ]

# Sort descending by energy
df_sorted <- df[order(df$energy, decreasing = TRUE), ]

# Cumulative proportion of total energy
df_sorted$cum_energy <- cumsum(df_sorted$energy) / sum(df_sorted$energy)

# Keep cells contributing to 90% of total energy
df_10 <- df_sorted[df_sorted$cum_energy <= percent, ]

return(df_10)

}

## Function for plotting top % of data as contour lines ##

topEnergy_contour<-function(df, percent, columnName) {

# Change df into a xyz data frame for transforming into a raster
plot1<-df[,c("x", "y")] 
plot1$z<-df[[columnName]]

# Turn into raster
rast1<-raster::rasterFromXYZ(plot1)
v <- values(rast1)

# Rank cells by energy
ord <- order(v, decreasing = TRUE)
cum <- numeric(length(v))
cum[ord] <- cumsum(v[ord]) / sum(v)

# Put cumulative surface back into raster
r_cum <- rast1
values(r_cum) <- NA
values(r_cum)[!is.na(values(rast1))] <- cum
r_terra <- rast(r_cum)
contours <- as.contour(
  r_terra,
  levels = c(percent)
)

# Turn back into data frame
r_df<-as.data.frame(r_terra, xy=TRUE)

return(r_df)

}

### Function for extracting activity budget for a given model colony ###
extractBudget<-function(modelColonySub, speciesSub, monthSub) {

# List table with Ids & colony names #
idCatalogue<-read.csv("./results/tables/main/table1_idcatalogue.csv")

# Translate the model colony name from the text in Per's rasters to the ones in mine #
modelcoloniesTranslate<-data.frame(modelColony_Per=c("Kap_Hoegh","Bjoernoeya","Hornsund","Isfjorden", "Franz_Josef_Land", "Witless_Bay","Isle_of_May","Faroe_Islands", "Vestmannaeyjar","Papey","Grimsey","Runde_and_Aalesund"
,"Sklinna", "Roest", "Anda", "Hjelmsoeya","Hornoeya", "Skellig_Michael", "Little_Saltee" ,"Eynhallow", "Jan_Mayen","Breidafjordur" ,"Skjalfandi","Jarsteinen", "Alkefjellet",
"Sermilinnguaq", "Kippaku", "Isle_of_Canna", "Nuuk", "Langanes", "Kongsfjorden", "Anda", "Cape_Krutik", "Kara_Gate", "Russkaya_Gavan",
"Latrabjarg", "Cape_Gorodetskiy", "Coats_Island",  "Oranskie_Islands"),
modelColony_saved=c("kaphoegh","bjornoya","hornsund","isfjorden", "franzjosefland", "witlessbay","isleofmay","faroeislands", "vestmannaeyjar","papey","grimsey","rundeandaalesund"
,"sklinna", "rost", "anda", "hjelmsoya","hornoya", "skelligmichael", "littlesaltee" ,"eynhallow", "janmayen","breidafjordur" ,"skjalfandi","jarsteinen", "alkefjellet",
"sermilinnguaq", "kippaku", "isleofcanna", "nuuk", "langanes", "kongsfjorden", "anda", "capekrutik", "karagate", "russkayagavan",
"latrabjarg", "capegorodetskiy", "coatsisland",  "oranskieislands"),  
modelColony_database=c("Kap Höegh", "Bjørnøya", "Hornsund", "Isfjorden", "Franz Josef Land", "Witless Bay", "Isle of May", "Faroe Island", "Vestmannaeyjar", "Papey", "Grimsey", "Runde and Ålesund", "Sklinna", "Røst", "Anda",
"Hjemlsøya", "Hornøya", "Skellig Michael", "Little Saltee", "Eynhallow", "Jan Mayen", "Breidafjordur", "Skjalfandi", "Jarsteinen", "Alkefjellet",
"Sermilinnguaq", "Kippaku", "Isle of Canna", "Nuuk", "Langanes", "Kongsfjorden", "Anda", "Cape Krutik", "Kara Gate", "Russkaya Gavan",
"Latrabjarg", "Cape Gorodetskiy", "Coats Island",  "Oranskie Islands"))

# Find correct model colony translation
modelColMatch<-subset(modelcoloniesTranslate, modelColony_Per==modelColonySub)$modelColony_database[1]
modelSave<-subset(modelcoloniesTranslate, modelColony_Per==modelColonySub)$modelColony_saved[1]

if(speciesSub=="Northern fulmar" & modelSave=="jarsteinen") {
modelSave<-"karmoy"
}

# Make sure we were able to find a colony match #
if (length(modelColMatch)<1) stop(print("Error: no match with model colony name"))

# Transform species sub into how it's written in the saved files
speciesMatch<-gsub(" ", "", speciesSub)
speciesMatch<-gsub("-", "", speciesMatch)
speciesMatch<-gsub("ü", "u", speciesMatch)
speciesMatch<-gsub("'", "", speciesMatch)

# But then also filter for files
allBudgets<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp3/", full.names=TRUE)
allBudgets_colony<-allBudgets[grep(paste0("_", modelSave, "_"), allBudgets)] # Subset to colony of interest
allBudgets_monthly<-allBudgets_colony[grep("monthly", allBudgets_colony)] # Subset to monthly data
allBudgets_species<-allBudgets_monthly[grep(speciesMatch, allBudgets_monthly)]

# make sure we identified only one budget file #
if (!length(allBudgets_species) %in% c(1)) stop (print("Error: wrong number of budget files"))

# Now we open it #
budgets<-read.csv(allBudgets_species[1])

# Filter to month of interest #
budgets_sub<-subset(budgets, month==monthSub)

# Filter out birds with not enough days & sum per month #
birdMonth<-budgets_sub %>%
ungroup() %>%
dplyr::group_by(species, colony, month) %>%
dplyr::summarise(FlightHrs_mean=mean(tFlightMean), LandHrs_mean=mean(tLandMean), RestWaterHrs_mean=mean(tRestWaterMean), ActiveHrs_mean=mean(tActiveMean),
ForageHrs_mean=mean(tForageMean), .groups="drop")

return(birdMonth)

}

### Function for making randomized polygons. They are made in a new location and rotated. It's so I can compare values in these fake ones to real ones ###

# Function made by chat gpt so let's see how this goes

st_rotate <- function(x, angle_deg, center = st_centroid(st_union(x))) {
  angle <- angle_deg * pi/180
  M <- matrix(c(cos(angle), -sin(angle),
                sin(angle),  cos(angle)), nrow = 2, byrow = TRUE)

  ctr <- st_coordinates(center)[1, 1:2]

  g <- st_geometry(x)
  g_rot <- (g - ctr) * M + ctr

  st_set_geometry(x, g_rot)
}

make_random_poly <- function(patchNo, numberPolys, patchProj, availTrans) {

print("Generating random patches...")

# PatchProj is projected patch #

for (i in 1:numberPolys) {

repeat{

# Rotate polygon at random #
patchProj <- st_make_valid(patchProj)
#rotateRandom <- rotate.polygon(patchProj, angle=sample(0:360, 1))
rotateRandom <- st_rotate(patchProj, sample(0:360, 1))

# Find out extent of terra raster
e <- ext(availTrans)

xmin <- e$xmin
xmax <- e$xmax
ymin <- e$ymin
ymax <- e$ymax

# Determine maximum distances this patch can move by #
distancex<-xmax-xmin
distancey<-ymax-ymin

# Choose random amount to move patch by
dx <-sample(0:distancex, 1)
dy <-sample(0:distancey, 1)   # meters south

# Move polygon
poly_moved <- rotateRandom
st_geometry(poly_moved) <- st_geometry(rotateRandom) + c(dx, dy)
st_crs(poly_moved)<-projection_NA # Specify projection

# Extract values from our availability raster
vpoly <- vect(poly_moved)   # convert sf -> terra vector
vals <- terra::extract(availTrans, vpoly)
vals[is.na(vals)] <- 0 # to make my ifelse statement below easier

# all(vals == TRUE)

if (all(vals == 1)) {
break
} 

}

# Create an sf to return?

if (i==1) {

ids<-c(paste0(patchNo, "_observed"), paste0(patchNo, "_control_1"))
polys <- rbind(patchProj, poly_moved)

} else {

newId<-paste0(patchNo, "_control_", i)
ids<-append(ids, newId)
polys<-rbind(polys, poly_moved)

}

}

# Give the polygons new names #
polysFinal<-polys %>%
dplyr::mutate(patches=ids)

return(polysFinal)

}

### Extract LCT_air ####

# Here we extract air temp for some locations #

extract_air_temp<-function(data) {

# Determine location of air temp files
airTemp<-data.frame(Files=list.files("../data/temp/", full.names=TRUE))
airTemp$year<-sub("^[^_]+_([^_]+)_.*", "\\1", airTemp$Files)
airTemp$month<-sub(".*_(\\d+)\\.nc$", "\\1", airTemp$Files)

# First we determine unique months/year #
data$month<-as.numeric(substr(data$date, 6, 7))
data$year<-as.numeric(substr(data$date, 1, 4))
uniqueDates<-data %>%
dplyr::group_by(year) %>%
dplyr::count(month)

# now we extract this information row by row

dataSaveAll<-list()

print("Extracting data for date...")

for (j in 1:nrow(uniqueDates)) {

print(paste0(j, "/", nrow(uniqueDates)))

# subset the correct raster file #
airTempSub<-subset(airTemp, month==uniqueDates[j,]$month & year==uniqueDates[j,]$year)

# open the raster
#print(airTempSub$Files[1])
airTempRast<-terra::rast(airTempSub$Files[1])
#print("Opened raster")

# Change temp from kelvin to degrees
airTempRastCelcius <- airTempRast - 273.15

# Subset these months/years in the dataset
dataExtract<-subset(data, month==uniqueDates[j,]$month & year==uniqueDates[j,]$year)
coordinates<-data.frame(lon=dataExtract$mean.lon, lat=dataExtract$mean.lat)

# Extract corresponding air temp values
airTempValues<-terra::extract(airTempRastCelcius, coordinates)
colnames(airTempValues)<-c("row", "airTemp")

# Add this information to the main data frame
dataSave<-dataExtract %>%
dplyr::mutate(airTemp=airTempValues$airTemp)

dataSaveAll<-rbind(dataSaveAll, dataSave)

}

return(dataSaveAll)

}

#### LOOVC FUNCTION #####

conduct_loovc<-function(model_formula, data, responseVar, testType) { # Conduct loovc {