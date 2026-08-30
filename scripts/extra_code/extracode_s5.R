### Next step ####


# Begin by assigning season #
energySpeciesDf$season<-ifelse(energySpeciesDf$month %in% c(9, 10, 11), "Autumn", "Winter")
energySpeciesDf$season<-ifelse(energySpeciesDf$month %in% c(3, 4), "Spring", energySpeciesDf$season)

# Summarize the data by season #
energySeason<-energySpeciesDf %>%
 dplyr::group_by(season, index) %>%
 dplyr::summarise(months=n_distinct(month), energyTot=sum(energyPopkJ_mean), energyPerMonth=energyTot/months, birdsTot=sum(NoBirds_mean), ForageFish_tot=sum(ForageFish_g), ForageFishPerMonth=ForageFish_tot/months, 
 HerbZoo_tot=sum(HerbZoo_g), HerbZooperMonth=HerbZoo_tot/months, MacroZoo_tot=sum(MacroZoo_g), MacroZooPerMonth=MacroZoo_tot/months)
 
# Re-attach x & y coordinates (i am doing this as the decimals don't match very well in R)
index<-energySpeciesDf %>%
dplyr::group_by(index) %>%
dplyr::slice(1) %>%
dplyr::select(index, x, y)

energySeasonLox<-energySeason %>%
dplyr::inner_join(index, by=c("index")) %>%
dplyr::mutate(season=factor(season, levels=c("Spring", "Winter", "Autumn")))

#### Step 2: Highlight where hotspots occur ####

print("Step 2: Extracting areas of top energy expenditure...")

# Here I highlight where the top energy-demanding locations are in comparison to where birds are #

# Scaling functions
scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# First we scale values #
energySeasonLox$energyDemand_scaled<-scale01(energySeasonLox$energyTot)
energySeasonLox$fishDemand_scaled<-scale01(energySeasonLox$ForageFish_tot)
energySeasonLox$zooDemand_scaled<-scale01(energySeasonLox$MacroZoo_tot)

# Now we make a new data frame which will make it easier to compare them #
energy1<-energySeasonLox %>%
dplyr::select(x, y, season, energyDemand_scaled) %>%
dplyr::rename(scaledVal=energyDemand_scaled) %>%
dplyr::mutate(metric="Energy_demand")

energy2<-energySeasonLox %>%
dplyr::select(x, y, season, fishDemand_scaled) %>%
dplyr::rename(scaledVal=fishDemand_scaled) %>%
dplyr::mutate(metric="Fish_demand")

energy3<-energySeasonLox %>%
dplyr::select(x, y, season, zooDemand_scaled) %>%
dplyr::rename(scaledVal=zooDemand_scaled) %>%
dplyr::mutate(metric="Zoo_demand")

energyPlot<-rbind(energy1, energy2, energy3)

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# make a plot for viewing
energyDemandPlots<-ggplot() +
 geom_tile(data=energyPlot, aes(x=x, y=y, fill=scaledVal)) +
 scale_fill_gradient2(
  low = "#2c7bb6",   # blue
  mid = "white",
  high = "#d7191c",  # red
  midpoint = 0.5
) + 
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(energyPlot$x), max(energyPlot$x)), ylim=c(min(energyPlot$y) + 20, max(energyPlot$y))) +
 theme_minimal() +
 labs(fill="") +
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")) + metric, nrow=3)
 
 pdf("./results/figures/speciesMaps/energy/demands.pdf", width=10)
 plot(energyDemandPlots)
 dev.off()
 
 # make a plot for viewing
energyDemandPlots2<-ggplot() +
 geom_tile(data=energySeasonLox, aes(x=x, y=y, fill=energyPerMonth)) +
 #scale_fill_gradientn(
 # colours = c("#0d0887", "#2c7bb6", "#7ad151", "#fdae61", "#d73027")
#) +
 scale_fill_gradientn(
  colours = c("#08306b", "#2171b5", "#ffff66", "#fdae61", "#d7191c")
) +
 #scale_fill_gradient2(
 # low = "#2c7bb6",   # blue
 # mid = "#FD9040",
 # high = "#d7191c",  # red
 # midpoint = 279440275/2
#) + 
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(energyPlot$x), max(energyPlot$x)), ylim=c(min(energyPlot$y) + 20, max(energyPlot$y))) +
 theme_minimal() +
 labs(fill="", x="", y="") +
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1)
 
 pdf("./results/figures/speciesMaps/energy/demands2.pdf", width=10, height=8)
 plot(energyDemandPlots2)
 dev.off()

# Plot top 25% of energy in both cases #

seasonsUnique<-unique(energySeasonLox$season)

# Make some empty lists to save information in #
topenergy<-list()
topbirds<-list()
topfish<-list()
topmacrozoo<-list()

for (i in 1:length(seasonsUnique)) {

# Subset to season i
seasonSub<-seasonsUnique[i]

# Calculate top % of data
top25_energy<-topEnergy(subset(energySeasonLox, season==seasonSub), 0.25, 'energyTot')
top25_birds<-topEnergy(subset(energySeasonLox, season==seasonSub), 0.10, 'birdsTot')
top25_fish<-topEnergy(subset(energySeasonLox, season==seasonSub), 0.10, 'ForageFish_tot')
top25_macrozoo<-topEnergy(subset(energySeasonLox, season==seasonSub), 0.10, 'MacroZoo_tot')

# Add seasonal info
top25_energy$season<-seasonSub
top25_birds$season<-seasonSub
top25_fish$season<-seasonSub
top25_macrozoo$season<-seasonSub

# Save results
topenergy<-rbind(topenergy, top25_energy)
topbirds<-rbind(topbirds, top25_birds)
topfish<-rbind(topfish, top25_fish)
topmacrozoo<-rbind(topmacrozoo, top25_macrozoo)

}

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# make a plot for viewing
hotspots<-ggplot() +
 #geom_tile(data=topenergy, aes(x=x, y=y, fill="10%_energy")) +
 #geom_tile(data=topbirds, aes(x=x, y=y, fill="10%_birds"), alpha=0.5) +
 geom_tile(data=topfish, aes(x=x, y=y, fill="10%_fish"), alpha=0.5) +
 geom_tile(data=topmacrozoo, aes(x=x, y=y, fill="10%_zoo"), alpha=0.5) +
 scale_fill_manual(values=c("#0072B2", "#E69F00", "darkred", "darkgreen")) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(topenergy$x), max(topenergy$x) + 25), ylim=c(min(topenergy$y), max(topenergy$y) + 10)) +
 theme_minimal() +
 ggtitle("Top 10% of data") +
 labs(fill="") +
 xlab("") +
 ylab("") + 
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1)
 
 pdf("./results/figures/speciesMaps/energy/hotspots.pdf", width=10)
 plot(hotspots)
 dev.off()
 
hotspots2<-ggplot() +
 geom_tile(data=topenergy, aes(x=x, y=y, fill="10%_energy")) +
 #geom_tile(data=topbirds, aes(x=x, y=y, fill="10%_birds"), alpha=0.5) +
 #geom_tile(data=topfish, aes(x=x, y=y, fill="10%_fish"), alpha=0.5) +
 #geom_tile(data=topmacrozoo, aes(x=x, y=y, fill="10%_zoo"), alpha=0.5) +
 scale_fill_manual(values=c("#0072B2", "#E69F00", "darkred", "darkgreen")) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(topenergy$x), max(topenergy$x) + 25), ylim=c(min(topenergy$y), max(topenergy$y) + 10)) +
 theme_minimal() +
 ggtitle("Top 25% of data") +
 labs(fill="") +
 xlab("") +
 ylab("") + 
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1)
 
 pdf("./results/figures/speciesMaps/energy/hotspots2.pdf", width=10)
 plot(hotspots2)
 dev.off()
 
 # North Atlantic extent
bathy <- getNOAA.bathy(
  lon1 = -80, lon2 =80,
  lat1 = 40,  lat2 = 85,
  resolution = 50
)

bathy_df <- fortify.bathy(bathy)

land <- ne_countries(scale = "medium", returnclass = "sf")

depthvshotspot<-ggplot() +
  geom_raster(
    data = bathy_df,
    aes(x = x, y = y, fill = z)
  ) +
  geom_contour(
    data = bathy_df,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -500, -1000, -3000),
    colour = "grey30",
    linewidth = 0.2
  ) +
  geom_sf(data = land, fill = "grey80", colour = "grey40", linewidth = 0.2) +
  geom_tile(data=topenergy, aes(x=x, y=y), fill="darkgreen", alpha=.4) +
  coord_sf(crs=projection_84, xlim=c(min(topenergy$x), max(topenergy$x) + 25), ylim=c(min(topenergy$y), max(topenergy$y) + 10))  +
  scale_fill_gradientn(
    colours = c("navy", "steelblue", "lightblue", "white"),
    name = "Depth / elevation (m)"
  ) +
  theme_minimal() +
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1) +
 theme(legend.position="bottom") +
 xlab("") +
 ylab("")
 
 pdf("./results/figures/speciesMaps/energy/depthvshotspot.pdf", width=10)
 plot(depthvshotspot)
 dev.off()
 
# Now we will make a plot with SST fronts #

# Determine location of SST data #
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

# Determine location of chla data #
chlMonthly<-list.files("./data/chl", full.names=TRUE)
chlMap<-raster::stack(chlMonthly) # Stack all layers (lol)
chlKey<-data.frame(layerNo=1:nlayers(chlMap)) # Make a key with dates so sstMap is easy to access
# Set start and end dates
start_date <- as.Date("2012-01-01")
end_date   <- as.Date("2022-12-01")
date_seq <- seq(from = start_date, to = end_date, by = "month") # Generate sequence for the 1st of each month
# Annotate key
chlKey$dates<-date_seq
chlKey$Year<-substr(chlKey$dates, 1, 4)
chlKey$Month<-as.numeric(substr(chlKey$dates, 6, 7))

# Now we do the same but for NEMO forecasts
chlNEMO<-list.files("./data/chl_copernicus/", full.names=TRUE)
chlMap_nemo<-raster::stack(chlNEMO) # Stack all layers (lol)
chlKey2<-data.frame(layerNo=1:nlayers(chlMap_nemo)) # Make a key with dates so sstMap is easy to access
# Set start and end dates
start_date <- as.Date("2012-01-01")
end_date   <- as.Date("2022-12-01")
date_seq <- seq(from = start_date, to = end_date, by = "month") # Generate sequence for the 1st of each month
# Annotate key
chlKey2$dates<-date_seq
chlKey2$Year<-substr(chlKey2$dates, 1, 4)
chlKey2$Month<-as.numeric(substr(chlKey2$dates, 6, 7))

# Make a list to save result in
seasonalMeanAll_SST<-list()
seasonalMeanAll_Chl<-list()
seasonalMeanAll_Chl2<-list()

for (i in 1:length(seasonsUnique)) {

print(paste0("Averaging sst for season ", i, "/", length(seasonsUnique)))

# subset to season i
seasonSub<-seasonsUnique[[i]]

# Figure out which months that is #
seasonSubLox<-subset(energySpeciesDf, season %in% c(seasonSub))
monthsSelect<-unique(seasonSubLox$month)

monthlyMeanSST<-list()
monthlyMeanChl<-list()
monthlyMeanChl2<-list()

for (j in 1:length(monthsSelect)) {

print(paste0("Month", j))

# Subset sst from months... J
sstSub<-subset(sstKey, Month==monthsSelect[j])
chlSub<-subset(chlKey, Month==monthsSelect[j])
chlSub2<-subset(chlKey2, Month==monthsSelect[j])

# Stack & summarise
sstStack<-subset(sstMap, sstSub$layerNo)
sstMean<-overlay(sstStack, fun="mean", na.rm=TRUE)

chlStack<-subset(chlMap, chlSub$layerNo)
chlProj<-projectRaster(chlStack, sstMean) # make coarser res
chlMean<-overlay(chlProj, fun="mean", na.rm=TRUE)

chlStack2<-subset(chlMap_nemo, chlSub2$layerNo)
chlProj2<-projectRaster(chlStack2, sstMean) # make coarser res
chlMean2<-overlay(chlProj2, fun="mean", na.rm=TRUE)

# Save in a list 
monthlyMeanSST<-append(monthlyMeanSST, sstMean)
monthlyMeanChl<-append(monthlyMeanChl, chlMean)
monthlyMeanChl2<-append(monthlyMeanChl2, chlMean2)

}

# Stack & average seasonal Mean #
seasonalStackSST<-stack(monthlyMeanSST)
seasonalMeanSST<-overlay(seasonalStackSST, fun="mean", na.rm=TRUE)

seasonalStackChl<-stack(monthlyMeanChl)
seasonalMeanChl<-overlay(seasonalStackChl ,fun="mean", na.rm=TRUE)

seasonalStackChl2<-stack(monthlyMeanChl2)
seasonalMeanChl2<-overlay(seasonalStackChl2 ,fun="mean", na.rm=TRUE)

# Turn into data frame
seasonalMeanDf<-as.data.frame(seasonalMeanSST, xy=TRUE)
seasonalMeanDf$season<-seasonSub
seasonalMeanAll_SST<-rbind(seasonalMeanAll_SST, seasonalMeanDf)

seasonalMeanDf2<-as.data.frame(seasonalMeanChl, xy=TRUE)
seasonalMeanDf2$season<-seasonSub
seasonalMeanAll_Chl<-rbind(seasonalMeanAll_Chl, seasonalMeanDf2)

seasonalMeanDf3<-as.data.frame(seasonalMeanChl2, xy=TRUE)
seasonalMeanDf3$season<-seasonSub
seasonalMeanAll_Chl2<-rbind(seasonalMeanAll_Chl2, seasonalMeanDf3)

}

fronts<-ggplot() +
  geom_tile(data=seasonalMeanAll_SST, aes(x=x, y=y, fill=layer)) + 
  geom_sf(data = land, fill = "grey80", colour = "grey40", linewidth = 0.2) +
  geom_tile(data=topenergy, aes(x=x, y=y), fill="darkgreen", alpha=.4) +
  coord_sf(crs=projection_84, xlim=c(min(topenergy$x), max(topenergy$x) + 25), ylim=c(min(topenergy$y), max(topenergy$y) + 10))  +
  scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), limits=c(-2.5, 26)) +
  theme_minimal() +
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1) +
 theme(legend.position="bottom") +
 xlab("") +
 ylab("")
 
 pdf("./results/figures/speciesMaps/energy/frontsvshotspot.pdf", width=10)
 plot(fronts)
 dev.off()
 
 chlMap<-ggplot() +
  geom_tile(data=seasonalMeanAll_Chl, aes(x=x, y=y, fill=layer)) + 
  geom_sf(data = land, fill = "grey80", colour = "grey40", linewidth = 0.2) +
  geom_tile(data=topenergy, aes(x=x, y=y), fill="darkgreen", alpha=.4) +
  coord_sf(crs=projection_84, xlim=c(min(topenergy$x), max(topenergy$x) + 25), ylim=c(min(topenergy$y), max(topenergy$y) + 10))  +
  scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  theme_minimal() +
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1) +
 theme(legend.position="bottom") +
 xlab("") +
 ylab("")
 
 pdf("./results/figures/speciesMaps/energy/chlvshotspot.pdf", width=10)
 plot(chlMap)
 dev.off()
 
 chlMap2<-ggplot() +
  geom_tile(data=seasonalMeanAll_Chl2, aes(x=x, y=y, fill=layer^(1/4))) + 
  geom_sf(data = land, fill = "grey80", colour = "grey40", linewidth = 0.2) +
  #geom_tile(data=topenergy, aes(x=x, y=y), fill="darkgreen", alpha=.2) +
  coord_sf(crs=projection_84, xlim=c(min(topenergy$x), max(topenergy$x) + 25), ylim=c(min(topenergy$y), max(topenergy$y) + 10))  +
  scale_fill_gradientn(
  colours = c("#2c7bb6", "#00a6ca", "#00ccbc", "#90eb9d",
              "#ffff8c", "#f9d057", "#f29e2e", "#e76818", "#d7191c"),
  name = "NPP (mg/m3/day)"
) +
  #scale_fill_gradientn('NPP (mg/m3/day)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  theme_minimal() +
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1) +
 theme(legend.position="bottom") +
 xlab("") +
 ylab("")
 
 pdf("./results/figures/speciesMaps/energy/chlvshotspot2.pdf", width=10)
 plot(chlMap2)
 dev.off()

### Step 3: Convert into polygons ####

print("Step 3: Converting into polygons...")

# Need to turn these files into polygons somehow (e.g. polygon 1, polygon 2 etc.) #

# First we turn into a raster #

seasonalHotspots<-list() # for energy demand
seasonalHotspots_fish<-list()
seasonalHotspots_zoo<-list()

seasons<-unique(energySeasonLox$season)

for (i in 1:length(seasons)) {

print(i)
print(seasons[[i]])

# Subset to season i
seasonalMap<-energySeasonLox %>%
dplyr::filter(season==seasons[[i]]) %>%
ungroup() %>%
dplyr::select(index, x, y)

# Trasnform dataset for making into raster
topenergy_season<-topenergy %>%
dplyr::filter(season==seasons[[i]]) %>%
dplyr::mutate(topenergy=1) %>%
dplyr::select(index, topenergy, season) %>%
dplyr::full_join(seasonalMap, by=c("index")) %>%
replace_na(list(topenergy=0)) %>%
arrange(index) %>%
ungroup() %>%
dplyr::select(x, y, topenergy) %>%
dplyr::rename(z=topenergy)

topfish_season<-topfish %>%
dplyr::filter(season==seasons[[i]]) %>%
dplyr::mutate(topfish=1) %>%
dplyr::select(index, topfish, season) %>%
dplyr::full_join(seasonalMap, by=c("index")) %>%
replace_na(list(topfish=0)) %>%
arrange(index) %>%
ungroup() %>%
dplyr::select(x, y, topfish) %>%
dplyr::rename(z=topfish)

topzoo_season<-topmacrozoo %>%
dplyr::filter(season==seasons[[i]]) %>%
dplyr::mutate(topzoo=1) %>%
dplyr::select(index, topzoo, season) %>%
dplyr::full_join(seasonalMap, by=c("index")) %>%
replace_na(list(topzoo=0)) %>%
arrange(index) %>%
ungroup() %>%
dplyr::select(x, y, topzoo) %>%
dplyr::rename(z=topzoo)

# Turn into raster
energyRast<-rasterFromXYZ(topenergy_season) 
energyRast2<-rasterFromXYZ(topfish_season) 
energyRast3<-rasterFromXYZ(topzoo_season) 

# Convert to terra object
r_terra <- rast(energyRast)
r_terra2 <- rast(energyRast2)
r_terra3 <- rast(energyRast3)

# Use a distance threshold to seperate 'hotspots' 
d <- 1 

# Convert hotspot cells into 1, others into NA
r1 <- r_terra == 1
r1 <- mask(r1, r1, maskvalues = 0)   
r2 <- r_terra2 == 1
r2 <- mask(r2, r2, maskvalues = 0) 
r3 <- r_terra3 == 1
r3 <- mask(r3, r3, maskvalues = 0) 

# distance from every cell to the nearest '1' cell
dist_to_1 <- distance(r1)
dist_to_2 <- distance(r2)
dist_to_3 <- distance(r3)

# expand the blobs by distance d (a morphological "buffer" in raster space)
expanded <- dist_to_1 <= d
expanded <- mask(expanded, expanded, maskvalues = 0)

expanded2 <- dist_to_2 <= d
expanded2 <- mask(expanded2, expanded2, maskvalues = 0)

expanded3 <- dist_to_3 <= d
expanded3 <- mask(expanded3, expanded3, maskvalues = 0)

# label connected components in expanded space
grp <- patches(expanded, directions = 8)
grp2 <- patches(expanded2, directions = 8)
grp3 <- patches(expanded3, directions = 8)

# assign each original '1' cell to a group (mask keeps only original 1-cells)
grp_on_ones <- mask(grp, r1)
grp_on_ones2 <- mask(grp2, r2)
grp_on_ones3 <- mask(grp3, r3)

# polygonize groups
pol <- as.polygons(grp_on_ones, dissolve = TRUE)
pol <- pol[!is.na(values(pol)[,1]), ]

pol2 <- as.polygons(grp_on_ones2, dissolve = TRUE)
pol2 <- pol2[!is.na(values(pol2)[,1]), ]

pol3 <- as.polygons(grp_on_ones3, dissolve = TRUE)
pol3 <- pol3[!is.na(values(pol3)[,1]), ]

# Convert to sf
pol_sf <- st_as_sf(pol)
pol_sf2 <- st_as_sf(pol2)
pol_sf3 <- st_as_sf(pol3)

# Save results
seasonalHotspots[[i]]<-pol_sf
seasonalHotspots_fish[[i]]<-pol_sf2
seasonalHotspots_zoo[[i]]<-pol_sf3

}

# Now we extract the amount of energy in each polygon, which we will turn into kJ of food per km2 and then g of fish and g of zooplankton #

#### Step 4b: convert rest of the ocean into a polygon ####

print("Step 4: making control polygons...")

# I want to place them randomly in other sea-areas & rotate them randomly too #

# First we will make an 'avaiable' area raster based on a random SST layer #
sstMonthly<-list.files("./data/sstPopMaps", full.names=TRUE)
sstMaps<-sstMonthly[grepl("tif", sstMonthly)] # Choose the tif file
sstMap<-raster::stack(sstMaps[1]) # Stack all layers (lol)
sstMapSub<-subset(sstMap, 1)

# availability raster: cells that exist in your seasonal raster grid
avail_r <- !is.na(sstMapSub)          # TRUE where cell exists
avail_r <- mask(avail_r, avail_r, maskvalues = 0)
avail_r<-rast(avail_r)

# Now we will loop through my seasons & patches to make null patches #

seasonalHotspots_control<-list()

for (i in 1:length(seasonalHotspots)) {

print(paste0("Making control patches for season ", i))

# Subset to season i
patchesSub<-seasonalHotspots[[i]]

# Now we will loop through the seperate patches? #
patches<-nrow(patchesSub)

sf_all<-c()

for (j in 1:patches) {

# Subset to patch j
poly1 <- patchesSub[j, ]

sf_all <- rbind(sf_all, poly1)

}

# remove this polygon from available space for control polygons
poly1_v  <- vect(sf_all)  # Turn into a vector             
patch_r <- rasterize(poly1_v, avail_r, field = 1, background = NA) # Rasterize patch to the same grid as avail_r
avail_no_patch <- mask(avail_r, patch_r, inverse = TRUE) # Update availability: set cells covered by patch to NA

# Replace seasonal hotspot patches with these real ones with the control ones
seasonalHotspots_control[[i]]<-avail_no_patch

}

#### Step 5: Extract characteristics of polygons ###

print("Step 5: Extracting characteristics of polygons...")

# First we start with... species/population composition #

# I will make a key identifying the months that need to be extracted for specific seasons #
seasonsMonth<-data.frame(month=c(9, 10, 11, 12, 1, 2, 3, 4), season=c("Autumn", "Autumn", "Autumn", "Winter", "Winter", "Winter", "Spring", "Spring"))

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

# We will have to re-calculate diets here #
seabirdDiet<-data.frame(functionalGroups=c("Forage fish", "Macrozooplankton", "Herbivorous zooplankton", "Other"))
seabirdDiet$AE<-0.75 # Assimilation efficiency -> assume same for all seabird species & prey groups for now
seabirdDiet$energyDensity<-c(4.5, 4, 3.1, 2) # (kJ per g wet mass) #
# Energy density for forage fish based on this paper for Capelin: https://d1wqtxts1xzle7.cloudfront.net/46124708/Assimilation_efficiency_of_adult_Kittiwa20160601-4336-10nyk6f-libre.pdf?1464771742=&response-content-disposition=inline%3B+filename%3DAssimilation_efficiency_of_adult_Kittiwa.pdf&Expires=1776686451&Signature=FiV260FL8sA6XfT611m3V764hctijC01w2o6Cx9yHcrlMSq-hw3PRJQT-duIKO3~V9lh63q3HFgDZO2kv2tTPPH54yr9~xc7cSAmUGZTisvbtVaO5i5Dd3UUGIHqehy~dG4AyzRfAQYGZNFpN~s2ZgO--emdvI~GM88W3jbDydiJf-2EEbz3bcBav70q9l~4AdDjWEXb-2UfqHWx9-tVWGDvJTtJd4hx6hpYOtO1nDklqsIXrmcqDydleX1lijaaFGHDoUJENqx3U12SVPfzjL4mF5-9IwM9wDc7A4GMRiFe2UDN0vQHQX0rP4u97wOOYWq-IyvckPnI49hu9pVL5w__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA
# Energy density for herb. zoo. based on this paper for Calanus: https://www.sciencedirect.com/science/article/pii/S0079661116302105#s0010
# Energy density for macro. zoo. based on this paper for amphipods: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190345

# Determine proportion of each class in diet #
seabirdSpecies<-data.frame(species=rep(c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot"), each=4))
seabirdSpecies$functionalGroups<-rep(c("Forage fish", "Macrozooplankton", "Herbivorous zooplankton", "Other"), 6)
seabirdSpecies$propGroup<-c(0.94, 0.06, 0, 0, # BLK https://www.int-res.com/articles/meps2007/349/m349p269.pdf (I made an average of the Other category)
0.93, 0, 0, 0.07, # NoFu # file:///C:/Users/caitlin.frankish/Downloads/s002270050613.pdf
0.66, 0.02, 0, 0.41, # AtPu : https://onlinelibrary.wiley.com/doi/pdf/10.1111/ibi.12272
0.12, 0.87, 0.01, 0,# LiAu: https://link.springer.com/article/10.1007/s00300-013-1379-4/tables/2
1, 0, 0, 0,# CoGu # https://www.tandfonline.com/doi/epdf/10.1080/17451000802279636?needAccess=true
0.39, 0.55, 0, 0.06) #BrGu: https://cdnsciencepub.com/doi/epdf/10.1139/cjz-2021-0120 # But there are 6% other species so I will jus attribute this to the two classes  

# Attach energy density values #
seabirdSpeciesDiet<-seabirdSpecies %>%
dplyr::left_join(seabirdDiet, by=c("functionalGroups"))

# Loop through each season
seasonsRes<-list()
seasonsRes2<-list()
seasonsRes3<-list()

# Determine unique list of seasns
seasons<-unique(seasonsMonth$season)

for (i in 1:length(seasons)) {

# Subset to season i 
seasonSub<-seasons[i]

# Determine where pop energy demand rasters are located #
energyMaps<-list.files("tmp5/", full.names=TRUE)
energyMaps_v1<-energyMaps[grep("v1", energyMaps)] 

# Join polygons together #
seasonalHotspots_sub<-rbind(seasonalHotspots[[i]])
seasonalHotspots_fish_sub<-rbind(seasonalHotspots_fish[[i]])
seasonalHotspots_zoo_sub<-rbind(seasonalHotspots_zoo[[i]])
#remainingLandscape<-seasonalHotspots_control[[i]]

# list to save values in
birdComp_season<-list()
birdComp_season2<-list()
birdComp_season3<-list()

for (j in 1:length(energyMaps_v1)) {

#for (j in 1:10) {

print(paste0("Extracting from map ", j, "/", length(energyMaps_v1)))

# open energy demand map J
energyMap_j<-fread(energyMaps_v1[j])

# Subset to correct months #
monthSub<-subset(seasonsMonth, season==seasonSub)
energyMap_j_season<-subset(energyMap_j, month %in% unique(monthSub$month))

# Subset to species of interest
seabirdSpeciesDiet_sub<-subset(seabirdSpeciesDiet, species==energyMap_j_season$species[1])

# Find out area size of every cell
energySpatialSub<-subset(energyMap_j, month==9) # subset to a random month
energyRast<-energySpatialSub %>%
dplyr::select(x, y, energyPopkJ_mean) %>%
rename(z=energyPopkJ_mean)
energyRast<-rasterFromXYZ(energyRast) # Turn into a raster
area_raster <- area(energyRast) # calculate area of every cell
area_raster_df<-as.data.frame(area_raster, xy=TRUE) # Turn into a data frame
colnames(area_raster_df)<-c("x", "y", "areaKm2") # Change col names

# Join this to main dataset? 
energyMap_j_season<-energyMap_j_season %>%
dplyr::left_join(area_raster_df, by=c("x", "y"))

# Determine monthly prey requirements #
fishInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Forage fish")
energyMap_j_season$ForageFish_g<-((energyMap_j_season$energyPopkJ_mean/(fishInfo$AE[1]*fishInfo$energyDensity))*fishInfo$propGroup)/energyMap_j_season$areaKm2

herbZooInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Herbivorous zooplankton")
energyMap_j_season$HerbZoo_g<-((energyMap_j_season$energyPopkJ_mean/(herbZooInfo$AE[1]*herbZooInfo$energyDensity))*herbZooInfo$propGroup)/energyMap_j_season$areaKm2

macroZooInfo<-subset(seabirdSpeciesDiet_sub, functionalGroups=="Macrozooplankton")
energyMap_j_season$MacroZoo_g<-((energyMap_j_season$energyPopkJ_mean/(macroZooInfo$AE[1]*macroZooInfo$energyDensity))*macroZooInfo$propGroup)/energyMap_j_season$areaKm2

# Add up energy by cell
energyMap_j_sum<-energyMap_j_season %>%
dplyr::group_by(species, colony, x, y) %>%
dplyr::summarise(totBirds=sum(NoBirds_mean, na.rm=TRUE), totEnergy=sum(energyPopkJ_mean, na.rm=TRUE), meanSST=mean(sst, na.rm=TRUE), totFish=sum(ForageFish_g, na.rm=TRUE), totMacroZoo=sum(MacroZoo_g, na.rm=TRUE))

# Make a density raster #
density<-energyMap_j_sum %>%
ungroup() %>%
dplyr::rename(z=totBirds) %>%
dplyr::select(x, y, z)

densityRast<-rasterFromXYZ(density)

# Make an energy raster #
energy<-energyMap_j_sum %>%
ungroup() %>%
dplyr::rename(z=totEnergy) %>%
dplyr::select(x, y, z)

energyRast<-rasterFromXYZ(energy)

# Make a fish raster
fish<-energyMap_j_sum %>%
ungroup() %>%
dplyr::rename(z=totFish) %>%
dplyr::select(x, y, z)

fishRast<-rasterFromXYZ(fish)

# Make a zoo raster 
zoo<-energyMap_j_sum %>%
ungroup() %>%
dplyr::rename(z=totMacroZoo) %>%
dplyr::select(x, y, z)

zooRast<-rasterFromXYZ(zoo)

# Make an sst raster #
sst<-energyMap_j_sum %>%
ungroup() %>%
dplyr::rename(z=meanSST) %>%
dplyr::select(x, y, z)

sstRast<-rasterFromXYZ(sst)

# Extract values from all of these #
poly_sp <- as(seasonalHotspots_sub, "Spatial")
valsDensity <- raster::extract(densityRast, poly_sp)
valsDensity_sumbyid<-sapply(valsDensity, sum, na.rm = TRUE)

valsEnergy <- raster::extract(energyRast, poly_sp)
valsEnergy_sumbyid<-sapply(valsEnergy, sum, na.rm = TRUE)

poly_sp2 <- as(seasonalHotspots_fish_sub, "Spatial")
valsFish <- raster::extract(fishRast, poly_sp2)
valsFish_sumbyid<-sapply(valsFish, sum, na.rm = TRUE)

poly_sp3 <- as(seasonalHotspots_zoo_sub, "Spatial")
valsZoo<- raster::extract(zooRast, poly_sp3)
valsZoo_sumbyid<-sapply(valsZoo, sum, na.rm = TRUE)

valsSST <- raster::extract(sstRast, poly_sp)
valsSSTsumbyid<-sapply(valsSST, mean, na.rm = TRUE)

# Save these results

polygonComp<-data.frame(polygonID=1:nrow(seasonalHotspots[[i]]), totBirds=valsDensity_sumbyid, totEnergy=valsEnergy_sumbyid)
polygonComp$species<-energyMap_j_sum$species[1]
polygonComp$colony<-energyMap_j_sum$colony[1]
polygonComp$season<-seasons[[i]]

polygonComp2<-data.frame(polygonID=1:nrow(seasonalHotspots_fish[[i]]), totFish=valsFish_sumbyid)
polygonComp2$species<-energyMap_j_sum$species[1]
polygonComp2$colony<-energyMap_j_sum$colony[1]
polygonComp2$season<-seasons[[i]]

polygonComp3<-data.frame(polygonID=1:nrow(seasonalHotspots_zoo[[i]]), totZoo=valsZoo_sumbyid)
polygonComp3$species<-energyMap_j_sum$species[1]
polygonComp3$colony<-energyMap_j_sum$colony[1]
polygonComp3$season<-seasons[[i]]

# Save results
birdComp_season<-rbind(birdComp_season, polygonComp)
birdComp_season2<-rbind(birdComp_season2, polygonComp2)
birdComp_season3<-rbind(birdComp_season3, polygonComp3)

}

PropEnergyContribution<-birdComp_season %>%
ungroup() %>%
dplyr::group_by(season, polygonID) %>%
dplyr::mutate(totEnergyAll=sum(totEnergy)) %>%
dplyr::ungroup() %>%
dplyr::mutate(propPolygon=totEnergy/totEnergyAll) %>%
dplyr::select(season, polygonID, species, colony, propPolygon) %>%
arrange(polygonID, desc(propPolygon)) 

PropFishContribution<-birdComp_season2 %>%
ungroup() %>%
dplyr::group_by(season, polygonID) %>%
dplyr::mutate(totFishAll=sum(totFish)) %>%
dplyr::ungroup() %>%
dplyr::mutate(propPolygon=totFish/totFishAll) %>%
dplyr::select(season, polygonID, species, colony, propPolygon) %>%
arrange(polygonID, desc(propPolygon))

PropZooContribution<-birdComp_season3 %>%
ungroup() %>%
dplyr::group_by(season, polygonID) %>%
dplyr::mutate(totZooAll=sum(totZoo)) %>%
dplyr::ungroup() %>%
dplyr::mutate(propPolygon=totZoo/totZooAll) %>%
dplyr::select(season, polygonID, species, colony, propPolygon) %>%
arrange(polygonID, desc(propPolygon))

# Save results
seasonsRes<-rbind(seasonsRes, PropEnergyContribution)
seasonsRes2<-rbind(seasonsRes2, PropFishContribution)
seasonsRes3<-rbind(seasonsRes3, PropZooContribution)

}

#### Plot pie circles of species proportions at different hotspots ####

seasonRes<-read.csv("./results/tables/energyhotspots.csv")

# Summarise proportion of energy demand by species at different polygons #

propHotspots<-seasonRes %>%
ungroup() %>%
dplyr::group_by(season, polygonID, species) %>%
dplyr::summarise(propPolygonSum=sum(propPolygon)) %>%
arrange(season, polygonID, desc(propPolygonSum)) %>%
rename(patches=polygonID)

seasonalHotspots[[i]]<-pol_sf

# Now we extract hotspots and plot by polygonID? #
seasonalHotspots_autumn<-rbind(seasonalHotspots[[1]])
seasonalHotspots_autumn$patches<-1:nrow(seasonalHotspots_autumn)
seasonalHotspots_autumn<-subset(seasonalHotspots_autumn, patches %in% c(1, 4, 5, 7))
#seasonalHotspots_autumn<-subset(rbind(seasonalHotspots[[1]]), patches %in% c(1, 4, 5, 10)) # 10 becomes 7

seasonalHotspots_spring<-rbind(seasonalHotspots[[2]])
seasonalHotspots_spring$patches<-1:nrow(seasonalHotspots_spring)
seasonalHotspots_spring<-subset(rbind(seasonalHotspots[[2]]), patches %in% c(1, 3, 5))

seasonalHotspots_winter<-rbind(seasonalHotspots[[3]])
seasonalHotspots_winter$patches<-1:nrow(seasonalHotspots_winter)

seasonalHotspots_autumn$season<-"Autumn"
seasonalHotspots_spring$season<-"Spring"
seasonalHotspots_winter$season<-"Winter"

allSeasonalHotspots<-rbind(seasonalHotspots_autumn, seasonalHotspots_spring, seasonalHotspots_winter)

st_crs(allSeasonalHotspots)<-projection_84

# Now I join the centroid of these polygons with my other dataset for plotting #

# 1. Get centroid of each polygon
patch_centroids <- allSeasonalHotspots %>%
group_by(season, patches) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  mutate(geometry = st_point_on_surface(geometry)) %>%   # or st_centroid()
  mutate(
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(season, patches, X, Y)
  
# 2. Join centroid coordinates to your bird proportions data
bird_pies <- propHotspots %>%
  left_join(patch_centroids, by = c("season", "patches"))
  
# Trasnform into correct format for scatterpie
pie_df <- bird_pies %>%
  tidyr::pivot_wider(
    id_cols = c(season, patches, X, Y),
    names_from = species,
    values_from = propPolygonSum,
    values_fill = 0
  ) %>%
  dplyr::filter(!is.na(X)) %>%
  ungroup() %>%
  droplevels()

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

names(pie_df) <- c(
  "season", "patches", "x", "y",
  "atlantic_puffin",
  "brunnichs_guillemot",
  "common_guillemot",
  "black_legged_kittiwake",
  "northern_fulmar",
  "little_auk"
)

# make a plot for viewing
energyComposition<-ggplot() +
 geom_sf(data=allSeasonalHotspots, fill="lightblue", alpha=0.5) +
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1) + 
 scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
 geom_sf(data=coast) +
 geom_scatterpie(data=pie_df, aes(x=x ,y=y, r=4), cols=c("black_legged_kittiwake", "northern_fulmar", "atlantic_puffin", "little_auk", "common_guillemot", "brunnichs_guillemot")) +
 coord_sf(crs=projection_84, xlim=c(min(topenergy$x), max(topenergy$x) + 25), ylim=c(min(topenergy$y), max(topenergy$y) + 10)) +
 theme_minimal() +
 ggtitle("Top 25% of data") +
 labs(fill="") +
 xlab("") +
 ylab("") + 
 facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1) +
 theme(legend.position="bottom")
 
 pdf("./results/figures/speciesMaps/energy/hotspots3.pdf", width=10)
 plot(energyComposition)
 dev.off()
 
#### How much food does this represent? ####

# Main functional groups, AE and energy content in Kj.g of wet weight #
seabirdDiet<-data.frame(functionalGroups=c("Forage fish", "Macrozooplankton", "Herbivorous zooplankton", "Other"))
seabirdDiet$AE<-0.75 # Assimilation efficiency -> assume same for all seabird species & prey groups for now
seabirdDiet$energyDensity<-c(4.5, 4, 3.1, 2) # (kJ per g wet mass) #
# Energy density for forage fish based on this paper for Capelin: https://d1wqtxts1xzle7.cloudfront.net/46124708/Assimilation_efficiency_of_adult_Kittiwa20160601-4336-10nyk6f-libre.pdf?1464771742=&response-content-disposition=inline%3B+filename%3DAssimilation_efficiency_of_adult_Kittiwa.pdf&Expires=1776686451&Signature=FiV260FL8sA6XfT611m3V764hctijC01w2o6Cx9yHcrlMSq-hw3PRJQT-duIKO3~V9lh63q3HFgDZO2kv2tTPPH54yr9~xc7cSAmUGZTisvbtVaO5i5Dd3UUGIHqehy~dG4AyzRfAQYGZNFpN~s2ZgO--emdvI~GM88W3jbDydiJf-2EEbz3bcBav70q9l~4AdDjWEXb-2UfqHWx9-tVWGDvJTtJd4hx6hpYOtO1nDklqsIXrmcqDydleX1lijaaFGHDoUJENqx3U12SVPfzjL4mF5-9IwM9wDc7A4GMRiFe2UDN0vQHQX0rP4u97wOOYWq-IyvckPnI49hu9pVL5w__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA
# Energy density for herb. zoo. based on this paper for Calanus: https://www.sciencedirect.com/science/article/pii/S0079661116302105#s0010
# Energy density for macro. zoo. based on this paper for amphipods: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190345

# Determine proportion of each class in diet #
seabirdSpecies<-data.frame(species=rep(c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot"), each=4))
seabirdSpecies$functionalGroups<-rep(c("Forage fish", "Macrozooplankton", "Herbivorous zooplankton", "Other"), 6)
seabirdSpecies$propGroup<-c(0.94, 0.06, 0, 0, # BLK https://www.int-res.com/articles/meps2007/349/m349p269.pdf (I made an average of the Other category)
0.93, 0, 0, 0.07, # NoFu # file:///C:/Users/caitlin.frankish/Downloads/s002270050613.pdf
0.66, 0.02, 0, 0.41, # AtPu : https://onlinelibrary.wiley.com/doi/pdf/10.1111/ibi.12272
0.12, 0.87, 0.01, 0,# LiAu: https://link.springer.com/article/10.1007/s00300-013-1379-4/tables/2
1, 0, 0, 0,# CoGu # https://www.tandfonline.com/doi/epdf/10.1080/17451000802279636?needAccess=true
0.39, 0.55, 0, 0.06) #BrGu: https://cdnsciencepub.com/doi/epdf/10.1139/cjz-2021-0120 # But there are 6% other species so I will jus attribute this to the two classes  

# Attach energy density values #
seabirdSpeciesDiet<-seabirdSpecies %>%
dplyr::left_join(seabirdDiet, by=c("functionalGroups"))

foodDemandSeason<-list()

for (i in 1:length(seasonalHotspots)) {

# Subset to season i 
seasonalHotsposts_i<-seasonalHotspots[[i]]

# Extract how much energy from a raster #
seasonalMap<-energySeasonLox %>%
dplyr::filter(season==seasons[[i]]) %>%
ungroup() %>%
dplyr::select(x, y, energyTot) %>%
rename(z=energyTot) 

# Turn into raster
energyRast<-rasterFromXYZ(seasonalMap) # This shows kJ per km2

# mean raster value inside each polygon
polys_sp <- as(seasonalHotsposts_i, "Spatial")
vals_sum <- raster::extract(energyRast, polys_sp, fun=sum, na.rm=TRUE) # better to sum first?

# Calculate area of polygons (but project to equal-area first)
polys <- st_set_crs(seasonalHotsposts_i, projection_84)  # WGS84 (lat/lon)
polys_trans<-st_transform(polys, projection_NA)
polys$area_m2 <- as.numeric(st_area(polys_trans))
polys$area_km2 <- polys$area_m2 / 1e6

# Save results (centre coordinates, sum of energy, area km2)
energyPerHotspot<-data.frame(season=seasons[[i]], patches=polys$patches, energyTotkJ=vals_sum, area_km2=polys$area_km2)
energyPerHotspot$patches<-1:nrow(energyPerHotspot)

# Subset bird pie dataset
bird_pies_sub<-subset(bird_pies, season==seasons[[i]])
bird_pies_join<-bird_pies_sub %>%
dplyr::left_join(energyPerHotspot, by=c("season", "patches")) %>%
dplyr::filter(!is.na(X))

# Fish First
fishDemand<-subset(seabirdSpeciesDiet, functionalGroups=="Forage fish")
bird_pies_join_fish<-bird_pies_join %>%
dplyr::left_join(fishDemand, by=c("species")) %>%
dplyr::mutate(ForageFish_g=(energyTotkJ/(AE*energyDensity))*propGroup)

# MacroZoo First
MacroZooDemand<-subset(seabirdSpeciesDiet, functionalGroups=="Macrozooplankton")
bird_pies_join_macrozoo<-bird_pies_join %>%
dplyr::left_join(MacroZooDemand, by=c("species")) %>%
dplyr::mutate(MacroZoo_g=(energyTotkJ/(AE*energyDensity))*propGroup)

# HerbZoo First
HerbZooDemand<-subset(seabirdSpeciesDiet, functionalGroups=="Herbivorous zooplankton")
bird_pies_join_herbzoo<-bird_pies_join %>%
dplyr::left_join(HerbZooDemand, by=c("species")) %>%
dplyr::mutate(HerbZoo_g=(energyTotkJ/(AE*energyDensity))*propGroup)

# Join all results together and summarize... #
foodDemand<-bird_pies_join_fish %>%
ungroup() %>%
dplyr::mutate(MacroZoo_g=bird_pies_join_macrozoo$MacroZoo_g) %>%
dplyr::mutate(HerbZoo_g=bird_pies_join_herbzoo$HerbZoo_g) %>%
dplyr::group_by(season, patches) %>%
dplyr::summarise(x=mean(X), y=mean(Y), area_km2=mean(area_km2), totFish_g=sum(ForageFish_g), Fish_kg_km2=(totFish_g/area_km2)/1000, totMacroZoo_g=sum(MacroZoo_g), MacroZoo_kg_km2=(totMacroZoo_g/area_km2)/1000,
totHerbZoo_g=sum(HerbZoo_g), HerbZoo_kg_km2=(totHerbZoo_g/area_km2)/1000) %>%
ungroup() %>%
dplyr::group_by(season) %>%
dplyr::mutate(patches=row_number())

# Join results #
foodDemandSeason<-rbind(foodDemandSeason, foodDemand)

}

preyperhotspot<-ggplot() +
geom_col(data=foodDemandSeason, aes(x=patches, y=Fish_kg_km2, fill="Forage fish")) +
geom_col(data=foodDemandSeason, aes(x=patches, y=MacroZoo_kg_km2, fill="MacroZoo")) +
geom_col(data=foodDemandSeason, aes(x=patches, y=HerbZoo_kg_km2, fill="HerbZoo")) +
ylab("Kg per km2") +
xlab("Patch #") +
theme_bw() +
facet_wrap(~factor(season, levels=c("Autumn", "Winter", "Spring")), nrow=1, scales="free_x") +
theme_bw() +
scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73")) +
labs(fill="Prey type")
 
 pdf("./results/figures/speciesMaps/energy/hotspots4.pdf", width=10)
 plot(preyperhotspot)
 dev.off()

 
 pdf("./results/figures/speciesMaps/energy/hotspots3.pdf", width=10)
 plot(energyComposition)
 dev.off()

# Save output files
print("Saving output files")

# Number 2
output_file1 <- args[1]
print("Saving output file 1")
write.csv(seasonsRes, file = output_file1, row.names = FALSE) # characteristics of hotspots and control hotspots

output_file2 <- args[2]
print("Saving output file 1")
write.csv(seasonsRes2, file = output_file2, row.names = FALSE) # characteristics of hotspots and control hotspots

output_file3 <- args[3]
print("Saving output file 1")
write.csv(seasonsRes3, file = output_file3, row.names = FALSE) # characteristics of hotspots and control hotspots	