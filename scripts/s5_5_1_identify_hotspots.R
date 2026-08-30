# This script will be used to identify where hotspots in energy 'demand' occur #

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
library(dbscan)
library(terra)
library(spatialEco)
library(scatterpie)
library(marmap)
library(smoothr)

#### Step 0: setting up basic conditions ####

# Set-up some projections
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # Flat projection
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"

# Set-up number of iterations...
overall.iterations<-5 # how many times this is calculated per individual
print(paste0("Determining activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
#input_file1 <- args[1] # This will read in a species-name

#print(input_file1) # input_file1<-"Littleauk"

# Open-up species-specific monthly energy expenditure #
energyFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables", full.names=TRUE)
energySpeciesAll<-energyFiles[grepl(paste0("allenergy_v1"), energyFiles)]
energySpeciesDf<-fread(energySpeciesAll[1])

#### Step 1a: Look into SST sensitivity ####
energySpeciesDf$season<-ifelse(energySpeciesDf$month %in% c(4, 5, 6, 7, 8, 9), "Breeding", "Non-breeding")

### First we will just look at temporal differences ###

# Calculate percent increase per cell
energySpeciesDf$percent_increase<-(energySpeciesDf$energyPopkJ_uppersst-energySpeciesDf$energyPopkJ_mean)/energySpeciesDf$energyPopkJ_mean
energySpeciesDf$percent_increase_sst<-(energySpeciesDf$sst_ci_upper-energySpeciesDf$sst)

# Calculate percent decrease per cell
energySpeciesDf$percent_decrease<-(energySpeciesDf$energyPopkJ_mean - energySpeciesDf$energyPopkJ_lowersst)/energySpeciesDf$energyPopkJ_mean
energySpeciesDf$percent_decrease_sst<-(energySpeciesDf$sst - energySpeciesDf$sst_ci_lower)

# Calculate mean increase and decrease per month accross cells
monthlyEnergyCell<-energySpeciesDf %>%
dplyr::group_by(season, month) %>%
dplyr::summarise(meanIncrease=mean(percent_increase, na.rm=TRUE), sdIncrease=sd(percent_increase, na.rm=TRUE),
meanDecrease=mean(percent_decrease, na.rm=TRUE), sdDecrease=sd(percent_decrease, na.rm=TRUE),
meanIncrease2=mean(percent_increase_sst, na.rm=TRUE), sdIncrease2=sd(percent_increase_sst, na.rm=TRUE),
meanDecrease2=mean(percent_decrease_sst, na.rm=TRUE), sdDecrease2=sd(percent_decrease_sst, na.rm=TRUE)
)

# Make a plot to show this sensitivity:

sst_temporal<-monthlyEnergyCell %>%
dplyr::filter(season=="Non-breeding") %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3))) %>%
  ggplot(aes(y = factor(month))) +
  geom_col(aes(x = meanIncrease * 100),
           fill = "firebrick") +
  geom_col(aes(x = -meanDecrease * 100),
           fill = "steelblue") +
  geom_errorbarh(aes(
    xmin = (meanIncrease - sdIncrease) * 100,
    xmax = (meanIncrease + sdIncrease) * 100
  ),
  height = 0.2) +
  geom_errorbarh(aes(
    xmin = -(meanDecrease + sdDecrease) * 100,
    xmax = -(meanDecrease - sdDecrease) * 100
  ),
  height = 0.2) +
  labs(x = "Percent change (%)",
       y = "Month") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Month") +
  xlab("Percent change in energy expenditure +/- SD") +
  theme_bw()
  
  sst_temporal2<-monthlyEnergyCell %>%
dplyr::filter(season=="Non-breeding") %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3))) %>%
  ggplot(aes(y = factor(month))) +
  geom_col(aes(x = meanIncrease2 ),
           fill = "firebrick") +
  geom_col(aes(x = -meanDecrease2 ),
           fill = "steelblue") +
  geom_errorbarh(aes(
    xmin = (meanIncrease2 - sdIncrease2) ,
    xmax = (meanIncrease2 + sdIncrease2) 
  ),
  height = 0.2) +
  geom_errorbarh(aes(
    xmin = -(meanDecrease2 + sdDecrease2) ,
    xmax = -(meanDecrease2 - sdDecrease2) 
  ),
  height = 0.2) +
  labs(x = "Absolute change",
       y = "Month") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Month") +
  xlab("Absolute change in sst +/- SD") +
  theme_bw()
  
 # Try and calculate percent change energy with degree celcius
 
 energySpeciesDf <- energySpeciesDf %>%
  mutate(
    energy_percent_increase =
      100 * (energyPopkJ_uppersst - energyPopkJ_mean) /
      energyPopkJ_mean,

    energy_percent_decrease =
      100 * (energyPopkJ_mean - energyPopkJ_lowersst) /
      energyPopkJ_mean,

    energy_percent_increase_per_C =
      energy_percent_increase / percent_increase_sst,

    energy_percent_decrease_per_C =
      energy_percent_decrease / percent_decrease_sst
  )
  
 # Calculate mean increase and decrease per month accross cells
monthlyEnergyCell2<-energySpeciesDf %>%
#dplyr::filter(percent_increase_sst>0.1) %>%
dplyr::group_by(season, month) %>%
dplyr::summarise(meanIncrease=mean(energy_percent_increase_per_C , na.rm=TRUE), sdIncrease=sd(energy_percent_increase_per_C , na.rm=TRUE),
meanDecrease=mean(energy_percent_decrease_per_C , na.rm=TRUE), sdDecrease=sd(energy_percent_decrease_per_C , na.rm=TRUE)
)

sst_temporal3<-monthlyEnergyCell2 %>%
dplyr::filter(season=="Non-breeding") %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3))) %>%
  ggplot(aes(y = factor(month))) +
  geom_col(aes(x = meanIncrease ),
           fill = "firebrick") +
  geom_col(aes(x = -meanDecrease ),
           fill = "steelblue") +
  geom_errorbarh(aes(
    xmin = (meanIncrease - sdIncrease) ,
    xmax = (meanIncrease + sdIncrease) 
  ),
  height = 0.2) +
  geom_errorbarh(aes(
    xmin = -(meanDecrease + sdDecrease) ,
    xmax = -(meanDecrease - sdDecrease) 
  ),
  height = 0.2) +
  labs(x = "Percent change (%)",
       y = "Month") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Month") +
  xlab("Percent change in energy expenditure per deg C +/- SD") +
  theme_bw()
  
# Now we map spatial variation in these changes #

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# make a plot for viewing

nonbreeding<-subset(energySpeciesDf, season=="Non-breeding")

sst_spatial<-ggplot() +
 geom_tile(data=filter(nonbreeding, month==12), aes(x=x, y=y, fill=energy_percent_increase_per_C)) +
 scale_fill_gradientn(
  colours = c("#08306b", "#2171b5", "#ffff66", "#fdae61", "#d7191c")
) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(nonbreeding$x), max(nonbreeding$x)), ylim=c(min(nonbreeding$y) + 20, max(nonbreeding$y))) +
 theme_minimal() +
 labs(fill="") +
 facet_wrap(~month) +
 theme(legend.position="bottom") +
 xlab("") +
 ylab("") +
 labs(fill = "Percent increase in energy exp per degree C")
 
 pdf("./results/figures/speciesMaps/energy/Figure1.pdf")
 plot(Figure1)
 dev.off()

#### Step 1: Divide the data into general seasons and plot ####

scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

energySpeciesDf$season<-ifelse(energySpeciesDf$month %in% c(4, 5, 6, 7, 8, 9), "Breeding", "Non-breeding")

seasonalEnergy<-energySpeciesDf %>%
dplyr::group_by(season, x, y) %>%
dplyr::summarise(totEnergy=sum(energyPopkJ_mean, na.rm=TRUE)) %>%
ungroup()

breeding<-seasonalEnergy %>%
dplyr::filter(season=="Breeding") %>%
dplyr::mutate(Period="Breeding") %>%
#dplyr::filter(totEnergy>1) %>%
ungroup() %>%
dplyr::mutate(Energy_demand=scale01(totEnergy)) %>%
dplyr::select(Period, x, y, Energy_demand, totEnergy)

Nonbreeding<-seasonalEnergy %>%
dplyr::filter(season=="Non-breeding") %>%
dplyr::mutate(Period="Non-breeding") %>%
#dplyr::filter(totEnergy>1) %>%
ungroup() %>%
dplyr::mutate(Energy_demand=scale01(totEnergy)) %>%
dplyr::select(Period, x, y, Energy_demand, totEnergy)

yearlyEnergy<-energySpeciesDf %>%
dplyr::group_by(x, y) %>%
dplyr::summarise(totEnergy=sum(energyPopkJ_mean, na.rm=TRUE)) %>%
dplyr::mutate(Period="Yearly") %>%
ungroup() %>%
dplyr::mutate(Energy_demand=scale01(totEnergy)) %>%
dplyr::select(Period, x, y, Energy_demand, totEnergy)

energyDemandPlots<-rbind(breeding, Nonbreeding)

# Find out area size of every cell
energySpatialSub<-subset(energySpeciesDf, month==9) # subset to a random month
energyRast<-energySpatialSub %>%
dplyr::select(x, y, energyPopkJ_mean) %>%
rename(z=energyPopkJ_mean)
energyRast<-rasterFromXYZ(energyRast) # Turn into a raster
area_raster <- area(energyRast) # calculate area of every cell
area_raster_df<-as.data.frame(area_raster, xy=TRUE) # Turn into a data frame
colnames(area_raster_df)<-c("x", "y", "areaKm2") # Change col names

# Attach back to dataset
energyDemandPlots2<-energyDemandPlots %>%
dplyr::left_join(area_raster_df, by=c("x", "y")) %>%
dplyr::mutate(EnergyPerKm2=totEnergy/areaKm2)

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# make a plot for viewing
Figure1<-ggplot() +
 geom_tile(data=energyDemandPlots, aes(x=x, y=y, fill=Energy_demand)) +
 scale_fill_gradientn(
  colours = c("#08306b", "#2171b5", "#ffff66", "#fdae61", "#d7191c")
) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(energyDemandPlots$x), max(energyDemandPlots$x)), ylim=c(min(energyDemandPlots$y) + 20, max(energyDemandPlots$y))) +
 theme_minimal() +
 labs(fill="") +
 facet_wrap(~factor(Period, levels=c("Yearly", "Breeding", "Non-breeding")), nrow=3) +
 theme(legend.position="bottom") +
 xlab("") +
 ylab("") +
 labs(fill = "Scaled energy demand (0-1)", tag="B)")
 
 pdf("./results/figures/speciesMaps/energy/Figure1.pdf")
 plot(Figure1)
 dev.off()
 
#### Step 2: Make plots of top energy per cell etc. ####

Nonbreeding_patterns <- energySpeciesDf %>%
  dplyr::filter(season == "Non-breeding") %>%
  dplyr::mutate(Period = "Non-breeding") %>%
  dplyr::group_by(Period, x, y) %>%
  dplyr::slice_max(
    order_by = energyPopkJ_mean,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Period, x, y, month, energyPopkJ_mean, speciesNo, NoBirds_mean)
  
plot_max_val_per_cell<-ggplot() +
 geom_tile(data=Nonbreeding_patterns, aes(x=x, y=y, fill=energyPopkJ_mean/10e06)) +
 scale_fill_gradientn(
  colours = c("#08306b", "#2171b5", "#ffff66", "#fdae61", "#d7191c")
) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(energyDemandPlots$x), max(energyDemandPlots$x)), ylim=c(min(energyDemandPlots$y) + 20, max(energyDemandPlots$y))) +
 theme_minimal() +
 labs(fill="") +
 #facet_wrap(~factor(Period, levels=c("Yearly", "Breeding", "Non-breeding")), nrow=3) +
 theme(legend.position="bottom") +
 xlab("") +
 ylab("") +
 labs(fill = "Max energy demand per cell (kJ x 10e06)")
 
 cb_palette <- c(
  "#0072B2",  # blue
  "#E69F00",  # orange
  "#009E73",  # bluish green
  "#CC79A7",  # reddish purple
  "#D55E00",  # vermillion
  "#56B4E9"   # sky blue
)

Nonbreeding_patterns <- Nonbreeding_patterns %>%
  dplyr::mutate(
    month_max = factor(
      month,
      levels = c(10, 11, 12, 1, 2, 3),
      labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")
    )
  )
  
# Determine what is top 70% of max data values or something
aboveZeroValues<-subset(Nonbreeding_patterns, energyPopkJ_mean>0)
topValues<-quantile(aboveZeroValues$energyPopkJ_mean, 0.75)

plot_max_val_per_cell_month<-ggplot() +
  geom_tile(
    data = filter(Nonbreeding_patterns, energyPopkJ_mean > topValues),
    aes(x = x, y = y, fill = month_max)
  ) +
  scale_fill_manual(
    values = cb_palette,
    drop = FALSE
  ) +
  geom_sf(data = world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data = coast) +
  coord_sf(
    crs = projection_84,
    xlim = c(min(energyDemandPlots$x), max(energyDemandPlots$x)),
    ylim = c(min(energyDemandPlots$y) + 20, max(energyDemandPlots$y))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "",
    y = "",
    fill = "Month of max value"
  )
  
# Calculate minimum number of cells which sum to 30% of total energy demands #

# Fraction of total energy to capture
target_fraction <- 0.30

top30_cells <- Nonbreeding_patterns %>%
  filter(energyPopkJ_mean > 0) %>%
  arrange(desc(energyPopkJ_mean)) %>%
  mutate(
    cumulative_energy = cumsum(energyPopkJ_mean),
    cumulative_fraction = cumulative_energy / sum(energyPopkJ_mean)
  ) %>%
  filter(cumulative_fraction <= target_fraction)

# Add the first cell that exceeds 30%
if (max(top30_cells$cumulative_fraction) < target_fraction) {
  top30_cells <- Nonbreeding_patterns %>%
    filter(energyPopkJ_mean > 0) %>%
    arrange(desc(energyPopkJ_mean)) %>%
    mutate(
      cumulative_energy = cumsum(energyPopkJ_mean),
      cumulative_fraction = cumulative_energy / sum(energyPopkJ_mean)
    ) %>%
    slice(1:(nrow(top30_cells) + 1))
}

nrow(top30_cells)
sum(top30_cells$energyPopkJ_mean) / sum(Nonbreeding_patterns$energyPopkJ_mean)

topThirty<-ggplot() +
  geom_tile(
    data = top30_cells,
    aes(x = x, y = y, fill = energyPopkJ_mean / 1e7)
  ) +
  scale_fill_gradientn(
    colours = c("#08306b", "#2171b5", "#ffff66", "#fdae61", "#d7191c")
  ) +
  geom_sf(data = world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data = coast) +
  coord_sf(
    crs = projection_84,
    xlim = c(min(energyDemandPlots$x), max(energyDemandPlots$x)),
    ylim = c(min(energyDemandPlots$y) + 20, max(energyDemandPlots$y))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(fill = "Energy demand (kJ × 10⁷)") +
  ylab("") +
  xlab("")
  
# Turn these into polygons #

# Estimate grid-cell dimensions from coordinate spacing
cell_width  <- min(diff(sort(unique(Nonbreeding_patterns$x))))
cell_height <- min(diff(sort(unique(Nonbreeding_patterns$y))))

# Convert each hotspot cell centre into a rectangular polygon
hotspot_cells_sf <- top30_cells %>%
  rowwise() %>%
  mutate(
    geometry = list(
      st_polygon(list(matrix(
        c(
          x - cell_width / 2, y - cell_height / 2,
          x + cell_width / 2, y - cell_height / 2,
          x + cell_width / 2, y + cell_height / 2,
          x - cell_width / 2, y + cell_height / 2,
          x - cell_width / 2, y - cell_height / 2
        ),
        ncol = 2,
        byrow = TRUE
      )))
    )
  ) %>%
  ungroup() %>%
  st_as_sf(crs = projection_84)
  
# project
hotspot_cells_projected <- st_transform(
  hotspot_cells_sf,
  crs = projection_NA
)

# Merge touching cells
hotspot_union <- hotspot_cells_projected %>%
  st_union() %>%
  st_make_valid()

# Split disconnected regions into individual polygons
hotspot_polygons <- hotspot_union %>%
  st_cast("POLYGON") %>%
  st_sf() %>%
  mutate(
    hotspot_id = row_number(),
    area_km2 = as.numeric(st_area(geometry)) / 1e6
  )
  
minimum_area_km2 <- 5000  # adjust this

hotspot_polygons_filtered <- hotspot_polygons %>%
  filter(area_km2 >= minimum_area_km2)
  
smooth_distance <- 250000  # metres

hotspot_polygons_smoothed <- hotspot_polygons_filtered %>%
  st_buffer(smooth_distance) %>%
  st_union() %>%
  st_buffer(-smooth_distance) %>%
  st_make_valid() %>%
  st_cast("POLYGON") %>%
  st_sf() %>%
  mutate(
    hotspot_id = row_number(),
    area_km2 = as.numeric(st_area(geometry)) / 1e6
  ) %>%
  st_transform(projection_84)
 
# Map this around my original map #
topThirty <- ggplot() +
  geom_tile(
    data = top30_cells,
    aes(x = x, y = y, fill = energyPopkJ_mean / 1e7)
  ) +
  geom_sf(
    data = hotspot_polygons_smoothed,
    fill = NA,
    aes(colour = factor(hotspot_id)),
    linewidth = 0.8
  ) +
  scale_fill_gradientn(
    colours = c("#08306b", "#2171b5", "#ffff66", "#fdae61", "#d7191c")
  ) +
  geom_sf(
    data = world,
    color = "#E5E5E5",
    fill = "#E5E5E5"
  ) +
  geom_sf(data = coast) +
  coord_sf(
    crs = projection_84,
    xlim = c(min(energyDemandPlots$x), max(energyDemandPlots$x)),
    ylim = c(min(energyDemandPlots$y) + 20, max(energyDemandPlots$y))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(fill = "Energy demand (kJ × 10⁷)") +
  xlab("") +
  ylab("")
  
# Make the cutter in longitude/latitude.
# Adjust these coordinates until the rectangle crosses the unwanted neck.
greenland_cut_ll <- st_as_sfc(
  st_bbox(
    c(
      xmin = -60,
      xmax = -40,
      ymin = 55.6,
      ymax = 57.8
    ),
    crs = st_crs(projection_84)
  )
)

# Transform it to the projected CRS used for your hotspot polygons
greenland_cut <- st_transform(
  greenland_cut_ll,
  st_crs(hotspot_polygons_smoothed)
)

# Remove the cutter area from the hotspot geometry
hotspot_polygons_separated <- hotspot_polygons_smoothed %>%
  st_difference(greenland_cut) %>%
  st_make_valid() %>%
  st_cast("POLYGON") %>%
  filter(!st_is_empty(geometry)) %>%
  mutate(
    hotspot_id = row_number(),
    area_km2 = as.numeric(st_area(geometry)) / 1e6
  )
  

hotspot_polygons_separated <- hotspot_polygons_separated %>%
  filter(area_km2 >= minimum_area_km2)
  
ggplot() +

  # Energy cells
  geom_tile(
    data = top30_cells,
    aes(
      x = x,
      y = y,
      fill = energyPopkJ_mean / 1e7
    ),
    alpha = 0.8
  ) +

  scale_fill_gradientn(
    colours = c(
      "#08306b",
      "#2171b5",
      "#ffff66",
      "#fdae61",
      "#d7191c"
    )
  ) +

  # White halo
  geom_sf(
    data = hotspot_polygons_separated,
    fill = NA,
    colour = "white",
    linewidth = 2.5
  ) +

  # Coloured hotspot outline
  geom_sf(
    data = hotspot_polygons_separated,
    aes(colour = factor(hotspot_id)),
    fill = NA,
    linewidth = 1.2
  ) +
  
  # Basemap first
  geom_sf(
    data = world,
    colour = "#E5E5E5",
    fill = "#E5E5E5"
  ) +

  geom_sf(
    data = coast,
    fill = NA
  ) +

  coord_sf(
    crs = projection_84,
    xlim = c(-110, 60),
    ylim = c(38, 85),
    expand = FALSE
  ) +

  theme_minimal() +
  theme(
    legend.position = "bottom"
  ) +

  labs(
    fill = "Energy demand (kJ × 10⁷)",
    colour = "Hotspot #",
    x = NULL,
    y = NULL
  )
 
#### Step 2: Figure out number of pops, birds and energy per species ####

# Number of pops per species #
allPopFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/tmp5", full.names=TRUE) # Determine where all population files are kept
allPopFilesv1<-allPopFiles[grep("v1", allPopFiles)] # Subset to temporal files for now

# Make a data frame out of This
speciesInfo<-data.frame(popfiles=allPopFilesv1, species=sub(".*/tmp5/([^_]+)_.*", "\\1", allPopFilesv1))
speciesInfoMeta<-speciesInfo %>%
dplyr::group_by(species) %>%
dplyr::summarise(popNo=n_distinct(popfiles))

# Open up species files one by one & attach information on total number of birds & tot energy demand during the winter months # (or all year-round?)
allSpeciesFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/", full.names=TRUE)
allSpeciesFiles<-allSpeciesFiles[grepl("monthly_map_v1.csv", allSpeciesFiles)]
allSpeciesFilesDf<-data.frame(files=allSpeciesFiles, species=sub(".*tables//([^_]+)_.*", "\\1", allSpeciesFiles))

popInfo<-list()

for (i in 1:length(allSpeciesFiles)) {

print(paste0("Step 2 : Summarizing info from species ", i)) 

# Open up file i 
mapSub<-fread(allSpeciesFiles[i])

# Extract pop number & tot energy demand
energyPopInfo<-mapSub %>%
dplyr::group_by(species) %>%
dplyr::summarise(totBirds=sum(NoBirds_mean, na.rm=TRUE), totBirds_lower=sum(NoBirds_ci_lower, na.rm=TRUE), totBirds_upper=sum(NoBirds_ci_upper, na.rm=TRUE), 
totEnergy=sum(energyPopkJ_mean), totEnergy_lower=sum(energyPopkJ_lowerci, na.rm=TRUE), totEnergy_upper=sum(energyPopkJ_upperci, na.rm=TRUE)) %>%
dplyr::ungroup() 

energyPopInfo$species<-allSpeciesFilesDf$species[i]

# Make a summary data frame to join with the main one #
dfSpecies<-data.frame(species=allSpeciesFilesDf$species[i], totBirds=max(energyPopInfo$totBirds), totEnergykJ=sum(energyPopInfo$totEnergy))

# Save results
popInfo<-rbind(popInfo, energyPopInfo)

}

# Join information back again #
speciesInfoMeta2<-speciesInfoMeta %>%
dplyr::left_join(popInfo, by=c("species"))

# Re-arrange this dataset for plotting #
popNo<-speciesInfoMeta2 %>%
dplyr::select(species, popNo) %>%
dplyr::rename(value=popNo) %>%
dplyr::mutate(metric="Tot pops")

birdNo<-speciesInfoMeta2 %>%
dplyr::select(species, totBirds, totBirds_lower, totBirds_upper) %>%
dplyr::rename(value=totBirds, value_lower=totBirds_lower, value_upper=totBirds_upper) %>%
dplyr::mutate(metric="Tot birds (x 10e06)") %>%
dplyr::mutate(value=value/10^6, value_lower=value_lower/10^6, value_upper=value_upper/10^6)

energyNo<-speciesInfoMeta2 %>%
dplyr::select(species, totEnergy, totEnergy_lower, totEnergy_upper) %>%
dplyr::rename(value=totEnergy, value_lower=totEnergy_lower, value_upper=totEnergy_upper) %>%
dplyr::mutate(value=value/10^9, value_lower=value_lower/10^9, value_upper=value_upper/10^9) %>%
dplyr::mutate(metric="Tot energy demand (kJ x 10e09)")

metrics<-rbind(birdNo, energyNo)

# Determine factor level based on highest to lowest number of birds #
birdNoArranged<-birdNo %>%
ungroup() %>%
arrange(desc(value))

levelsBirds<-unique(birdNoArranged$species)

# plot metrics #
Figure1b<-metrics %>%
dplyr::mutate(species=factor(species, levels=c(
      levelsBirds
    ))) %>%
dplyr::mutate(metric=factor(metric, levels=c("Tot pops", "Tot birds (x 10e06)", "Tot energy demand (kJ x 10e09)"))) %>%
dplyr::filter(!metric %in% c("Tot pops")) %>%
ungroup() %>%
ggplot(aes(x=species, y=value)) +
scale_fill_manual(
  values = c(
    "#008856",
	"#C3A600",
    "#E25822",
    "#0072b2",
    "#875692",
	"#BE0032"
  )
) +
geom_col(aes(fill=species)) + 
geom_errorbar(
  aes(ymin = value - value_lower,
      ymax = value + value_upper),
  width = 0.2
)+
facet_wrap(~metric, scales="free_y", nrow=3) +
theme_bw() +
#coord_flip() +
theme(legend.position="none", axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
xlab("") +
ylab("") +
labs(tag="A)")


pdf("./results/figures/speciesMaps/energy/Figure1b.pdf", width=3, height=5)
plot(Figure1b)
dev.off()

 
#### Step 3: Extract monthly polygons for NB season only ####

energySpeciesDf$season<-ifelse(energySpeciesDf$month %in% c(4, 5, 6, 7, 8, 9), "Breeding", "Non-breeding")

# Plot top 25% of energy in both cases #

nonBreedingEnergy<-subset(energySpeciesDf, season=="Non-breeding")
nonBreedingEnergy$month<-factor(nonBreedingEnergy$month, levels=c(10, 11, 12, 1, 2, 3))
monthsUnique<-levels(nonBreedingEnergy$month)

topenergy<-list()

for (i in 1:length(monthsUnique)) {

# Subset to season i
seasonSub<-monthsUnique[i]

print(seasonSub)

# Calculate top % of data
top25_energy<-topEnergy(subset(nonBreedingEnergy, month==seasonSub), 0.10, 'energyPopkJ_mean')

# Add seasonal info
top25_energy$month<-seasonSub

# Save results
topenergy<-rbind(topenergy, top25_energy)

}

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# make a plot for viewing
hotspots<-ggplot() +
 geom_tile(data=topenergy, aes(x=x, y=y, fill="10%_energy")) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(topenergy$x), max(topenergy$x) + 25), ylim=c(min(topenergy$y), max(topenergy$y) + 10)) +
 theme_minimal() +
 ggtitle("Top 10% of data") +
 labs(fill="") +
 xlab("") +
 ylab("") + 
 facet_wrap(~month)
 
### Step 3: Convert into polygons ####

print("Step 3: Converting into polygons...")

# Need to turn these files into polygons somehow (e.g. polygon 1, polygon 2 etc.) #

# First we turn into a raster #

seasonalHotspots<-list() # for energy demand

monthsUnique<-levels(nonBreedingEnergy$month)

for (i in 1:length(monthsUnique)) {

print(monthsUnique[i])

# Subset to season i
seasonalMap<-nonBreedingEnergy %>%
dplyr::filter(month==monthsUnique[[i]]) %>%
ungroup() %>%
dplyr::select(index, x, y)

# Trasnform dataset for making into raster
topenergy_season<-topenergy %>%
dplyr::filter(month==monthsUnique[[i]]) %>%
dplyr::mutate(topenergy=1) %>%
dplyr::select(index, topenergy, season) %>%
dplyr::full_join(seasonalMap, by=c("index")) %>%
replace_na(list(topenergy=0)) %>%
arrange(index) %>%
ungroup() %>%
dplyr::select(x, y, topenergy) %>%
dplyr::rename(z=topenergy)

# Turn into raster
energyRast<-rasterFromXYZ(topenergy_season) 

# Convert to terra object
r_terra <- rast(energyRast)

# Use a distance threshold to seperate 'hotspots' 
d <- 1 

# Convert hotspot cells into 1, others into NA
r1 <- r_terra == 1
r1 <- mask(r1, r1, maskvalues = 0)   

# distance from every cell to the nearest '1' cell
dist_to_1 <- distance(r1)

# expand the blobs by distance d (a morphological "buffer" in raster space)
expanded <- dist_to_1 <= d
expanded <- mask(expanded, expanded, maskvalues = 0)

# label connected components in expanded space
grp <- patches(expanded, directions = 8)

# polygonize expanded groups, not just original 1-cells
pol <- as.polygons(grp, dissolve = TRUE)
pol <- pol[!is.na(values(pol)[,1]), ]

# Convert to sf
pol_sf <- st_as_sf(pol)

# Remove small polygons and renumber remaining polygons
min_area <- 0.1  # choose your threshold in map units^2

pol_sf <- pol_sf %>%
  dplyr::mutate(area = sf::st_area(.)) %>%
  dplyr::filter(area >= min_area) %>%
  dplyr::mutate(patches = dplyr::row_number()) %>%
  dplyr::mutate(month=monthsUnique[i])

# Save results
seasonalHotspots[[i]]<-pol_sf

}

# Bind these all together
seasonalHotspotsAll <- dplyr::bind_rows(seasonalHotspots) %>%
  dplyr::mutate(
    patches = factor(patches),
    month = factor(month, levels = monthsUnique)
  )
  
# Set CRS
seasonalHotspotsAll <- st_set_crs(seasonalHotspotsAll, 4326)

# Smooth polygons
seasonalHotspotsAll_smooth <- smoothr::smooth(
  seasonalHotspotsAll,
  method = "ksmooth",
  smoothness = 2
)

# Plot results #
cb_cols <- c(
  "#000000", # black
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7"  # reddish purple
)

Figure2a<-ggplot() +
  geom_sf(
    data = seasonalHotspotsAll_smooth,
    aes(color = month),
    fill = NA,
    alpha = 0.7, linewidth = 0.6
  ) +
  geom_sf(data = world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data = coast) +
  scale_color_manual(values = cb_cols) +
  coord_sf(
    crs = projection_84,
    xlim = c(min(topenergy$x), max(topenergy$x) + 25),
    ylim = c(min(topenergy$y), max(topenergy$y) + 10)
  ) +
  theme_minimal() +
  ggtitle("Top 10% of data") +
  labs(color= "Month") +
  xlab("") +
  ylab("")
  
pdf("./results/figures/speciesMaps/energy/Figure2a.pdf")
plot(Figure2a)
dev.off()

### Step 4: Extract proportion of energy demand by species for different hotspots ####

print("Step 4: Extracting characteristics of polygons...")

# First here I am grouping them together based on what I can see visually #

hotspotMetadata<-data.frame(month=c(10, 10, 10, 10, 10, 10, 
11, 11, 11, 11,
12, 12, 12, 12, 12,
1, 1, 1, 1, 1, 1,
2, 2, 2,
3, 3, 3
), 
patches=c(1, 2, 3, 4, 5, 6,
1, 2, 3, 4,
1, 2, 3, 4, 5,
1, 2, 3, 4, 5, 6,
1, 2, 3,
1, 2, 3
),
hotspotNo=c(1, 2, 3, 4, 5, 5,
1, 2, 2, 3,
1, 2, 3, 4, 4,
1, 2, 3, 4, 4, 4,
1, 2, 2,
1, 2, 3
),
name=c("Barents_Sea", "NW_Iceland", "SW_Greenland", "UK", "Grand_banks", "Grand_banks",
"NW_Iceland", "SW_Greenland", "SW_Greenland", "Grand_banks",
"Barents_Sea", "NW_Iceland", "SW_Greenland", "Grand_banks", "Grand_banks",
"Barents_Sea", "NW_Iceland", "SW_Greenland", "Grand_banks", "Grand_banks", "Grand_banks",
"NW_Iceland", "Grand_banks", "Grand_banks",
"Barents_Sea", "NW_Iceland", "UK"))

hotspotMetadata$month<-factor(hotspotMetadata$month)
hotspotMetadata$patches<-factor(hotspotMetadata$patches)

# Loop through each population & then month
hotspotComposition<-list() # for energy demand

# Determine unique months
monthsUnique<-levels(nonBreedingEnergy$month)

# Determine where pop energy demand rasters are located #
energyMaps<-list.files("tmp5/", full.names=TRUE)
energyMaps_v1<-energyMaps[grep("v1", energyMaps)] 

for (i in 1:length(energyMaps_v1)) {

print(paste0("Extracting from map ", i, "/", length(energyMaps_v1)))

# Open population i
energyMap_i<-fread(energyMaps_v1[i])

# list to save values in
birdComp_season<-list()

# Now we loop through the months #

for (j in 1:length(monthsUnique)) {

#for (j in 1:10) {

print(paste0("Extracting from month ", j, "/", length(monthsUnique)))

# Determine month j
seasonSub<-monthsUnique[j]

# Subset to correct months #
energyMap_i_season<-subset(energyMap_i, month %in% unique(seasonSub))

# Find out area size of every cell
energySpatialSub<-subset(energyMap_i, month==9) # subset to a random month
energyRast<-energySpatialSub %>%
dplyr::select(x, y, energyPopkJ_mean) %>%
rename(z=energyPopkJ_mean)
energyRast<-rasterFromXYZ(energyRast) # Turn into a raster
area_raster <- area(energyRast) # calculate area of every cell
area_raster_df<-as.data.frame(area_raster, xy=TRUE) # Turn into a data frame
colnames(area_raster_df)<-c("x", "y", "areaKm2") # Change col names

# Join this to main dataset? 
energyMap_i_season<-energyMap_i_season %>%
dplyr::left_join(area_raster_df, by=c("x", "y"))

# Make an energy raster #
energy<-energyMap_i_season %>%
ungroup() %>%
dplyr::rename(z=energyPopkJ_mean) %>%
dplyr::select(x, y, z)

energyRast<-rasterFromXYZ(energy)

# Extract values from all of these #
#polySub<-seasonalHotspotsAll_smooth %>%
#dplyr::filter(month==seasonSub) %>%
#dplyr::left_join(hotspotMetadata, by=c("month", "patches"))

polySub<-hotspot_polygons_separated 

# Calculate area of polygons to join later
polySub$area_m2 <- st_area(polySub)

valsEnergy <- raster::extract(energyRast, polySub)
valsEnergy_sumbyid<-sapply(valsEnergy, sum, na.rm = TRUE)

# Save these results

polygonComp<-data.frame(patches=1:nrow(polySub), totEnergy=valsEnergy_sumbyid)
polygonComp$species<-energyMap_i_season$species[1]
polygonComp$colony<-energyMap_i_season$colony[1]
polygonComp$month<-seasonSub
polygonComp<-polygonComp %>%
dplyr::mutate(patches=factor(patches), month=monthsUnique[j]) %>%
#dplyr::left_join(hotspotMetadata, by=c("month", "patches")) %>%
dplyr::mutate(area_m2=polySub$area_m2)

# Save results
birdComp_season<-rbind(birdComp_season, polygonComp)

}

# Save results
hotspotComposition<-rbind(hotspotComposition, birdComp_season)

}

energyMaps_v1 <- list.files(
  "tmp5/",
  pattern = "v1",
  full.names = TRUE
)

#--------------------------------------------------
# 1. Prepare polygons once
#--------------------------------------------------

poly_sf <- hotspot_polygons_separated %>%
  dplyr::mutate(
    patches = seq_len(dplyr::n())
  )

# Calculate polygon area only once.
# Replace projection_projected with an appropriate projected CRS.
polygon_area <- data.table(
  patches = poly_sf$patches,
  area_m2 = as.numeric(
    sf::st_area(
      sf::st_transform(poly_sf, projection_projected)
    )
  )
)

poly_vect <- terra::vect(poly_sf)

#--------------------------------------------------
# 2. Create a template raster once
#--------------------------------------------------

template_data <- fread(
  energyMaps_v1[1],
  select = c("x", "y", "month", "energyPopkJ_mean")
)

template_xyz <- template_data[
  month == 9,
  .(
    x,
    y,
    z = 0
  )
]

template_raster <- terra::rast(
  template_xyz,
  type = "xyz",
  crs = "EPSG:4326"
)

# Lookup table linking coordinates to raster cell numbers
cell_lookup <- as.data.table(
  terra::xyFromCell(
    template_raster,
    seq_len(terra::ncell(template_raster))
  )
)

cell_lookup[, cell := seq_len(.N)]

setkey(cell_lookup, x, y)

#--------------------------------------------------
# 3. Determine which raster cells belong to polygons
#    Do this only once
#--------------------------------------------------

cell_membership <- terra::extract(
  template_raster,
  poly_vect,
  cells = TRUE
)

cell_membership <- as.data.table(cell_membership)

# terra usually calls the polygon identifier "ID"
cell_membership <- unique(
  cell_membership[, .(
    patches = ID,
    cell
  )]
)

setkey(cell_membership, cell)

#--------------------------------------------------
# 4. Process each file
#--------------------------------------------------

results_list <- lapply(
  seq_along(energyMaps_v1),
  function(i) {

    message(
      "Extracting map ",
      i,
      "/",
      length(energyMaps_v1)
    )

    energy_dt <- fread(
      energyMaps_v1[i],
      select = c(
        "x",
        "y",
        "month",
        "energyPopkJ_mean",
        "species",
        "colony"
      )
    )

    # Keep only required months
    energy_dt <- energy_dt[
      month %in% monthsUnique
    ]

    # Add raster cell numbers using x/y coordinates
    energy_dt <- cell_lookup[
      energy_dt,
      on = .(x, y),
      nomatch = 0
    ]

    # Attach polygon membership to each raster cell
    polygon_values <- cell_membership[
      energy_dt,
      on = "cell",
      nomatch = 0,
      allow.cartesian = TRUE
    ]

    # Sum energy by polygon and month
    output <- polygon_values[
      ,
      .(
        totEnergy = sum(
          energyPopkJ_mean,
          na.rm = TRUE
        )
      ),
      by = .(
        patches,
        month,
        species,
        colony
      )
    ]

    # Add polygon area
    output <- polygon_area[
      output,
      on = "patches"
    ]

    output
  }
)

#--------------------------------------------------
# 5. Combine everything once
#--------------------------------------------------

hotspotComposition <- rbindlist(
  results_list,
  use.names = TRUE,
  fill = TRUE
)

hotspotComposition[
  ,
  `:=`(
    patches = factor(patches),
    month = factor(month)
  )
]

# Plot these results #

Figure2b<-hotspotComposition %>%
  ungroup() %>%
  group_by(patches, month) %>%
  mutate(
    energyAllSpecies = sum(totEnergy, na.rm = TRUE),
    areaAll = sum(area_m2, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(patches, month, species) %>%
  summarise(
    EnergyOverall = sum(totEnergy, na.rm = TRUE),
    propEnergy = EnergyOverall / max(energyAllSpecies, na.rm = TRUE),
    totEnergy = max(energyAllSpecies, na.rm = TRUE),
    totArea = max(areaAll, na.rm = TRUE),
	totArea_km2 =totArea / 1e6,
	EnergyPerkm2=EnergyOverall/totArea_km2,
    .groups = "drop"
  ) %>%
  mutate(
    month_plot = match(month, c(10, 11, 12, 1, 2, 3)),
    species = factor(
      species,
      levels = c(
        "Black-legged kittiwake", "Northern fulmar",
        "Atlantic puffin", "Little auk",
        "Common guillemot", "Brünnich's guillemot"
      )
    )
  ) %>%
  #dplyr::mutate(name=factor(name, levels=c("Barents_Sea", "NW_Iceland", "SW_Greenland", "UK", "Grand_banks"))) %>%
  ggplot(aes(x = month_plot, y = EnergyPerkm2, fill = species)) +
  geom_col() +
  facet_wrap(~patches, nrow = 4) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")
  ) +
  scale_fill_manual(values = c(
    "#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"
  )) +
  theme_bw(base_size = 14) +
  theme(legend.position="bottom") +
 ylab(expression("Total energy demand (kJ km"^-2*" month"^-1*")")) +
  xlab("") 
  
pdf("./results/figures/speciesMaps/energy/Figure2b.pdf", width=9, height=7)
plot(Figure2b)
dev.off()

# Save output files
print("Saving output files")

# Number 2
output_file1 <- args[1]
print("Saving output file 1")
write.csv(hotspotComposition, file = output_file1, row.names = FALSE) # characteristics of hotspots and control hotspots