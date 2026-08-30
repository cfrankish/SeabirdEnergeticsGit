library(data.table)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(gganimate)
library(raster)
library(gifski)

# Make an animated map for WSC #

scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

energySpeciesDf<-fread("tablex_allenergy_v1.csv")

energySpeciesDf$season<-ifelse(energySpeciesDf$month %in% c(4, 5, 6, 7, 8, 9), "Breeding", "Non-breeding")

energySpeciesDf$season<-ifelse(energySpeciesDf$month %in% c(4, 5, 6, 7, 8, 9), "Breeding", "Non-breeding")

seasonalEnergy<-energySpeciesDf %>%
  dplyr::group_by(season, month, x, y) %>%
  #dplyr::summarise(totEnergy=sum(energyPopkJ_mean, na.rm=TRUE)) %>%
  ungroup()

Nonbreeding<-seasonalEnergy %>%
  dplyr::filter(season=="Non-breeding") %>%
  dplyr::mutate(Period="Non-breeding") %>%
  #dplyr::filter(totEnergy>1) %>%
  ungroup() %>%
  dplyr::mutate(Energy_demand=scale01(energyPopkJ_mean)) %>%
  dplyr::select(month, Period, x, y, Energy_demand, energyPopkJ_mean)

# Now we map spatial variation in these changes #

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# First we need to project the maps #

# Set-up some projections
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # Flat projection
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"

monthsUnique<-c(10, 11, 12, 1, 2, 3)

monthlyrasters<-list()

for (i in 1:length(monthsUnique)) {

  print(monthsUnique[i])
  
nonbreeding_rast<-Nonbreeding %>%
  dplyr::filter(month==monthsUnique[i]) %>%
  dplyr::select(x, y, Energy_demand) %>%
  rename(z=Energy_demand)

NB_demand_rast<-rasterFromXYZ(nonbreeding_rast, crs=projection_84)

NB_proj <- raster::projectRaster(
  NB_demand_rast,
  crs = projection_NA,
  method = "bilinear"
)

NB_proj_df<-as.data.frame(NB_proj, xy=TRUE)

NB_proj_df$month<-monthsUnique[i]

monthlyrasters<-rbind(monthlyrasters, NB_proj_df)

}

library(gganimate)

month1<-monthlyrasters
month1$month<-factor(month1$month)

plot_conference<-ggplot() +
  geom_raster(
    data = month1,
    aes(x = x, y = y, fill = z)
  ) +
  scale_fill_gradientn(
    colours = c("#08306b", "#2171b5", "#ffff66", "#fdae61", "#d7191c"),
    limits = c(0, 1)
  ) +
  geom_sf(
    data = world,
    color = "#E5E5E5",
    fill = "#E5E5E5"
  ) +
  geom_sf(data = coast) +
  coord_sf(
    crs = projection_NA,
    xlim = c(-4300000, 1900000),
    ylim = c(-3000000, 2700000)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlab("") +
  ylab("") +
  labs(fill="Energy demand scaled 0-1")

animation = plot_conference +
  ggtitle("Month: {closest_state}") +
  transition_states(states = month)

gganimate::animate(animation, 
        renderer = gifski_renderer("map.gif"))
