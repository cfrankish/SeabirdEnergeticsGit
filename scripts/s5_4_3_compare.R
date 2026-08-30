## Here we map energetics using static activity (just changes by month but spatial SST ##

library(dplyr)
library(data.table)
library(ncdf4)
library(raster)
library(sp)
library(sf)
library(lubridate)
library(data.table)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggplot2)
library(gridExtra)
library(MASS)
library(terra)
#library(biscale)
#library(cowplot)
#library(classInt)
#library(spdep)

#### Step 0: setting up basic conditions ####

# Set-up number of iterations...
print(paste0("Mapping energy..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file1 <- args[1] # This will read in a species name
input_file2<-args[2]

# open these files
map_v1<-fread(input_file1) # map_v1<-fread("results/tables/tablex_allenergy_v1.csv")
map_v2<-fread(input_file2) # map_v2<-fread("results/tables/tablex_allenergy_v2.csv")
#sum_v1<-fread(input_file3) # sum_v1<-fread("results/tables/Littleauk_energyMap_monthly_sum_v1.csv")
#sum_v2<-fread(input_file4) # sum_v2<-fread("results/tables/Littleauk_energyMap_monthly_sum_v2.csv")

#### Step 1: make some comparison plots ####

#print("Prepping map info")
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Scaling functions
scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

monthLoop<-c(9, 10, 11, 12, 1, 2, 3, 4)

pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/energy/AllBirds_energy_compare.pdf"), width=10, height=7)

for (m in 1:length(monthLoop)) {

print(m)

# Subset maps to month i
monthSub_v1<-subset(map_v1, month==monthLoop[m])
monthSub_v2<-subset(map_v2, month==monthLoop[m])

# Subset sums to month i
#monthSub_v1_sum<-subset(sum_v1, month==monthLoop[m])
#monthSub_v2_sum<-subset(sum_v2, month==monthLoop[m])

# Make scaled datasets to compare

energySpecies1<-monthSub_v1 %>%
#dplyr::filter(energyPopkJ_mean>0) %>%
dplyr::ungroup() %>%
dplyr::select(x, y, month, energyPopkJ_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric="Energy_expenditure_v1") %>%
dplyr::select(-energyPopkJ_mean)

energySpecies2<-monthSub_v2 %>%
#dplyr::filter(energyPopkJ_mean>0) %>%
dplyr::ungroup() %>%
dplyr::select( x, y, month, energyPopkJ_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric="Energy_expenditure_v2") %>%
dplyr::select(-energyPopkJ_mean)

energySpecies1$metric_scaled2<-energySpecies2$metric_scaled
energySpecies1$diff<-energySpecies1$metric_scaled - energySpecies1$metric_scaled2

# Figuring out a common extent #
energyBirds<-subset(energySpecies1, metric_scaled>0)
min1<-min(energyBirds$metric_scaled)
min2<-min(energyBirds$metric_scaled2)
min<-ifelse(min1<min2, min1, min2)
max1<-max(energyBirds$metric_scaled)
max2<-max(energyBirds$metric_scaled2)
max<-ifelse(max1>max2, max1, max2)

plot1_energy1<-ggplot() +
geom_tile(data=subset(energySpecies1, metric_scaled>0), aes(x=x, y=y, fill=metric_scaled)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(min, max)) +
#facet_wrap(~metric) +
theme_bw() +
ggtitle(paste0("V1: Month ", monthLoop[m]))

plot1_energy2<-ggplot() +
geom_tile(data=subset(energySpecies1, metric_scaled>0), aes(x=x, y=y, fill=metric_scaled2)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(min, max)) +
#facet_wrap(~metric) +
theme_bw() +
ggtitle(paste0("V2: Month ", monthLoop[m]))

plot1_energy_diff<-ggplot() +
geom_tile(data=subset(energySpecies1, metric_scaled>0), aes(x=x, y=y, fill=diff)) +
scale_fill_gradientn(
  colours = c(
    "#364B9A", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#A50026"
  ),
  limits = c(min(energySpecies1$diff, na.rm=TRUE), max(energySpecies1$diff, na.rm=TRUE)),
  values = scales::rescale(c(min(energySpecies1$diff, na.rm=TRUE), 0, max(energySpecies1$diff, na.rm=TRUE)))
) +
theme_bw() +
ggtitle(paste0("v1-v2; Month ", monthLoop[m]))

# Add map of where hotspots occur # (energy, bird density)
topEnergy_v1<-topEnergy(monthSub_v1, 0.1, "energyPopkJ_mean") 
topEnergy_v2<-topEnergy(monthSub_v2, 0.1, "energyPopkJ_mean") 
topBirds<-topEnergy(monthSub_v1, 0.1, "NoBirds_mean") 

topData<-ggplot() +
geom_tile(data=subset(monthSub_v1, NoBirds_mean>0), aes(x=x, y=y), fill="grey", alpha=0.5)+
geom_tile(data=topBirds, aes(x=x, y=y, fill="BirdNo"),  alpha=0.6) +
geom_tile(data=topEnergy_v1, aes(x=x, y=y, fill="MEE_v1"),  alpha=0.3) +
#geom_tile(data=topEnergy_v2, aes(x=x, y=y, fill="MEE_v2"),  alpha=0.3) +
scale_fill_manual(values=c("#0072B2", "#E69F00", "#009E73")) +
theme_minimal() +
  ggtitle("top 5% data contour polygons")

# other option
topBird_contour_v1<-topEnergy_contour(monthSub_v1, 0.05, "NoBirds_mean")
topEnergy_contour_v1<-topEnergy_contour(monthSub_v1, 0.05, "energyPopkJ_mean")
topEnergy_contour_v2<-topEnergy_contour(monthSub_v2, 0.05, "energyPopkJ_mean")

contours<-ggplot() +
geom_tile(data=subset(monthSub_v1, NoBirds_mean>0), aes(x=x, y=y), fill="grey", alpha=0.5)+
geom_contour(data = topBird_contour_v1,
               aes(x = x, y = y, z = z, color="BirdNo"),
               breaks = c(0.05),
               linewidth = 1) +
 geom_contour(data = topEnergy_contour_v1,
               aes(x = x, y = y, z = z, color="MEE_v1"),
               breaks = c(0.05),
               linewidth = 1, alpha=0.3) +
#geom_contour(data = topEnergy_contour_v2,
             #  aes(x = x, y = y, z = z, color="MEE_v2"),
              # breaks = c(0.05),
             #  linewidth = 1, alpha=0.3) +
scale_color_manual(values=c("#0072B2", "#E69F00", "#009E73")) +
  theme_minimal() +
  ggtitle("top 5% data contour lines")


grid.arrange(plot1_energy1, plot1_energy2, plot1_energy_diff, topData, nrow=2)

}

#dev.off()

### Step 2: Add division into LMEs because it looks fun ###

# Open LME data # (modified by Benny) 
lme_polygons <- st_read("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/data/LME/modified_LME.shp")

# Name the three extra areas that Benny added according to this paper (https://doi.org/10.3354/meps13580)
lme_polygons$LME_NAME[is.na(lme_polygons$LME_NAME)] <- c("Labrador_Sea","Central_North_Atlantic",  "Mid-Atlantic")
lme_polygons$ID <- 1:nrow(lme_polygons)

# Intersect this with my data somehow so I can an overview of how much energy is expended per month in every LME
lme_polygons <- st_make_valid(lme_polygons) # To avoid some kind of problem
lme_polygons <- lme_polygons[st_is_valid(lme_polygons), ] # Similar to above
points_sf <- st_as_sf(map_v1, coords = c("x", "y"), crs = st_crs(lme_polygons)) # Look for intersections

# Turn back into a data frame

points_with_LME <- st_join(points_sf, lme_polygons, left = TRUE)

points_with_LME_df <- cbind(
  st_coordinates(points_with_LME),  # adds X and Y columns
  st_drop_geometry(points_with_LME) # adds all other attributes
)

# Summarize energy per LME per month #

energyPerLME<-points_with_LME_df %>%
ungroup() %>%
dplyr::group_by(month, LME_NAME) %>%
dplyr::summarise(totEnergy=sum(energyPopkJ_mean, na.rm=TRUE), totBirds=sum(NoBirds_mean, na.rm=TRUE)) %>%
dplyr::ungroup() %>%
dplyr::group_by(month) %>%
dplyr::mutate(totEnergyAll=sum(totEnergy), totBirdsAll=sum(totBirds)) %>%
dplyr::mutate(propEnergy=totEnergy/totEnergyAll, propBirds=totBirds/totBirdsAll) %>%
dplyr::mutate(LME_NAME2=ifelse(propEnergy<0.1, "Other", LME_NAME)) 

# Plot the results #
stat1<-energyPerLME %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))) %>%
ggplot(aes(x=month, y=propEnergy)) +
geom_col(aes(fill=LME_NAME2)) +
scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#D55E00")) +
labs(fill="LME") +
xlab("Month") +
ylab("Prop MEE") +
theme(legend.position="bottom")

stat2<-energyPerLME %>%
dplyr::mutate(month=factor(month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))) %>%
ggplot(aes(x=month, y=propBirds)) +
geom_col(aes(fill=LME_NAME2)) +
scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#D55E00")) +
labs(fill="LME") +
xlab("Month") +
ylab("Prop tot birds") +
theme(legend.position="bottom")

# Now I add my nice maps of contours to the LME-annotated section #

# I want to summarize by 'winter' #

winterVals_v1<-map_v1 %>%
dplyr::filter(month %in% c(11, 12, 1, 2)) %>%
dplyr::group_by(x, y) %>%
dplyr::summarise(totEnergy=sum(energyPopkJ_mean, na.rm=TRUE), totBirds=sum(NoBirds_mean))

winterVals_v2<-map_v2 %>%
dplyr::filter(month %in% c(11, 12, 1, 2)) %>%
dplyr::group_by(x, y) %>%
dplyr::summarise(totEnergy=sum(energyPopkJ_mean, na.rm=TRUE), totBirds=sum(NoBirds_mean))

# Add map of where hotspots occur # (energy, bird density)

topEnergy_v1<-topEnergy(winterVals_v1, 0.1, "totEnergy") 
topEnergy_v2<-topEnergy(winterVals_v2, 0.1, "totEnergy") 
topBirds<-topEnergy(winterVals_v1, 0.1, "totBirds") 

topData<-ggplot() +
geom_sf(data=lme_polygons, color="black") +
geom_tile(data=subset(monthSub_v1, NoBirds_mean>0), aes(x=x, y=y), fill="grey", alpha=0.5)+
geom_tile(data=topBirds, aes(x=x, y=y, fill="BirdNo"),  alpha=0.6) +
geom_tile(data=topEnergy_v1, aes(x=x, y=y, fill="MEE_v1"),  alpha=0.3) +
#geom_tile(data=topEnergy_v2, aes(x=x, y=y, fill="MEE_v2"),  alpha=0.3) +
scale_fill_manual(values=c("#0072B2", "#E69F00", "#009E73")) +
theme_minimal() +
ggtitle("top 10% data polygons (winter)") +
labs(fill="") 
 

#pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/energy/AllBirds_energy_compare_LME.pdf"), width=10, height=7)

grid.arrange(stat1, stat2, topData, nrow=2)

dev.off()
