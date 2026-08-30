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
input_file3<-args[3]
input_file4<-args[4]

# open these files
map_v1<-fread(input_file1) # map_v1<-fread("results/tables/Commonguillemot_energyMap_monthly_map_v1.csv")
map_v2<-fread(input_file2) # map_v2<-fread("results/tables/Commonguillemot_energyMap_monthly_map_v2.csv")
sum_v1<-fread(input_file3) # sum_v1<-fread("results/tables/Commonguillemot_energyMap_monthly_sum_v1.csv")
sum_v2<-fread(input_file4) # sum_v2<-fread("results/tables/Commonguillemot_energyMap_monthly_sum_v2.csv")

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

pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/energy/", map_v1$species[1], "_energy_compare.pdf"), width=10, height=7)

for (m in 1:length(monthLoop)) {

print(m)

# Subset maps to month i
monthSub_v1<-subset(map_v1, month==monthLoop[m])
monthSub_v2<-subset(map_v2, month==monthLoop[m])

# Subset sums to month i
monthSub_v1_sum<-subset(sum_v1, month==monthLoop[m])
monthSub_v2_sum<-subset(sum_v2, month==monthLoop[m])

# Make scaled datasets to compare

energySpecies1<-monthSub_v1 %>%
#dplyr::filter(energyPopkJ_mean>0) %>%
dplyr::ungroup() %>%
dplyr::select(species, x, y, month, energyPopkJ_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric="Energy_expenditure_v1") %>%
dplyr::select(-energyPopkJ_mean)

energySpecies2<-monthSub_v2 %>%
#dplyr::filter(energyPopkJ_mean>0) %>%
dplyr::ungroup() %>%
dplyr::select(species, x, y, month, energyPopkJ_mean) %>%
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

# Make a plot with the sums in it
monthSub_v1_sum$method<-"V1"
monthSub_v2_sum$method<-"V2"
sumsAll<-rbind(monthSub_v1_sum, monthSub_v2_sum)

totalEnergy<-ggplot() +
geom_pointrange(data=sumsAll, aes(x=method, y=energyMean/1e09, ymin=(energyMean-energySD)/1e09, ymax=(energyMean + energySD)/1e09)) +
xlab("Energy method") +
ylab("Energy expenditure kJ.month.10e9") +
theme_bw() +
ggtitle("Monthly energy expenditure")

# Plot top 25% of energy in both cases #
top25_v1<-topEnergy(monthSub_v1, 0.25, 'energyPopkJ_mean')
top25_v2<-topEnergy(monthSub_v2, 0.25, 'energyPopkJ_mean')

plot_compare<-ggplot() +
 geom_tile(data=top25_v1, aes(x=x, y=y, fill="25%_v1")) +
 geom_tile(data=top25_v2, aes(x=x, y=y, fill="25%_v2")) +
 scale_fill_manual(values=c("#0072B2", "#E69F00")) +
 theme_minimal() +
 ggtitle("Top 25% of MEE") +
 labs(fill="")

grid.arrange(plot1_energy1, plot1_energy2, plot1_energy_diff, plot_compare, nrow=2)

}

dev.off()

### Now we do a second loop to compare DEE not multiplied by number of birds ###
# The aim here is to see if there are any differences in energy expenditure when using spatial vs non-spatial activity #

monthLoop<-c(9, 10, 11, 12, 1, 2, 3, 4)

pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/energy/", map_v1$species[1], "_energy_compare_MEE.pdf"), width=10, height=7)

for (m in 1:length(monthLoop)) {

print(m)

# Subset maps to month i
monthSub_v1<-subset(map_v1, month==monthLoop[m])
monthSub_v2<-subset(map_v2, month==monthLoop[m])

# Subset sums to month i
monthSub_v1_sum<-subset(sum_v1, month==monthLoop[m])
monthSub_v2_sum<-subset(sum_v2, month==monthLoop[m])

# Make scaled datasets to compare

energySpecies1<-monthSub_v1 %>%
#dplyr::filter(energyPopkJ_mean>0) %>%
dplyr::ungroup() %>%
dplyr::select(species, x, y, month, DEEkJ_mean, NoBirds_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(DEEkJ_mean)) %>%
dplyr::mutate(metric="Energy_expenditure_v1") %>%
dplyr::select(-DEEkJ_mean)

energySpecies2<-monthSub_v2 %>%
#dplyr::filter(energyPopkJ_mean>0) %>%
dplyr::ungroup() %>%
dplyr::select(species, x, y, month, DEEkJ_mean, NoBirds_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(DEEkJ_mean)) %>%
dplyr::mutate(metric="Energy_expenditure_v2") %>%
dplyr::select(-DEEkJ_mean)

energySpecies1$metric_scaled2<-energySpecies2$metric_scaled
energySpecies1$diff<-energySpecies1$metric_scaled - energySpecies1$metric_scaled2

# Figuring out a common extent #
energyBirds<-subset(energySpecies1, NoBirds_mean>0)
min1<-min(energyBirds$metric_scaled)
min2<-min(energyBirds$metric_scaled2)
min<-ifelse(min1<min2, min1, min2)
max1<-max(energyBirds$metric_scaled)
max2<-max(energyBirds$metric_scaled2)
max<-ifelse(max1>max2, max1, max2)

plot1_energy1<-ggplot() +
geom_tile(data=subset(energySpecies1, NoBirds_mean>0), aes(x=x, y=y, fill=metric_scaled)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limit=c(min, max)) +
#facet_wrap(~metric) +
theme_bw() +
ggtitle(paste0("V1: Month ", monthLoop[m]))

plot1_energy2<-ggplot() +
geom_tile(data=subset(energySpecies1, NoBirds_mean>0), aes(x=x, y=y, fill=metric_scaled2)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limit=c(min, max)) +
#facet_wrap(~metric) +
theme_bw() +
ggtitle(paste0("V2: Month ", monthLoop[m]))

plot1_energy_diff<-ggplot() +
geom_tile(data=subset(energySpecies1, NoBirds_mean>0), aes(x=x, y=y, fill=diff)) +
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

# Make a plot with the sums in it
monthSub_v1_sum$method<-"V1"
monthSub_v2_sum$method<-"V2"
sumsAll<-rbind(monthSub_v1_sum, monthSub_v2_sum)

totalEnergy<-ggplot() +
geom_pointrange(data=sumsAll, aes(x=method, y=energyMean/1e09, ymin=(energyMean-energySD)/1e09, ymax=(energyMean + energySD)/1e09)) +
xlab("Energy method") +
ylab("Energy expenditure kJ.month.10e9") +
theme_bw() +
ggtitle("Monthly energy expenditure")

grid.arrange(plot1_energy1, plot1_energy2, plot1_energy_diff, totalEnergy, nrow=2)

}

dev.off()

# Save something #
sum_v1$method<-"V1"
sum_v2$method<-"V2"
sumsFinal<-rbind(sum_v1, sum_v2)
sumsFinal$month<-factor(sumsFinal$month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))


temporalChange<-ggplot() +
geom_pointrange(data=sumsFinal, aes(x=factor(month), group=method, color=method, y=energyMean/1e09, ymin=(energyMean-energySD)/1e09, ymax=(energyMean + energySD)/1e09), position=position_dodge2(width=0.2)) +
geom_ribbon(data=sumsFinal, aes(x=month, group=method, fill=method, y=energyMean/1e09, ymin=(energyMean-energySD)/1e09, ymax=(energyMean + energySD)/1e09), position=position_dodge2(width=0.2), alpha=0.2) +
geom_smooth(data=sumsFinal, aes(x=month, y=energyMean/1e09, group=method)) +
xlab("Energy method") +
ylab("Energy expenditure kJ.month.10e9") +
theme_bw() +
ggtitle("Monthly energy expenditure")

print("Saving results...")

# Save result
output_file1 <- args[5]

# Save metadata as main output file
write.csv(sumsFinal, file=output_file1, row.names=FALSE)



