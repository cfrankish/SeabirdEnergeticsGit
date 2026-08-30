## Here we aggregate energy maps to make estimates accross all species ##

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
library(biscale)
library(cowplot)
library(classInt)
library(spdep)

#### Step 0: setting up basic conditions ####

# Set-up number of iterations...
print(paste0("Mapping energy..."))

# Source all necessary functions
source("./scripts/functions.R")

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
#input_file1 <- args[1] # This will read in a species name
#print(input_file1)
#input_file1<-"Commonguillemot"

# Now we determine the list of activity/energy data to loop through
allFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables/", full.names=TRUE)
energyFiles<-allFiles[grep("energyMap_monthly", allFiles)] # extract energy files

# Create a world map for plotting
print("Prepping map info")
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#### Step 1: Loop through species files, open & add ###

for (i in 1:length(energyFiles)) {

# Print status #
print(paste0("Opening file", i, "/", length(energyFiles)))

# open file i
energyFile_i<-fread(energyFiles[i])

# if it's the first file, we create a master file #
if (i==1) {

energyAll<-energyFile_i %>%
dplyr::select(month, index, energyPopkJ_mean, NoBirds_mean)

} else {

# Otherwise they are summed together #

energyFile_i<-energyFile_i %>%
dplyr::select(month, index, energyPopkJ_mean, NoBirds_mean)

energyAll<-rbind(energyAll, energyFile_i)

energyAll<-energyAll %>%
ungroup() %>%
dplyr::group_by(month, index) %>%
dplyr::summarise(energyPopkJ_mean=sum(energyPopkJ_mean, na.rm=TRUE), NoBirds_mean=sum(NoBirds_mean, na.rm=TRUE))

}

}

### Step 2: plot outputs ####

print("Step 2: plotting outputs...")

scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Open up random file to get x & y coordinates 
energyFile_i<-fread(energyFiles[1])

coords<-energyFile_i %>%
ungroup() %>%
dplyr::group_by(index) %>%
dplyr::slice(1) %>%
dplyr::select(index, x, y)

# So I will make seperate datasets for energy expenditure, no of birds and the difference between these two

energySpecies1<-energyAll %>%
dplyr::ungroup() %>%
dplyr::left_join(coords, by=c("index")) %>%
dplyr::select(x, y, month, energyPopkJ_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric_scaled=ifelse(energyPopkJ_mean==0, NA, metric_scaled)) %>%
dplyr::mutate(metric="Energy_expenditure") %>%
dplyr::select(-energyPopkJ_mean)

energySpecies2<-energyAll %>%
dplyr::ungroup() %>%
dplyr::left_join(coords, by=c("index")) %>%
dplyr::select(x, y, month, NoBirds_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled=scale01(NoBirds_mean)) %>%
dplyr::mutate(metric_scaled=ifelse(NoBirds_mean==0, NA, metric_scaled)) %>%
dplyr::mutate(metric="No_Birds") %>%
dplyr::select(-NoBirds_mean)

energySpeciesAll<-rbind(energySpecies1, energySpecies2)

energySpecies3<-energyAll %>%
dplyr::ungroup() %>%
dplyr::left_join(coords, by=c("index")) %>%
dplyr::select(x, y, month, energyPopkJ_mean, NoBirds_mean) %>%
dplyr::group_by(month) %>%
dplyr::mutate(metric_scaled1=scale01(energyPopkJ_mean)) %>%
dplyr::mutate(metric_scaled2=scale01(NoBirds_mean)) %>%
dplyr::mutate(metric_diff=metric_scaled1 - metric_scaled2)%>%
dplyr::mutate(metric_diff=ifelse(energyPopkJ_mean==0, NA, metric_diff)) %>%
dplyr::mutate(metric="Energy_expenditure - No_Birds") 

# Now we plot them side by side by month

# Plot these together #

pdf(paste0("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/figures/speciesMaps/energy/allBirds_energy_v1.pdf"), width=10, height=7)

print("Saving PDF")

monthsLoop<-c(9, 10, 11, 12, 3, 4)

for (m in 1:length(monthsLoop)) {

print(paste0("Month", m))

plot1_energy<-ggplot() +
geom_tile(data=filter(energySpeciesAll, metric=="Energy_expenditure" & month==monthsLoop[m]), aes(x=x, y=y, fill=metric_scaled)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, 1), na.value = "grey80") +
facet_wrap(~metric) +
theme_bw() +
ggtitle(paste0("Month ", monthsLoop[m]))

plot1_density<-ggplot() +
geom_tile(data=filter(energySpeciesAll, metric=="No_Birds"& month==monthsLoop[m]), aes(x=x, y=y, fill=metric_scaled)) +
scale_fill_gradientn('scaled 0-1', colors=c("#364B9A", "#4A7BB7","#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"), limits=c(0, 1)) +
facet_wrap(~metric) +
theme_bw()

# Now I make a plot showing the difference between the two

plot2_diff<-ggplot() +
geom_tile(data=filter(energySpecies3, month==monthsLoop[m]), aes(x=x, y=y, fill=metric_diff)) +
scale_fill_gradientn(
  colours = c(
    "#364B9A", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#A50026"
  ),
  limits = c(min(energySpecies3$metric_diff, na.rm=TRUE), max(energySpecies3$metric_diff, na.rm=TRUE)),
  values = scales::rescale(c(min(energySpecies3$metric_diff, na.rm=TRUE), 0, max(energySpecies3$metric_diff, na.rm=TRUE)))
) +
facet_wrap(~metric) +
theme_bw()

grid.arrange(plot1_energy, plot1_density, plot2_diff, nrow=2)

}

dev.off()

# Save outputs

print("Saving outputs...")

output_file1 <- args[2]
write.csv(energyAll, file=output_file1, row.names=FALSE)