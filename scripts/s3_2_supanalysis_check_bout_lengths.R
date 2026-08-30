### This script does some checks for the supplementary section on dry bout allocation ###
## Input file is the id catalogue to know which IDs to loop through ('table1_idcatalogue.csv') ###
### It looks into the proportion of dry bouts randomly allocated to land vs. flight (which happens in some cases) -> "./results/tables/supplementary/table2_timeRandom.csv" ###
### It also measures the duration of flight bouts during darkness to make sure these are realistic -> "./results/tables/supplementary/table3_timeDark.csv" ###
### This script also outputs three figures: Figures S5-S7 showing the results of the two analyses above ###

library(ggplot2)
library(dplyr)
library(fields)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(fasterize)
library(gdistance)
library(sf)
library(sp)
library(tidyr)
library(suncalc)
library(lubridate)
library(data.table)
library(terra)
library(ncdf4)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

### Step s3_2_0: set-up sample size & iteration number ###

print("Step s3_2_0: setting up initial parameters")

# Set up minimum sample size & number of iterations
minSampleSize<-5
print(paste0("min sample size per colony is: ", minSampleSize))
reps<-50
print(paste0("min iteration number is: ", reps))

### Step s3_2_1: read in id catalogue ###

# This is so we can loop through them later #

print("Step s3_2_1: Load id catalogue")
input_file <- args[1]
energyAll <- read.csv(input_file)

### Step 1: Estimate some important metrics ###

print("Step 2: Estimate prop bouts attributed at random & at night...")

# This first one estimates proportion of dry bouts allocated at random per species and Month. 
# This script also calculates flight time during darkness for kittiwakes & fulmars. 

allResults<-list.files("./tmp", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# All birds: prop bouts in day allocated at random #
timeRandom<-list()
timeDark<-list()

speciesList<-unique(energyAll$species)

for (i in 1:length(speciesList)) {

# subset to species i
speciesSub<-speciesList[i]
speciesSub_df<-subset(energyAll, species==speciesSub)

# Subset relevant Ids to loop through
ids<-unique(speciesSub_df$individ_id)

# Save species-specific results in a list
timeRandom_species<-list()
timeDark_species<-list()

for (j in 1:length(ids)) {

print(paste0("Step s3_2_2: Prop dry bouts random, Species ", i, ", Bird ", j))

# Subset to id j
idSub<-ids[j]

# Open id j
birdSub<-read.csv(energyRes_day[grepl(ids[j], energyRes_day)])
      
# subset to correct months
birdSub_months<-birdSub %>%
dplyr::filter(Month %in% c(9, 10, 11, 12, 1, 2, 3, 4))

if(nrow(birdSub_months)<1) {
next
}

# Add 'boutsRandom' as a column if it doesn't exist (this is because BLki have no bouts attributed at random...)
if ("boutsRandom" %in% colnames(birdSub_months) == FALSE) {
      
      birdSub_months$boutsRandom<-0
    }
	
if ("propflight_dark" %in% colnames(birdSub_months) == FALSE) {
      
      birdSub_months$propflight_dark<-0
    }	
	
if ("flightTimeMins_dark" %in% colnames(birdSub_months) == FALSE) {
      
      birdSub_months$flightTimeMins_dark<-0
    }	
      
# Summarise average prop bout random & range
birdSub_months_random<-birdSub_months %>%
ungroup() %>%
dplyr::group_by(rep, species, Month) %>%
dplyr::summarise(meanRandom=mean(boutsRandom), minRandom=min(boutsRandom), maxRandom=max(boutsRandom))

# Summarise time spent flying during darkness
birdSub_months_random_dark<-birdSub_months %>%
dplyr::filter(tDarkness>0) %>%
ungroup() %>%
dplyr::group_by(rep, species, Month) %>%
replace_na(list(propflight_dark=0, flightTimeMins_dark=0)) %>%
dplyr::summarise(meanpropflight_dark=mean(propflight_dark), 
meanflightTimeMins_dark=mean(flightTimeMins_dark))

# Save resutls
timeRandom_species<-rbind(timeRandom_species, birdSub_months_random)
timeDark_species<-rbind(timeDark_species, birdSub_months_random_dark)

}

# Summarize so it isn't massive
timeRandom_species_mean<-timeRandom_species %>%
ungroup() %>%
dplyr::group_by(species, Month) %>%
dplyr::summarise(meanRandomOverall=mean(meanRandom), sdRandom=sd(meanRandom))

timeDark_species_mean<-timeDark_species %>%
ungroup() %>%
dplyr::group_by(species, Month) %>%
dplyr::summarise(meanpropflightOverall=mean(meanpropflight_dark), sdpropflight=sd(meanpropflight_dark),
meanflighttimeOverall=mean(meanflightTimeMins_dark), sdflighttime=sd(meanflightTimeMins_dark))

# Save results for all birds #
timeRandom<-rbind(timeRandom, timeRandom_species_mean)
timeDark<-rbind(timeDark, timeDark_species_mean)

}

output_file1 <- args[2]
print("Saving output file 1")
write.csv(timeRandom, file = output_file1, row.names = FALSE) # Monthly activity data

# Save output as plot (proportion of dry bouts allocated at random)
FigureS7<-timeRandom %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(Month=factor(Month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))) %>%
ggplot() +
geom_pointrange(aes(x=Month, y=meanRandomOverall, ymin=meanRandomOverall - sdRandom, ymax=meanRandomOverall + sdRandom, color=species, group=interaction(species, Month)), position=position_dodge2(width=0.6)) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822")) +
ylim(-0.1, 1) +
theme_bw()+
ylab("Proportion dry bouts allocated at random (mean +/- SD)") +
xlab("Month")+
labs(color="Species")

pdf("./results/figures/supplementary/FigureS7.pdf", width=8, height=5)
grid.arrange(FigureS7)
dev.off()

# Save output as plot (flightTime)
FigureS5a<-timeDark %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(Month=factor(Month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))) %>%
dplyr::filter(species %in% c("Black-legged kittiwake", "Northern fulmar")) %>%
ggplot() +
geom_pointrange(aes(x=Month, y=meanpropflightOverall, ymin=meanpropflightOverall - sdpropflight, ymax=meanpropflightOverall + sdpropflight, color=species, group=interaction(species, Month)), position=position_dodge2(width=0.6)) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822")) +
#ylim(0, 1) +
theme_bw()+
ylab("Proportion darkness spent in flight") +
xlab("Month")+
labs(color="Species", tag="A)")

FigureS5b<-timeDark %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(Month=factor(Month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))) %>%
dplyr::filter(species %in% c("Black-legged kittiwake", "Northern fulmar")) %>%
ggplot() +
geom_pointrange(aes(x=Month, y=meanflighttimeOverall, ymin=meanflighttimeOverall - sdflighttime, ymax=meanflighttimeOverall + sdflighttime, color=species, group=interaction(species, Month)), position=position_dodge2(width=0.6)) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822")) +
#ylim(0, 1) +
theme_bw()+
ylab("Flight time during darkness (mins)") +
xlab("Month")+
labs(color="Species", tag="B)")

pdf("./results/figures/supplementary/FigureS5.pdf")
grid.arrange(FigureS5a, FigureS5b, nrow=2)
dev.off()

output_file2 <- args[3]
print("Saving output file 2")
write.csv(timeDark, file = output_file2, row.names = FALSE) # Monthly activity data

### Step s3_2_3: Check duration of fulmar dry bouts ####

print("Step s3_2_3: Checking duration of fulmar dry bouts...")

# Determine location of fulmar dry bouts #
loxBouts<-list.files("./results/tables/supplementary/fulmarBoutLengths/", full.name=TRUE)

# Create data frame with individual id & rep # no
loxBoutsDf<-data.frame(fileNames=loxBouts)
loxBoutsDf$individ_id<-sub(".*flightbouts_(.*?)_rep.*", "\\1", loxBoutsDf$fileNames) # extract ID name
loxBoutsDf$rep<-as.numeric(sub(".*rep(.*?).csv.*", "\\1", loxBoutsDf$fileNames)) # extract repetition

# Subset to 50 reps per individual
loxBoutsDf_50<-subset(loxBoutsDf, rep<=50)
loxBoutsDf_50<-loxBoutsDf_50 %>%
dplyr::mutate(birdNo=as.numeric(factor(individ_id))) %>%
arrange(birdNo, rep)

# Create empty list to save results in
flightLengths<-list()

# Create a loop to go through these #
for (i in 1:nrow(loxBoutsDf_50)) {

print(paste0("Step s3_2_3: fulmar bout, Bird ", loxBoutsDf_50$birdNo[i], " /", n_distinct(loxBoutsDf_50$birdNo), " , rep ", loxBoutsDf_50$rep[i], "/", reps))

# Open bout from individual i
loxBoutsSub<-fread(loxBoutsDf_50$fileNames[i])

# Choose only necessary information
loxBoutsSub_info<-loxBoutsSub %>%
dplyr::filter(month %in% c(9, 10, 11, 12, 1, 2, 3, 4)) %>%
dplyr::filter(Activity=="Flight")%>%
dplyr::group_by(BoutNo) %>%
dplyr::summarise(times=n_distinct(date_time))%>%
dplyr::mutate(flightLengthMins=times*10)

print(paste("Bird: ", loxBoutsDf_50$individ_id[i], " Max flight: ", max(loxBoutsSub_info$flightLengthMins), " mins"))

# Save all results
flightLengths<-rbind(flightLengths, loxBoutsSub_info)

}

# Save as output # 
saveRDS(flightLengths, file="./results/tables/supplementary/flightBoutLengths_fulmars.rds")

# Make a ggplot to save
FigureS6<-flightLengths %>%
ggplot() +
geom_histogram(aes(x=flightLengthMins)) +
theme_bw() +
xlab("Duration of flight bouts (mins)")

pdf("./results/figures/supplementary/FigureS6.pdf")
plot(FigureS6)
dev.off()