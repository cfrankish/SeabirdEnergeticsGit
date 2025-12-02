### Here I plot how each parameter affects TEE_NB (total energy expended) & WEE_COV_NB (variation in weekly energy expenditure) ###
### There are official input files but they are individual daily energy files (tmp2/id_energyDay.csv) produced by the previous script ###
### Output file is a csv showing how change in each parameter affects TEE_NB & WEE_COV_NB: ("./results/tables/supplementary/table0_sensitivityResults.csv") ###
### This step also outputs two supplementary figures (XXX) ###

# load functions
library(ggplot2)
library(dplyr)
library(fields)
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
library(tibble)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

### Step 1: Determine location of sensitivity files... ###

print(paste0("Step 1: determine location of sensitivity files..."))

lox.results<-data.frame(files=list.files("tmp2/", full.names=TRUE))
lox.results$size<-file.info(lox.results$files)$size

# remove files with nothing inside
lox.results.loop<-subset(lox.results, size>3)

### Step 2: loop through results ###

# Here I go through every individual file, subset to study period, and calculate total time spent activities + energy expended ###
# I do this for every parameter - run combination #

# Determine study period
dates<-data.frame(dateKeep=seq(as.Date("2021-09-15"), as.Date("2022-04-15"), 1))
dates$doy<-1:nrow(dates)
dates$month<-as.numeric(substr(dates$date, 6, 7))
dates$day<-as.numeric(substr(dates$date, 9, 10))

# Add week number for summarizing information
dates_weekly<-dates %>%
  dplyr::mutate(weekNo=ceiling(doy/7)) %>%
  dplyr::group_by(weekNo) %>%
  dplyr::mutate(days=n_distinct(dateKeep)) %>%
  dplyr::filter(days==7) %>%
  dplyr::select(-days)
 
day1<-as.Date(min(dates_weekly$dateKeep)-1)
day2<-as.Date(max(dates_weekly$dateKeep)+1)
date1<-data.frame(dateKeep=day1, doy=0, month=as.numeric(substr(day1, 6, 7)), day=as.numeric(substr(day1, 9, 10)), weekNo=0)
date2<-data.frame(dateKeep=day2, doy=0, month=as.numeric(substr(day2, 6, 7)), day=as.numeric(substr(day2, 9, 10)), weekNo=0)
dates_weekly2<-rbind(dates_weekly)
  
# divide DEE by weight to get in kJ.g
species<-data.frame(species=c("Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot", "Little auk", "Atlantic puffin"))
species$allometryCoef<-c(0.765, 0.717, 0.689, 0.689, 0.689, 0.689)

# Make a list to save results in
sensitivityRes<-list()

print("Step 2: calculate total time in activity + energy expended...")

for (i in 1:nrow(lox.results.loop)) {

print(paste0("Bird", i, "/", nrow(lox.results.loop)))

# Open up file #
print("Opening file")

sensRes<-read.csv(lox.results.loop[i,]$files[1])
print(sensRes$species[1])

# Subset to dates of interest
sensRes$day<-as.numeric(substr(sensRes$date, 9, 10))
sensRes$month<-as.numeric(substr(sensRes$date, 6, 7))
sensResStudy<-sensRes %>%
ungroup() %>%
dplyr::inner_join(dates_weekly2, by=c("month", "day")) %>%
dplyr::group_by(Parameter, Run) %>%
dplyr::mutate(totalDays=n_distinct(date)) 

# Make sure the number of days is correct
print("Subsetting to study period")
studyPeriod<-unique(sensResStudy$totalDays)
if (length(studyPeriod)>1 || ! studyPeriod %in% c(nrow(dates_weekly2))) {
 (print("Error: number of days off"))
 next
}

# add tactive as a column
if (!any(colnames(sensResStudy)%in% c("tActive"))) {
sensResStudy$tActive<-0
}

# Calculate output metrics for sensitivity analysis
print("Calculating output metrics")
sensResOutputs<-sensResStudy %>%
dplyr::filter(weekNo>0) %>% # remove first & last buffer days
ungroup() %>%
dplyr::left_join(species, by=c("species")) %>%
dplyr::mutate(DEE=DEEkJ/weight^allometryCoef) %>%
dplyr::group_by(species, colony, individ_id, Parameter, Run, weekNo, L1, Th1, Th2, L1_colony_min, L1_colony_max, dist_colony, pLand_prob, c, RMR, 
c1, c2, c3, c4, TC, Beta_active, Beta_rest, year) %>%
dplyr::summarize(weeklyDEE=sum(DEE), tFlight=sum(tFlight), tForage=sum(tForage), tActive=sum(tActive), tLand=sum(tLand), tRestWater=sum(tRestWater)) %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony, individ_id, Parameter, Run, L1, Th1, Th2, L1_colony_min, L1_colony_max,dist_colony, pLand_prob, c, RMR, 
c1, c2, c3, c4, TC, Beta_active, Beta_rest, year) %>%
dplyr::mutate(weeklyDEE_nb=mean(weeklyDEE), sdweeklyDEE_nb=sd(weeklyDEE), devianceDEE=abs(weeklyDEE-weeklyDEE_nb)/weeklyDEE_nb) %>%
dplyr::summarise(devianceDEE_nb=sum(devianceDEE), TEE_nb=sum(weeklyDEE), tFlight_nb=sum(tFlight),
tForage_nb=sum(tForage), tActive_nb=sum(tActive), tLand_nb=sum(tLand), tRestWater_nb=sum(tRestWater), DEE_cov_nb=(sd(weeklyDEE)/mean(weeklyDEE))*100)

# Calculate difference compared to baseline scenario
sensResCompare<-sensResOutputs %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id, Parameter) %>%
dplyr::mutate(devianceDEE_nb_diff=(devianceDEE_nb-first(devianceDEE_nb))/first(devianceDEE_nb), TEE_nb_diff=(TEE_nb-first(TEE_nb))/first(TEE_nb),
tFlight_nb_diff=(tFlight_nb-first(tFlight_nb))/first(tFlight_nb), tForage_nb_diff=(tForage_nb-first(tForage_nb))/first(tForage_nb), tLand_nb_diff=(tLand_nb-first(tLand_nb))/first(tLand_nb),
tActive_nb_diff=(tActive_nb-first(tActive_nb))/first(tActive_nb), tRest_nb_diff=(tRestWater_nb-first(tRestWater_nb))/first(tRestWater_nb), DEE_cov_nb_diff=((DEE_cov_nb)-first(DEE_cov_nb))/(first(DEE_cov_nb))) %>%
dplyr::mutate(TEE_nb_diff=ifelse(Parameter=="year", abs(TEE_nb_diff), TEE_nb_diff))

# Save results
print("Saving result")
sensitivityRes<-rbind(sensitivityRes, sensResCompare)

}

saveRDS(sensitivityRes, file="birdsRound2.rds")

# Determine order of parameters for plotting #
paramOrder<-c("year", "L1", "Th1", "Th2", "L1_colony", "ice", "dist_colony",
"pLand_prob", "c", "sst", "RMR", "c1", "c2", "c3", "c4", "TC", "Beta_active", "Beta_rest")

# Plot results #
print("Plotting results")

sensitivityRes$Run2<-ifelse(sensitivityRes$Parameter =="year" & sensitivityRes$Run > 1, 2, sensitivityRes$Run)

sensitivityResSum<-sensitivityRes %>%
ungroup() %>%
dplyr::group_by(species, Parameter, Run2) %>%
summarise(across(devianceDEE_nb_diff:DEE_cov_nb_diff, list(
                     mean = ~ mean(.x, na.rm = TRUE),
                     sd   = ~ sd(.x, na.rm = TRUE),
                     n    = ~ sum(!is.na(.x)),
                     se   = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
                     ci_lower = ~ mean(.x, na.rm = TRUE) - 
                                 1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
                     ci_upper = ~ mean(.x, na.rm = TRUE) + 
                                 1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
                   ))) %>%
dplyr::mutate(Parameter=factor(Parameter, levels=rev(paramOrder)))	%>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))	   
												   
# Calculate range of differences #
rangeDiff<-sensitivityResSum %>%
ungroup() %>%
dplyr::filter(Run2>1) %>%
dplyr::group_by(species, Parameter) %>%
dplyr::summarise(mean_change1=mean(abs(devianceDEE_nb_diff_mean)), mean_change2=mean(abs(TEE_nb_diff_mean))) %>%
dplyr::filter(mean_change1>0)

rangeDiffMean<-rangeDiff %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::summarise(mean_change1_all=mean(mean_change1), minci1=min(mean_change1), maxci1=max(mean_change1), mean_change2_all=mean(mean_change2), minci2=min(mean_change2), maxci2=max(mean_change2))

# Seperate these further so everything can be facet-wrapped?
output0<-sensitivityResSum %>%
dplyr::select(species, Parameter, Run2, devianceDEE_nb_diff_mean, devianceDEE_nb_diff_ci_lower, devianceDEE_nb_diff_ci_upper) %>%
dplyr::rename(mean=devianceDEE_nb_diff_mean, ci_lower=devianceDEE_nb_diff_ci_lower, ci_upper=devianceDEE_nb_diff_ci_upper) %>%
dplyr::mutate(outputMetric="WEE_deviance_nb")

output1<-sensitivityResSum %>%
dplyr::select(species, Parameter, Run2, DEE_cov_nb_diff_mean, DEE_cov_nb_diff_ci_lower, DEE_cov_nb_diff_ci_upper) %>%
dplyr::rename(mean=DEE_cov_nb_diff_mean, ci_lower=DEE_cov_nb_diff_ci_lower, ci_upper=DEE_cov_nb_diff_ci_upper) %>%
dplyr::mutate(outputMetric="WEE_cov_nb")

output2<-sensitivityResSum %>%
dplyr::select(species, Parameter, Run2, TEE_nb_diff_mean, TEE_nb_diff_ci_lower, TEE_nb_diff_ci_upper) %>%
dplyr::rename(mean=TEE_nb_diff_mean, ci_lower=TEE_nb_diff_ci_lower, ci_upper=TEE_nb_diff_ci_upper) %>%
dplyr::mutate(outputMetric="TEE_nb")

#output3<-sensitivityResSum %>%
#dplyr::select(species, Parameter, Run, tFlight_nb_diff_mean, tFlight_nb_diff_ci_lower, tFlight_nb_diff_ci_upper) %>%
#dplyr::rename(mean=tFlight_nb_diff_mean, ci_lower=tFlight_nb_diff_ci_lower, ci_upper=tFlight_nb_diff_ci_upper) %>%
#dplyr::mutate(outputMetric="TFlight_nb")

#output4<-sensitivityResSum %>%
#dplyr::select(species, Parameter, Run, tForage_nb_diff_mean, tForage_nb_diff_ci_lower, tForage_nb_diff_ci_upper) %>%
#dplyr::rename(mean=tForage_nb_diff_mean, ci_lower=tForage_nb_diff_ci_lower, ci_upper=tForage_nb_diff_ci_upper) %>%
#dplyr::mutate(outputMetric="TForage_nb")

#output5<-sensitivityResSum %>%
#dplyr::select(species, Parameter, Run, tLand_nb_diff_mean, tLand_nb_diff_ci_lower, tLand_nb_diff_ci_upper) %>%
#dplyr::rename(mean=tLand_nb_diff_mean, ci_lower=tLand_nb_diff_ci_lower, ci_upper=tLand_nb_diff_ci_upper) %>%
#dplyr::mutate(outputMetric="TLand_nb")

#output6<-sensitivityResSum %>%
#dplyr::select(species, Parameter, Run, tActive_nb_diff_mean, tActive_nb_diff_ci_lower, tActive_nb_diff_ci_upper) %>%
#dplyr::rename(mean=tActive_nb_diff_mean, ci_lower=tActive_nb_diff_ci_lower, ci_upper=tActive_nb_diff_ci_upper) %>%
#dplyr::mutate(outputMetric="TActive_nb")

#output7<-sensitivityResSum %>%
#dplyr::select(species, Parameter, Run, tRest_nb_diff_mean, tRest_nb_diff_ci_lower, tRest_nb_diff_ci_upper) %>%
#dplyr::rename(mean=tRest_nb_diff_mean, ci_lower=tRest_nb_diff_ci_lower, ci_upper=tRest_nb_diff_ci_upper) %>%
#dplyr::mutate(outputMetric="TWater_nb")

sensitivityPlot<-rbind(output2)

# Seperate runs	
sens1<-subset(sensitivityPlot, Run2==2 & !Parameter %in% c("year"))
sens2<-subset(sensitivityPlot, Run2==3 & !Parameter %in% c("year"))
sens3<-subset(sensitivityPlot, !Run2 %in% c(1) & Parameter %in% c("year"))

sens1<-sens1 %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))	
sens2<-sens2 %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))	
sens3<-sens3 %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))													   

# Plot 1

# Position dodge object
pd <- position_dodge(width = 0.2)

sensitivitySave1<-ggplot() +
geom_col(data=sens1, aes(x=Parameter, y=mean, fill=species, group=interaction(species, Parameter, outputMetric)), position=position_dodge2(width=0.2)) +
geom_col(data=sens2, aes(x=Parameter, y=mean, fill=species, group=interaction(species, Parameter, outputMetric)), position=position_dodge2(width=0.2)) +
geom_col(data=sens3, aes(x=Parameter, y=mean, fill=species, group=interaction(species, Parameter, outputMetric)), position=position_dodge2(width=0.2)) +
geom_errorbar(data=sens1, aes(x=Parameter, ymin=ci_lower, ymax=ci_upper, group=species), width=0.2, position=position_dodge2(width=0.2)) +
geom_errorbar(data=sens2, aes(x=Parameter, ymin=ci_lower, ymax=ci_upper, group=species), width=0.2, position=position_dodge2(width=0.2)) +
geom_errorbar(data=sens3, aes(x=Parameter, ymin=ci_lower, ymax=ci_upper,  group=species), width=0.2, position=position_dodge2(width=0.2)) +
coord_flip() +
ylim(-1.2, 1.2) +
geom_hline(yintercept=0, linetype="dashed") +
#scale_fill_manual(values=c("grey")) +
theme_bw() +
facet_wrap(~outputMetric + species) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab("")

pdf("results/figures/supplementary/FigureS10.pdf")
plot(sensitivitySave1)
dev.off()

# Save the other plot #

sensitivityPlot2<-rbind(output1)

# Seperate runs	
sens1<-subset(sensitivityPlot2, Run2==2 & !Parameter %in% c("year"))
sens2<-subset(sensitivityPlot2, Run2==3 & !Parameter %in% c("year"))
sens3<-subset(sensitivityPlot2, !Run2 %in% c(1) & Parameter %in% c("year"))

sens1<-sens1 %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))	
sens2<-sens2 %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))	
sens3<-sens3 %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                   "Little auk", "Common guillemot", "Brünnich's guillemot")))													   

# Plot 1

# Position dodge object
pd <- position_dodge(width = 0.2)

sensitivitySave2<-ggplot() +
geom_col(data=sens1, aes(x=Parameter, y=mean, fill=species, group=interaction(species, Parameter, outputMetric)), position=position_dodge2(width=0.2)) +
geom_col(data=sens2, aes(x=Parameter, y=mean, fill=species, group=interaction(species, Parameter, outputMetric)), position=position_dodge2(width=0.2)) +
geom_col(data=sens3, aes(x=Parameter, y=mean, fill=species, group=interaction(species, Parameter, outputMetric)), position=position_dodge2(width=0.2)) +
geom_errorbar(data=sens1, aes(x=Parameter, ymin=ci_lower, ymax=ci_upper, group=species), width=0.2, position=position_dodge2(width=0.2)) +
geom_errorbar(data=sens2, aes(x=Parameter, ymin=ci_lower, ymax=ci_upper, group=species), width=0.2, position=position_dodge2(width=0.2)) +
geom_errorbar(data=sens3, aes(x=Parameter, ymin=ci_lower, ymax=ci_upper,  group=species), width=0.2, position=position_dodge2(width=0.2)) +
coord_flip() +
ylim(-1.2, 1.2) +
geom_hline(yintercept=0, linetype="dashed") +
#scale_fill_manual(values=c("grey")) +
theme_bw() +
facet_wrap(~outputMetric + species) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab("")

pdf("results/figures/supplementary/FigureS9.pdf")
plot(sensitivitySave2)
dev.off()

# Save intermediate output file #
print("Saving all results")

output_file1 <- args[1]
write.csv(sensitivityRes, file = output_file1, row.names = FALSE) # Daily energy data 






