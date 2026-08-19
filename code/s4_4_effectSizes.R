# This code looks into how many reps are needed to accurately estimate effect sizes #
# It is an additional figure I added at a later point in response to a reviewer comment #
# This file has as input files all the stats I conduct in the previous step:
#input_files1 = f"./results/tables/main/table6_totalNBCosts.csv"
#input_files2 = f"./results/tables/main/table7_stats_WEE_vs_TEE.csv"
#input_files3 = f"./results/tables/main/table8_stats_WEE_vs_pred.csv"
#input_files4 = f"./results/tables/main/table9_stats_TEE_vs_pred.csv"
# And outputs Figure S13

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
library(lme4)
#library(MuMIn)
#library(pracma)
#library(lmerTest)
library(cluster)
library(ggrepel)
library(tibble)
library(rstatix)

args <- commandArgs(trailingOnly = TRUE) # This allows R to read in arguments written in the workflow file

### Step s4_4_0: set-up sample size & iteration number ###

print("Step s4_4_0: setting up initial parameters")

# Set up minimum sample size & number of iterations
minSampleSize<-5
print(paste0("min sample size per colony is: ", minSampleSize))
reps<-50
print(paste0("min iteration number is: ", reps))

### Step s4_4_1: open up all input files ###

# Assign arguments to input files 
input_file1 <- args[3]
input_file2 <- args[4]
input_file3 <- args[5]
input_file4 <- args[6]

# open these
totalCosts <- read.csv(input_file1)
lmRes_tot<-read.csv(input_file2)# WEE_cov_nb vs TEE_nb
lmRes_tot2<-read.csv(input_file3) # WEE_cov_nb vs other behaviors 
lmRes_tot3<-read.csv(input_file4) # TEE_nb vs other predictors

### Step s4_4_2: Calculate effect sizes for an increasing number of reps ###

results_h1_h2_all<-list()
results_h3_all<-list()
results_h4_all<-list()

for (i in 2:reps) {

print(i)

# Summarize effect sizes for H1 & H2 #
results_h1_h2<-totalCosts %>%
dplyr::filter(rep<=i) %>%
dplyr::group_by(rep, species) %>%
dplyr::slice(1) %>%
dplyr::select(rep, species, popTest_nbcosts_effectSize, popTest_cov_effectSize) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::summarise(repsTot=n_distinct(reps), eff_TEE_nb_mean=mean(popTest_nbcosts_effectSize), eff_TEE_nb_median=median(popTest_nbcosts_effectSize), eff_TEE_nb_sd=sd(popTest_nbcosts_effectSize), 
eff_TEE_nb_se=eff_TEE_nb_sd/sqrt(repsTot), eff_TEE_nb_low=eff_TEE_nb_median - 1.96*eff_TEE_nb_se, eff_TEE_nb_high=eff_TEE_nb_median + 1.96*eff_TEE_nb_se,
eff_WEE_cov_nb_mean=mean(popTest_cov_effectSize), eff_WEE_cov_nb_median=median(popTest_cov_effectSize), eff_WEE_cov_nb_sd=sd(popTest_cov_effectSize), 
eff_WEE_cov_nb_se=eff_WEE_cov_nb_sd/sqrt(repsTot), eff_WEE_cov_nb_low=eff_WEE_cov_nb_median - 1.96*eff_WEE_cov_nb_se, eff_WEE_cov_nb_high=eff_WEE_cov_nb_median + 1.96*eff_WEE_cov_nb_se) %>%
ungroup() %>%
dplyr::mutate(reps=i)

# Summarize effect sizes for H3 #
results_h3<-lmRes_tot %>%
  dplyr::filter(rep<=i) %>%
  dplyr::group_by(species, type, predictor, predictor_val) %>%
  dplyr::summarise(meanEstimate=mean(fit), meanSE=mean(sevals), reps=n_distinct(rep), 
                   meanr2=mean(r2), sdr2=sd(r2), ser2=sdr2/sqrt(reps), meanCoef=mean(coefficient), sdCoef=sd(coefficient), seCoef=sdCoef/sqrt(reps), within_model_var=mean(se^2), between_model_var=var(coefficient),
				   totalVar=within_model_var + (1+1/reps)*between_model_var, se_pooled=sqrt(totalVar)) %>%
  dplyr::mutate(low=meanCoef-1.96*se_pooled, high=meanCoef + 1.96*se_pooled, sig=ifelse(meanCoef - 1.96*se_pooled <0 & meanCoef + 1.96*se_pooled>0, 0, 1)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",                                               "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::filter(type=="Un_weighted") %>%
  dplyr::mutate(reps=i)
  
# Summarize environmental effects #
results_h4<-lmRes_tot2 %>%
dplyr::bind_rows(lmRes_tot3) %>%
dplyr::filter(rep<=i) %>%
dplyr::filter(!predictors %in% c("(Intercept)")) %>%
ungroup() %>%
dplyr::mutate(upper=Estimate + 1.96*Std..Error, lower=Estimate - 1.96*Std..Error) %>%
dplyr::group_by(species, type, test, predictors) %>%
dplyr::summarise(repsTot=n_distinct(rep), mean.fit=mean(Estimate), sd.fit=sd(Estimate), se.fit=sd.fit/sqrt(reps), within_model_var=mean(Std..Error^2), between_model_var=var(Estimate), totalVar=within_model_var + (1+1/repsTot)*between_model_var,
se_pooled=sqrt(totalVar), mean.se=mean(Std..Error), max.se=max(Std..Error), min.se=min(Std..Error), meanCI=1.96*mean.se, mean.r2=mean(r2), sd.r2=sd(r2), se.r2=1.96*sd.r2/sqrt(reps), mean.p=mean(Pr...t..), sdp=sd(Pr...t..), sep=sdp/sqrt(reps), minlower=min(lower), maxupper=max(upper), meanlower=mean(lower), meanupper=mean(upper),
CI_low=mean.fit-1.96*se_pooled, CI_high=mean.fit+1.96*se_pooled, r2low=mean.r2-1.96*se.r2, r2high=mean.r2+1.96*se.r2) %>% 
dplyr::mutate(sig=ifelse(mean.fit-1.96*se_pooled< 0 & mean.fit + 1.96*se_pooled > 0, 0, 1)) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(test=ifelse(test=="WEE_cov_nb_vs_predictors", "WEE_COV_NB", "TEE_NB")) %>%
dplyr::mutate(predictors=ifelse(predictors=="Migr_scale", "MigratoryDist", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="sst_gain_scale", "SST_Gain", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="sst_start_scale", "SST_Start", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="Migr_scale:sst_gain_scale", "MigratoryDist:SST_Gain", predictors)) %>%
dplyr::mutate(predictor=factor(predictors, levels=c("MigratoryDist", "SST_Start", "SST_Gain", "MigratoryDist:SST_Gain"))) %>%
dplyr::mutate(test=factor(test, levels=c("WEE_COV_NB", "TEE_NB"))) %>%
dplyr::filter(type=="Un-weighted") %>%
dplyr::mutate(reps=i)

# Save all results #
results_h1_h2_all<-rbind(results_h1_h2_all, results_h1_h2)
results_h3_all<-rbind(results_h3_all, results_h3)
results_h4_all<-rbind(results_h4_all, results_h4)

}

### Step s4_4_4: plot results ###

# Factorize species in the same order #
results_h1_h2_all<-results_h1_h2_all %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 
												 
results_h3_all<-results_h3_all %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 


# Summary plots #

plot1<-ggplot() +
geom_line(data=results_h1_h2_all, aes(x=reps, y=eff_WEE_cov_nb_mean, color=species)) +
geom_ribbon(data=results_h1_h2_all, aes(x=reps, y=eff_WEE_cov_nb_mean, ymin=eff_WEE_cov_nb_mean - 1.96*eff_WEE_cov_nb_se, ymax=eff_WEE_cov_nb_mean + 1.96*eff_WEE_cov_nb_se, fill=species), alpha=0.2) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab("KW #1: Effect size (mean +/- 95% CI)") +
xlab("Number of iterations") +
theme_bw() +
labs(tag="A)")

plot2<-ggplot() +
geom_line(data=results_h1_h2_all, aes(x=reps, y=eff_TEE_nb_mean, color=species)) +
geom_ribbon(data=results_h1_h2_all, aes(x=reps, y=eff_TEE_nb_mean, ymin=eff_TEE_nb_mean - 1.96*eff_TEE_nb_se, ymax=eff_TEE_nb_mean + 1.96*eff_TEE_nb_se, fill=species), alpha=0.2) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab("KW #2: Effect size (mean +/- 95% CI)") +
xlab("Number of iterations") +
theme_bw()+
labs(tag="B)")

plot3<-ggplot() +
geom_line(data=results_h3_all, aes(x=reps, y=meanCoef  , color=species)) +
geom_ribbon(data=results_h3_all, aes(x=reps, y=meanCoef  , ymin=meanCoef - 1.96*se_pooled , ymax=meanCoef + 1.96*se_pooled , fill=species), alpha=0.2) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab(expression(WEE[COV_NB] ~ "vs." ~ TEE[NB] ~ ": Effect size (mean +/- 95% CI)")) +
xlab("Number of iterations") +
theme_bw() +
labs(tag="C)")

plot4<-ggplot() +
geom_line(data=filter(results_h4_all, test=="WEE_COV_NB"), aes(x=reps, y=mean.fit  , color=species)) +
geom_ribbon(data=filter(results_h4_all, test=="WEE_COV_NB"), aes(x=reps, y=mean.fit , ymin=mean.fit - 1.96*mean.se , ymax=mean.fit + 1.96*mean.se , fill=species, group=interaction(species, predictors)), alpha=0.2) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab(expression(WEE[COV_NB] ~ "vs. other predictors: Effect size (mean +/- 95% CI)")) +
xlab("Number of iterations") +
theme_bw() +
facet_wrap(~predictors)+
labs(tag="D)")

plot5<-ggplot() +
geom_line(data=filter(results_h4_all, test=="TEE_NB"), aes(x=reps, y=mean.fit  , color=species)) +
geom_ribbon(data=filter(results_h4_all, test=="TEE_NB"), aes(x=reps, y=mean.fit , ymin=mean.fit - 1.96*mean.se , ymax=mean.fit + 1.96*mean.se , fill=species, group=interaction(species, predictors)), alpha=0.2) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
ylab(expression(TEE[NB] ~ "vs. other predictors: Effect size (mean +/- 95% CI)")) +
xlab("Number of iterations") +
theme_bw() +
facet_wrap(~predictors)+
labs(tag="E)")

pdf("./results/figures/supplementary/FigureS13.pdf")
grid.arrange(plot1, plot2)
grid.arrange(plot1, plot3)
plot(plot4)
plot(plot5)
dev.off()