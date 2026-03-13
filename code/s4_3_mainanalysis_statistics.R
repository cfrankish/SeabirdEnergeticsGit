# This script conducts the main statistics in my manuscript #
# Input files are start & end of the study period #
# There are many output files:
# output_files1 = f"./results/tables/main/table9_totalNBCosts.csv" -> has colony & species-level patterns in WEEnb & TEEnb
# output_files2 = f"./results/tables/main/table10_stats3_dredge_nbCosts.csv"
# output_files3 = f"./results/tables/main/table11_cluster_pop.csv"
# output_files4 = f"./results/tables/main/table12_stats4_dredge_nbCosts_species.csv"
# output_files5 = f"./results/tables/main/table13_stats5_dredge_nbCosts_deviance.csv" -> # result of statistics (link between deviance, nonbreeding costs & possible explanatory variables)
# output_files6 = f"./results/tables/main/table14_stats6_dredge_nbCosts_cov.csv" -> result of statistics (link between cov, nonbreeding costs)

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
library(MuMIn)
library(pracma)
library(lmerTest)
library(cluster)
library(ggrepel)
library(tibble)
library(rstatix)

args <- commandArgs(trailingOnly = TRUE) # This allows R to read in arguments written in the workflow file

### Step 00: define start & end of study period ###

startDate<-args[1] # Read-in start of study period
endDate<-args[2] # Read-in end date of study period

### Step 0: set-up sample size & iteration number ###

print("Step 0: setting up initial parameters")

# Set up minimum sample size & number of iterations
minSampleSize<-5
print(paste0("min sample size per colony is: ", minSampleSize))
reps<-5
print(paste0("min iteration number is: ", reps))

### Step 1: Estimate total non-breeding costs ###

print("Step 1: Estimate total NB costs ..:")

## First we make a dataset listing the study period dates so we can subset our individual energy expenditure datasets ##

# First we create a sequence of dates going from start to end date # 
dates<-data.frame(dateKeep=seq(as.Date(startDate), as.Date(endDate), 1))
dates$doy<-1:nrow(dates)
dates$month<-as.numeric(substr(dates$date, 6, 7))
dates$day<-as.numeric(substr(dates$date, 9, 10))

# Add week number for summarizing information #
dates_weekly<-dates %>%
  dplyr::mutate(weekNo=ceiling(doy/7)) %>%
  dplyr::group_by(weekNo) %>%
  dplyr::mutate(days=n_distinct(dateKeep)) %>%
  dplyr::filter(days==7) %>%
  dplyr::select(-days)

## Now we open up a datset with WEE_cov_nb and migratory distance per bird & rep that we will conduct stats on ##

# Determine list of unique ids
allResults_deviation<-list.files("./results/tables/main/", full.names=TRUE)
deviation_base<-allResults_deviation[grepl("/energy_metrics_", allResults_deviation)]

# Create list to save results in #
devianceMean<-list()

# Loop through species-specific datasets

for (m in 1:length(deviation_base)) {

deviance_sub<-readRDS(deviation_base[m])

devianceMean<-rbind(devianceMean, deviance_sub)

}

# Save unique ids
ids<-unique(devianceMean$individ_id)

# Define number of species so we can loop through these 
speciesList<-unique(devianceMean$species)

# divide DEE by weight to get in kJ.g
species<-data.frame(species=c("Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot", "Little auk", "Atlantic puffin"))
species$allometryCoef<-c(0.765, 0.717, 0.689, 0.689, 0.689, 0.689)
species$LCT<-c(9, 12.5, 14.18, 14.18, 14.18, 14.18)

# Define location of other information (sst?)
deviation_weekly<-allResults_deviation[grepl("/weeklydeviance_", allResults_deviation)]

# Set options for dredging
options(na.action = "na.fail") 

# open colony information & join so we have longitude & latitude & project into correct format
colony.summary<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")
colonyCoords<-colony.summary %>%
dplyr::ungroup() %>%
dplyr::group_by(colony) %>%
dplyr::slice(1) %>%
dplyr::select(colony, col_lat)

# Make list to save results in
totalCosts<-list() # raw data
lmRes_tot<-list() # WEE_cov_nb vs TEE_nb
lmRes_tot2<-list() # WEE_cov_nb vs other behaviors 
lmRes_tot3<-list() # TEE_nb vs other predictors

for (i in 1:reps) {
  
  repSub<-i  
  
  # Loop through species
  
  totalCosts_reps<-list() # Raw data for plotting results agaisn't
  lmRes_reps<-list()
  lmRes_reps2<-list()
  lmRes_reps3<-list()
  
  for (j in 1:length(speciesList)) {
    
    #print(paste0("Calculating amplitude for rep ", i, " species ", j ))
    print(paste0("Rep ", i, " Species ", j, ": Estimating NB costs...")) 
    
    # Subset to species j
	speciesSub<-speciesList[j]
	devianceWeekly<-readRDS(deviation_weekly[grepl(speciesSub, deviation_weekly)])
    
    # Subset species data frame to rep i
    speciesDaily<-subset(devianceWeekly, rep==i)
    
    # Summarise SST (at NB locations and at colony) so I can join as predictor variables
	# I also look at proportion of time spent below LCT for investigating other questions #
	
    speciesNB_stats<-speciesDaily %>%
	dplyr::left_join(species, by=c("species")) %>%
      ungroup() %>%
	  dplyr::group_by(rep, species, colony, individ_id) %>%
	  dplyr::mutate(meanNBSST=mean(meanSST)) %>%
	  ungroup() %>%
	  dplyr::group_by(rep, species, colony, individ_id, weekNo) %>%
	  dplyr::mutate(underLCT=ifelse(meanSST<LCT, 1, 0), underLCT_col=ifelse(meanSST_colony < LCT, 1, 0), weeks=n_distinct(weekNo)) %>%
	  ungroup() %>%
      dplyr::group_by(rep, species, colony, individ_id) %>%
      dplyr::summarise(SST=mean(meanSST), SST_col=mean(meanSST_colony), SSTdiff=SST-SST_col, 
	  propWarm=1-sum(underLCT)/n_distinct(weekNo), propCold_col=sum(underLCT_col)/n_distinct(weekNo))   
	
    # scale variables for conducting stats # 
    speciesNB_stats5<-speciesNB_stats %>%
      dplyr::left_join(species, by=c("species")) %>%
	  dplyr::left_join(devianceMean, by=c("rep", "species", "colony", "individ_id")) %>%
      ungroup() %>%
      dplyr::mutate(TEE_nb_scale=scale(TEE_nb), sst_scale=scale(SST), Migr_scale=scale(MigratoryDistKm), WEE_cov_nb_scale=scale(WEE_cov_nb), 
	  sst_start_scale=scale(SST_col), sst_gain_scale=scale(SSTdiff))  
	  
    # Remove colonies with small sample sizes
    speciesNB_stats6<-speciesNB_stats5 %>%
      ungroup() %>%
      dplyr::group_by(rep, species, colony) %>%
      dplyr::mutate(birdsTot=n_distinct(individ_id)) %>%
      dplyr::filter(birdsTot >= minSampleSize) %>%
	  dplyr::mutate(colonyWeights=1/birdsTot) %>%
	  dplyr::left_join(colonyCoords, by=c("colony")) %>%
	  ungroup() %>%
	  droplevels()
	
    # Turn colony into a factor 	
	speciesNB_stats6$colony<-factor(speciesNB_stats6$colony)
	  
	#### Stats # 1: Does total TEE_nb differ between populations? (kruskal-wallis) ####
	normalityTest<-shapiro.test(speciesNB_stats6$TEE_nb)
	pval<-normalityTest$p.value
	speciesNB_stats6$nbcosts_normality<-pval

	# data is not normally distirbuted and we conduct a kruskal-wallis test #
	test1<-kruskal.test(TEE_nb ~ colony, data = speciesNB_stats6)
	pvaltest<-test1$p.value
	df<-test1$parameter
	H<-test1$statistic
	speciesNB_stats6$popTest_nbcosts_pval<-pvaltest # Would be good to conduct a post-hoc test later to determine which populations are different
	speciesNB_stats6$popTest_nbcosts_df<-df 
	speciesNB_stats6$popTest_nbcosts_H<-H 
	effSize<-rstatix::kruskal_effsize(data=speciesNB_stats6, formula=TEE_nb ~ colony) # Calculate effect size
	speciesNB_stats6$popTest_nbcosts_effectSize<-effSize$effsize # Now we extract effect size 
	
    #### Stats # 2: Does WEE_cov_nb differ between populations ####
	normalityTest2<-shapiro.test(speciesNB_stats6$WEE_cov_nb)
	pval2<-normalityTest2$p.value
	speciesNB_stats6$cov_normality<-pval2

	# data is not normally distributed and we conduct a kruskal-wallis test #
	test2<-kruskal.test(WEE_cov_nb ~ colony, data = speciesNB_stats6)
	pvaltest2<-test2$p.value
	df2<-test2$parameter
	H2<-test2$statistic
	speciesNB_stats6$popTest_cov_pval<-pvaltest2 # Would be good to conduct a post-hoc test later to determine which populations are different
    speciesNB_stats6$popTest_cov_df<-df2
    speciesNB_stats6$popTest_cov_H<-H2
	effSize2<-kruskal_effsize(data=speciesNB_stats6, formula=WEE_cov_nb ~ colony) # Calculate effect size
	speciesNB_stats6$popTest_cov_effectSize<-effSize2$effsize 

    #### Stats # 3: How are these two linked to each other & environmental variables? ####
	
    # Fit model #1 # (weighted)
	lm1<-lmer(TEE_nb ~ WEE_cov_nb + (1|colony), weights=colonyWeights, data=speciesNB_stats6) # do the weights change anything? 
	
	lm1_sum<-summary(lm1) # Extract summary
    lm1_sum_df<-data.frame(coefficients(lm1_sum)) # Turn into a data frame
    pred1<-data.frame(WEE_cov_nb=seq(min(devianceWeekly$covDEE), max(devianceWeekly$covDEE), 0.5)) # Make a data frame to predict fitted values from this model (we use the range of values accross all reps)
	pred1$colony<-speciesNB_stats6$colony[1] # Add colony 
    pred1$fit<-predict(lm1, newdata=pred1, re.form=NA) # Make predictions (fit)
	se<-predict(lm1, newdata=pred1, re.form=NA, se.fit=TRUE) # Calculate SE
	pred1$se<-se$se.fit # Add SE
    r2_values <- r.squaredGLMM(lm1) # Calculate R2
    pred1$r2<-r2_values[1] # Add r2
    pred1$rep<-speciesNB_stats6$rep[1] # add iteration #
    pred1$species<-speciesNB_stats6$species[1] # Add species
    pred1$predictor<-"WEE_cov_nb" # Add predictor name
    pred1$pvalue<-lm1_sum$Pr...t..[2] # Add p-value
    pred1<-pred1%>%
      rename(predictor_val=WEE_cov_nb)
	pred1$intercept<-coefficients(lm1_sum)[1,1]
	pred1$coefficient<-coefficients(lm1_sum)[2,1]
	pred1$type<-"Weighted"
	
	# Model 1: un-weighted #
	lm1_2<-lmer(TEE_nb ~ WEE_cov_nb + (1|colony), data=speciesNB_stats6) 
	
	lm1_2_sum<-summary(lm1_2) # Extract summary
    lm1_2_sum_df<-data.frame(coefficients(lm1_2_sum)) # Turn into a data frame
    pred1_2<-data.frame(WEE_cov_nb=seq(min(devianceWeekly$covDEE), max(devianceWeekly$covDEE), 0.5)) # Make a data frame to predict fitted values from this model (we use the range of values accross all reps)
	pred1_2$colony<-speciesNB_stats6$colony[1] # Add colony 
    pred1_2$fit<-predict(lm1_2, newdata=pred1_2, re.form=NA) # Make predictions (fit)
	se<-predict(lm1_2, newdata=pred1_2, re.form=NA, se.fit=TRUE) # Calculate SE
	pred1_2$se<-se$se.fit # Add SE
    r2_values <- r.squaredGLMM(lm1_2) # Calculate R2
    pred1_2$r2<-r2_values[1] # Add r2
    pred1_2$rep<-speciesNB_stats6$rep[1] # add iteration #
    pred1_2$species<-speciesNB_stats6$species[1] # Add species
    pred1_2$predictor<-"WEE_cov_nb" # Add predictor name
    pred1_2$pvalue<-lm1_2_sum$Pr...t..[2] # Add p-value
    pred1_2<-pred1_2%>%
      rename(predictor_val=WEE_cov_nb)
	pred1_2$intercept<-coefficients(lm1_2_sum)[1,1]
	pred1_2$coefficient<-coefficients(lm1_2_sum)[2,1]
	pred1_2$type<-"Un_weighted"
	
	# WEE_cov_nb vs- environmental variables #
	
	# Weighted #
    lm2<-lmer(WEE_cov_nb_scale ~ Migr_scale*sst_gain_scale + sst_start_scale +  (1 |colony), weights=colonyWeights, data=speciesNB_stats6)
	
	# Extract coefficients
    lm2_sum<-summary(lm2)
    lm2_sum_df<-data.frame(coefficients(lm2_sum))
	sum1<-lm2_sum_df
	sum1$predictors<-row.names(sum1)
    r2_values <- r.squaredGLMM(lm2)
    sum1$r2<-r2_values[1]
    sum1$rep<-speciesNB_stats6$rep[1]
    sum1$species<-speciesNB_stats6$species[1]
    sum1$test<-"WEE_cov_nb_vs_predictors"
	sum1$type<-"Weighted"
	
	# Un-weighted #
    lm2_2<-lmer(WEE_cov_nb_scale ~ Migr_scale*sst_gain_scale + sst_start_scale +  (1 |colony), data=speciesNB_stats6)
	
	# Extract coefficients
    lm2_2_sum<-summary(lm2_2)
    lm2_2_sum_df<-data.frame(coefficients(lm2_2_sum))
	sum1_2<-lm2_2_sum_df
	sum1_2$predictors<-row.names(sum1_2)
    r2_values <- r.squaredGLMM(lm2_2)
    sum1_2$r2<-r2_values[1]
    sum1_2$rep<-speciesNB_stats6$rep[1]
    sum1_2$species<-speciesNB_stats6$species[1]
    sum1_2$test<-"WEE_cov_nb_vs_predictors"
	sum1_2$type<-"Un-weighted"
    
    # TEE_nb vs- environmental variables #
	
	# Weighted #
    lm3<-lmer(TEE_nb_scale ~ Migr_scale*sst_gain_scale + sst_start_scale +  (1 |colony), weights=colonyWeights, data=speciesNB_stats6)
	
	# Extract coefficients
	lm3_sum<-summary(lm3)
    lm3_sum_df<-data.frame(coefficients(lm3_sum))
	sum3<-lm3_sum_df
	sum3$predictors<-row.names(sum3)
    r2_values <- r.squaredGLMM(lm3)
    sum3$r2<-r2_values[1]
    sum3$rep<-speciesNB_stats6$rep[1]
    sum3$species<-speciesNB_stats6$species[1]
    sum3$test<-"TEE_nb_vs_predictors"
	sum3$type<-"Weighted"
	
	# Un-weighted #
    lm3_2<-lmer(TEE_nb_scale ~ Migr_scale*sst_gain_scale + sst_start_scale +  (1 |colony), data=speciesNB_stats6)
	
	# Extract coefficients
	lm3_2_sum<-summary(lm3_2)
    lm3_2_sum_df<-data.frame(coefficients(lm3_2_sum))
	sum3_2<-lm3_2_sum_df
	sum3_2$predictors<-row.names(sum3_2)
    r2_values <- r.squaredGLMM(lm3_2)
    sum3_2$r2<-r2_values[1]
    sum3_2$rep<-speciesNB_stats6$rep[1]
    sum3_2$species<-speciesNB_stats6$species[1]
    sum3_2$test<-"TEE_nb_vs_predictors"
	sum3_2$type<-"Un-weighted"
    
    # Collate results
    lmRes<-rbind(pred1, pred1_2) # TEE_nb vs WEE_cov_nb
	lmRes2<-rbind(sum1, sum1_2) # WEE_cov_nb vs other predictions
	lmRes3<-rbind(sum3, sum3_2) # TEE_nb vs other predictors
	
    # Save results of models
    lmRes_reps<-rbind(lmRes_reps, lmRes) # TEE_nb vs WEE_cov_nb
	lmRes_reps2<-rbind(lmRes_reps2, lmRes2) # WEE_cov_nb vs other predictions
	lmRes_reps3<-rbind(lmRes_reps3, lmRes3) # TEE_nb vs other predictors
    
    # Save results : raw data for plotting 
	speciesNB_stats6$rep<-i
    totalCosts_reps<-rbind(totalCosts_reps, speciesNB_stats6)
    
  }
  
  # Save all results
  totalCosts<-rbind(totalCosts, totalCosts_reps) # individual-level scaled data for plotting
  lmRes_tot<-rbind(lmRes_tot, lmRes_reps) # species seperatey, # TEE_nb vs WEE_cov_nb
  lmRes_tot2<-rbind(lmRes_tot2, lmRes_reps2) # species seperatey, WEE_cov_nb vs other predictions
  lmRes_tot3<-rbind(lmRes_tot3, lmRes_reps3) # species seperatey, TEE_nb vs other predictions
  
} 

print("Preparing datasets for plotting...")

#### Main results ##### (Figure 2A, B & C)

# Summarize effect sizes from tests H1 & H2:

results_h1_h2<-totalCosts %>%
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
dplyr::group_by(species) %>%
dplyr::mutate(sig_h1=ifelse(eff_WEE_cov_nb_low>0 & eff_WEE_cov_nb_high >0, 1, 0)) %>%
dplyr::mutate(sig_h2=ifelse(eff_TEE_nb_low >0 & eff_TEE_nb_high > 0, 1, 0)) 

# Make species mean

speciesRes<-totalCosts %>%
   ungroup() %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::mutate(birds=n_distinct(individ_id)) %>%
  dplyr::filter(birds>=minSampleSize)  %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, colony, individ_id) %>%
  dplyr::summarise(TEE_nb_id=mean(TEE_nb), WEE_cov_nb_id=mean(WEE_cov_nb)) %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>% 
  dplyr::summarize(TEE_nb_pop=mean(TEE_nb_id), WEE_cov_nb_pop=mean(WEE_cov_nb_id)) %>%
  ungroup() %>%
  dplyr::group_by(species) %>% 
  dplyr::summarize(pops=n_distinct(colony), TEE_nb_sp_mean=mean(TEE_nb_pop), TEE_nb_sp_sd=sd(TEE_nb_pop), TEE_nb_sp_se=TEE_nb_sp_sd/sqrt(pops), 
  WEE_cov_nb_sp_mean=mean(WEE_cov_nb_pop), WEE_cov_nb_sp_sd=sd(WEE_cov_nb_pop), WEE_cov_nb_sp_se=WEE_cov_nb_sp_sd/sqrt(pops)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::left_join(results_h1_h2, by=c("species"))
  
# Make colony mean 

colonyRes<-totalCosts %>%
   ungroup() %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::mutate(birds=n_distinct(individ_id)) %>%
  dplyr::filter(birds>=minSampleSize)  %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, colony, individ_id) %>%
  dplyr::summarise(TEE_nb_id=mean(TEE_nb), WEE_cov_nb_id=mean(WEE_cov_nb), sst_mean_id=mean(SST), sst_col_mean_id=mean(SST_col), sst_diff_mean_id=mean(SSTdiff), migratory_dist_mean_id=mean(MigratoryDistKm)) %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>% 
  dplyr::summarize(birds=n_distinct(individ_id), TEE_nb_pop_mean=mean(TEE_nb_id), TEE_nb_pop_sd=sd(TEE_nb_id), TEE_nb_pop_se=TEE_nb_pop_sd/sqrt(birds), 
  WEE_cov_nb_pop_mean=mean(WEE_cov_nb_id), WEE_cov_nb_pop_sd=sd(WEE_cov_nb_id), WEE_cov_nb_pop_se=WEE_cov_nb_pop_sd/sqrt(birds), 
  sst_mean_pop=mean(sst_mean_id),  sst_col_mean_pop=mean(sst_col_mean_id), sst_diff_mean_pop=mean(sst_diff_mean_id), migratory_dist_mean_pop=mean(migratory_dist_mean_id)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  ungroup() %>%
  dplyr::group_by(species) %>%
  dplyr::arrange(TEE_nb_pop_mean)%>% # Order by TEEnb
  dplyr::mutate(colonySp=paste0(colony, "_", species))
  
 colonyRes2<-totalCosts %>%
   ungroup() %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::mutate(birds=n_distinct(individ_id)) %>%
  dplyr::filter(birds>=minSampleSize)  %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, colony, individ_id) %>%
  dplyr::summarise(TEE_nb_id=mean(TEE_nb), WEE_cov_nb_id=mean(WEE_cov_nb)) %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>% 
  dplyr::summarize(birds=n_distinct(individ_id), TEE_nb_pop_mean=mean(TEE_nb_id), TEE_nb_pop_sd=sd(TEE_nb_id), TEE_nb_pop_se=TEE_nb_pop_sd/sqrt(birds), 
  WEE_cov_nb_pop_mean=mean(WEE_cov_nb_id), WEE_cov_nb_pop_sd=sd(WEE_cov_nb_id), WEE_cov_nb_pop_se=WEE_cov_nb_pop_sd/sqrt(birds)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  ungroup() %>%
  dplyr::group_by(species) %>%
  dplyr::arrange(WEE_cov_nb_pop_mean)%>% # Order by WEE-cov-nb
  dplyr::mutate(colonySp=paste0(colony, "_", species))								 

# Order colonyRes populations according to value
colonyOrder<-unique(colonyRes$colonySp)
colonyRes$colonySp<-factor(colonyRes$colonySp, levels=colonyOrder)	

# Add stars to denote how significant #
speciesRes$stars<-ifelse(speciesRes$sig_h2==1, "*", "")

print("Saving Figures 2A & B")
		
Figure2B<-ggplot() +
   geom_pointrange(data=colonyRes, aes(x=TEE_nb_pop_mean, y=species, xmin=TEE_nb_pop_mean - 1.96*TEE_nb_pop_se, xmax=TEE_nb_pop_mean + 1.96*TEE_nb_pop_se, color=species), position = position_dodge2(width=0.6), cex=0.2, alpha=0.2) +
  geom_pointrange(data=speciesRes, aes(x=TEE_nb_sp_mean, y=species, xmin=TEE_nb_sp_mean - 1.96*TEE_nb_sp_se, xmax=TEE_nb_sp_mean + 1.96*TEE_nb_sp_se, color=species), position = position_dodge2(width=0.6), cex=0.4) +
  scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+
  geom_text(data=filter(speciesRes, sig_h2==1), aes(x=5900, y=species, label=stars)) +
  theme_bw() +
  xlab(expression("TEE"[NB]*" (kJ.g"^-1*")")) +
  ylab("") +
  theme(legend.position = "none") +
  ggtitle("B)") + 
  theme(
       panel.grid = element_blank(),
		axis.text=element_text(size=12),
		axis.title=element_text(size=12)) 

pdf("./results/figures/main/Figure2B.pdf", width=7, height=5)
plot(Figure2B)
dev.off()

# Coefficient of variation version

# Add stars to denote how significant #
speciesRes$stars3<-ifelse(speciesRes$sig_h1==1, "*", "")
		
Figure2A<-ggplot() +
 geom_pointrange(data=colonyRes2, aes(x=WEE_cov_nb_pop_mean, y=species, xmin=WEE_cov_nb_pop_mean - 1.96*WEE_cov_nb_pop_se, xmax=WEE_cov_nb_pop_mean + 1.96*WEE_cov_nb_pop_se, color=species), position = position_dodge2(width=0.6), cex=0.2, alpha=0.2) +
 geom_pointrange(data=speciesRes, aes(x=WEE_cov_nb_sp_mean, y=species, xmin=WEE_cov_nb_sp_mean - 1.96*WEE_cov_nb_sp_se, xmax=WEE_cov_nb_sp_mean + 1.96*WEE_cov_nb_sp_se, color=species), position = position_dodge2(width=0.6), cex=0.4) +
 scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+
 geom_text(data=filter(speciesRes, sig_h1==1), aes(x=34, y=species, label=stars3)) +
 theme_bw() +
 xlab(expression("WEE"[COV_NB]*"")) +
 ylab("") +
 theme(legend.position = "none") +
 ggtitle("A)") + 
  theme(
		panel.grid = element_blank(),
		axis.text=element_text(size=12),
		axis.title=element_text(size=12),
		ylim(0, 35))

pdf("./results/figures/main/Figure2A.pdf", width=7, height=5)
plot(Figure2A)
dev.off()

### Figure 2C... ###

# Summarize results of WEE_cov_nb vs. TEE_nb

lmRes_sum<-lmRes_tot %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::group_by(species, predictor, predictor_val) %>%
  dplyr::summarise(meanEstimate=mean(fit), meanSE=mean(se), reps=n_distinct(rep), 
                   meanr2=mean(r2), sdr2=sd(r2), ser2=sdr2/sqrt(reps), meanCoef=mean(coefficient), sdCoef=sd(coefficient), seCoef=sdCoef/sqrt(reps)) %>%
  dplyr::mutate(sig=ifelse(meanCoef - 1.96*seCoef <0 & meanCoef + 1.96*seCoef>0, 0, 1)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 
  
# Create some raw data to plot below these predictions #

totalCosts_sum<-totalCosts %>%
dplyr::ungroup() %>%
dplyr::group_by(species, individ_id) %>%
dplyr::summarise(TEE_nb_mean=mean(TEE_nb), TEE_nb_sd=sd(TEE_nb), WEE_cov_nb_mean=mean(WEE_cov_nb), WEE_cov_nb_sd=sd(WEE_cov_nb)) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 

totalCosts_colony_sum<-totalCosts %>%
dplyr::ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(TEE_nb_mean=mean(TEE_nb), WEE_cov_nb_mean=mean(WEE_cov_nb)) %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(TEE_nb_mean2=mean(TEE_nb_mean), sd1=sd(TEE_nb_mean), WEE_cov_nb_mean2=mean(WEE_cov_nb_mean), sd2=sd(WEE_cov_nb_mean)) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 

# Shorten estimate lines to range of data

minTEE<-totalCosts_sum %>%
dplyr::group_by(species) %>%
arrange(TEE_nb_mean) %>%
dplyr::slice(1) %>%
dplyr::select(species, TEE_nb_mean) %>%
rename(minTEE=TEE_nb_mean)

maxTEE<-totalCosts_sum %>%
dplyr::group_by(species) %>%
arrange(TEE_nb_mean) %>%
dplyr::slice(n()) %>%
dplyr::select(species, TEE_nb_mean) %>%
rename(maxTEE=TEE_nb_mean)

minWEE<-totalCosts_sum %>%
dplyr::group_by(species) %>%
arrange(WEE_cov_nb_mean) %>%
dplyr::slice(1) %>%
dplyr::select(species, WEE_cov_nb_mean) %>%
rename(minWEE=WEE_cov_nb_mean)

maxWEE<-totalCosts_sum %>%
dplyr::group_by(species) %>%
arrange(WEE_cov_nb_mean) %>%
dplyr::slice(n()) %>%
dplyr::select(species, WEE_cov_nb_mean) %>%
rename(maxWEE=WEE_cov_nb_mean)

lmRes_short<-lmRes_sum %>%
dplyr::left_join(minTEE, by=c("species")) %>%
dplyr::left_join(maxTEE, by=c("species")) %>%
dplyr::left_join(minWEE, by=c("species")) %>%
dplyr::left_join(maxWEE, by=c("species")) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::filter(meanEstimate >= minTEE & meanEstimate <= maxTEE) %>%
dplyr::filter(predictor_val >= minWEE & predictor_val <= maxWEE) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 

# plot these

Figure2C<-ggplot() +
geom_point(data=totalCosts_sum, aes(x=WEE_cov_nb_mean, y=TEE_nb_mean, color=species), alpha=0.1) +
geom_line(data=lmRes_short, aes(y=meanEstimate, x=predictor_val, color=species, group=interaction(species, predictor))) +
geom_ribbon(data=lmRes_short, aes(x=predictor_val, y=meanEstimate, ymin=meanEstimate-meanSE*1.96, ymax=meanEstimate + meanSE*1.96, group=species, fill=species), alpha=0.2) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
theme_bw() +
xlab(expression("WEE"[COV_NB]*"")) +
ylab(expression("TEE"[NB]*" (kJ.g"^-1*")")) +
theme_bw() +
ggtitle("C)") + 
theme(legend.position = "none",
panel.grid = element_blank(),
axis.text=element_text(size=12),
		axis.title=element_text(size=12)) +
facet_wrap(~species)  
  
pdf("./results/figures/main/Figure2C.pdf", width=6, height=5)
plot(Figure2C)
dev.off()

### Figure 3A ####

Figure3A<-lmRes_tot2 %>%
dplyr::bind_rows(lmRes_tot3) %>%
dplyr::filter(!predictors %in% c("(Intercept)")) %>%
ungroup() %>%
dplyr::mutate(upper=Estimate + 1.96*Std..Error, lower=Estimate - 1.96*Std..Error) %>%
dplyr::group_by(species, test, predictors) %>%
dplyr::summarise(repsTot=n_distinct(rep), mean.fit=mean(Estimate), sd.fit=sd(Estimate), se.fit=sd.fit/sqrt(reps), mean.se=mean(Std..Error), max.se=max(Std..Error), min.se=min(Std..Error), meanCI=1.96*mean.se, mean.r2=mean(r2), sd.r2=sd(r2), se.r2=1.96*sd.r2/sqrt(reps), mean.p=mean(Pr...t..), sdp=sd(Pr...t..), sep=sdp/sqrt(reps), minlower=min(lower), maxupper=max(upper), meanlower=mean(lower), meanupper=mean(upper)) %>% 
dplyr::mutate(sig=ifelse(mean.fit-1.96*mean.se < 0 & mean.fit + 1.96*mean.se > 0, 0, 1)) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(test=ifelse(test=="WEE_cov_nb_vs_predictors", "WEE_COV_NB", "TEE_NB")) %>%
dplyr::mutate(predictors=ifelse(predictors=="Migr_scale", "MigratoryDist", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="sst_gain_scale", "SST_Gain", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="sst_start_scale", "SST_Start", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="Migr_scale:sst_gain_scale", "MigratoryDist:SST_Gain", predictors)) %>%
dplyr::mutate(predictor=factor(predictors, levels=c("MigratoryDist", "SST_Start", "SST_Gain", "MigratoryDist:SST_Gain"))) %>%
dplyr::mutate(test=factor(test, levels=c("WEE_COV_NB", "TEE_NB"))) %>%
ggplot(aes(x=mean.fit, y=predictors)) +
geom_pointrange(aes(y=predictors, xmin=mean.fit - 1.96*mean.se, xmax=mean.fit + 1.96*mean.se, group=interaction(species, test), color=species, alpha=factor(sig)), position=position_dodge(width=0.9)) +
geom_segment(aes(x=minlower, y=predictors, xend=maxupper, group=interaction(species, test), color=species, alpha=factor(sig)), position=position_dodge2(width=0.9), linetype="dotted") +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+ 
theme_bw() +
theme(legend.position="bottom", axis.title.y = element_text(size = 16)) +
scale_shape_manual(values=c(1, 16)) +
geom_vline(xintercept=0, linetype="dashed") +
ylab("Predictor variables") +
xlab("Estimate +/- 95% confidence intervals") +
facet_wrap(~test, nrow=2) +
guides(color="none") +
labs(tag="A)", alpha="Sig", linetype="Start temp")

pdf("./results/figures/main/Figure3A.pdf", width=5, height=7.5)
grid.arrange(Figure3A)
dev.off()

#### Preparing data for supplementary plot making ####

# Make a supplementary figure showing Min & max values per population #

# open colony information & join so we have longitude & latitude & project into correct format
colony.summary<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")

speciesMatch<-data.frame(speciesLatin=c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle"), 
                         species=c("Brünnich's guillemot", "Black-legged kittiwake", "Common guillemot", "Atlantic puffin", "Northern fulmar", "Little auk"))

colonyMatch<-colony.summary %>%
 dplyr::filter(species %in% c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle")) %>%
  dplyr::group_by(colony) %>%
  dplyr::slice(1) %>%
  dplyr::select(species, colony) %>%
  rename(speciesLatin=species) %>%
  dplyr::left_join(speciesMatch, by=c("speciesLatin"))

colonyNames<-colony.summary %>%
  dplyr::filter(species %in% c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle")) %>%
  ungroup() %>%
  dplyr::group_by(colony) %>%
  dplyr::slice(1) %>%
  dplyr::select(colony, col_lon, col_lat) %>%
  arrange(desc(col_lat)) %>%
  ungroup() %>%
  dplyr::mutate(country=c("RU", "NO", "NO", "NO", "GR", "RU", "NO", "GR", "GR", "RU", "NO", "CA", "GR", "CA", "NO", "NO", "GR",
                          "RU", "NO", "GR", "RU", "RU", "NO", "RU", "NO", "CA", "IC", "IC", "IC", "IC", "IC", "IC", 
                          "GR", "NO", "IC", "IC", "GR", "IC", "IC", "NO", "IC", "CA", "CA", "NO", "FA", "GR", "UK", "NO", "UK", "UK", "DK",
                          "UK", "UK", "UK", "IR", "IR", "CA", "IR", "IR", "IR", "UK", "CA", "CA")) %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(colonyNo=row_number()) %>%
  dplyr::mutate(colonyName=paste0(country, colonyNo)) # Here i make a special naming system with country first and colony second. It is ordered from North to South
  
 colonies_lox<-colonyMatch %>%
 dplyr::left_join(colonyNames, by=c("colony")) %>%
 ungroup() %>%
 arrange(colonyName) %>%
 dplyr::group_by(colonyName) %>%
 dplyr::slice(1) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
												 dplyr::select(species, colony, col_lon, col_lat)
												 
#### Figures S8 & S12 Pop-specific plots of variation in TEEnb & WEEnb ####

colonies_lox2<-colonies_lox %>%
dplyr::select(-species)

colonyRes_lox<-colonyRes %>%
dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ungroup() %>%
dplyr::arrange(col_lat) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

colonyOrder<-unique(colonyRes_lox$colony)
colonyOrder2<-unique(colonyRes_lox$colonyName)

colonyRes_lox$colony<-factor(colonyRes_lox$colony, levels=colonyOrder)
colonyRes_lox$colonyName<-factor(colonyRes_lox$colonyName, levels=colonyOrder2)

# Add DEE estimated using information from shaffer et al for comparison
weightsdf<-data.frame(species=c("Black-legged kittiwake", "Atlantic puffin","Common guillemot", "Northern fulmar", "Brünnich's guillemot" , "Little auk"), 
weight=c(392, 395, 940, 728, 980, 149) , allometry=c(0.717, 0.689, 0.689, 0.765, 0.689, 0.689))

# Here we estimate average DEE to get an idea of whether calculations are any good
colonyRes_lox2<-colonyRes_lox %>%
ungroup() %>%
dplyr::left_join(weightsdf, by=c("species")) %>%
dplyr::mutate(totalNBcosts=TEE_nb_pop_mean*weight^allometry, costs_DEE=totalNBcosts/n_distinct(dates_weekly$dateKeep), 
sdDEE_tot=TEE_nb_pop_sd*weight^allometry, sdDay=sdDEE_tot/n_distinct(dates_weekly$dateKeep)) %>%
dplyr::select(species, colonyName, costs_DEE, sdDay) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 
												 
FMR_litt<-data.frame(species=c("Black-legged kittiwake", "Atlantic puffin","Common guillemot", "Northern fulmar", "Brünnich's guillemot" , "Little auk"), 
costs_DEE=c(995, 874, 1789, 1444, 2036, 609.9), sdDay=c(290, 151, 265, 720.6, 552, 26), colonyName="A")

colonyRes_lox3<-rbind(colonyRes_lox2, FMR_litt) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 


print("Saving plots S8 & S12...")

FiguresS8<-ggplot() +
 geom_pointrange(data=filter(colonyRes_lox3, !colonyName %in% c("A")), aes(y=costs_DEE, x=colonyName, ymin=costs_DEE - sdDay, ymax=costs_DEE + sdDay, color=species)) +
 geom_pointrange(data=filter(colonyRes_lox3, colonyName %in% c("A")), aes(y=costs_DEE, x=colonyName, ymin=costs_DEE - sdDay, ymax=costs_DEE + sdDay), color="black") +
 facet_wrap(~species, scales="free_y", nrow=3) +
 scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+
  theme_bw() +
  ylab(expression("DEE kJ day"^{-1})) +
  xlab("") +
 # geom_text(data=colonyRes_lox2, aes(x=colonyName, y=meanDEEg + 1.96*seDEE +  50, label=round(costs_DEE)), size=1.7) +
  theme_bw() +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("./results/figures/supplementary/FigureS8.pdf", width=15, height=9)
plot(FiguresS8)
dev.off()

FigureS12<-ggplot() +
 geom_pointrange(data=colonyRes_lox, aes(y=WEE_cov_nb_pop_mean, x=colonyName, ymin=WEE_cov_nb_pop_mean - 1.96*WEE_cov_nb_pop_se, ymax=WEE_cov_nb_pop_mean + 1.96*WEE_cov_nb_pop_se, color=species)) +
 facet_wrap(~species, scales="free_y", nrow=3) +
 scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+
  theme_bw() +
  ylab(expression("WEE"[COV_NB]*"")) +
  xlab("") +
  theme_bw() +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("./results/figures/supplementary/FigureS12.pdf", width=15, height=9)
plot(FigureS12)
dev.off()
 
### Figures S22, 24, 26, 28, 32, 34: Plot spatial distribution of energy strategies  ###

# First we scale energy from 0-1 for all species & then we plot by colony I guess #

scale_01 <- function(x) {
  if (!is.numeric(x)) stop("Input must be numeric")
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Subset to just each colony
colonies<-colonies_lox %>%
ungroup() %>%
dplyr::group_by(colony) %>%
dplyr::slice(1) %>%
dplyr::select(colony, col_lon, col_lat)

# Join with scaled dataset
colonyRes_scale<-colonyRes %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::mutate(TEE_nb_scale=scale_01(TEE_nb_pop_mean), WEE_cov_nb_scaled=scale_01(WEE_cov_nb_pop_mean)) %>%
dplyr::left_join(colonies, by=c( "colony")) 

# Open country layers
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Project coordinates 
coordinates(colonyRes_scale)<-~col_lon + col_lat
proj4string(colonyRes_scale)<-projection_84
colonies_lox_trans<-data.frame(spTransform(colonyRes_scale, projection_NA))
colonies_lox_trans<-colonies_lox_trans %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))


# Try and summarize overall #
coloniesSum<-colonies_lox_trans %>%
ungroup() %>%
dplyr::group_by(colony, species) %>%
dplyr::group_by(species) %>%
dplyr::group_by(colony) %>%
dplyr::summarise(coords.x1=mean(coords.x1), coords.x2=mean(coords.x2))
							
# In the supplementary - > we run this for every species individually #

print("Saving maps of energy strategies per species...")

speciesList<-unique(colonies_lox_trans$species)

for (i in 1:length(speciesList)) {

# Subset to species i
energyplot<-subset(colonies_lox_trans, species==speciesList[i]) 

# Map
Figure4a_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=TEE_nb_pop_mean), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(expression("TEE"[NB]*" (kJ.g"^-1*")"), colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), breaks = c(min(energyplot$TEE_nb_pop_mean), max(energyplot$TEE_nb_pop_mean))) +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  xlab("") +
  ylab("") +
  theme(legend.position="bottom") +
  labs(tag="A)")
  
  Figure4b_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=WEE_cov_nb_pop_mean), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(name=expression("WEE"[COV_NB]*""), colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  xlab("") +
  ylab("") +
  theme(legend.position="bottom") +
  labs(tag="B)")
  
  Figure4c_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=sst_col_mean_pop), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(name=expression("SST"[Start]*" (°C)"), colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  xlab("") +
  ylab("") +
  theme(legend.position="bottom") +
  labs(tag="C)")
  
  Figure4d_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=sst_diff_mean_pop), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(name=expression("SST"[Gain]*" (°C)"), colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  #labs(shape=expression("WEE"[devianceNB]*"")) +
  theme(legend.position="bottom") +
  labs(tag="D)")
  
  Figure4e_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=migratory_dist_mean_pop), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(name="MigratoryDist (km)", colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  #labs(shape=expression("WEE"[devianceNB]*"")) +
  theme(legend.position="bottom") +
  labs(tag="E)")
  
  
pdf(paste0("./results/figures/supplementary/energyStrategyPlots_", speciesList[i], ".pdf"), width=9)
grid.arrange(Figure4a_sup, Figure4b_sup, Figure4c_sup, Figure4d_sup, Figure4e_sup, nrow=2)
dev.off()

}


# Figure S13: Distribution of raw data -> WEE vs TEE #

FigureS13<-ggplot() +
 geom_pointrange(data=colonyRes_lox, aes(x=WEE_cov_nb_pop_mean, y=TEE_nb_pop_mean, xmin=WEE_cov_nb_pop_mean - 1.96*WEE_cov_nb_pop_se, xmax=WEE_cov_nb_pop_mean + 1.96*WEE_cov_nb_pop_se, color=species), alpha=0.1) +
 geom_pointrange(data=colonyRes_lox, aes(x=WEE_cov_nb_pop_mean, y=TEE_nb_pop_mean, ymin=TEE_nb_pop_mean - 1.96*TEE_nb_pop_se, ymax=TEE_nb_pop_mean + 1.96*TEE_nb_pop_se, color=species), alpha=0.1) +
 geom_text_repel(data=colonyRes_lox, aes(x=WEE_cov_nb_pop_mean, y=TEE_nb_pop_mean, label=colonyName), max.overlaps=30) + 
  geom_line(data=lmRes_short, aes(y=meanEstimate, x=predictor_val, color=species, group=interaction(species, predictor))) +
geom_ribbon(data=lmRes_short, aes(x=predictor_val, y=meanEstimate, ymin=meanEstimate-meanSE*1.96, ymax=meanEstimate + meanSE*1.96, group=species, fill=species), alpha=0.2) +
 scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
 scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
  theme_bw() +
  xlab(expression("WEE"[COV_NB]*"")) +
  ylab(expression("TEE"[NB]*" (kJ.g"^-1*")")) +
  theme_bw() +
  #ggtitle("C) Energetic strategies") + 
  theme(legend.position = "none",
  axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
	axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16))  +
  facet_wrap(~species, scales="free")
  
pdf("./results/figures/supplementary/FigureS13.pdf", width=18, height=10)
plot(FigureS13)
dev.off()

### Figure S14: Make plots showing distribution of data for supplementary (and to understand what data looks like) ###

DistributionData<-totalCosts %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(meanSST=mean(SST),
meanMigr=mean(MigratoryDistKm), meanSST_col=mean(SST_col), meanSST_diff=mean(SSTdiff))

dist1<-DistributionData %>%
dplyr::mutate(Metric=c("SST_NB"), value=meanSST)

dist2<-DistributionData %>%
dplyr::mutate(Metric=c("SST_COL"), value=meanSST_col)

dist3<-DistributionData %>%
dplyr::mutate(Metric=c("MIGRATORY_DIST"), value=meanMigr)

dist4<-DistributionData %>%
dplyr::mutate(Metric=c("SST_DIFF"), value=meanSST_diff)

alldist<-rbind(dist1, dist2, dist3, dist4) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) 
FigureS14a<-alldist %>%
dplyr::mutate(Metric=factor(Metric, levels=c("SST_COL", "SST_NB", "SST_DIFF", "MIGRATORY_DIST"))) %>%
dplyr::filter(Metric %in% c("SST_NB", "SST_COL")) %>%
ggplot() +
geom_density(aes(x=value, group=species, color=species)) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+ 
theme_bw() +
theme(legend.position="bottom") +
facet_wrap(~Metric, nrow=2) +
xlab("")

FigureS14b<-alldist %>%
dplyr::mutate(Metric=factor(Metric, levels=c("SST_COL", "SST_NB", "SST_DIFF", "MIGRATORY_DIST"))) %>%
dplyr::filter(Metric %in% c("SST_DIFF", "MIGRATORY_DIST")) %>%
ggplot() +
geom_density(aes(x=value, group=species, color=species)) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+ 
theme_bw() +
theme(legend.position="bottom") +
facet_wrap(~Metric, nrow=2, scales="free") +
xlab("")

pdf("./results/figures/supplementary/FigureS14a.pdf", height=10)
grid.arrange(FigureS14a)
dev.off()

pdf("./results/figures/supplementary/FigureS14b.pdf", height=10)
grid.arrange(FigureS14a)
dev.off()

### Figures S21, 23, 25, 27, 31 & 33: Plotting examples of how different pops behave (high & low deviance) ###

print("Plotting extreme cases of seabirds...pop-level by deviance")

# Show case some extreme birds: high energy cost, low deviance
DistributionData2<-totalCosts %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(meanSST=mean(SST), meanSST_dev=mean(totDevianceSST),
meanFlight=mean(totFlightHrs), meanFlight_dev=mean(totDevianceFlight),
meanForage=mean(totForageHrs), meanForage_dev=mean(totDevianceForage),
meanMigr=mean(MigratoryDistKm), meanDEE=mean(totalDEE), meanDev=mean(totDeviance),
meanSST_cov=mean(covSST), meanCov=mean(totCov)) %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::summarise(birds=n_distinct(individ_id), meanDEECol=mean(meanDEE), sdDEE=sd(meanDEE), seDEE=sdDEE/sqrt(birds), meanDevCol=mean(meanDev), meanCovCol=mean(meanCov)) %>%
dplyr::filter(birds>=minSampleSize) %>%
arrange(desc(meanCovCol)) %>%
dplyr::slice(c(1, n())) %>%
dplyr::select(species, colony)

DistributionData1<-totalCosts %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(meanSST=mean(SST), meanSST_dev=mean(totDevianceSST),
meanFlight=mean(totFlightHrs), meanFlight_dev=mean(totDevianceFlight),
meanForage=mean(totForageHrs), meanForage_dev=mean(totDevianceForage),
meanMigr=mean(MigratoryDistKm), meanDEE=mean(totalDEE), meanDev=mean(totDeviance),
meanSST_cov=mean(covSST), meanCov=mean(totCov)) %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::mutate(birds=n_distinct(individ_id), meanDEECol=mean(meanDEE), sdDEE=sd(meanDEE), seDEE=sdDEE/sqrt(birds), meanDevCol=mean(meanDev), meanCovCol=mean(meanCov)) %>%
arrange(meanCovCol) %>%
dplyr::inner_join(DistributionData2, by=c("species", "colony"), relationship = "many-to-many")

# I want to plot: weekly deviance in SST, energy expenditure, + activity budgets (same plot?), total energy expenditure, locations (gridded points?) + location of colonies with colony names 

# Determine location of processed locations (IRMA data) 
irma.files<-list.files("./data/positionsIRMA/", full.names=TRUE)
metadata<-load(irma.files[8])
irma.files<-irma.files[grepl("IRMAlocs", irma.files)]
irma.files.df<-data.frame(irma.files)
colnames(irma.files.df)<-c("FileName")
irma.files.df$species<-c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot")

	speciesList<-unique(DistributionData2$species) 

	for (i in 1:length(speciesList)) {

	# Subset species i
	speciesSub<-speciesList[i]
	speciesdf<-subset(DistributionData1, species==speciesSub)

	# Create plot 1: general deviance vs. DEE
	energyMetrics<-speciesdf %>%
	dplyr::mutate(metric="WEE_COV_NB", value=meanCov) 
	energyMetrics2<-speciesdf %>%
	dplyr::mutate(metric="TEE_NB", value=meanDEE) 
	energyMetricsAll<-rbind(energyMetrics, energyMetrics2)

	plota<-energyMetricsAll %>%
	dplyr::left_join(colonies_lox2, by=c("colony")) %>%
	ggplot(aes(x=colonyName, y=value)) +
	geom_boxplot() +
	facet_wrap(~metric, scales="free_y") +
	theme_bw() +
	ylab("Value") +
	xlab("Population code") +
	labs(tag="A)")

	print("Plotting deviation in behavior...")

	# Create plot 2 = showing deviance in behavior & sst # (one plot)

	devianceWeekly<-readRDS(deviation_weekly[grepl(speciesSub, deviation_weekly)])

	# Summarize mean deviation
	devianceColonyRep<-devianceWeekly %>%
	dplyr::inner_join(DistributionData2, by=c("species", "colony")) %>%
	  ungroup() %>%
	  dplyr::group_by(species, colony, individ_id, weekNo) %>%
	  dplyr::summarise(propFlight=mean(propFlight), weeklyDEE=mean(weeklyDEE), propActive=mean(propActive), propRest=mean(propRest), propForage=mean(propForage), meanSST=mean(meanSST), propLand=mean(propLand)) %>%
	  ungroup() %>%
	  dplyr::group_by(species, colony, individ_id) %>%
	  dplyr::mutate(meanPropFlight=mean(propFlight), meanPropActive=mean(propActive), meanPropRest=mean(propRest), meanPropLand=mean(propLand), meanPropForage=mean(propForage), annualSST=mean(meanSST), meanDEE=mean(weeklyDEE), sdDEE=sd(weeklyDEE)) %>%
	  ungroup() %>%
	  dplyr::group_by(species, colony, individ_id, weekNo) %>%
	  dplyr::mutate(devianceFlight=(propFlight - meanPropFlight)/meanPropFlight, devianceActive=(propActive - meanPropActive)/meanPropActive, devianceRest=(propRest - meanPropRest)/meanPropRest, devianceLand=(propLand - meanPropLand)/meanPropLand, devianceForage=(propForage - meanPropForage)/meanPropForage, devianceSST=(meanSST - annualSST)/annualSST,
	  devianceDEE=(weeklyDEE-meanDEE)/meanDEE, covDEE=(sdDEE/meanDEE)*100)  %>%
	  ungroup() 
	  
	 behv1<-devianceColonyRep %>%
	 dplyr::group_by(species, colony, weekNo) %>%
	 dplyr::summarise(birds=n_distinct(individ_id), mean=mean(devianceDEE), sd=sd(devianceDEE), se=sd/sqrt(birds)) %>%
	 dplyr::mutate(metric="Energy") 
	 
	 behv2<-devianceColonyRep %>%
	 dplyr::group_by(species, colony, weekNo) %>%
	 dplyr::summarise(birds=n_distinct(individ_id), mean=mean(meanSST), sd=sd(meanSST), se=sd/sqrt(birds)) %>%
	 dplyr::mutate(metric="SST")
	 
	 behv3<-devianceColonyRep %>%
	 dplyr::group_by(species, colony, weekNo) %>%
	 dplyr::summarise(birds=n_distinct(individ_id), mean=mean(propFlight), sd=sd(propFlight), se=sd/sqrt(birds)) %>%
	 dplyr::mutate(metric="Flight")
	 
	 behv4<-devianceColonyRep %>%
	 dplyr::group_by(species, colony, weekNo) %>%
	 dplyr::summarise(birds=n_distinct(individ_id), mean=mean(propForage), sd=sd(propForage), se=sd/sqrt(birds)) %>%
	 dplyr::mutate(metric="Forage")
 
 behv5<-devianceColonyRep %>%
 dplyr::group_by(species, colony, weekNo) %>%
 dplyr::summarise(birds=n_distinct(individ_id), mean=mean(devianceRest), sd=sd(devianceRest), se=sd/sqrt(birds)) %>%
 dplyr::mutate(metric="Rest")
 
 behv6<-devianceColonyRep %>%
 dplyr::group_by(species, colony, weekNo) %>%
 dplyr::summarise(birds=n_distinct(individ_id), mean=mean(devianceActive), sd=sd(devianceActive), se=sd/sqrt(birds)) %>%
 dplyr::mutate(metric="Active")
 
 behvAll<-rbind(behv1, behv2, behv3, behv4, behv5, behv6) %>%
 dplyr::left_join(colonies_lox2, by=c("colony")) 

plotb<-behvAll %>%
filter(metric %in% c("Energy")) %>%
#dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ggplot(aes(x=weekNo, y=mean)) +
geom_line(aes(group=interaction(metric,colony), color=colonyName)) +
geom_ribbon(aes(group=interaction(metric, colony), ymin=mean - 1.96*se, ymax=mean + 1.96*se, fill=colonyName), alpha=0.2) +
ylab(expression("WEE"[deviance]*"")) +
#facet_wrap(~colonyName) +
theme_bw() +
xlab("Week #") +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
labs(tag="B)")
  
# Now we plot energy vs. sst

meanSST<-median(subset(behvAll, metric =="SST")$mean)
meanDev<-mean(subset(behvAll, metric =="Energy")$mean)
coef1<-meanSST/meanDev

plotc<-speciesdf %>%
dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ggplot(aes(x=colonyName, y=meanSST)) +
geom_boxplot() +
theme_bw() +
ylab("Encountered sst (°C)") +
xlab("Population code") +
labs(tag="C)")

plotd<-behvAll %>%
filter(metric %in% c("SST")) %>%
#dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ggplot(aes(x=weekNo, y=mean)) +
geom_line(aes(group=interaction(metric,colony), color=colonyName)) +
geom_ribbon(aes(group=interaction(metric, colony), ymin=mean - 1.96*se, ymax=mean + 1.96*se, fill=colonyName), alpha=0.2) +
#facet_wrap(~colonyName) +
theme_bw() +
xlab("Week #") +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
labs(tag="D)") +
ylab("SST (°C)")

# Coef 1: between energy & flight
meanFlight<-median(subset(behvAll, metric =="Flight")$mean)
meanDev<-median(subset(behvAll, metric =="Energy")$mean)
coef2<-meanFlight/meanDev

# tot hours
totHours<-n_distinct(dates_weekly$dateKeep)*24

plote<-speciesdf %>%
dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ungroup() %>%
dplyr::mutate(propFlight=meanFlight/totHours) %>%
ggplot(aes(x=colonyName, y=propFlight)) +
geom_boxplot() +
theme_bw() +
ylab("Proportion time in flight") +
xlab("Population code") +
labs(tag="E)")

plotf<-behvAll %>%
filter(metric %in% c("Flight")) %>%
#dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ggplot(aes(x=weekNo, y=mean)) +
geom_line(aes(group=interaction(metric,colony), color=colonyName)) +
geom_ribbon(aes(group=interaction(metric, colony), ymin=mean - 1.96*se, ymax=mean + 1.96*se, fill=colonyName), alpha=0.2) +
#facet_wrap(~colonyName) +
theme_bw() +
xlab("Week #") +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
labs(tag="F)") +
ylab("Proportion of time in flight")

# Foraging if it's a surface-forager
# Coef 3: between energy & forage
meanForage<-median(subset(behvAll, metric =="Forage")$mean)
meanDev<-median(subset(behvAll, metric =="Energy")$mean)
coef3<-meanForage/meanDev

ploth<-behvAll %>%
filter(metric %in% c("Forage")) %>%
#dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ggplot(aes(x=weekNo, y=mean)) +
geom_line(aes(group=interaction(metric,colony), color=colonyName)) +
geom_ribbon(aes(group=interaction(metric, colony), ymin=mean - 1.96*se, ymax=mean + 1.96*se, fill=colonyName), alpha=0.2) +
#facet_wrap(~colonyName) +
theme_bw() +
xlab("Week #") +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
labs(tag="H)") +
ylab("Proportion time foraging")

# tot hours
totHours<-n_distinct(dates_weekly$dateKeep)*24

plotg<-speciesdf %>%
dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ungroup() %>%
dplyr::mutate(propForage=meanForage/totHours) %>%
ggplot(aes(x=colonyName, y=propForage)) +
geom_boxplot() +
theme_bw() +
ylab("Proportion time foraging") +
xlab("Population code") +
labs(tag="G)")

# Plot where they are: find IRMA locations
ids<-unique(devianceColonyRep$individ_id)

# Determine where locations are stored
allResults<-list.files("./tmp/", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# Make a list to save locations in
allLox<-list()

print("Finding locations...")

for (j in 1:length(ids)) {

energyOpen<-fread(energyRes_day[grepl(ids[j], energyRes_day)])

# Subset to relevant months #
energyMonth<-subset(energyOpen, Month %in% unique(dates_weekly$month))
energyMean<-energyMonth %>%
ungroup() %>%
dplyr::mutate(day=as.numeric(substr(date, 9, 10))) %>%
dplyr::mutate(month=Month) %>%
dplyr::inner_join(dates_weekly, by=c("day", "month")) %>%
dplyr::group_by(species, colony, individ_id, Month, day) %>%
dplyr::summarise(meanLon=mean(mean.lon), meanLat=mean(mean.lat), meanForage=mean(tForage), meanFlight=mean(tFlight), meanSST=mean(sst_random))

# Save results
allLox<-rbind(allLox, energyMean)

}

# project coordinates & save within some grid squares & note proportion of days in each square # 

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Set-up extent 
fake.rast<-raster(xmn=(-180), xmx=180, ymn=(-60), ymx=90, res=5)
values(fake.rast)<-0
proj4string(fake.rast)<-projection_84
WorldProject<-projectRaster(fake.rast, crs=projection_NA)
grid<-as.data.frame(WorldProject, xy=TRUE)
grid$loxCol1<-0	
grid$loxCol2<-0	
grid$timeFlight1<-0 # time spent in flight in hours at colony 1
grid$timeFlight2<-0 # time spent in flight in hours at colony 2
grid$timeForage1<-0 # time spent foraging in hours at colony 1
grid$timeForage2<-0 # time spent foraging in hours at colony 2
grid$meanSST1<-NA # mean SST encountered at colony 1
grid$meanSST2<-NA # mean SST at colony 2
grid<-grid %>%
arrange(x, y) 

# project locations
coordinates(allLox)<-~ meanLon + meanLat
proj4string(allLox)<-projection_84
allLoxTrans<-data.frame(spTransform(allLox, projection_NA))
colonies<-unique(allLoxTrans$colony)

# Attach all individual ids to accurately estimate behaviors
birdIds<-allLoxTrans %>%
ungroup() %>%
dplyr::group_by(species, colony) %>%
dplyr::count(individ_id)

print("Gridding locations...")

for (m in 1:nrow(grid)) {

# Subset grid x
gridSub<-grid[m,]

resx<-res(WorldProject)[1]
resy<-res(WorldProject)[2]

# Subset coordinates which fit #
loxSub1<-subset(allLoxTrans, colony==colonies[1] & coords.x1 > gridSub$x & coords.x1< gridSub$x + resx & coords.x2 > gridSub$y & coords.x2 < gridSub$y + resy)
loxSub2<-subset(allLoxTrans, colony==colonies[2] & coords.x1 > gridSub$x & coords.x1< gridSub$x + resx & coords.x2 > gridSub$y & coords.x2 < gridSub$y + resy)

if (nrow(loxSub1)>0) {

birdIdsSub1<-subset(birdIds, colony==colonies[1])

# Calculate time spent in flight
timeFlight<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(meanFlight)) %>%
dplyr::full_join(birdIdsSub1, by=c("species", "colony", "individ_id")) %>%
replace_na(list(totTime=0))

# Calculate time spent foraging
timeForage<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(meanForage)) %>%
dplyr::full_join(birdIdsSub1, by=c("species", "colony", "individ_id")) %>%
replace_na(list(totTime=0))

# Calculate SST
timeSST<-loxSub1 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=mean(meanSST)) 

# Save results in grid for plotting
grid$timeFlight1[m]<-mean(timeFlight$totTime)
grid$timeForage1[m]<-mean(timeForage$totTime)
grid$meanSST1[m]<-mean(timeSST$totTime)

}

if (nrow(loxSub2)>0) {

birdIdsSub2<-subset(birdIds, colony==colonies[1])

# Calculate time spent in flight
timeFlight2<-loxSub2 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(meanFlight)) %>%
dplyr::full_join(birdIdsSub2, by=c("species", "colony", "individ_id")) %>%
replace_na(list(totTime=0))

# Calculate time spent foraging
timeForage2<-loxSub2 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=sum(meanForage)) %>%
dplyr::full_join(birdIdsSub2, by=c("species", "colony", "individ_id")) %>%
replace_na(list(totTime=0))

# Calculate SST
timeSST2<-loxSub2 %>%
ungroup() %>%
dplyr::group_by(species, colony, individ_id) %>%
dplyr::summarise(totTime=mean(meanSST)) 

# save results in grid
grid$timeFlight2[m]<-mean(timeFlight2$totTime)
grid$timeForage2[m]<-mean(timeForage2$totTime)
grid$meanSST2[m]<-mean(timeSST2$totTime)

}

# Add number of locations
grid$loxCol1[m]<-nrow(loxSub1)
grid$loxCol2[m]<-nrow(loxSub2)

}

# Open country layers
print("Making map S1...")
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

lox1<-grid %>%
dplyr::mutate(LoxNo=loxCol1, timeFlight=timeFlight1, timeForage=timeForage1, sst=meanSST1) %>%
dplyr::mutate(colony=colonies[1])

lox2<-grid %>%
dplyr::mutate(LoxNo=loxCol2, timeFlight=timeFlight2, timeForage=timeForage2, sst=meanSST2) %>%
dplyr::mutate(colony=colonies[2])

loxAll<-rbind(lox1, lox2)

loxAll2<-loxAll %>%
dplyr::group_by(colony) %>%
dplyr::mutate(propLox=LoxNo/sum(LoxNo)) %>%
dplyr::left_join(colonies_lox2, by=c("colony"))  

# open colony information & join so we have longitude & latitude & project into correct format
colony.summary<-readRDS("./data/positionsIRMA/SEATRACK_export_20241120_ringInfo.rds")

speciesMatch<-data.frame(speciesLatin=c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle"), 
                         species=c("Brünnich's guillemot", "Black-legged kittiwake", "Common guillemot", "Atlantic puffin", "Northern fulmar", "Litte auk"))

colonyMatch<-colony.summary %>%
  dplyr::filter(species %in% c("Uria_lomvia", "Rissa_tridactyla", "Uria_aalge", "Fratercula_arctica", "Fulmarus_glacialis", "Alle_alle")) %>%
  dplyr::group_by(species, colony) %>%
  dplyr::slice(1) %>%
  dplyr::select(species, colony) %>%
  rename(speciesLatin=species) %>%
  dplyr::left_join(speciesMatch, by=c("speciesLatin"))
  
  colonies_lox<-colonyMatch %>%
  dplyr::left_join(colonyNames, by=c("colony")) %>%
  ungroup() %>%
  arrange(colonyName) %>%
  dplyr::group_by(colonyName) %>%
  dplyr::slice(1)
  
coordinates(colonies_lox)<-~col_lon + col_lat
proj4string(colonies_lox)<-projection_84
colonies_lox_trans<-data.frame(spTransform(colonies_lox, projection_NA))

loxAll3<-loxAll2 %>%
dplyr::left_join(colonies_lox_trans, by=c("colony")) %>%
dplyr::left_join(colonies_lox2, by=c("colony"))

positions<-filter(loxAll2, LoxNo>0)

if (speciesSub %in% c("Atlantic puffin", "Little auk", "Common guillemot", "Brünnich's guillemot")) {

plotj<-ggplot() +
  geom_tile(data=filter(loxAll2, LoxNo>0), aes(x=x, y=y, fill=sst)) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Non-breeding SST (°C)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(positions$x)-100000, max(positions$x) + 100000), ylim=c(min(positions$y)-100000, max(positions$y) + 100000)) +
  geom_point(data=loxAll3, aes(x=coords.x1, coords.x2), color="yellow") + 
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  facet_wrap(~colonyName) +
  labs(tag="J)") +
  theme(legend.position="bottom") 
  
  } 
  
 if (speciesSub %in% c("Black-legged kittiwake")) {

plotj<-ggplot() +
  geom_tile(data=filter(loxAll2, LoxNo>0), aes(x=x, y=y, fill=timeFlight)) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent in flight (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(positions$x)-100000, max(positions$x) + 100000), ylim=c(min(positions$y)-100000, max(positions$y) + 100000)) +
  geom_point(data=loxAll3, aes(x=coords.x1, coords.x2), color="yellow") + 
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  facet_wrap(~colonyName) +
  labs(tag="J)") +
  theme(legend.position="bottom") 

}  

if (speciesSub %in% c("Northern fulmar")) {

plotj<-ggplot() +
  geom_tile(data=filter(loxAll2, LoxNo>0), aes(x=x, y=y, fill=timeForage)) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent foraging (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(min(positions$x)-100000, max(positions$x) + 100000), ylim=c(min(positions$y)-100000, max(positions$y) + 100000)) +
  geom_point(data=loxAll3, aes(x=coords.x1, coords.x2), color="yellow") + 
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  facet_wrap(~colonyName) +
  labs(tag="J)") +
  theme(legend.position="bottom") 

}  
  
 
  
# Make a migratory distance plot (d)

ploti<-speciesdf %>%
dplyr::left_join(colonies_lox2, by=c("colony")) %>%
ggplot(aes(x=colonyName, y=meanMigr)) +
geom_boxplot() +
theme_bw() +
ylab("Migratory distance (km)") +
xlab("Population code") +
labs(tag="I)")
  
if (speciesSub %in% c("Northern fulmar", "Black-legged kittiwake")) { 
 
pdf(paste0("./results/figures/supplementary/pop_variation_", speciesSub, ".pdf"), height=10, width=8)
grid.arrange(plota, plotb, plote, plotf, plotg, ploth, ploti, plotj, ncol=2)
dev.off()

} else {

pdf(paste0("./results/figures/supplementary/pop_variation_", speciesSub, ".pdf"), height=10, width=8)
grid.arrange(plota, plotb, plotc, plotd, plote, plotf, ploti, plotj, ncol=2)
dev.off()
  
  }
  
  }
  
###  Figure S29: Time spent on water vs. dry ####

propBehaviors<-totalCosts %>%
   ungroup() %>%
   dplyr::group_by(rep, species, individ_id) %>%
   dplyr::mutate(propLandNB=totLandHrs/totHrs, propRestNB=totRestHrs/totHrs, propWet=sum(propRestNB, propActiveNB), propDry=sum(propFlightNB, propForageNB, propLandNB)) %>%
   ungroup() %>%
   dplyr::group_by(rep, species) %>%
   dplyr::summarise(meanWet=mean(propWet), sdWet=sd(propWet), meanDry=mean(propDry), sdDry=sd(propDry)) %>%
   ungroup() %>%
   dplyr::group_by(species) %>%
   dplyr::summarise(meanWet=mean(meanWet), sdWet=mean(sdWet), meanDry=mean(meanDry), sdDry=mean(sdDry)) 
   
propBehaviors<-propBehaviors %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

prop1<-propBehaviors %>%
dplyr::select(species, meanWet, sdWet) %>%
dplyr::rename(meanBehavior=meanWet, sdBehavior=sdWet) %>%
dplyr::mutate(Behavior="Wet") 

prop2<-propBehaviors %>%
dplyr::select(species, meanDry, sdDry) %>%
dplyr::rename(meanBehavior=meanDry, sdBehavior=sdDry) %>%
dplyr::mutate(Behavior="Dry") 

propAll<-rbind(prop1, prop2)
												 
FigureS29<-ggplot() +
  geom_col(
    data = propAll,
    aes(x = Behavior, y = meanBehavior, fill = species),
    position = position_dodge(width = 0.8)  # wider dodge
  ) +
  geom_errorbar(
    data = propAll,
    aes(
      x = Behavior,
      y = meanBehavior,
      ymin = meanBehavior - sdBehavior,
      ymax = meanBehavior + sdBehavior,
      group = species
    ),
    position = position_dodge(width = 0.8),  # match dodge exactly
    width = 0.2
  ) +
  scale_fill_manual(values = c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
  theme_bw() +
  xlab("") +
  ylab("Prop time spent in behaviour (mean +/- SD)")

pdf("./results/figures/supplementary/FigureS29.pdf", width=6.5, height=5)
plot(FigureS29)
dev.off()

###  Figure S30: Plot time spent below LCT ####

propWarm<-totalCosts %>%
   ungroup() %>%
   dplyr::group_by(rep, species) %>%
   dplyr::summarise(meanpropWarm=mean(propWarm), sdpropWarm=sd(propWarm))%>%  
   ungroup() %>%
   dplyr::group_by(species) %>%
   dplyr::summarise(meanpropWarm=mean(meanpropWarm), sdpropWarm=mean(sdpropWarm))  
   
propWarm<-propWarm %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))
												 
FigureS30<-ggplot() +
  geom_col(
    data = propWarm,
    aes(x = species, y = meanpropWarm, fill = species),
    position = position_dodge(width = 0.8)  # wider dodge
  ) +
  geom_errorbar(
    data = propWarm,
    aes(
      x = species,
      y = meanpropWarm,
      ymin = meanpropWarm - sdpropWarm,
      ymax = meanpropWarm + sdpropWarm,
      group = species
    ),
    position = position_dodge(width = 0.8),  # match dodge exactly
    width = 0.2
  ) +
  scale_fill_manual(values = c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
  theme_bw() +
  xlab("") +
  ylab("Prop time spent > LCT") +
 theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

pdf("./results/figures/supplementary/FigureS30.pdf", width=6.5, height=5)
plot(FigureS30)
dev.off()
  
# Save output files
print("Saving output files...")

# Number 1
output_file1 <- args[3]
print("Saving output file 1")
write.csv(totalCosts, file = output_file1, row.names = FALSE) # Total non-breeding costs for plotting + kruskal-wallis test + clustering at individual-level

# Number 2
output_file2 <- args[4]
print("Saving output file 2")
write.csv(lmRes_tot, file = output_file2, row.names = FALSE) # result of statistics (link between deviance & tot nonbreeding costs)

# Number 3
output_file3 <- args[5]
print("Saving output file 3")
write.csv(totalCosts_colony, file = output_file3, row.names = FALSE) # result of statistics (clustering at population level)

# Number 4
output_file4<- args[6]
print("Saving output file 4")
write.csv(lmRes_allspecies, file = output_file4, row.names = FALSE) # result of statistics (link between deviance & tot nonbreeding costs - species level)

# Number 5
output_file5<- args[7]
print("Saving output file 5")
write.csv(lmRes_tot2, file = output_file5, row.names = FALSE) # result of statistics (link between deviance, nonbreeding costs & possible explanatory variables)

# Number 6
output_file5<- args[8]
print("Saving output file 6")
write.csv(lmRes_tot3, file = output_file5, row.names = FALSE) # result of statistics (link between cov, nonbreeding costs)
