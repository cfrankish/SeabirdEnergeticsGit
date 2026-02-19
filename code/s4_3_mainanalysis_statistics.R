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
reps<-50
print(paste0("min iteration number is: ", reps))

### Step 1: Estimate total non-breeding costs ###

print("Step 1: Estimate total NB costs, weighted for population size..:")

# Create list of weeks to roll through
dates<-data.frame(dateKeep=seq(as.Date(startDate), as.Date(endDate), 1))
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
dates_weekly2<-rbind(date1, dates_weekly, date2)

# Within sessions, subset to 'trackyears' if possible
day1_month<-as.numeric(substr(day1, 6, 7))
day2_month<-as.numeric(substr(day2, 6, 7))

# Determine list of unique ids
allResults_deviation<-list.files("./results/tables/main/", full.names=TRUE)
deviation_base<-allResults_deviation[grepl("/deviance_", allResults_deviation)]

devianceMean<-list()

for (m in 1:length(deviation_base)) {

deviance_sub<-readRDS(deviation_base[m])

devianceMean<-rbind(devianceMean, deviance_sub)

}
ids<-unique(devianceMean$individ_id)

# Define number of species
speciesList<-unique(devianceMean$species)

# divide DEE by weight to get in kJ.g
species<-data.frame(species=c("Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot", "Little auk", "Atlantic puffin"))
species$allometryCoef<-c(0.765, 0.717, 0.689, 0.689, 0.689, 0.689)
species$LCT<-c(9, 12.5, 14.18, 14.18, 14.18, 14.18)

# Define location of weekly energy information
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
totalCosts<-list()
lmRes_tot<-list() # deviance vs DEE
lmRes_tot2<-list() # deviance & DEE vs other behaviors 
lmRes_tot3<-list() # cov vs. DEE
lmRes_allspecies<-list()
totalCosts_colony<-list()
randomEffects_all<-list()

for (i in 1:reps) {
  
  repSub<-i  
  
  # Loop through species
  
  totalCosts_reps<-list()
  lmRes_reps<-list()
  lmRes_reps2<-list()
  lmRes_reps3<-list()
  randomEffects_species<-list()
  
  for (j in 1:length(speciesList)) {
    
    #print(paste0("Calculating amplitude for rep ", i, " species ", j ))
    print(paste0("Rep ", i, " Species ", j, ": Estimating NB costs...")) 
    
    # Subset to species j
	speciesSub<-speciesList[j]
	devianceWeekly<-readRDS(deviation_weekly[grepl(speciesSub, deviation_weekly)])
    
    # Subset species data frame to rep i
    speciesDaily<-subset(devianceWeekly, rep==i)
    
    # Summarise total costs as well as mean SST
    speciesNB_stats<-speciesDaily %>%
	dplyr::left_join(species, by=c("species")) %>%
      ungroup() %>%
	  dplyr::group_by(rep, species, colony, individ_id) %>%
	  dplyr::mutate(meanpropRest=mean(propRest), meanpropLand=mean(propLand), meanPropFlight=mean(propFlight), meanPropForage=mean(propForage), meanNBSST=mean(meanSST), sdPropFlight=sd(propFlight), sdPropForage=sd(propForage),
	  sdNBSST=sd(meanSST)) %>%
	  ungroup() %>%
	  dplyr::group_by(rep, species, colony, individ_id, weekNo) %>%
	  dplyr::mutate(devianceFlight=abs(propFlight-meanPropFlight), devianceForage=abs(propForage-meanPropForage), devianceSST=abs(meanSST-meanNBSST), 
	  underLCT=ifelse(meanSST<LCT, 1, 0), underLCT_col=ifelse(meanSST_colony < LCT, 1, 0), weeks=n_distinct(weekNo)) %>%
	  ungroup() %>%
      dplyr::group_by(rep, species, colony, individ_id) %>%
      dplyr::summarise(totalDEE=sum(weeklyDEE), totalDEE_col=sum(weeklyDEE_col), SST=mean(meanSST), SSTsd=sd(meanSST), SST_col=mean(meanSST_colony), SSTdiff=SST-SST_col, 
	  totDeviance=max(totDeviance), totLandHrs=sum(propLand*7*24), totRestHrs=sum(propRest*7*24), totFlightHrs=sum(propFlight*7*24), totForageHrs=sum(propForage*7*24), totActiveHrs=sum(propActive*7*24),
	  totHrs=n_distinct(weekNo)*168, totDevianceFlight=sum(devianceFlight), totDevianceForage=sum(devianceForage), totDevianceSST=sum(devianceSST), 
	  meanWeeklyDEE=mean(weeklyDEE), sdWeeklyDEE=sd(weeklyDEE), totCov=(sd(weeklyDEE)/mean(weeklyDEE))*100, covFlight=(sd(propFlight)/mean(propFlight))*100, 
	  covForage=(sd(propForage)/mean(propForage))*100, covActive=(sd(propActive)/mean(propActive))*100, covSST=(sd(meanSST)/mean(meanSST))*100, 
	  propFlightNB=totFlightHrs/totHrs, propForageNB=totForageHrs/totHrs, propActiveNB=totActiveHrs/totHrs, propWarm=1-sum(underLCT)/n_distinct(weekNo), propCold_col=sum(underLCT_col)/n_distinct(weekNo))   
      #ungroup() %>%
      #dplyr::mutate(type=ifelse(propCold_col < 1, "warm", "cold")) 	
	
    # turn DEE from kj to kj.g & scale variables
	
	if (speciesNB_stats$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {
	
    speciesNB_stats5<-speciesNB_stats %>%
      dplyr::left_join(species, by=c("species")) %>%
	  dplyr::left_join(devianceMean, by=c("rep", "species", "colony", "individ_id")) %>%
      ungroup() %>%
	  dplyr::filter(propFlightNB>0) %>% # Otherwise there are nas...
	  #dplyr::filter(propActiveNB>0) %>% 
      dplyr::mutate(DEE_scale=scale(totalDEE), sst_scale=scale(SST), deviance_scale=scale(totDeviance), Migr_scale=scale(MigratoryDistKm), Flight_scale=scale(propFlightNB), 
	  Forage_scale=scale(propForageNB) , Active_scale=scale(propActiveNB), Flight_scale2=scale(covFlight), 
	  Forage_scale2=scale(covForage), sst_scale2=scale(SSTsd), cov_scale=scale(totCov), Active_scale2=scale(covActive), sst_start_scale=scale(SST_col), sst_diff_scale=scale(SSTdiff), sst_prop_scale=scale(propWarm))  
	  
	  } else {
	  
	  speciesNB_stats5<-speciesNB_stats %>%
      dplyr::left_join(species, by=c("species")) %>%
	  dplyr::left_join(devianceMean, by=c("rep", "species", "colony", "individ_id")) %>%
      ungroup() %>%
	  dplyr::filter(propFlightNB>0) %>% # Otherwise there are nas...
	  dplyr::filter(propActiveNB>0) %>% 
      dplyr::mutate(DEE_scale=scale(totalDEE), sst_scale=scale(SST), deviance_scale=scale(totDeviance), Migr_scale=scale(MigratoryDistKm), Flight_scale=scale(propFlightNB), 
	  Forage_scale=scale(propForageNB) , Active_scale=scale(propActiveNB), Flight_scale2=scale(covFlight), 
	  Forage_scale2=scale(covForage), sst_scale2=scale(SSTsd), cov_scale=scale(totCov), Active_scale2=scale(covActive),
	  sst_start_scale=scale(SST_col), sst_diff_scale=scale(SSTdiff), sst_prop_scale=scale(propWarm)) 
	  
	  
	  }
    
    # Remove colonies with small sample sizes
    speciesNB_stats6<-speciesNB_stats5 %>%
      ungroup() %>%
      dplyr::group_by(rep, species, colony) %>%
      dplyr::mutate(birdsTot=n_distinct(individ_id)) %>%
      dplyr::filter(birdsTot >= minSampleSize) %>%
	  dplyr::mutate(colonyWeights=1/birdsTot) %>%
	  dplyr::left_join(colonyCoords, by=c("colony")) %>%
	  dplyr::ungroup() %>%
	  dplyr::group_by(species) %>%
	  dplyr::mutate(lat_scale=scale(col_lat))
	  
	speciesNB_stats6$colony<-factor(speciesNB_stats6$colony)
	  
	# Stats # 1: Does total NB differ between populations? (kruskal-wallis)
	normalityTest<-shapiro.test(speciesNB_stats6$totalDEE)
	pval<-normalityTest$p.value
	speciesNB_stats6$nbcosts_normality<-pval

	# data is not normally distirbuted and we conduct a kruskal-wallis test #
	
	test1<-kruskal.test(totalDEE ~ colony, data = speciesNB_stats6)
	pvaltest<-test1$p.value
	df<-test1$parameter
	H<-test1$statistic
	speciesNB_stats6$popTest_nbcosts_pval<-pvaltest # Would be good to conduct a post-hoc test later to determine which populations are different
	speciesNB_stats6$popTest_nbcosts_df<-df 
	speciesNB_stats6$popTest_nbcosts_H<-H 
	effSize<-kruskal_effsize(data=speciesNB_stats6, formula=totalDEE ~ colony) # Calculate effect size
	speciesNB_stats6$popTest_nbcosts_effectSize<-effSize$effsize
	
    # Stats # 2: Does deviance differ between populations
	
	normalityTest2<-shapiro.test(speciesNB_stats6$totDeviance)
	pval2<-normalityTest2$p.value
	speciesNB_stats6$deviance_normality<-pval2

	# data is not normally distributed and we conduct a kruskal-wallis test #
	
	test2<-kruskal.test(totDeviance ~ colony, data = speciesNB_stats6)
	pvaltest2<-test2$p.value
	df2<-test2$parameter
	H2<-test2$statistic
	speciesNB_stats6$popTest_deviance_pval<-pvaltest2 # Would be good to conduct a post-hoc test later to determine which populations are different
    speciesNB_stats6$popTest_deviance_df<-df2
    speciesNB_stats6$popTest_deviance_H<-H2
	
	# Stats # 3: Does COV (coefficient of variation) differ between populations
	
	normalityTest3<-shapiro.test(speciesNB_stats6$totCov)
	pval3<-normalityTest3$p.value
	speciesNB_stats6$cov_normality<-pval3

	# data is not normally distributed and we conduct a kruskal-wallis test #
	
	test3<-kruskal.test(totCov ~ colony, data = speciesNB_stats6)
	pvaltest3<-test3$p.value
	df3<-test3$parameter
	H3<-test3$statistic
	speciesNB_stats6$popTest_cov_pval<-pvaltest3 # Would be good to conduct a post-hoc test later to determine which populations are different
    speciesNB_stats6$popTest_cov_df<-df3
    speciesNB_stats6$popTest_cov_H<-H3
	effSize2<-kruskal_effsize(data=speciesNB_stats6, formula=totCov ~ colony) # Calculate effect size
	speciesNB_stats6$popTest_cov_effectSize<-effSize2$effsize

    # Stats # 3: How are these two linked? 	
    
    # Fit some stats - link between DEE & various variables 
	
	# Check correlation first (if it's significant then I remove time spent in flight, otherwise I keep #
	corTravel<-cor(speciesNB_stats6$Migr_scale, speciesNB_stats6$sst_diff_scale)
	
    # Fit model

  #lm1<-lmer(DEE_scale ~ Migr_scale + sst_scale + Flight_scale + Forage_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  #lm1<-lmer(DEE_scale ~ Migr_scale + sst_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm1<-lmer(DEE_scale ~ Migr_scale*sst_diff_scale + lat_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm2<-lmer(DEE_scale ~ Migr_scale*sst_diff_scale +  (1 |colony), weights=colonyWeights, data=speciesNB_stats6)
  lm7<-lmer(DEE_scale ~ Migr_scale*sst_diff_scale +  (1 + sst_diff_scale|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm9<-lmer(DEE_scale ~ Migr_scale + sst_diff_scale +  (1 + sst_diff_scale|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm11<-lmer(DEE_scale ~ Migr_scale*sst_diff_scale + sst_start_scale +  (1 |colony), weights=colonyWeights, data=speciesNB_stats6)
  #lm1<-lmer(DEE_scale ~ Migr_scale + sst_diff_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  #lm1<-lmer(DEE_scale ~ Migr_scale + sst_diff_scale*sst_start_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)

  #lm2<-lmer(deviance_scale ~ Migr_scale + sst_scale + Flight_scale + Forage_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  #lm3<-lmer(deviance_scale ~ Migr_scale + sst_scale2 + Flight_scale2 + Forage_scale2 +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  #lm5<-lmer(cov_scale ~ Migr_scale + sst_scale2 + Flight_scale2 + Forage_scale2 +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  #lm5<-lmer(cov_scale ~ Migr_scale + sst_scale2  +   (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm3<-lmer(cov_scale ~ Migr_scale*sst_diff_scale + lat_scale  +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm5<-lmer(cov_scale ~ Migr_scale*sst_diff_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm8<-lmer(cov_scale ~ Migr_scale*sst_diff_scale +  (1 + sst_diff_scale|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm10<-lmer(cov_scale ~ Migr_scale + sst_diff_scale +  (1 + sst_diff_scale|colony), weights=colonyWeights, data=speciesNB_stats6)
  lm12<-lmer(cov_scale ~ Migr_scale*sst_diff_scale + sst_start_scale +  (1 |colony), weights=colonyWeights, data=speciesNB_stats6)
  #lm5<-lmer(cov_scale ~ Migr_scale + sst_diff_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)
  #lm5<-lmer(cov_scale ~ Migr_scale + sst_diff_scale*sst_start_scale +  (1|colony), weights=colonyWeights, data=speciesNB_stats6)

    # Link between DEE & cov
	
	lm4<-lmer(totalDEE ~ totDeviance + (1|colony), weights=colonyWeights, data=speciesNB_stats6)
	lm6<-lmer(totalDEE ~ totCov + (1|colony), weights=colonyWeights, data=speciesNB_stats6)
	lm6_2<-lmer(totalDEE ~ totCov + (1|colony), data=speciesNB_stats6)
    
    # Extract coefficients
    lm1_sum<-summary(lm1)
    lm1_sum2<-data.frame(coefficients(lm1_sum))
	sum1<-lm1_sum2
	sum1$predictors<-row.names(sum1)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm1)
    sum1$r2<-r2_values[1]
    sum1$rep<-speciesNB_stats6$rep[1]
    sum1$species<-speciesNB_stats6$species[1]
    sum1$test<-"DEE_vs_behavior"
    
    lm2_sum<-summary(lm2)
    lm2_sum<-data.frame(coefficients(lm2_sum))
	sum2<-lm2_sum
	sum2$predictors<-row.names(sum2)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm2)
    sum2$r2<-r2_values[1]
    sum2$rep<-speciesNB_stats6$rep[1]
    sum2$species<-speciesNB_stats6$species[1]
    sum2$test<-"DEE_vs_behavior2"
	
	lm3_sum<-summary(lm3)
    lm3_sum<-data.frame(coefficients(lm3_sum))
	sum3<-lm3_sum
	sum3$predictors<-row.names(sum3)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm3)
    sum3$r2<-r2_values[1]
    sum3$rep<-speciesNB_stats6$rep[1]
    sum3$species<-speciesNB_stats6$species[1]
    sum3$test<-"Deviance_vs_behavior2"
	
	lm5_sum<-summary(lm5)
    lm5_sum<-data.frame(coefficients(lm5_sum))
	sum5<-lm5_sum
	sum5$predictors<-row.names(sum5)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm5)
    sum5$r2<-r2_values[1]
    sum5$rep<-speciesNB_stats6$rep[1]
    sum5$species<-speciesNB_stats6$species[1]
    sum5$test<-"COV_vs_behavior"
	
	lm7_sum<-summary(lm7)
    lm7_sum<-data.frame(coefficients(lm7_sum))
	sum7<-lm7_sum
	sum7$predictors<-row.names(sum7)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm7)
    sum7$r2<-r2_values[1]
    sum7$rep<-speciesNB_stats6$rep[1]
    sum7$species<-speciesNB_stats6$species[1]
    sum7$test<-"DEE_vs_behavior3"
	
	# Extract random effects #
	re_col1 <- ranef(lm7)$colony %>% as.data.frame() %>% rownames_to_column("colony")
	slopes1 <- re_col1 %>%
    mutate(
    sst_diff_scale = fixef(lm7)["sst_diff_scale"] + sst_diff_scale
    ) %>%
    dplyr::select(colony, sst_diff_scale)
	slopes1$test<-"DEE_vs_behavior3"
	slopes1$species<-speciesSub
	slopes1$rep<-i
	slopes1_2<-slopes1 %>%
	dplyr::select(rep, species, test, colony, sst_diff_scale) %>%
	rename(randomEffect=sst_diff_scale) %>%
	dplyr::mutate(predictor="sst_diff_scale")
	slopes1All<-rbind(slopes1_2)

	lm8_sum<-summary(lm8)
    lm8_sum<-data.frame(coefficients(lm8_sum))
	sum8<-lm8_sum
	sum8$predictors<-row.names(sum8)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm8)
    sum8$r2<-r2_values[1]
    sum8$rep<-speciesNB_stats6$rep[1]
    sum8$species<-speciesNB_stats6$species[1]
    sum8$test<-"COV_vs_behavior2"
	
	# Extract random effects #
	re_col2 <- ranef(lm8)$colony %>% as.data.frame() %>% rownames_to_column("colony")
	slopes2 <- re_col2 %>%
    mutate(
    sst_diff_scale = fixef(lm8)["sst_diff_scale"] + sst_diff_scale
    ) %>%
    dplyr::select(colony, sst_diff_scale)
	slopes2$test<-"COV_vs_behavior2"
	slopes2$species<-speciesSub
	slopes2$rep<-i
	slopes2_2<-slopes2 %>%
	dplyr::select(rep, species, test, colony, sst_diff_scale) %>%
	rename(randomEffect=sst_diff_scale) %>%
	dplyr::mutate(predictor="sst_diff_scale")
	slopes2All<-rbind(slopes2_2)
	
	## Model 9 ##
	lm9_sum<-summary(lm9)
    lm9_sum<-data.frame(coefficients(lm9_sum))
	sum9<-lm9_sum
	sum9$predictors<-row.names(sum9)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm9)
    sum9$r2<-r2_values[1]
    sum9$rep<-speciesNB_stats6$rep[1]
    sum9$species<-speciesNB_stats6$species[1]
    sum9$test<-"DEE_vs_behavior4"
	
	# Extract random effects #
	re_col3 <- ranef(lm9)$colony %>% as.data.frame() %>% rownames_to_column("colony")
	slopes3 <- re_col3 %>%
    mutate(
    sst_diff_scale = fixef(lm9)["sst_diff_scale"] + sst_diff_scale
    ) %>%
    dplyr::select(colony, sst_diff_scale)
	slopes3$test<-"DEE_vs_behavior4"
	slopes3$species<-speciesSub
	slopes3$rep<-i
	slopes3_2<-slopes3 %>%
	dplyr::select(rep, species, test, colony, sst_diff_scale) %>%
	rename(randomEffect=sst_diff_scale) %>%
	dplyr::mutate(predictor="sst_diff_scale")
	slopes3All<-rbind(slopes1_2)

	lm10_sum<-summary(lm10)
    lm10_sum<-data.frame(coefficients(lm10_sum))
	sum10<-lm10_sum
	sum10$predictors<-row.names(sum10)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm10)
    sum10$r2<-r2_values[1]
    sum10$rep<-speciesNB_stats6$rep[1]
    sum10$species<-speciesNB_stats6$species[1]
    sum10$test<-"COV_vs_behavior3"
	
	# Extract random effects #
	re_col4 <- ranef(lm10)$colony %>% as.data.frame() %>% rownames_to_column("colony")
	slopes4 <- re_col4 %>%
    mutate(
    sst_diff_scale = fixef(lm10)["sst_diff_scale"] + sst_diff_scale
    ) %>%
    dplyr::select(colony, sst_diff_scale)
	slopes4$test<-"COV_vs_behavior3"
	slopes4$species<-speciesSub
	slopes4$rep<-i
	slopes4_2<-slopes4 %>%
	dplyr::select(rep, species, test, colony, sst_diff_scale) %>%
	rename(randomEffect=sst_diff_scale) %>%
	dplyr::mutate(predictor="sst_diff_scale")
	slopes4All<-rbind(slopes4_2)
	
	## Model 9 ##
	lm11_sum<-summary(lm11)
    lm11_sum<-data.frame(coefficients(lm11_sum))
	sum11<-lm11_sum
	sum11$predictors<-row.names(sum11)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm11)
    sum11$r2<-r2_values[1]
    sum11$rep<-speciesNB_stats6$rep[1]
    sum11$species<-speciesNB_stats6$species[1]
    sum11$test<-"DEE_vs_behavior5"
	
	## Model 9 ##
	lm12_sum<-summary(lm12)
    lm12_sum<-data.frame(coefficients(lm12_sum))
	sum12<-lm12_sum
	sum12$predictors<-row.names(sum12)
    #sum1<-data.frame(Migr_scale=seq(0, 1, 0.1))
	#sum1$colony<-speciesNB_stats6$colony[1]
    #sum1$fit<-predict(lm1, newdata=sum1, re.form=NA)
    r2_values <- r.squaredGLMM(lm12)
    sum12$r2<-r2_values[1]
    sum12$rep<-speciesNB_stats6$rep[1]
    sum12$species<-speciesNB_stats6$species[1]
    sum12$test<-"COV_vs_behavior4"
	 
	lm4_sum<-summary(lm4)
    lm4_sum4<-data.frame(coefficients(lm4_sum))
    sum4<-data.frame(totDeviance=seq(min(devianceWeekly$totDeviance), max(devianceWeekly$totDeviance), 0.5))
	sum4$colony<-speciesNB_stats6$colony[1]
    sum4$fit<-predict(lm4, newdata=sum4, re.form=NA)
	se<-predict(lm4, newdata=sum4, re.form=NA, se.fit=TRUE)
	sum4$se<-se$se.fit
    sum4$r2<-lm4_sum$r.squared
    r2_values <- r.squaredGLMM(lm4)
    sum4$r2<-r2_values[1]
    sum4$rep<-speciesNB_stats6$rep[1]
    sum4$species<-speciesNB_stats6$species[1]
    sum4$predictor<-"totDeviance"
    sum4$pvalue<-lm4_sum4$Pr...t..[2]
    sum4<-sum4%>%
      rename(predictor_val=totDeviance)
	sum4$intercept<-coefficients(lm4_sum)[1,1]
	sum4$coefficient<-coefficients(lm4_sum)[2,1]
	
	devianceWeekly$totCov<-devianceWeekly$covDEE
	lm6_sum<-summary(lm6)
    lm6_sum6<-data.frame(coefficients(lm6_sum))
    sum6<-data.frame(totCov=seq(min(devianceWeekly$covDEE), max(devianceWeekly$covDEE), 0.5))
	sum6$colony<-speciesNB_stats6$colony[1]
    sum6$fit<-predict(lm6, newdata=sum6, re.form=NA)
	se<-predict(lm6, newdata=sum6, re.form=NA, se.fit=TRUE)
	sum6$se<-se$se.fit
    sum6$r2<-lm6_sum$r.squared
    r2_values <- r.squaredGLMM(lm6)
    sum6$r2<-r2_values[1]
    sum6$rep<-speciesNB_stats6$rep[1]
    sum6$species<-speciesNB_stats6$species[1]
    sum6$predictor<-"totCov"
    sum6$pvalue<-lm6_sum6$Pr...t..[2]
    sum6<-sum6%>%
      rename(predictor_val=totCov)
	sum6$intercept<-coefficients(lm6_sum)[1,1]
	sum6$coefficient<-coefficients(lm6_sum)[2,1]
	sum6$type<-"Weights"
	
	devianceWeekly$totCov<-devianceWeekly$covDEE
	lm6_sum2<-summary(lm6_2)
    lm6_sum6_2<-data.frame(coefficients(lm6_sum2))
    sum6_2<-data.frame(totCov=seq(min(devianceWeekly$covDEE), max(devianceWeekly$covDEE), 0.5))
	sum6_2$colony<-speciesNB_stats6$colony[1]
    sum6_2$fit<-predict(lm6_2, newdata=sum6_2, re.form=NA)
	se<-predict(lm6_2, newdata=sum6_2, re.form=NA, se.fit=TRUE)
	sum6_2$se<-se$se.fit
    sum6_2$r2<-lm6_sum2$r.squared
    r2_values <- r.squaredGLMM(lm6_2)
    sum6_2$r2<-r2_values[1]
    sum6_2$rep<-speciesNB_stats6$rep[1]
    sum6_2$species<-speciesNB_stats6$species[1]
    sum6_2$predictor<-"totCov"
    sum6_2$pvalue<-lm6_sum6_2$Pr...t..[2]
    sum6_2<-sum6_2%>%
      rename(predictor_val=totCov)
	sum6_2$intercept<-coefficients(lm6_sum2)[1,1]
	sum6_2$coefficient<-coefficients(lm6_sum2)[2,1]
	sum6_2$type<-"NoWeights"
    
    # Collate results
    lmRes<-rbind(sum4)
	lmRes2<-rbind(sum1, sum2, sum3, sum5, sum7, sum8, sum9, sum10, sum11, sum12)
	lmRes3<-rbind(sum6, sum6_2) # COV option (testing)
	
    # Save results of models
    lmRes_reps<-rbind(lmRes_reps, lmRes)
	lmRes_reps2<-rbind(lmRes_reps2, lmRes2)
	lmRes_reps3<-rbind(lmRes_reps3, lmRes3)
    
    # Save results : raw data for plotting 
	speciesNB_stats6$rep<-i
	speciesNB_stats6$corTravel<-corTravel
    totalCosts_reps<-rbind(totalCosts_reps, speciesNB_stats6)
	
	# Save random effect results
	 randomEffects_species<-rbind(slopes1All, slopes2All, slopes3All, slopes4All)
    
  }
  
  # Conduct stats on all species (effect of deviance on non-breeding costs - but maybe should have been a simple correlation?)
  
  totalCosts_reps_sum<-totalCosts_reps %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>%
  dplyr::mutate(birds=n_distinct(individ_id), weights=1/birds) %>%
  dplyr::mutate(colonySp=paste0(species, "_", colony))
  
  lm4<-lmer(DEE_scale ~ deviance_scale*species + (1|colonySp), weights=weights, data=totalCosts_reps_sum)
  
  d1<-dredge(lm4)
  d2<-subset(d1, !nested(.))
  d3<-subset(d2, delta<2)
  d3<-d3[1,]
   
  # Extract significant parameters
  PredictorsModel<-data.frame(t(d3))
  PredictorsModel$Predictors<-row.names(PredictorsModel)
  PredictorsModelFinal<-subset(PredictorsModel, !Predictors %in% c("(Intercept)", "df", "logLik", "AICc", "delta", "weight"))
  colnames(PredictorsModelFinal)<-c("Estimate", "Predictors")
  PredictorsModelFinal<-subset(PredictorsModelFinal, !is.na(Estimate))
    
  fixed_effects <- paste(PredictorsModelFinal$Predictors, collapse = " + ")
  random_effects <- " + (1|colonySp)"
  full_formula <- as.formula(paste0("DEE_scale ~", fixed_effects, random_effects))
    
  # Create best model
  lmFinal <- lmer(full_formula,   
                   weights=weights, data=totalCosts_reps_sum)
  
  # Make some predictions
  predictions<-data.frame(deviance_scale=rep(seq(0, 1, 0.1), 6), species=rep(unique(totalCosts_reps_sum$species), each=11), totalCosts_reps_sum$colonySp[1])
  predictions$fit<-predict(lmFinal, newdata=predictions, re.form=NA)
  predictions$rep<-i
  r2_values <- r.squaredGLMM(lmFinal)
  predictions$r2<-r2_values[1]
  
  # Summarize raw data for plotting later
  
  colonyData<-totalCosts_reps_sum %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>%
  dplyr::summarise(meanDEE=mean(DEE_scale), meanDev=mean(deviance_scale), meanCov=mean(cov_scale)) %>%
  ungroup() %>%
  dplyr::select(meanDEE, meanDev, meanCov) 
  
  individualData<-totalCosts_reps_sum %>%
  ungroup() %>%
  dplyr::select(DEE_scale, deviance_scale, cov_scale) 
  
  colonyDataSum<-totalCosts_reps_sum %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>%
  dplyr::summarise(meanDEE=mean(DEE_scale), meanDev=mean(deviance_scale), meanCov=mean(cov_scale)) %>%
  ungroup() %>%
  dplyr::mutate(rep=i)

individualDataSum<-totalCosts_reps_sum %>%
  ungroup() %>%
  dplyr::mutate(rep=i)
  
  # Save all results
  totalCosts_colony<-rbind(totalCosts_colony, colonyDataSum) # population-level clustering results
  totalCosts<-rbind(totalCosts, individualDataSum) # has raw data & a lot of other information (individual clustering results + kruskal-wallis tests)
  lmRes_tot<-rbind(lmRes_tot, lmRes_reps) # species seperatey, deviance vs. DEE
  lmRes_tot2<-rbind(lmRes_tot2, lmRes_reps2) # species seperatey, deviance & DEE vs other behaviors
  lmRes_tot3<-rbind(lmRes_tot3, lmRes_reps3) # species seperatey, cov vs. DEE
  lmRes_allspecies<-rbind(lmRes_allspecies, predictions) # species together 
  randomEffects_all<-rbind(randomEffects_all, randomEffects_species) # save random effects
  
} 

print("Preparing datasets for plotting...")

#### Main results ##### (Figure 2A, B & C)

# Make species mean
speciesRes<-totalCosts %>%
   ungroup() %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::summarise(birds=1/n_distinct(individ_id), meanDEE=mean(totalDEE), meanDev=mean(totDeviance), meanCov=mean(totCov), sdDEE=sd(totalDEE), sdDev=sd(totDeviance), sdCov=sd(totCov), totBirds=n_distinct(individ_id), pval1=mean(popTest_nbcosts_pval), pval2=mean(popTest_deviance_pval), pval3=mean(popTest_cov_pval)) %>%
  dplyr::filter(totBirds>=minSampleSize) %>%
  ungroup() %>%
  dplyr::group_by(rep, species) %>%
  dplyr::summarise(weightedMean=sum(meanDEE*birds)/sum(birds), minDEE=min(meanDEE), maxDEE=max(meanDEE), meanpval1=mean(pval1), sdpval1=sd(pval1), sepval1=sdpval1/sqrt(reps), weightedMean2=sum(meanDev*birds)/sum(birds), minDev=min(meanDev), maxDev=max(meanDev), meanpval2=mean(pval2), 
  sdpval2=sd(pval2), sepval2=sdpval2/sqrt(reps), weightedMean3=sum(meanCov*birds)/sum(birds), minCov=min(meanCov), maxCov=max(meanCov), meanpval3=mean(pval3), sdpval3=sd(pval3), sepval3=sdpval3/sqrt(reps),
  n_effective=sum(birds)^2/(sum(birds^2)), total_var1=sum(birds*(meanDEE - weightedMean)^2)/sum(n_effective-1), weightedSD_DEE=sqrt(total_var1), weightedSE_DEE=weightedSD_DEE/sqrt(sum(birds)),
  total_var2=sum(birds*(meanDev - weightedMean2)^2)/sum(n_effective-1), weightedSD_Dev=sqrt(total_var2), weightedSE_Dev=weightedSD_Dev/sqrt(sum(birds)),
  total_var3=sum(birds*(meanCov - weightedMean3)^2)/sum(n_effective-1), weightedSD_Cov=sqrt(total_var3), weightedSE_Cov=weightedSD_Cov/sqrt(sum(birds))) %>%
  arrange(species, rep) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(meanDEEg=mean(weightedMean), sd=mean(weightedSD_DEE), reps=n_distinct(rep), se=mean(weightedSE_DEE), minDEE=mean(minDEE), maxDEE=mean(maxDEE), pval1=mean(meanpval1), sdpval1=sd(meanpval1),sepval1=sdpval1/sqrt(reps),
  meanDev=mean(weightedMean2), sd2=mean(weightedSD_Dev), reps=n_distinct(rep), se2=mean(weightedSE_Dev), minDev=mean(minDev), maxDev=mean(maxDev), pval2=mean(meanpval2), sdpval2=sd(meanpval2), sepval2=sdpval2/sqrt(reps),
  meanCov=mean(weightedMean3), sd3=mean(weightedSD_Cov), se3=mean(weightedSE_Cov), minCov=mean(minCov), maxCov=mean(maxCov), pval3=mean(meanpval3), sdpval3=sd(meanpval3), sepval3=sdpval3/sqrt(reps)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::mutate(sig1=ifelse(pval1 + 1.96*sepval1<0.05, 1, 0)) %>%
  dplyr::mutate(sig2=ifelse(pval2 + 1.96*sepval2 <0.05, 1, 0)) %>%
  dplyr::mutate(sig3=ifelse(pval3 + 1.96*sepval3 <0.05, 1, 0))
  
# Make colony mean 
colonyRes<-totalCosts %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::summarise(birds=n_distinct(individ_id), meanDEE=mean(totalDEE), minDEE=min(totalDEE), maxDEE=max(totalDEE), sdDEE=sd(totalDEE), seDEE=sdDEE/sqrt(birds),
  meanSSTdiff=mean(SSTdiff), minSSTdiff=min(SSTdiff), maxSSTdiff=max(SSTdiff), sdSSTdiff=sd(SSTdiff), seSSTdiff=sdSSTdiff/sqrt(birds),
  meanSST_col=mean(SST_col), minSST_col=min(SST_col), maxSST_col=max(SST_col), sdSST_col=sd(SST_col), seSST_col=sdSST_col/sqrt(birds),
  meanSSTsd=mean(SSTsd), minSSTsd=min(SSTsd), maxSSTsd=max(SSTsd), sdSSTsd=sd(SSTsd), seSSTsd=sdSSTsd/sqrt(birds),
  meanMigratoryDistKm=mean(MigratoryDistKm), minMigratoryDistKm=min(MigratoryDistKm), maxMigratoryDistKm=max(MigratoryDistKm), sdMigratoryDistKm=sd(MigratoryDistKm), seMigratoryDistKm=sdMigratoryDistKm/sqrt(birds),
  meanPropFlight=mean(propFlightNB), meanPropForage=mean(propForageNB),
  meanDev=mean(totDeviance), minDev=min(totDeviance), maxDev=max(totDeviance), sdDev=sd(totDeviance), seDev=sdDev/sqrt(birds),
  meanCov=mean(totCov), minCov=min(totCov), maxCov=max(totCov), sdCov=sd(totCov), seCov=sdCov/sqrt(birds),
  meanDEE2=mean(DEE_scale), minDEE2=min(DEE_scale), maxDEE2=max(DEE_scale), sdDEE2=sd(DEE_scale), seDEE2=sdDEE2/sqrt(birds),
  meanDev2=mean(deviance_scale), minDev2=min(deviance_scale), maxDev2=max(deviance_scale), sdDev2=sd(deviance_scale), seDev2=sdDev2/sqrt(birds),
  meanCov2=mean(cov_scale), minCov2=min(cov_scale), maxCov2=max(cov_scale), sdCov2=sd(cov_scale), seCov2=sdCov2/sqrt(birds), ) %>% 
  dplyr::filter(birds>=minSampleSize) %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>%
  dplyr::summarise(meanDEEg=mean(meanDEE), sd=mean(sdDEE), seDEE=mean(seDEE), minDEE=mean(minDEE), maxDEE=mean(maxDEE),
  meanSSTcol=mean(meanSST_col), meanSSTdiff=mean(meanSSTdiff), meanMigratoryDist=mean(meanMigratoryDistKm), meanSSTsd=mean(meanSSTsd), meanPropFlight=mean(meanPropFlight),meanPropForage=mean(meanPropForage),
  meanDev=mean(meanDev), sd2=mean(sdDev), seDev=mean(seDev), minDev=mean(minDev), maxDev=mean(maxDev),
  meanCov=mean(meanCov), sd3=mean(sdCov), seCov=mean(seCov), minCov=mean(minCov), maxCov=mean(maxCov),
  meanDEE2=mean(meanDEE2), seDEE2=mean(seDEE2),
  meanDev2=mean(meanDev2), seDev2=mean(seDev2), 
  meanCov2=mean(meanCov2), seCov2=mean(seCov2), meanBirds=mean(birds)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  ungroup() %>%
  dplyr::group_by(species) %>%
  dplyr::arrange(meanDEEg)%>%
  dplyr::mutate(colonySp=paste0(colony, "_", species))
  
  colonyRes2<-totalCosts %>%
  dplyr::group_by(rep, species, colony) %>%
  dplyr::summarise(birds=n_distinct(individ_id), meanDEE=mean(totalDEE), minDEE=min(totalDEE), maxDEE=max(totalDEE), sdDEE=sd(totalDEE), seDEE=sdDEE/sqrt(birds),
  meanDev=mean(totDeviance), minDev=min(totDeviance), maxDev=max(totDeviance), sdDev=sd(totDeviance), seDev=sdDev/sqrt(birds),
  meanCov=mean(totCov), minCov=min(totCov), maxCov=max(totCov), sdCov=sd(totCov), seCov=sdCov/sqrt(birds), 
  meanDEE2=mean(DEE_scale), minDEE2=min(DEE_scale), maxDEE2=max(DEE_scale), sdDEE2=sd(DEE_scale), seDEE2=sdDEE2/sqrt(birds),
  meanDev2=mean(deviance_scale), minDev2=min(deviance_scale), maxDev2=max(deviance_scale), sdDev2=sd(deviance_scale), seDev2=sdDev2/sqrt(birds),
  meanCov2=mean(cov_scale), minCov2=min(cov_scale), maxCov2=max(cov_scale), sdCov2=sd(cov_scale), seCov2=sdCov2/sqrt(birds)) %>% 
  dplyr::filter(birds>=minSampleSize) %>%
  ungroup() %>%
  dplyr::group_by(species, colony) %>%
  dplyr::summarise(meanDEEg=mean(meanDEE), sd=mean(sdDEE), seDEE=mean(seDEE), minDEE=mean(minDEE), maxDEE=mean(maxDEE),
  meanDev=mean(meanDev), sd2=mean(sdDev), seDev=mean(seDev), minDev=mean(minDev), maxDev=mean(maxDev),
  meanCov=mean(meanCov), sd3=mean(sdCov), seCov=mean(seCov), minCov=mean(minCov), maxCov=mean(maxCov),
  meanDEE2=mean(meanDEE2), seDEE2=mean(seDEE2),
  meanDev2=mean(meanDev2), seDev2=mean(seDev2),
  meanCov2=mean(meanCov2), seCov2=mean(seCov2)) %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  ungroup() %>%
  dplyr::group_by(species) %>%
  dplyr::arrange(meanDev)%>%
  dplyr::mutate(colonySp=paste0(colony, "_", species))										 

# Order colonyRes populations according to value
colonyOrder<-unique(colonyRes$colonySp)
colonyRes$colonySp<-factor(colonyRes$colonySp, levels=colonyOrder)	

# Add stars to denote how significant #
speciesRes$stars<-ifelse(speciesRes$pval1 + 1.96*speciesRes$sepval1<0.05, "*", "")
speciesRes$stars<-ifelse(speciesRes$pval1 + 1.96*speciesRes$sepval1<0.01, "**", speciesRes$stars)
speciesRes$stars<-ifelse(speciesRes$pval1 + 1.96*speciesRes$sepval1<0.001, "***", speciesRes$stars)

print("Saving Figures 2A & B")
		
Figure2B<-ggplot() +
   geom_pointrange(data=colonyRes, aes(x=meanDEEg, y=species, xmin=meanDEEg - 1.96*seDEE, xmax=meanDEEg + 1.96*seDEE, color=species), position = position_dodge2(width=0.6), cex=0.2, alpha=0.2) +
  geom_pointrange(data=speciesRes, aes(x=meanDEEg, y=species, xmin=meanDEEg - 1.96*se, xmax=meanDEEg + 1.96*se, color=species), position = position_dodge2(width=0.6), cex=0.4) +
  scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+
  geom_text(data=filter(speciesRes, sig1==1), aes(x=5900, y=species, label=stars)) +
  theme_bw() +
  xlab(expression("TEE"[NB]*" (kJ.g"^-1*")")) +
  ylab("") +
  theme(legend.position = "none") +
  ggtitle("B)") + 
  theme(
       panel.grid = element_blank(),
		#axis.text.x = element_text(angle = 45, hjust = 1),
		axis.text=element_text(size=12),
		axis.title=element_text(size=12)) 
		#ylim(1000, 6300)
		
		#axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),

pdf("./results/figures/main/Figure2B.pdf", width=7, height=5)
plot(Figure2B)
dev.off()

# Coefficient of variation version

# Add stars to denote how significant #
speciesRes$stars3<-ifelse(speciesRes$pval3 + 1.96*speciesRes$sepval3<0.05, "*", "")
speciesRes$stars3<-ifelse(speciesRes$pval3 + 1.96*speciesRes$sepval3<0.01, "**", speciesRes$stars3)
speciesRes$stars3<-ifelse(speciesRes$pval3 + 1.96*speciesRes$sepval3<0.001, "***", speciesRes$stars3)
		
Figure2A<-ggplot() +
 geom_pointrange(data=colonyRes2, aes(x=meanCov, y=species, xmin=meanCov - 1.96*seCov, xmax=meanCov + 1.96*seCov, color=species), position = position_dodge2(width=0.6), cex=0.2, alpha=0.2) +
 geom_pointrange(data=speciesRes, aes(x=meanCov, y=species, xmin=meanCov - se3, xmax=meanCov + se3, color=species), position = position_dodge2(width=0.6), cex=0.4) +
 scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+
 geom_text(data=filter(speciesRes, sig3==1), aes(x=34, y=species, label=stars3)) +
 theme_bw() +
 xlab(expression("WEE"[COV_NB]*"")) +
 ylab("") +
 theme(legend.position = "none") +
 ggtitle("A)") + 
  theme(
		panel.grid = element_blank(),
		axis.text=element_text(size=12),
		axis.title=element_text(size=12),
		#axis.text.x = element_text(angle = 45, hjust = 1)) +
		ylim(0, 35))

pdf("./results/figures/main/Figure2A.pdf", width=7, height=5)
plot(Figure2A)
dev.off()

### Figure 2C... ###

# Summarize results of deviance vs total costs
lmRes_sum<-lmRes_tot %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::group_by(species, predictor, predictor_val) %>%
  dplyr::summarise(meanEstimate=mean(fit), meanSE=mean(se), mean.p=mean(pvalue), sd.p=sd(pvalue), reps=n_distinct(rep), se.p=sd.p/sqrt(reps),
                   meanr2=mean(r2), sdr2=sd(r2), ser2=sdr2/sqrt(reps)) %>%
  dplyr::mutate(sig=ifelse(mean.p + 1.96*se.p<0.05, 1, 0)) 

# Shorten estimate lines to range of data

minAEE<-colonyRes_lox %>%
dplyr::group_by(species) %>%
arrange(meanDEEg) %>%
dplyr::slice(1) %>%
dplyr::select(species, meanDEEg, seDEE) %>%
rename(minAEE=meanDEEg, seDEEmin=seDEE)

maxAEE<-colonyRes_lox %>%
dplyr::group_by(species) %>%
arrange(meanDEEg) %>%
dplyr::slice(n()) %>%
dplyr::select(species, meanDEEg, seDEE) %>%
rename(maxAEE=meanDEEg, seDEEmax=seDEE)

minDev<-colonyRes_lox %>%
dplyr::group_by(species) %>%
arrange(meanDev) %>%
dplyr::slice(1) %>%
dplyr::select(species, meanDev) %>%
rename(minDev=meanDev)

maxDev<-colonyRes_lox %>%
dplyr::group_by(species) %>%
arrange(meanDev) %>%
dplyr::slice(n()) %>%
dplyr::select(species, meanDev) %>%
rename(maxDev=meanDev)

minCov<-colonyRes_lox %>%
dplyr::group_by(species) %>%
arrange(meanCov) %>%
dplyr::slice(1) %>%
dplyr::select(species, meanCov, seCov) %>%
rename(minCov=meanCov, minCovSE=seCov)

maxCov<-colonyRes_lox %>%
dplyr::group_by(species) %>%
arrange(meanCov) %>%
dplyr::slice(n()) %>%
dplyr::select(species, meanCov, seCov) %>%
rename(maxCov=meanCov, maxCovSE=seCov)

lmRes_sum3<-lmRes_tot3 %>%
  dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
  dplyr::group_by(species, predictor, predictor_val) %>%
  dplyr::summarise(meanEstimate=mean(fit), meanSE=mean(se), mean.p=mean(pvalue), sd.p=sd(pvalue), reps=n_distinct(rep), se.p=sd.p/sqrt(reps),
                   meanr2=mean(r2), sdr2=sd(r2), ser2=sdr2/sqrt(reps)) %>%
  dplyr::mutate(sig=ifelse(mean.p + 1.96*se.p<0.05, 1, 0)) 

lmRes_short2<-lmRes_sum3 %>%
dplyr::left_join(minAEE, by=c("species")) %>%
dplyr::left_join(maxAEE, by=c("species")) %>%
dplyr::left_join(minDev, by=c("species")) %>%
dplyr::left_join(maxDev, by=c("species")) %>%
dplyr::left_join(minCov, by=c("species")) %>%
dplyr::left_join(maxCov, by=c("species")) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::filter(meanEstimate >= minAEE - 1.96*seDEEmin & meanEstimate <= maxAEE + 1.96*seDEEmax) %>%
dplyr::filter(predictor_val >= minCov - 1.96*minCovSE & predictor_val <= maxCov + 1.96*maxCovSE)

# plot these
Figure2C<-ggplot() +
geom_pointrange(data=colonyRes_lox, aes(x=meanCov, y=meanDEEg, xmin=meanCov - 1.96*seCov, xmax=meanCov + 1.96*seCov, color=species), alpha=0.1) +
geom_pointrange(data=colonyRes_lox, aes(x=meanCov, y=meanDEEg, ymin=meanDEEg - 1.96*seDEE, ymax=meanDEEg + 1.96*seDEE, color=species), alpha=0.1) +
geom_line(data=filter(lmRes_short2, predictor=="totCov"), aes(y=meanEstimate, x=predictor_val, color=species, group=interaction(species, predictor))) +
geom_ribbon(data=filter(lmRes_short2, predictor=="totCov"), aes(x=predictor_val, y=meanEstimate, ymin=meanEstimate-meanSE*1.96, ymax=meanEstimate + meanSE*1.96, group=species, fill=species), alpha=0.2) +
scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
scale_fill_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822")) +
theme_bw() +
xlab(expression("WEE"[COV_NB]*"")) +
ylab(expression("TEE"[NB]*" (kJ.g"^-1*")")) +
#ylab("Sum of weekly deviance in energy expenditure") +
#xlab("NB energy expenditure (kJ.g)") +
theme_bw() +
ggtitle("C)") + 
#xlim(0.9, 8)+
#ylim(1000, 6300) +
theme(legend.position = "none",
panel.grid = element_blank(),
axis.text=element_text(size=12),
		axis.title=element_text(size=12))  
  
pdf("./results/figures/main/Figure2C.pdf", width=6, height=5)
plot(Figure2C)
dev.off()

### Figure 3A ####

Figure3A<-lmRes_tot2 %>%
dplyr::filter(test %in% c("DEE_vs_behavior5", "COV_vs_behavior4")) %>%
dplyr::filter(!predictors %in% c("(Intercept)", "typeWarm")) %>%
ungroup() %>%
#dplyr::left_join(warmCoefs, by=c("rep", "species", "test", "predictors")) %>%
#dplyr::mutate(Estimate=ifelse(predictors=="sst_diff_scale", Estimate + newEstimate, Estimate)) %>%
dplyr::mutate(upper=Estimate + 1.96*Std..Error, lower=Estimate - 1.96*Std..Error) %>%
dplyr::group_by(species, test, predictors) %>%
dplyr::summarise(repsTot=n_distinct(rep), mean.fit=mean(Estimate), sd.fit=sd(Estimate), se.fit=sd.fit/sqrt(reps), mean.se=mean(Std..Error), max.se=max(Std..Error), min.se=min(Std..Error), meanCI=1.96*mean.se, mean.r2=mean(r2), sd.r2=sd(r2), se.r2=1.96*sd.r2/sqrt(reps), mean.p=mean(Pr...t..), sdp=sd(Pr...t..), sep=sdp/sqrt(reps), minlower=min(lower), maxupper=max(upper), meanlower=mean(lower), meanupper=mean(upper)) %>% 
dplyr::mutate(sig=ifelse(mean.fit-1.96*mean.se < 0 & mean.fit + 1.96*mean.se > 0, 0, 1)) %>%
#dplyr::mutate(sig=ifelse(mean.fit-1.96*se.fit < 0 & mean.fit + 1.96*se.fit > 0, 0, 1)) %>%
#dplyr::mutate(sig=ifelse(meanlower < 0 & meanupper > 0, 0, 1)) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot"))) %>%
dplyr::mutate(test=ifelse(test=="DEE_vs_behavior5", "TEE_NB", "WEE_COV_NB")) %>%
dplyr::mutate(predictors=ifelse(predictors=="Migr_scale", "MigratoryDist", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="sst_diff_scale", "SST_Gain", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="sst_start_scale", "SST_Start", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="Migr_scale:sst_diff_scale", "MigratoryDist:SST_Gain", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="Forage_scale", "Forage", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="Flight_scale", "Flight", predictors)) %>%
#dplyr::mutate(predictors=ifelse(predictors=="Migr_scale", "MigratoryDist", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="sst_scale2", "SST", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="Forage_scale2", "Forage", predictors)) %>%
dplyr::mutate(predictors=ifelse(predictors=="Flight_scale2", "Flight", predictors)) %>%
dplyr::mutate(predictor=factor(predictors, levels=c("MigratoryDist", "SST_Start", "SST_Gain", "MigratoryDist:SST_Gain"))) %>%
ggplot(aes(x=mean.fit, y=predictors)) +
#geom_errorbarh(aes(xmin=mean.fit-1.96*se.fit, xmax=mean.fit + 1.96*se.fit, y=predictors, group=interaction(species, test), color=species, height = 0.2, alpha=factor(sig)), position=position_dodge(width=0.9))+
#geom_pointrange(aes(xmin=mean.fit-1.96*max.se, xmax=mean.fit + 1.96*max.se, group=interaction(species, test), color=species, alpha=factor(sig)), position=position_dodge2(width=0.9)) +
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
weights<-data.frame(species=c("Black-legged kittiwake", "Atlantic puffin","Common guillemot", "Northern fulmar", "Brünnich's guillemot" , "Little auk"), 
weight=c(386, 460, 1025, 728, 980, 149) , allometry=c(0.717, 0.689, 0.689, 0.765, 0.689, 0.689))

# Here we estimate average DEE to get an idea of whether calculations are any good
colonyRes_lox2<-colonyRes_lox %>%
ungroup() %>%
dplyr::left_join(weights, by=c("species")) %>%
dplyr::mutate(totalNBcosts=meanDEEg*weight^allometry, costs_DEE=totalNBcosts/n_distinct(dates_weekly$dateKeep), 
seDEE_tot=seDEE*weight^allometry, CI=seDEE_tot*1.96) %>%
dplyr::mutate(species=factor(species, levels=c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
                                                 "Little auk", "Common guillemot", "Brünnich's guillemot")))

print("Saving plots S8 & S12...")

FiguresS8<-ggplot() +
 geom_pointrange(data=colonyRes_lox2, aes(y=meanDEEg, x=colonyName, ymin=meanDEEg - 1.96*seDEE, ymax=meanDEEg + 1.96*seDEE, color=species)) +
 facet_wrap(~species, scales="free_y", nrow=3) +
 scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#C3A600", "#0072b2", "#E25822"))+
  theme_bw() +
  ylab(expression("TEE"[NB]*" (kJ)")) +
  xlab("") +
  geom_text(data=colonyRes_lox2, aes(x=colonyName, y=meanDEEg + 1.96*seDEE +  50, label=round(costs_DEE)), size=1.7) +
  theme_bw() +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("./results/figures/supplementary/FigureS8.pdf", width=15, height=9)
plot(FiguresS8)
dev.off()

FigureS12<-ggplot() +
 geom_pointrange(data=colonyRes_lox, aes(y=meanCov, x=colonyName, ymin=meanCov - 1.96*seCov, ymax=meanCov + 1.96*seCov, color=species)) +
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
dplyr::mutate(WEE_tot_scale=scale_01(meanDEEg), WEE_dev_scaled=scale_01(meanDev), WEE_cov_scaled=scale_01(meanCov)) %>%
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
#dplyr::filter(!species %in% c("Northern fulmar")) %>%
ungroup() %>%
dplyr::group_by(colony, species) %>%
dplyr::mutate(weightedMean1=meanDEEg*(1/meanBirds), weightedMean2=meanDev*(1/meanBirds), weightedMean3=meanCov*(1/meanBirds)) %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::mutate(mean1_scale=scale_01(weightedMean1), mean3_scale=scale_01(weightedMean3)) %>%
dplyr::group_by(colony) %>%
dplyr::summarise(coords.x1=mean(coords.x1), coords.x2=mean(coords.x2), WEE_tot_scale=mean(meanDEEg), WEE_dev_scale=mean(meanDev), 
WEE_cov_scale=mean(meanCov), mean1=sum(weightedMean1)/sum(1/meanBirds), mean2=sum(weightedMean2)/sum(1/meanBirds), mean3=sum(weightedMean3)/sum(1/meanBirds), mean1_2=mean(mean1_scale), mean3_2=mean(mean3_scale)) 
							
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
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=meanDEEg), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(expression("TEE"[NB]*" (kJ.g"^-1*")"), colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026"), breaks = c(min(energyplot$meanDEEg), max(energyplot$meanDEEg))) +
  #scale_fill_gradientn(expression("WEE"[totalNB]*" (kJ.g"^-1*")"), colors=viridis_colors) +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  #labs(shape=expression("WEE"[devianceNB]*"")) +
  theme(legend.position="bottom") +
  labs(tag="A)")
  
  Figure4b_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=meanCov), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(name=expression("WEE"[COV_NB]*""), colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  #labs(shape=expression("WEE"[devianceNB]*"")) +
  theme(legend.position="bottom") +
  labs(tag="B)")
  
  Figure4c_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=meanSSTcol), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(name=expression("SST"[Start]*" (°C)"), colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  #labs(shape=expression("WEE"[devianceNB]*"")) +
  theme(legend.position="bottom") +
  labs(tag="C)")
  
  Figure4d_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=meanSSTdiff), color="black", shape=21, cex=3) + 
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
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=meanMigratoryDist), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(name="MigratoryDist (km)", colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  #labs(shape=expression("WEE"[devianceNB]*"")) +
  theme(legend.position="bottom") +
  labs(tag="E)")
  
  Figure4f_sup<-ggplot() +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=energyplot, aes(x=coords.x1, y=coords.x2, fill=meanSSTsd), color="black", shape=21, cex=3) + 
  scale_fill_gradientn(name="SST sd (°C)", colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #scale_color_manual(values=c("#875692", "#BE0032", "#008856", "#F3CC00", "#0072b2", "#E25822"))  +
  coord_sf(crs=projection_NA, xlim=c(min(coloniesSum$coords.x1-100000),max(coloniesSum$coords.x1+100000)) , ylim=c(min(coloniesSum$coords.x2-100000),max(coloniesSum$coords.x2+100000))) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  #labs(shape=expression("WEE"[devianceNB]*"")) +
  theme(legend.position="bottom") +
  labs(tag="F)")
  
  
pdf(paste0("./results/figures/supplementary/energyStrategyPlots_", speciesList[i], ".pdf"), width=9)
grid.arrange(Figure4a_sup, Figure4b_sup, Figure4c_sup, Figure4d_sup, Figure4e_sup, nrow=2)
dev.off()

}


# Figure S13: Distribution of raw data -> WEE vs TEE #

FigureS13<-ggplot() +
 geom_pointrange(data=colonyRes_lox, aes(x=meanCov, y=meanDEEg, xmin=meanCov - 1.96*seCov, xmax=meanCov + 1.96*seCov, color=species), alpha=0.1) +
 geom_pointrange(data=colonyRes_lox, aes(x=meanCov, y=meanDEEg, ymin=meanDEEg - 1.96*seDEE, ymax=meanDEEg + 1.96*seDEE, color=species), alpha=0.1) +
 geom_text_repel(data=colonyRes_lox, aes(x=meanCov, y=meanDEEg, label=colonyName), max.overlaps=30) + 
  geom_line(data=filter(lmRes_short2, predictor=="totCov"), aes(y=meanEstimate, x=predictor_val, color=species, group=interaction(species, predictor))) +
  geom_ribbon(data=filter(lmRes_short2, predictor=="totCov"), aes(x=predictor_val, y=meanEstimate, ymin=meanEstimate-meanSE, ymax=meanEstimate + meanSE, group=species, fill=species), alpha=0.2) +
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
dplyr::summarise(meanSST=mean(SST), meanSST_dev=mean(covSSTWeekly),
meanFlight=mean(totFlightHrs), meanFlight_dev=mean(totDevianceFlight),
meanForage=mean(totForageHrs), meanForage_dev=mean(totDevianceForage),
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
