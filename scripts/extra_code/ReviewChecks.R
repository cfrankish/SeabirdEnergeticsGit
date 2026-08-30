# Checks for MS review #

### How many gaps is the IRMA dataset filling in? ###

# open up all birds & IDS used in the analysis #
birds<-read.csv("./results/tables/main/table1_idcatalogue.csv")

# Find ids #
allResults_deviation<-list.files("./results/tables/main/", full.names=TRUE)
deviation_weekly<-allResults_deviation[grepl("/weeklydeviance_", allResults_deviation)]

# Make a species-list
speciesList<-unique(birds$species)

# determine start & end date of study period #
startDate<-"2021-09-15" # Read-in start of study period
endDate<-"2022-04-15" # Read-in end date of study period

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

# Find matching IRMA dataset #
irma.files<-list.files("./data/positionsIRMA/", full.names=TRUE)
metadata<-load(irma.files[8])
irma.files<-irma.files[grepl("IRMAlocs", irma.files)]
irma.files.df<-data.frame(irma.files)
colnames(irma.files.df)<-c("FileName")
irma.files.df$species<-c("Little auk", "Atlantic puffin", "Northern fulmar", "Black-legged kittiwake", "Common guillemot", "Brünnich's guillemot")

propFillSumAll<-list()
timeFillSumAll<-list()

for (i in 1:length(speciesList)) {

print(speciesList[i])

# subset to species i #
speciesSub<-speciesList[i]

# open corresponding IRMA file
irma.files.sub<-subset(irma.files.df, species==speciesSub)
irmaSpecies<-readRDS(irma.files.sub$FileName[1])

# Add month & day & track_year
irmaSpecies$month<-as.numeric(substr(irmaSpecies$timestamp, 6, 7))
irmaSpecies$day<-as.numeric(substr(irmaSpecies$timestamp, 9, 10))
irmaSpecies$individ_id_og<-irmaSpecies$individ_id
irmaSpecies$individ_id<-gsub("-", "_", irmaSpecies$individ_id_og)
irmaSpecies$individ_id<-gsub("BLUE", "blue", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("PUFFIN", "puffin", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("6T_0405", "6T-0405", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("6T_", "6T-", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("MK14", "mk14", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("HORNØYA", "Hornøya", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("CAPEKRUTIK", "CapeKrutik", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("GBT_EA35285", "GBT_C2513_2018_IoM_BK", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("KAPH", "KapH", irmaSpecies$individ_id)
irmaSpecies$individ_id<-gsub("4F_3042", "4F-3042", irmaSpecies$individ_id)

# filter to relevant Ids & track-years? #
birdsSub<-readRDS(deviation_weekly[grepl(speciesSub, deviation_weekly)]) # 605 individuals
birdIdsMatch<-birdsSub %>%
dplyr::group_by(individ_id, track_year) %>%
dplyr::slice(1) %>%
dplyr::select(individ_id, track_year)
irmaSpeciesSub<-irmaSpecies %>%
dplyr::filter(individ_id %in% unique(birdsSub$individ_id))

# Determine gap stats #
irmaSpeciesFilt<-irmaSpeciesSub %>%
dplyr::ungroup() %>%
dplyr::mutate(date=substr(timestamp, 1, 10)) %>%
dplyr::inner_join(dates_weekly, by=c("month", "day")) 

# Add track-year #
irmaSpeciesFilt$track_year<-ifelse(irmaSpeciesFilt$month %in% c(1, 2, 3, 4), paste0(as.numeric(substr(irmaSpeciesFilt$timestamp, 1, 4))-1, "_", substr(irmaSpeciesFilt$timestamp, 3, 4)), 
paste0(substr(irmaSpeciesFilt$timestamp, 1, 4), "_", substr(as.numeric(substr(irmaSpeciesFilt$timestamp, 1, 4)) + 1, 3, 4)))

# Summarize gaps
propFill<-irmaSpeciesFilt %>%
dplyr::mutate(track_year=as.character(track_year)) %>%
dplyr::inner_join(birdIdsMatch, by=c("individ_id", "track_year")) %>%
dplyr::ungroup() %>%
dplyr::group_by(individ_id, track_year) %>%
dplyr::count(loc_type) %>%
dplyr::mutate(totLoc=sum(n)) %>%
dplyr::mutate(locProp=n/totLoc) %>%
dplyr::filter(loc_type=="IRMA")

# Summarize duration of gaps
timeFill <- irmaSpeciesFilt %>%
dplyr::inner_join(birdIdsMatch, by=c("individ_id", "track_year")) %>%
dplyr::ungroup() %>%
#dplyr::group_by(individ_id, track_year) 
  arrange(individ_id, timestamp) %>%               # drop id if you don't have it
  group_by(individ_id, track_year) %>%                         # drop if no id
  mutate(
    is_irma = loc_type == "IRMA",
    start_flag = is_irma & !lag(is_irma, default = FALSE),
    end_flag   = !is_irma &  lag(is_irma, default = FALSE),
    interval_id = cumsum(start_flag)
  ) %>%
  # keep rows that matter: starts (IRMA) and ends (first non-IRMA after run)
  filter(start_flag | end_flag) %>%
  group_by(individ_id, track_year, interval_id) %>%
  summarise(
    start_time = timestamp[start_flag][1],
    end_time   = timestamp[end_flag][1],         # NA if run ends at dataset end
    duration   = difftime(end_time, start_time, units = "days"),
    .groups = "drop"
  )

# Calculate mean in days  
timeFill_mean<-timeFill %>%
dplyr::mutate(species=speciesSub) %>%
dplyr::ungroup() %>%
dplyr::group_by(species) %>%
dplyr::mutate(duration=as.numeric(duration)) %>%
dplyr::summarise(meanGap=mean(duration, na.rm=TRUE), sdGap=sd(duration, na.rm=TRUE)) 

# Check number of birds matching
print(paste0("total birds dataset:" ,n_distinct(birdIdsMatch$individ_id)))
print(paste0("total birds IRMA:" ,n_distinct(propFill$individ_id)))

# print unique lox total to make sure nothing is weird #
print(paste0("Total locations:", unique(propFill$totLoc)))

# Summarise mean +/- sd
propFillSum<-propFill %>%
ungroup() %>%
dplyr::mutate(species=speciesSub) %>%
dplyr::group_by(species) %>%
dplyr::summarise(meanFill=mean(locProp), sdFill=sd(locProp))

# Save result #
propFillSumAll<-rbind(propFillSumAll, propFillSum)
timeFillSumAll<-rbind(timeFillSumAll, timeFill_mean)

}

#### How much time below LCT? or above? #### + degrees below

lcts<-data.frame(species=c("Black-legged kittiwake", "Northern fulmar", "Little auk", "Common guillemot", "Brünnich's guillemots", "Atlantic puffin"),
LCT=c(12.5, 9, 14.18, 14.18, 14.18, 14.18))

# Open each individual file cry
idCatalogue<-read.csv("./results/tables/main/table1_idcatalogue.csv")

# Determine where daily files are
allResults<-list.files("./tmp/", full.names=TRUE)
energyRes_day<-allResults[grepl("energyDay", allResults)]

# Determine study period
startDate<-"2021-09-15" # Read-in start of study period
endDate<-"2022-04-15" # Read-in end date of study period

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

# Make a species list to iterate through #
speciesUnique<-unique(idCatalogue$species)

# Make a list to save results
species_lct_stats<-list()

for (i in 1:length(speciesUnique)) {

# subset to species i
speciesSub<-speciesUnique[i]

# make a list of ids to iterate through
idSub<-subset(idCatalogue, species==speciesSub)
idList<-unique(idSub$individ_id)

# list to save results
bird_lct_stats<-list()

for (j in 1:length(idList)) {

# Subset bird j
    birdSub<-idList[j]
	birdSub<-gsub("-", "_", birdSub)
	birdSub<-gsub("ø", "o", birdSub)
    
    # Find csv file
    print(paste0("Opening file..Species ", i, " Bird ", j, "/", length(idList)))
    birdSub<-fread(energyRes_day[grep(birdSub, energyRes_day)])
	
	# Determine month and year
	birdSub$day<-as.numeric(substr(birdSub$date, 9, 10))
	birdSub$month<-as.numeric(substr(birdSub$date, 6, 7))
	birdSub$year<-as.numeric(substr(birdSub$date, 1, 4))
	birdSub$track_year<-ifelse(birdSub$month<day1_month, paste0(birdSub$year-1, "_", substr(birdSub$year, 3, 4)), paste0(birdSub$year, "_", substr(birdSub$year + 1, 3, 4)))
    birdSub$session_year<-paste0(birdSub$session_id, "_", birdSub$track_year)

    # Determine session_year combos with enough data	
    dates_sessions<-birdSub%>%
    dplyr::mutate(date=substr(date, 1, 10)) %>%
    dplyr::group_by(session_id, track_year) %>%
    dplyr::count(date) %>%
    dplyr::mutate(month=as.numeric(substr(date, 6, 7))) %>%
    dplyr::mutate(day=as.numeric(substr(date, 9, 10))) %>%
    dplyr::inner_join(dates_weekly2, by=c("month", "day")) %>%
    ungroup() %>%
    dplyr::group_by(session_id, track_year) %>%
    dplyr::summarise(daysTot=n_distinct(date)) %>%
    dplyr::filter(daysTot==nrow(dates_weekly2)) %>%
    dplyr::mutate(session_year=paste0(session_id, "_", track_year))
    
	if (nrow(dates_sessions) < 1) {
      next}
	
    # Subset to one location per date to reduce file size
    birdcsv_reducted<-birdSub %>%
      ungroup() %>%
	  dplyr::filter(session_year %in% c(dates_sessions$session_year)) %>%
      dplyr::inner_join(dates_weekly2, by=c("day", "month")) %>%
	  dplyr::filter(weekNo>0) %>%
	  dplyr::select(rep, species, colony, individ_id, track_year, date, session_year, mean.lon, mean.lat, sst_random) 

   # Attach LCT
   bird_lct<-birdcsv_reducted %>%
   dplyr::left_join(lcts, by=c("species")) %>%
   dplyr::mutate(below_lct=ifelse(sst_random<LCT, 1, 0)) %>%
   dplyr::mutate(degrees_below_lct=LCT-sst_random) %>%
   dplyr::group_by(rep, species, colony, individ_id, track_year) %>%
   dplyr::summarise(propBelow=sum(below_lct)/n_distinct(date), degrees_below_mean=mean(degrees_below_lct), degrees_below_min=min(degrees_below_lct), degrees_below_max=max(degrees_below_lct)) %>%
   ungroup() %>%
   dplyr::group_by(species, colony, individ_id) %>%
   dplyr::summarise(meanpropBelow=mean(propBelow), minpropBelow=min(propBelow), maxpropBelow=max(propBelow), mean_degrees_below_lct=mean(degrees_below_mean), min_degrees_below_lct=min(degrees_below_min), max_degrees_below_lct=max(degrees_below_max))

 bird_lct_stats<-rbind(bird_lct_stats, bird_lct)
 
 }
 
 species_lct_stats<-rbind(species_lct_stats, bird_lct_stats)
 
 }
 
### Check breeding colony range (in latitude per species) ###

# open up all birds & IDS used in the analysis #
birds<-read.csv("./results/tables/main/table1_idcatalogue.csv")

birds_colonies<-birds %>%
ungroup() %>%
dplyr::group_by(species) %>%
dplyr::count(colony)

# Find colony latitude & longitude #
irma.files<-list.files("./data/positionsIRMA/", full.names=TRUE)
metadata<-load(irma.files[8])

# Attach information #
birds_colonies_lox<-birds_colonies %>%
dplyr::left_join(colony.summary, by=c("colony")) %>%
dplyr::ungroup() %>%
dplyr::filter(!is.na(col_lat)) %>%
dplyr::group_by(species) %>%
dplyr::arrange(col_lat) %>%
dplyr::slice(c(1, n()))

# Attach seabird body mass #
seabirdBodyMass<-data.frame(species=c("Black-legged kittiwake", "Northern fulmar", "Common guillemot",
"Little auk", "Atlantic puffin", "Brünnich's guillemot"), mass=c(392, 728, 940, 149, 395, 980))

# Calculate range of LCT-air #
birds_colonies_lox %>%
dplyr::left_join(seabirdBodyMass, by=c("species")) %>%
dplyr::mutate(LCTair=43.15-6.58*log10(mass)-0.26*col_lat) %>%
dplyr::mutate(LCTwater_min=LCTair + 5, LCTwater_max=LCTair + 14)

