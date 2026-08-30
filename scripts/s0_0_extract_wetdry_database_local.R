### This code is for users of the SEATRACK database (as it requires access to the raw data) ###
### Here is just some code to extract raw wet-dry data & individual information ### 

library(seatrackR)

### S1_0_1: extract wet-dry data by species ###

# Determine species list to loop through #
speciesList<-c("Black-legged kittiwake", "Northern fulmar", "Atlantic puffin",
               "Common guillemot", "Brünnich's guillemot", "Little auk")

# Determine data of today
today <- Sys.Date()

# Connect to database
connectSeatrack(Username = "caitlin_frankish", Password = "Alca_torda",
                host = "seatrack.nina.no")

# Loop through wet-dry data and individual information data #
for (i in 1:length(speciesList)) {

# print an update message to follow status 
print(paste0("Extracting data for species ", i, "/", length(speciesList)))
    
# Subset to species i 
speciesName<-speciesList[i]

# Download wet-dry data from species of interest				
wetdry_recording<-getRecordings(type="activity", species=speciesName)

# Save this information
saveRDS(wetdry_recording, file=paste0("data/wetdry_database/wetdry_database_", speciesName, "_", today, ".rds"))
  
}

# Download individual metadata from all species
info <-getIndividInfo()

# Save this information
saveRDS(info, file=paste0("data/wetdry_database/metadata.rds"))
