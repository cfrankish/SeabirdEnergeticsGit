library(heatwaveR)
library(dplyr)
library(lubridate)
library(dplyr)
library(ncdf4)
library(raster)
library(sp)
library(sf)
library(lubridate)
library(data.table)
library(tidyr)
library(ggplot2)
library(rnaturalearthdata)
library(rnaturalearth)
library(dbscan)
library(terra)
library(spatialEco)

scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

### Step 0: Open December energy Map ###

# Set-up some projections
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # Flat projection
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"

# Set-up number of iterations...
overall.iterations<-5 # how many times this is calculated per individual
print(paste0("Determining activity distributions for ", overall.iterations, " iterations..."))

# Source all necessary functions
source("./scripts/functions.R")

# Open-up species-specific monthly energy expenditure #
energyFiles<-list.files("/cluster/projects/nn11080k/cfrank93/cbirdEnergy/results/tables", full.names=TRUE)
energySpeciesAll<-energyFiles[grepl(paste0("allenergy_v1"), energyFiles)]
energySpeciesDf<-fread(energySpeciesAll[1])
energyDec<-subset(energySpeciesDf, month=="12")
energyDecPlot<-energyDec %>%
dplyr::select(x, y, energyPopkJ_mean) %>%
dplyr::filter(energyPopkJ_mean>0) %>%
dplyr::mutate(scaled=scale01(energyPopkJ_mean)) %>%
dplyr::select(x, y, scaled) %>%
dplyr::mutate(metric="EnergyDemand")

# Turn into a raster #
energyPrep<-energyDec %>%
dplyr::select(x, y, energyPopkJ_mean) %>%
rename(z=energyPopkJ_mean)

energyRast<-rasterFromXYZ(energyPrep)
proj4string(energyRast)<-projection_84

### Step 1: Open Chl map / nutrient map ####

nutrientFiles<-list.files("data/chlNutrients/", full.names=TRUE)

# Stack monthly files #

chlList<-list()
nutrientList<-list()

for (i in 1:length(nutrientFiles)) {

chlSub<-raster(nutrientFiles[i], varname="chl")
nutrSub<-raster(nutrientFiles[i], varname="po4")

chlList<-append(chlList, chlSub)
nutrientList<-append(nutrientList, nutrSub)

}

# Calculate a mean for each variable #
chlMean<-overlay(stack(chlList), fun="mean")
pMean<-overlay(stack(nutrientList), fun="mean")

# re-sample so they match the other raster
chlMatch<-projectRaster(chlMean, energyRast)
pMatch<-projectRaster(pMean, energyRast)

# Mask 0 values in energy demand #
mask_r <- energyRast
mask_r[mask_r == 0] <- NA
r1_masked <- mask(chlMatch, mask_r)
r2_masked <- mask(pMatch, mask_r)

# Turn into data frames #
chlDf<-as.data.frame(r1_masked, xy=TRUE)
chlPlot<-chlDf %>%
dplyr::filter(layer>0) %>%
mutate(scaled=scale01(layer)) %>%
dplyr::select(x, y, scaled) %>%
dplyr::mutate(metric="Chla")

poDf<-as.data.frame(r2_masked, xy=TRUE)
poPlot<-poDf %>%
dplyr::filter(layer>0) %>%
mutate(scaled=scale01(layer)) %>%
dplyr::select(x, y, scaled) %>%
dplyr::mutate(metric="po4")

# Make some comparative plots

# Create a world map for plotting
coast <- ne_coastline(scale = "small", returnclass = "sf")
world <- ne_countries(scale = "small", returnclass = "sf")

# Set up projection values
projection_NA<-"+proj=laea +x_0=0 +y_0=0 +lon_0=-9 +lat_0=61"
projection_84<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# make a plot for viewing

variablesPlot<-rbind(energyDecPlot, poPlot)

examplePlots<-ggplot() +
geom_tile(data=variablesPlot, aes(x=x, y=y, fill=scaled)) +
 geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_84, xlim=c(min(energyDecPlot$x), max(energyDecPlot$x)), ylim=c(min(energyDecPlot$y), max(energyDecPlot$y))) +
  scale_fill_gradientn(
  colours = c(
    "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"
  )) +
 theme_minimal() +
 ggtitle("Top 10% of data") +
 labs(fill="") +
 facet_wrap(~ metric, nrow=1)
 
pdf("examplePlots.pdf", height=8, width=10)
plot(examplePlots)
dev.off()

### Step 2: analyse storm data ####

####Classification of winter cyclones between 2000-2016####
#Preparation of the dataset: Cyclone locations were extracted from the "Northern Hemisphere Cyclone Locations and Characteristics from NCEP/NCAR Reanalysis Data, Version 1" online dataset 
t5808<-read.table("data/storm/ncepstorms_1958_2008.txt",header=TRUE)
t0008<-subset(t5808,t5808$year<=8) #we select cyclones from 2000

t2009<-read.table("data/storm/ncepstorms_2009.txt",header=FALSE)
colnames(t2009)<-colnames(t5808)

t2010<-read.table("data/storm/ncepstorms_2010.txt",header=FALSE)
colnames(t2010)<-colnames(t5808)

t2011<-read.table("data/storm/ncepstorms_2011.txt",header=FALSE)
colnames(t2011)<-colnames(t5808)

t2012<-read.table("data/storm/ncepstorms_2012.txt",header=FALSE)
colnames(t2012)<-colnames(t5808)

t2013<-read.table("data/storm/ncepstorms_2013.txt",header=FALSE)
colnames(t2013)<-colnames(t5808)

t2014<-read.table("data/storm/ncepstorms_2014.txt",header=FALSE)
colnames(t2014)<-colnames(t5808)

t2015<-read.table("data/storm/ncepstorms_2015.txt",header=FALSE)
colnames(t2015)<-colnames(t5808)

t2016<-read.table("data/storm/ncepstorms_2016.txt",header=FALSE)
colnames(t2016)<-colnames(t5808)

#Merging and saving
total<-rbind(t0008,t2009,t2010,t2011,t2012,t2013,t2014,t2015,t2016)
#setwd("/Volumes/")
#write.csv(total,"cyc2000_2016.csv")

#Selection of winter cyclones (between October and February)
tothiv<-subset(total,total$month>=10|total$month<=2)

#Resetting longitudes between -180 and 180 instead of 0 and 360
for (i in 1:nrow(tothiv)){
  if(tothiv$lon[i]<180){
    tothiv$longbis[i]<-tothiv$lon[i]
  }else{
    tothiv$longbis[i]<-tothiv$lon[i]-360
  }
}
write.csv(tothiv,"cyc2000_2016_hiv.csv")
#Cropping on our studied area
tothivatl<-subset(tothiv,tothiv$longbis>=-100&tothiv$longbis<=100)

#Classification by intensity using Dvorak's scale: see STAR METHOD
#Class-1: air pressure >1009 hPa
#Class-2: 1009 hPa> air pressure>1005 hPa
#Class-3: 1005 hPa> air pressure>987 hPa
#Class-4: air pressure< 987 hPa

pt<-subset(tothivatl,tothivatl$cent_pressure>1009)
pt$classe<-1

dt<-subset(tothivatl,tothivatl$cent_pressure>1005 & tothivatl$cent_pressure<=1009)
dt$classe<-2

tt<-subset(tothivatl,tothivatl$cent_pressure>987 & tothivatl$cent_pressure<=1005)
tt$classe<-3

h1<-subset(tothivatl,tothivatl$cent_pressure<=987)
h1$classe<-4

fin<-rbind(pt,dt,tt,h1)
write.table(fin,"CycloneNAClasse.csv",sep=";",row.names=FALSE)

####Calculation of the average number of cyclone per pixel (250km) of the studied area (100°W-100°E, 0°N-90°N), for each month and cyclone category####
terre<-readOGR("/Volumes/terre84.shp")# shapefile contening polygones of each states. Projection WSG84
redc<-readOGR("/Volumes/redcotepoly.shp")#Shapefile delineating the continent
redcp<- spTransform(redc, CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))# re-projection with the North Pole Lambert Azimutal Equal Area.

#The position of the cyclones corresponds to the centre of the pixel where it was detected: these pixels are derived from a 250km*250km EASE grid. We therefore download the grid to use it.
Temp<-nc_open("/Volumes/EASE2_N25km.geolocation.v0.9.nc")
lon <- ncvar_get(Temp, "x")
lat <- ncvar_get(Temp, "y")
lonlat <- as.matrix(expand.grid(lon,lat))
Tempdf <- data.frame(lonlat)
names(Tempdf) <- c("lon","lat")
Tempdf$val<-0
rastease<-rasterFromXYZ(Tempdf,crs="+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs")
rastease250<-aggregate(rastease,10) #re scaling (250km*250km)

#Cropping on our studied area
rasteasecr<-crop(rastease250,extent(c(-6340823,6340823,-6362500,0)))

#Calculation of the average number of cyclone per pixel (250km) of the studied area (100°W-100°E, 0°N-90°N) for each month and each cyclone category
for (i in c(1,2,10:12)){ #for each month
  mois<-subset(fin,fin$month==i)
  for(c in 1:4){ # for each cyclone category
    cl<-subset(mois,mois$classe==c)
    for (y in 0:16){ # for each year
      annee<-subset(cl,cl$year==y)
      if(nrow(annee)!=0){
        un<-unique(annee[c("lat", "longbis","sys_num")])
        sp<-un
        sp$count<-1
        coordinates(sp) <- ~longbis+lat
        proj4string(sp)<-proj4string(terre)
        #conversion on the Ease grid (250km)
        spease<-spTransform(sp,CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))
        raster<-rasterize(spease,rasteasecr,spease$count,fun="count",background=0)
        #saving
        writeRaster(raster,paste("NBcyclone_classe",c,"_annee",y,"_mois",i,".tif"),overwrite=TRUE)
      }else{
        raster<-rasteasecr
        raster[]<-0
        writeRaster(raster,paste("NBcyclone_classe",c,"_annee",y,"_mois",i,".tif"),overwrite=TRUE)
      }
      #calculation of the average number of cyclone for the "cl" category and the "i" month
      if (y==0){
        st<-raster
      }else{
        st<-stack(st,raster) 
      }
      
    }
    moy<-mean(st)
    
    #Suppression of terrestrial pixel 
    moyterre<-mask(moy,redcp,inverse=TRUE)
    writeRaster(moyterre,paste("NBmoycyclone_classe",c,"_mois",i,".tif"),overwrite=TRUE) 
  } 
}

#Graphical representation of the average number of winter cyclones in the North Atlantic Ocean between 2000 and 2016 for each month and category (Figure S2)

nb<- data.frame(matrix(NA,ncol=4,nrow=1))
colnames(nb)<-c("Classe","Mois","Nombremoyen","SD") # Class, month, average number, standard devidation

moy<- data.frame(matrix(NA,ncol=1,nrow=1))
colnames(moy)<-c("Nombremoyen")

for(i in c(1,2,10:12)){
  for (j in 1:4){
    nb$Classe<-j
    nb$Mois<-i
    for (y in 0:16){
      mois<-subset(fin,fin$month==i)# selection of cyclones for the month "i" 
      cl<-subset(mois,mois$classe==j)# selection of category "j" cyclones 
      yr<-subset(cl,cl$year==y)
      # Count
      un<-unique(yr[c("sys_num")])
      moy$Nombremoyen<-nrow(un)
      if(y==0){
        moyannee<-moy  
      }else{
        moyannee<-rbind(moyannee,moy)  
      }
    }
    nb$Nombremoyen<-mean(moyannee$Nombremoyen)
    nb$SD<-sd(moyannee$Nombremoyen)
    
    if(j==1){
      moisfin<-nb  
    }else{
      moisfin<-rbind(moisfin,nb)  
    }
  }
  
  if(i==1){
    nbfinal<-moisfin  
  }else{
    nbfinal<-rbind(nbfinal,moisfin)  
  }
}
# Saving
write.table(nbfinal,"NombremoyencycloneAtlantiqueOcean.csv",sep=";",row.names=FALSE)

nbfinal$Month2 <- factor(nbfinal$Mois, levels=unique(nbfinal$Mois))
nbfinal$Month2 <- factor(nbfinal$Month2,levels(nbfinal$Month2)[c(3,4,5,1,2)])
nbfinal$Classe <- as.factor(nbfinal$Classe)

p<-ggplot(nbfinal, aes(x=Month2, y=Nombremoyen,color=Classe,group=Classe))+
  scale_colour_manual(name=c("Cyclones intensities"),values=c("#99CCFF","#3399CC","#0066CC","#003366"),labels=c("Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("1","2","3","4"))+ 
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Nombremoyen-SD, ymax=Nombremoyen+SD),colour="black", width=.1,alpha=0.4,position=position_dodge(.05))+
  xlab("Month") +
  ylab("Mean number of cyclones in the North Atlantic Ocean 
       between 2000 and 2016\n") +
  scale_x_discrete(labels=c("Oct","Nov","Dec","Jan","Feb"))+
  theme_classic() 


ggsave(filename ="NbmoyencycloneAtl.jpeg",width = 20, height = 20, units = "cm") 

# We also calculated the mean number of cyclone for each wintering month, considering all cyclone categories, between 2000 and 2016
maxv <- 4
minv <- 0
cols <-c("#FFFFFF",brewer.pal(n =9, name = "Reds"),"#330000")
pal <- colorRampPalette(cols)
brks <- seq(minv,maxv,by=0.2)
nbrks <- length(brks)-1


for (i in c(1,2,10:12)){
  cl<-subset(fin,fin$month==i)
  for(j in 0:16){
    annee<-subset(cl,cl$year==j)
    un<-unique(annee[c("lat", "longbis","sys_num")])
    sp<-un
    sp$count<-1
    coordinates(sp) <- ~longbis+lat
    proj4string(sp)<-proj4string(terre)
    spease<-spTransform(sp,CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))
    raster<-rasterize(spease,rasteasecr,spease$count,fun="count",background=0)
    
    if (j==0){
      st<-raster
    }else{
      st<-stack(st,raster) 
    }
    
  }
  moy<-mean(st)# average calculation
  #we deleted terrestrial pixels
  moyterre<-mask(moy,redcp,inverse=TRUE)
  writeRaster(moyterre,paste("HIVNBmoycyclone_mois",i,".tif"),overwrite=TRUE)
  
  #Graphical representation
  t<-paste("HivNBmoycylcone_mois",i,".jpeg",sep="")
  jpeg(filename =t, width = 3200, height = 2800, units = "px", res=300)
  plot(moyterre,breaks=brks,col=pal(nbrks),legend=T,useRaster=FALSE,box=FALSE,axes=F,main=paste("mois_",i,sep=""))
  plot(redcp,col="lightgrey",border=c("#999999"),add=T)
  plot(gratp,col="#FFFFFF",border=c("#FFFFFF"),cex=0.2,add=T)
  dev.off() 
  
}

####Focusing on the 1000km*1000km area off Newfounland: calculation of the total number of cyclones in the area between 2000-2016 for each cyclone category and month####
#Defining the area
rastdf<-as.data.frame(rasteasecr,xy=T)
terren<-subset(rastdf,rastdf$x>=-3625000 & rastdf$x<=-2875000& rastdf$y>=-3375000 & rastdf$y<=-2625000)
coordinates(terren) <- ~x+y
terreneuverast<-crop(rasteasecr,extent(c(-3675000,-2825000,-3375000,-2625000)))

#calculation of the total number of cyclones in the area between 2000-2016 for each category and month
nb<- data.frame(matrix(NA,ncol=3,nrow=1))
colnames(nb)<-c("Classe","Mois","Nombre") #Classn month, Number

for(i in c(1,2,10:12)){
  for (j in 1:4){
    nb$Classe<-j
    nb$Mois<-i
    
    mois<-subset(fin,fin$month==i)# selection of cyclones for the month "i"
    cl<-subset(mois,mois$classe==j)# selection of category "j" cyclones 
    coordinates(cl) <- ~longbis+lat
    proj4string(cl)<-proj4string(terre)
    cl<-spTransform(cl,CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))#changing the projection to North Pole Lambert Azimutal Equal Area
    cldt<-as.data.frame(cl)
    zone<-subset(cldt,cldt$longbis>=-3675000 & cl$longbis<=-2825000& cl$lat>=-3375000 & cl$lat<=-2625000)# selection of cyclone in the area off Newfoundland
    # Count
    un<-unique(zone[c("sys_num","year")])
    nb$Nombre<-nrow(un)
    
    if(j==1){
      moisfin<-nb  
    }else{
      moisfin<-rbind(moisfin,nb)  
    }
  }
  
  if(i==1){
    nbfinal<-moisfin  
  }else{
    nbfinal<-rbind(nbfinal,moisfin)  
  }
}
# Saving
write.table(nbfinal,"NombrecycloneTerreNeuve.csv",sep=";",row.names=FALSE)
#Graphical representation (see Figure S3)
for(i in 1:nrow(nbfinal)){
  if (nbfinal$Mois[i]== 1){
    nbfinal$ML[i]<-"Jan"
  }  else if (nbfinal$Mois[i]== 2){
    nbfinal$ML[i]<-"Feb"
  }  else if (nbfinal$Mois[i]== 10){
    nbfinal$ML[i]<-"Oct"
  }else if (nbfinal$Mois[i]== 11){
    nbfinal$ML[i]<-"Nov"
  }else{
    nbfinal$ML[i]<-"Dec"
  }
}
nbfinal$ML<-factor(nbfinal$ML)

nbfinal$Month2 <- factor(nbfinal$ML, levels=unique(nbfinal$ML))
nbfinal$Month2 <- factor(nbfinal$Month2,levels(nbfinal$Month2)[c(3,4,5,1,2)])
nbfinal$Classe <- as.factor(nbfinal$Classe)

p<-ggplot(nbfinal, aes(x=Month2, y=Nombre,color=Classe,group=Classe))+
  scale_colour_manual(name=c("Cyclones intensities"),values=c("#99CCFF","#3399CC","#0066CC","#003366"),labels=c("Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("1","2","3","4"))+ 
  geom_line(size=1)+
  geom_point(size=3)+
  xlab("Month") +
  ylab("Average number of cyclone\n")+ 
  scale_x_discrete(labels=c("Oct","Nov","Dec","Jan","Feb"))+
  theme_classic() 


ggsave(filename ="NbcycloneTerreNeuveFIGmatsup4.jpeg",width = 20, height = 20, units = "cm")

####Extraction of environmental characteristics (cyclonic and non-cyclonic) in this 1000km*1000km area off Newfoundland####
#Extraction of environmental characteristics for each cyclones within the area of Newfoundland
#Sea Surface Temperature: NOAA/OAR/ESRL PSL see STAR Method

for (i in min(fin$year): max(fin$year)){
  g<-i+2000
  Temp<-nc_open(paste("/Volumes/sst.day.mean.",g,".nc",sep="")) 
  lon <- ncvar_get(Temp, "lon")
  lat <- ncvar_get(Temp, "lat")
  time<-ncvar_get(Temp,"time")
  tunits <- ncatt_get(Temp,"time","units")
  tas_time <- nc.get.time.series(Temp)
  tas_time<-as.character(tas_time)
  
  yhiv<-subset(fin,fin$year==i) #selection of the year "i"
  for (j in c(1,2,10:12)){
    mhiv<-subset(yhiv,yhiv$month==j) # selection of the year "j"
    for (d in min(mhiv$day):max(mhiv$day)){
      dhiv<-subset(mhiv,mhiv$day==d) 
      if(nrow(dhiv)!=0) {
        mhivspace<-dhiv
        coordinates(mhivspace) <- ~longbis+lat
        proj4string(mhivspace)<-proj4string(terre) 
        mhivspaceease<-spTransform(mhivspace,CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))
        
        if(g==2000|g==2004|g==2008|g==2012|g==2016){ # leap year
          if (j==1){
            x<-d  
          }else if(j==2){
            x<-31+d
          }else if(j==10){
            x<-274+d
          }else if(j==11){
            x<-305+d
          }else{
            x<-335+d
          }
        }else{# regular year 
          if (j==1){
            x<-d 
          }else if(j==2){
            x<-31+d
          }else if(j==10){
            x<-273+d
          }else if(j==11){
            x<-304+d
          }else{
            x<-334+d
          } 
        }
        
        start <-c(1,1,x)
        count<- c(1440,720,1)	
        temp <-ncvar_get(Temp,"sst", start=start, count=count)
        temp_long <- as.vector(temp)
        Temp_mat <- matrix(temp_long, nrow=dim(lon)*dim(lat), ncol=1)
        lonlat <- as.matrix(expand.grid(lon,lat))
        Temp_df <- data.frame(cbind(lonlat,Temp_mat))              
        names(Temp_df) <- c("lon","lat","SST") 
        coordinates(Temp_df) <- ~lon + lat 
        bfrast<-rasterFromXYZ(Temp_df)
        bfrot<-rotate(bfrast)
        proj4string(bfrot)<-proj4string(terre)
        
        bag<-projectRaster(bfrot,rastease250) #set to the 250km resolution and North Pole Lambert Equal Area Projection
        names(bag)<-"SST"
        
        extr<-extract(bag,mhivspaceease,method="simple",sp=TRUE)
        extrtemp<-as.data.frame(extr)
        
      }else{ # if there is no cyclone for this day
        extrtemp<- data.frame(matrix(NA,ncol=23,nrow=1))
        colnames(extrtemp)<-c(colnames(dhiv),"SST") 
      } 
      
      if(d==1){ 
        mois<-extrtemp
      }else{
        mois<-rbind(mois,extrtemp)
      }
    }
    
    if(j==1){ 
      annee<-mois
    }else{
      annee<-rbind(annee,mois)
    } 
  }
  
  
  if(i==min(fin$year)){
    tempfin<-annee  
  }else{
    tempfin<-rbind(tempfin,annee)
  }   
}

tempfin<-tempfin[complete.cases(tempfin[,15]),] 
write.table(tempfin,"cyclone_SST.csv",sep=";",row.names=FALSE)#sauvegarde 

#Minimum and Maximum Air Temperature: NCEP/NCAR Reanalysis dataset 
#same principle than for SST
rast<- raster(ext=extent(c(0,360,-90,90)), resolution=2.5)
for (i in min(fin$year): max(fin$year)){
  g<-i+2000
  Temp<-nc_open(paste("/Volumes/air.sig995.",g,".nc",sep="")) 
  lon <- ncvar_get(Temp, "lon")
  lat <- ncvar_get(Temp, "lat")
  time<-ncvar_get(Temp,"time")
  tunits <- ncatt_get(Temp,"time","units")
  tas_time <- nc.get.time.series(Temp)
  tas_time<-as.character(tas_time)
  
  yhiv<-subset(fin,fin$year==i) 
  for (j in c(1,2,10:12)){
    mhiv<-subset(yhiv,yhiv$month==j) 
    for (d in min(mhiv$day):max(mhiv$day)){
      dhiv<-subset(mhiv,mhiv$day==d) 
      if(nrow(dhiv)!=0) {
        mhivspace<-dhiv
        coordinates(mhivspace) <- ~longbis+lat
        proj4string(mhivspace)<-proj4string(terre) 
        mhivspaceease<-spTransform(mhivspace,CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))
        
        start <- rep(1,3)	
        if(g==2000|g==2004|g==2008|g==2012|g==2016){ 
          if (j==1){
            x<-d*4-3  
          }else if(j==2){
            x<-31*4+d*4-3 
          }else if(j==10){
            x<-274*4+d*4-3
          }else if(j==11){
            x<-305*4+d*4-3
          }else{
            x<-335*4+d*4-3
          }
        }else{
          if (j==1){
            x<-d*4-3 
          }else if(j==2){
            x<-31*4+d*4-3
          }else if(j==10){
            x<-273*4+d*4-3
          }else if(j==11){
            x<-304*4+d*4-3
          }else{
            x<-334*4+d*4-3
          } 
        }
        
        start <-c(1,1,x)
        count<- c(144,73,4)	
        temp <-ncvar_get(Temp,"air", start=start, count=count)
        temp_long <- as.vector(temp)
        Temp_mat <- matrix(temp_long, nrow=dim(lon)*dim(lat), ncol=4)
        lonlat <- as.matrix(expand.grid(lon,lat))
        Temp_df <- data.frame(cbind(lonlat,Temp_mat))              
        names(Temp_df) <- c("lon","lat","T0","T6","T12","T18") 
        Temp_df$Tmin<-apply(Temp_df[,3:6],1,min)
        Temp_df$Tmax<-apply(Temp_df[,3:6],1,max)
        coordinates(Temp_df) <- ~lon + lat 
        bfrastmin<-rasterize(Temp_df,rast,Temp_df$Tmin)
        bfrastmax<-rasterize(Temp_df,rast,Temp_df$Tmax)
        bfrotmin<-rotate(bfrastmin)
        bfrotmax<-rotate(bfrastmax)
        proj4string(bfrotmin)<-proj4string(terre)
        proj4string(bfrotmax)<-proj4string(terre) 
        bfrotmin<-bfrotmin-273.15
        bfrotmax<-bfrotmax-273.15#Conversion to Kalvin
        names(bfrotmin)<-"Tmin"
        names(bfrotmax)<-"Tmax"
        bfcrmin<-crop(bfrotmin,extent(c(-180,180,0,90)))
        bfcrmax<-crop(bfrotmax,extent(c(-180,180,0,90)))
        bfrotminfin<-projectRaster(bfcrmin,rastease250)
        bfrotmaxfin<-projectRaster(bfcrmax,rastease250)
        
        
        extrmin<-extract(bfrotminfin,mhivspaceease,method="simple",sp=TRUE)
        extrtempmin<-as.data.frame(extrmin)
        
        extrmax<-extract(bfrotmaxfin,mhivspaceease,method="simple",sp=TRUE)
        extrtempmax<-as.data.frame(extrmax)
        
      }else{ 
        extrtempmin<- data.frame(matrix(NA,ncol=23,nrow=1))
        colnames(extrtempmin)<-c(colnames(dhiv),"Tmin")
        
        extrtempmax<- data.frame(matrix(NA,ncol=23,nrow=1))
        colnames(extrtempmax)<-c(colnames(dhiv),"Tmax") 
      } 
      
      if(d==1){ 
        moismin<-extrtempmin
        moismax<-extrtempmax
      }else{
        moismin<-rbind(moismin,extrtempmin)
        moismax<-rbind(moismax,extrtempmax)
      }
    }
    
    if(j==1){ 
      anneemin<-moismin
      anneemax<-moismax
    }else{
      anneemin<-rbind(anneemin,moismin)
      anneemax<-rbind(anneemax,moismax)
    } 
  }
  
  
  if(i==min(fin$year)){
    tempfinmin<-anneemin
    tempfinmax<-anneemax
  }else{
    tempfinmin<-rbind(tempfinmin,anneemin)
    tempfinmax<-rbind(tempfinmax,anneemax)
  }   
}

tempfinmin<-tempfinmin[complete.cases(tempfinmin[,15]),] 
tempfinmax<-tempfinmax[complete.cases(tempfinmax[,15]),] 
write.table(tempfinmin,"cyclone_Tmin.csv",sep=";",row.names=FALSE)#sauvegarde 
write.table(tempfinmax,"cyclone_Tmax.csv",sep=";",row.names=FALSE)#sauvegarde 

#Minimum and Maximum relative humidity:NCEP/NCAR Reanalysis dataset 

rast<- raster(ext=extent(c(0,360,-90,90)), resolution=2.5)
for (i in min(fin$year): max(fin$year)){
  g<-i+2000
  Temp<-nc_open(paste("/Volumes/rhum.sig995.",g,".nc",sep="")) 
  lon <- ncvar_get(Temp, "lon")
  lat <- ncvar_get(Temp, "lat")
  time<-ncvar_get(Temp,"time")
  tunits <- ncatt_get(Temp,"time","units")
  tas_time <- nc.get.time.series(Temp)
  tas_time<-as.character(tas_time)
  
  yhiv<-subset(fin,fin$year==i) 
  for (j in c(1,2,10:12)){
    mhiv<-subset(yhiv,yhiv$month==j) 
    for (d in min(mhiv$day):max(mhiv$day)){
      dhiv<-subset(mhiv,mhiv$day==d) 
      if(nrow(dhiv)!=0) {
        mhivspace<-dhiv
        coordinates(mhivspace) <- ~longbis+lat
        proj4string(mhivspace)<-proj4string(terre) 
        mhivspaceease<-spTransform(mhivspace,CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))
        
        start <- rep(1,3)	
        if(g==2000|i==2004|g==2008|g==2012|g==2016){
          if (j==1){
            x<-d*4-3  
          }else if(j==2){
            x<-31*4+d*4-3 
          }else if(j==10){
            x<-274*4+d*4-3
          }else if(j==11){
            x<-305*4+d*4-3
          }else{
            x<-335*4+d*4-3
          }
        }else{
          if (j==1){
            x<-d*4-3 
          }else if(j==2){
            x<-31*4+d*4-3
          }else if(j==10){
            x<-273*4+d*4-3
          }else if(j==11){
            x<-304*4+d*4-3
          }else{
            x<-334*4+d*4-3
          } 
        }
        
        start <-c(1,1,x)
        count<- c(144,73,4)	
        temp <-ncvar_get(Temp,"rhum", start=start, count=count)
        temp_long <- as.vector(temp)
        Temp_mat <- matrix(temp_long, nrow=dim(lon)*dim(lat), ncol=4)
        lonlat <- as.matrix(expand.grid(lon,lat))
        Temp_df <- data.frame(cbind(lonlat,Temp_mat))              
        names(Temp_df) <- c("lon","lat","Hum0","Hum6","Hum12","Hum18") 
        Temp_df$Hummin<-apply(Temp_df[,3:6],1,min)
        Temp_df$Hummax<-apply(Temp_df[,3:6],1,max)
        coordinates(Temp_df) <- ~lon + lat 
        bfrastmin<-rasterize(Temp_df,rast,Temp_df$Hummin)
        bfrastmax<-rasterize(Temp_df,rast,Temp_df$Hummax)
        bfrotmin<-rotate(bfrastmin)
        bfrotmax<-rotate(bfrastmax)
        proj4string(bfrotmin)<-proj4string(terre)
        proj4string(bfrotmax)<-proj4string(terre) 
        names(bfrotmin)<-"Hummin"
        names(bfrotmax)<-"Hummax"
        bfcrmin<-crop(bfrotmin,extent(c(-180,180,0,90)))
        bfcrmax<-crop(bfrotmax,extent(c(-180,180,0,90)))
        bfrotminfin<-projectRaster(bfcrmin,rastease250)
        bfrotmaxfin<-projectRaster(bfcrmax,rastease250)
        
        extrmin<-extract(bfrotminfin,mhivspaceease,method="simple",sp=TRUE)
        extrtempmin<-as.data.frame(extrmin)
        
        extrmax<-extract(bfrotmaxfin,mhivspaceease,method="simple",sp=TRUE)
        extrtempmax<-as.data.frame(extrmax)
        
      }else{ 
        extrtempmin<- data.frame(matrix(NA,ncol=23,nrow=1))
        colnames(extrtempmin)<-c(colnames(dhiv),"Hummin")
        
        extrtempmax<- data.frame(matrix(NA,ncol=23,nrow=1))
        colnames(extrtempmax)<-c(colnames(dhiv),"Hummax") 
      } 
      
      if(d==1){ 
        moismin<-extrtempmin
        moismax<-extrtempmax
      }else{
        moismin<-rbind(moismin,extrtempmin)
        moismax<-rbind(moismax,extrtempmax)
      }
    }
    
    if(j==1){ 
      anneemin<-moismin
      anneemax<-moismax
    }else{
      anneemin<-rbind(anneemin,moismin)
      anneemax<-rbind(anneemax,moismax)
    } 
  }
  
  
  if(i==min(fin$year)){
    tempfinmin<-anneemin
    tempfinmax<-anneemax
  }else{
    tempfinmin<-rbind(tempfinmin,anneemin)
    tempfinmax<-rbind(tempfinmax,anneemax)
  }   
}

tempfinmin<-tempfinmin[complete.cases(tempfinmin[,15]),] 
tempfinmax<-tempfinmax[complete.cases(tempfinmax[,15]),] 
write.table(tempfinmin,"cyclone_Hummin.csv",sep=";",row.names=FALSE)#sauvegarde 
write.table(tempfinmax,"cyclone_Hummax.csv",sep=";",row.names=FALSE)#sauvegarde 

#Calculation of average characteristics per month and per cyclone category
sst<-read.csv("/Volumes/cyclone_SST.csv",sep=";") 
tmin<-read.csv("/Volumes/cyclone_Tmin.csv",sep=";") 
tmax<-read.csv("/Volumes/cyclone/cyclone_Tmax.csv",sep=";") 
humin<-read.csv("/Volumes/cyclone_Hummin.csv",sep=";") 
humax<-read.csv("/Volumes/cyclone_Hummax.csv",sep=";") 

nb<- data.frame(matrix(NA,ncol=9,nrow=1))
colnames(nb)<-c("Classe","Mois","Annee","SST","Tmin","Tmax","Tmean","Hmin","Hmax") 


for(i in c(1,2,10:12)){#For each month
  for (j in 1:4){#For each cyclone category
    for (y in 0:16){ # For each year 
      nb$Classe<-j
      nb$Mois<-i
      nb$Annee<-y
      #sst
      mois<-subset(sst,sst$month==i)# selection of cyclone for the month "i"
      cl<-subset(mois,mois$classe==j)# selection of cyclone category "j"
      yr<-subset(cl,cl$year==y)
      zone<-subset(yr,yr$longbis>=-3675000 & yr$longbis<=-2825000& yr$lat>=-3375000 & yr$lat<=-2625000)# sélection des cyclones dans la zone
      un<-unique(zone[c("day","sys_num","SST","longbis","lat")])
      #case where a cyclone moves over several cells in the same day: Its daily characteristics are considered to be the average of those encountered
      #On the other hand, from one day to the next, cyclones are considered independent
      bis<-unique(un[c("day","sys_num")])
      if(nrow(bis)!=0){
        
        for (d in 1:nrow(bis)){
          sel<-subset(un,un$day==bis[d,1]&un$sys_num==bis[d,2])
          bis$SST[d]<-mean(sel$SST,na.rm=TRUE)
        } 
        nb$SST<-mean(bis$SST,na.rm=TRUE)
      }else{
        nb$SST<-"NA"
      }
      
      #Tmax,Tmin, Tmean
      mois<-subset(tmin,tmin$month==i)
      cl<-subset(mois,mois$classe==j)
      yr<-subset(cl,cl$year==y)
      zone<-subset(yr,yr$longbis>=-3675000 & yr$longbis<=-2825000& yr$lat>=-3375000 & yr$lat<=-2625000)# sélection des cyclones dans la zone
      un<-unique(zone[c("day","sys_num","Tmin","longbis","lat")])
      bis<-unique(un[c("day","sys_num")])
      if(nrow(bis)!=0){
        
        for (d in 1:nrow(bis)){
          sel<-subset(un,un$day==bis[d,1]&un$sys_num==bis[d,2])
          bis$Tmin[d]<-mean(sel$Tmin,na.rm=TRUE)
        } 
        nb$Tmin<-mean(bis$Tmin,na.rm=TRUE)
      }else{
        nb$Tmin<-"NA"
      }
      
      mois<-subset(tmax,tmax$month==i)
      cl<-subset(mois,mois$classe==j)
      yr<-subset(cl,cl$year==y)
      zone<-subset(yr,yr$longbis>=-3675000 & yr$longbis<=-2825000& yr$lat>=-3375000 & yr$lat<=-2625000)# sélection des cyclones dans la zone
      un<-unique(zone[c("day","sys_num","Tmax","longbis","lat")])
      bis<-unique(un[c("day","sys_num")])
      if(nrow(bis)!=0){
        
        for (d in 1:nrow(bis)){
          sel<-subset(un,un$day==bis[d,1]&un$sys_num==bis[d,2])
          bis$Tmax[d]<-mean(sel$Tmax,na.rm=TRUE)
        } 
        nb$Tmax<-mean(bis$Tmax,na.rm=TRUE)
      }else{
        nb$Tmax<-"NA"
      }
      
      
      #Hmin,Hmax
      mois<-subset(humin,humin$month==i)
      cl<-subset(mois,mois$classe==j)
      yr<-subset(cl,cl$year==y)
      zone<-subset(yr,yr$longbis>=-3675000 & yr$longbis<=-2825000& yr$lat>=-3375000 & yr$lat<=-2625000)# sélection des cyclones dans la zone
      un<-unique(zone[c("day","sys_num","Hummin","longbis","lat")])
      bis<-unique(un[c("day","sys_num")])
      if(nrow(bis)!=0){
        
        for (d in 1:nrow(bis)){
          sel<-subset(un,un$day==bis[d,1]&un$sys_num==bis[d,2])
          bis$Hummin[d]<-mean(sel$Hummin,na.rm=TRUE)
        } 
        nb$Hmin<-mean(bis$Hummin,na.rm=TRUE)
      }else{
        nb$Hmin<-"NA"
      }
      
      mois<-subset(humax,humax$month==i)
      cl<-subset(mois,mois$classe==j)
      yr<-subset(cl,cl$year==y)
      zone<-subset(yr,yr$longbis>=-3675000 & yr$longbis<=-2825000& yr$lat>=-3375000 & yr$lat<=-2625000)# sélection des cyclones dans la zone
      un<-unique(zone[c("day","sys_num","Hummax","longbis","lat")])
      bis<-unique(un[c("day","sys_num")])
      if(nrow(bis)!=0){
        
        for (d in 1:nrow(bis)){
          sel<-subset(un,un$day==bis[d,1]&un$sys_num==bis[d,2])
          bis$Hummax[d]<-mean(sel$Hummax,na.rm=TRUE)
        } 
        nb$Hmax<-mean(bis$Hummax,na.rm=TRUE)
      }else{
        nb$Hmax<-"NA"
      }
      
      if(y==0){
        annee<-nb 
      }else{
        annee<-rbind(annee,nb)  
      }
      
    }
    if(j==1){
      classe<-annee
    }else{
      classe<-rbind(classe,annee)  
    }
  }
  if(i==1){
    final<-classe
  }else{
    final<-rbind(final,classe)  
  }
}
final$Tmin<-as.numeric(final$Tmin)
final$Tmax<-as.numeric(final$Tmax)
final$Tmean<-apply(final[,5:6],1,mean,na.rm=TRUE)
write.table(final,"DescriptioncycloneTerreNeuve.csv",sep=";",row.names=FALSE)

# Extraction of average non-cyclonic conditions within the off Newfoundland areas for each month 
#Extraction for each pixel of the area and averaging 
#We excluded days for which cyclones occured

nb<- data.frame(matrix(NA,ncol=10,nrow=1))
colnames(nb)<-c("Mois","Annee","SST","Tmin","Tmax","Tmean","Hmin","Hmax","Windmin","Windmax") 

annee<-data.frame(matrix(NA,ncol=10,nrow=5))
colnames(annee)<-c("Mois","Annee","SST","Tmin","Tmax","Tmean","Hmin","Hmax","Windmin","Windmax")

mois<-data.frame(matrix(NA,ncol=10,nrow=85))
colnames(mois)<-c("Mois","Annee","SST","Tmin","Tmax","Tmean","Hmin","Hmax","Windmin","Windmax")

terreneuverast<-crop(rasteasecr,extent(c(-3675000,-2825000,-3375000,-2625000)))
TNdt<-as.data.frame(terreneuverast,xy=TRUE)
coordinates(TNdt) <- ~x+y
rast<- raster(ext=extent(c(0,360,-90,90)), resolution=2.5)

hum<-read.csv("/Volumes/cyclone_Hummax.csv",sep=";")

for(j in 0:16){
  print(j)
  g<-j+2000 
  
  #SST
  Temp<-nc_open(paste("/Volumes/sst.day.mean.",g,".nc",sep="")) 
  lon <- ncvar_get(Temp, "lon")
  lat <- ncvar_get(Temp, "lat")
  time<-ncvar_get(Temp,"time")
  tunits <- ncatt_get(Temp,"time","units")
  tas_time <- nc.get.time.series(Temp)
  tas_time<-as.character(tas_time)
  
  for(i in c(1,2,10:12)){
    print(i)
    nb$Mois<-i
    nb$Annee<-g
    humsel<-subset(hum,hum$month==i&hum$year==j)
    humzone<-subset(humsel,humsel$longbis>=-3675000 & humsel$longbis<=-2825000& humsel$lat>=-3375000 & humsel$lat<=-2625000)# sélection des cyclones dans la zone
    jourcyclo<-unique(humzone$day)+2
    
    start <- rep(1,3)	
    if(g==2000|g==2004|g==2008|g==2012|g==2016){ 
      if (i==1){
        x<-1 
        dur<-31
      }else if(i==2){
        x<-32
        dur<-29
      }else if(i==10){
        x<-275
        dur<-31
      }else if(i==11){
        x<-306
        dur<-30
      }else{
        x<-336
        dur<-31
      }
    }else{
      if (i==1){
        x<-1 
        dur<-31
      }else if(i==2){
        x<-32
        dur<-28
      }else if(i==10){
        x<-274
        dur<-31
      }else if(i==11){
        x<-305
        dur<-30
      }else{
        x<-335
        dur<-31
      } 
    }
    
    start <-c(1,1,x)
    count<- c(1440,720,dur)	
    temp <-ncvar_get(Temp,"sst", start=start, count=count)
    temp_long <- as.vector(temp)
    Temp_mat <- matrix(temp_long, nrow=dim(lon)*dim(lat), ncol=dur)
    lonlat <- as.matrix(expand.grid(lon,lat))
    Temp_df <- data.frame(cbind(lonlat,Temp_mat))              
    names(Temp_df) <- c("lon","lat",tas_time[x:(x+dur-1)]) 
    nc_close(Temp)
    bfrast<-rasterFromXYZ(Temp_df)
    bfrot<-rotate(bfrast)
    proj4string(bfrot)<-proj4string(terre)
    bag<-projectRaster(bfrot,rastease250)
    bgcr<-crop(bag,extent(c(-3675000,-2825000,-3375000,-2625000)))#crop for the area off Newfoundland
    bgfr<-as.data.frame(bgcr,xy=TRUE)#conversion in dataframe
    #delete days concerned by cyclones
    bgfrclean<-bgfr[,-jourcyclo] 
    #Calcul of the average on the area 
    bgfrclean$SSTmeanpix<-apply(bgfrclean[,3:ncol(bgfrclean)],1,mean,na.rm=T)
    Temp_dfmean<-bgfrclean[,c(1,2,ncol(bgfrclean))]
    nb$SST<-mean(Temp_dfmean$SSTmeanpix,na.rm=T)
    
    if(i==1){
      annee[1,]<-nb
    }else if (i==2){
      annee[2,]<-nb   
    }else if (i==10){
      annee[3,]<-nb   
    }else if (i==11){
      annee[4,]<-nb  
    }else if (i==12){
      annee[5,]<-nb  
    } 
    
  }
  
  mois[((j*5)+1):((j*5)+5),]<-annee
  
}

write.table(mois,"DescriptionBaselineTerreNeuveSST.csv",sep=";",row.names=FALSE)

for(j in 0:16){
  g<-j+2000 
  #Tmin,Tmax, Tmean
  Temp<-nc_open(paste("/Volumes/air.sig995.",g,".nc",sep="")) 
  lon <- ncvar_get(Temp, "lon")
  lat <- ncvar_get(Temp, "lat")
  time<-ncvar_get(Temp,"time")
  tunits <- ncatt_get(Temp,"time","units")
  tas_time <- nc.get.time.series(Temp)
  tas_time<-as.character(tas_time)
  lonlat <- as.matrix(expand.grid(lon,lat))
  
  for(i in c(1,2,10:12)){
    
    print(i)
    print(j)
    
    nb$Mois<-i
    nb$Annee<-g
    humsel<-subset(hum,hum$month==i&hum$year==j)
    humzone<-subset(humsel,humsel$longbis>=-3675000 & humsel$longbis<=-2825000& humsel$lat>=-3375000 & humsel$lat<=-2625000)
    jourcyclo<-unique(humzone$day)+2 
    
    if(g==2000|g==2004|g==2008|g==2012|g==2016){ 
      if (i==1){
        dur<-31
      }else if(i==2){
        dur<-29
      }else if(i==10){
        dur<-31
      }else if(i==11){
        dur<-30
      }else{
        dur<-31
      }
    }else{# année standard
      if (i==1){
        dur<-31
      }else if(i==2){
        dur<-28
      }else if(i==10){
        dur<-31
      }else if(i==11){
        dur<-30
      }else{
        dur<-31
      } 
    }
    
    for (d in 1:dur){
      print(d)
      start <- rep(1,3)	
      if(g==2000|g==2004|g==2008|g==2012|g==2016){ 
        if (i==1){
          x<-d*4-3 
        }else if(i==2){
          x<-31*4+d*4-3 
        }else if(i==10){
          x<-274*4+d*4-3
        }else if(i==11){
          x<-305*4+d*4-3
        }else{
          x<-335*4+d*4-3
        }
      }else{
        if (i==1){
          x<-d*4-3 
        }else if(i==2){
          x<-31*4+d*4-3
        }else if(i==10){
          x<-273*4+d*4-3
        }else if(i==11){
          x<-304*4+d*4-3
        }else{
          x<-334*4+d*4-3
        } 
      }
      
      start <-c(1,1,x)
      count<- c(144,73,4)	
      temp <-ncvar_get(Temp,"air", start=start, count=count)
      temp_long <- as.vector(temp)
      Temp_mat <- matrix(temp_long, nrow=dim(lon)*dim(lat), ncol=4)
      Temp_df <- data.frame(cbind(lonlat,Temp_mat))              
      names(Temp_df) <- c("lon","lat","T0","T6","T12","T18") 
      Temp_df$Tmin<-apply(Temp_df[,3:6],1,min)
      Temp_df$Tmax<-apply(Temp_df[,3:6],1,max)
      coordinates(Temp_df) <- ~lon + lat 
      bfrastmin<-rasterize(Temp_df,rast,Temp_df$Tmin)
      bfrastmax<-rasterize(Temp_df,rast,Temp_df$Tmax)
      bfrotmin<-rotate(bfrastmin)
      bfrotmax<-rotate(bfrastmax)
      proj4string(bfrotmin)<-proj4string(terre)
      proj4string(bfrotmax)<-proj4string(terre) 
      bfrotmin<-bfrotmin-273.15
      bfrotmax<-bfrotmax-273.15#Conversion Kalvin
      names(bfrotmin)<-paste("Tmin",d,sep="")
      names(bfrotmax)<-paste("Tmax",d,sep="")
      bfcrmin<-crop(bfrotmin,extent(c(-180,180,0,90)))
      bfcrmax<-crop(bfrotmax,extent(c(-180,180,0,90)))
      bfrotminfin<-projectRaster(bfcrmin,rastease250)
      bfrotmaxfin<-projectRaster(bfcrmax,rastease250)
      
      if(d==1){
        stmin<-bfrotminfin
        stmax<-bfrotmaxfin
      }else{
        stmin<-stack(stmin,bfrotminfin)
        stmax<-stack(stmax,bfrotmaxfin) 
      }
    }
    stmincr<-crop(stmin,extent(c(-3675000,-2825000,-3375000,-2625000)))
    stminfr<-as.data.frame(stmincr,xy=TRUE)
    #We deleted days for which cyclones happened
    stminfrclean<-stminfr[,-jourcyclo] 
    stmaxcr<-crop(stmax,extent(c(-3675000,-2825000,-3375000,-2625000)))
    stmaxfr<-as.data.frame(stmaxcr,xy=TRUE)
    stmaxfrclean<-stmaxfr[,-jourcyclo]
    
    #Average calculation for the area
    stminfrclean$Tminpix<-apply(stminfrclean[,3:ncol(stminfrclean)],1,mean,na.rm=T)
    Tempmin<-stminfrclean[,c(1,2,ncol(stminfrclean))]
    
    stmaxfrclean$Tmaxpix<-apply(stmaxfrclean[,3:ncol(stmaxfrclean)],1,mean,na.rm=T)
    Tempmax<-stmaxfrclean[,c(1,2,ncol(stmaxfrclean))]
    
    nb$Tmin<-mean(Tempmin$Tminpix,na.rm=T)
    nb$Tmax<-mean(Tempmax$Tmaxpix,na.rm=T)
    nb$Tmean<-(nb$Tmin+nb$Tmax)/2 
    
    if(i==1){
      annee[1,]<-nb
    }else if (i==2){
      annee[2,]<-nb   
    }else if (i==10){
      annee[3,]<-nb   
    }else if (i==11){
      annee[4,]<-nb  
    }else if (i==12){
      annee[5,]<-nb  
    } 
    
  }
  
  mois[((j*5)+1):((j*5)+5),]<-annee
  
}
write.table(mois,"DescriptionBaselineTerreNeuveTemp.csv",sep=";",row.names=FALSE) 

for(j in 0:16){
  g<-j+2000
  #Hmin,Hmax
  Temp<-nc_open(paste("/Volumes/rhum.sig995.",g,".nc",sep="")) 
  lon <- ncvar_get(Temp, "lon")
  lat <- ncvar_get(Temp, "lat")
  time<-ncvar_get(Temp,"time")
  tunits <- ncatt_get(Temp,"time","units")
  tas_time <- nc.get.time.series(Temp)
  tas_time<-as.character(tas_time)
  for(i in c(1,2,10:12)){
    
    print(i)
    print(j)
    
    nb$Mois<-i
    nb$Annee<-g
    humsel<-subset(hum,hum$month==i&hum$year==j)
    humzone<-subset(humsel,humsel$longbis>=-3675000 & humsel$longbis<=-2825000& humsel$lat>=-3375000 & humsel$lat<=-2625000)
    jourcyclo<-unique(humzone$day)+2        
    
    if(g==2000|g==2004|g==2008|g==2012|g==2016){ 
      if (i==1){
        dur<-31
      }else if(i==2){
        dur<-29
      }else if(i==10){
        dur<-31
      }else if(i==11){
        dur<-30
      }else{
        dur<-31
      }
    }else{
      if (i==1){
        dur<-31
      }else if(i==2){
        dur<-28
      }else if(i==10){
        dur<-31
      }else if(i==11){
        dur<-30
      }else{
        dur<-31
      } 
    }
    
    for (d in 1:dur){
      print(d)
      start <- rep(1,3)	
      if(g==2000|g==2004|g==2008|g==2012|g==2016){ 
        if (i==1){
          x<-d*4-3  
        }else if(i==2){
          x<-31*4+d*4-3 
        }else if(i==10){
          x<-274*4+d*4-3
        }else if(i==11){
          x<-305*4+d*4-3
        }else{
          x<-335*4+d*4-3
        }
      }else{
        if (i==1){
          x<-d*4-3 
        }else if(i==2){
          x<-31*4+d*4-3
        }else if(i==10){
          x<-273*4+d*4-3
        }else if(i==11){
          x<-304*4+d*4-3
        }else{
          x<-334*4+d*4-3
        } 
      }
      start <-c(1,1,x)
      count<- c(144,73,4)
      temp <-ncvar_get(Temp,"rhum", start=start, count=count)
      temp_long <- as.vector(temp)
      Temp_mat <- matrix(temp_long, nrow=dim(lon)*dim(lat), ncol=4)
      lonlat <- as.matrix(expand.grid(lon,lat))
      Temp_df <- data.frame(cbind(lonlat,Temp_mat))              
      names(Temp_df) <- c("lon","lat","Hum0","Hum6","Hum12","Hum18") 
      #nc_close(Temp)
      Temp_df$Hummin<-apply(Temp_df[,3:6],1,min)
      Temp_df$Hummax<-apply(Temp_df[,3:6],1,max)
      coordinates(Temp_df) <- ~lon + lat 
      bfrastmin<-rasterize(Temp_df,rast,Temp_df$Hummin)
      bfrastmax<-rasterize(Temp_df,rast,Temp_df$Hummax)
      bfrotmin<-rotate(bfrastmin)
      bfrotmax<-rotate(bfrastmax)
      proj4string(bfrotmin)<-proj4string(terre)
      proj4string(bfrotmax)<-proj4string(terre) 
      names(bfrotmin)<-paste("Hummin",d,sep="")
      names(bfrotmax)<-paste("Hummax",d,sep="")
      bfcrmin<-crop(bfrotmin,extent(c(-180,180,0,90)))
      bfcrmax<-crop(bfrotmax,extent(c(-180,180,0,90)))
      bfrotminfin<-projectRaster(bfcrmin,rastease250)
      bfrotmaxfin<-projectRaster(bfcrmax,rastease250)
      
      if(d==1){
        stmin<-bfrotminfin
        stmax<-bfrotmaxfin
      }else{
        stmin<-stack(stmin,bfrotminfin)
        stmax<-stack(stmax,bfrotmaxfin) 
      }
    }
    
    stmincr<-crop(stmin,extent(c(-3675000,-2825000,-3375000,-2625000)))
    stminfr<-as.data.frame(stmincr,xy=TRUE)
    stminfrclean<-stminfr[,-jourcyclo] 
    stmaxcr<-crop(stmax,extent(c(-3675000,-2825000,-3375000,-2625000)))
    stmaxfr<-as.data.frame(stmaxcr,xy=TRUE)
    stmaxfrclean<-stmaxfr[,-jourcyclo]
    
    #Averaging
    stminfrclean$Hminpix<-apply(stminfrclean[,3:ncol(stminfrclean)],1,mean,na.rm=T)
    Tempmin<-stminfrclean[,c(1,2,ncol(stminfrclean))]
    
    stmaxfrclean$Hmaxpix<-apply(stmaxfrclean[,3:ncol(stmaxfrclean)],1,mean,na.rm=T)
    Tempmax<-stmaxfrclean[,c(1,2,ncol(stmaxfrclean))]
    
    
    nb$Hmin<-mean(Tempmin$Hminpix,na.rm=T)
    nb$Hmax<-mean(Tempmax$Hmaxpix,na.rm=T)
    
    if(i==1){
      annee[1,]<-nb
    }else if (i==2){
      annee[2,]<-nb   
    }else if (i==10){
      annee[3,]<-nb   
    }else if (i==11){
      annee[4,]<-nb  
    }else if (i==12){
      annee[5,]<-nb  
    } 
    
  }
  
  mois[((j*5)+1):((j*5)+5),]<-annee
  
}

write.table(mois,"DescriptionBaselineTerreNeuveHum.csv",sep=";",row.names=FALSE)

for(j in 0:16){
  g<-j+2000 
  #Windmin, Windmax
  Temp<-nc_open(paste("/Volumes/uwnd.10m.gauss.",g,".nc",sep="")) 
  lon <- ncvar_get(Temp, "lon")
  lat <- ncvar_get(Temp, "lat")
  time<-ncvar_get(Temp,"time")
  tunits <- ncatt_get(Temp,"time","units")
  tas_time <- nc.get.time.series(Temp)
  tas_time<-as.character(tas_time)
  
  Tempv<-nc_open(paste("/Volumes/vwnd.10m.gauss.",g,".nc",sep="")) 
  lonv <- ncvar_get(Tempv, "lon")
  latv <- ncvar_get(Tempv, "lat")
  for(i in c(1,2,10:12)){
    print(i)
    print(j)
    
    nb$Mois<-i
    nb$Annee<-g
    humsel<-subset(hum,hum$month==i&hum$year==j)
    humzone<-subset(humsel,humsel$longbis>=-3675000 & humsel$longbis<=-2825000& humsel$lat>=-3375000 & humsel$lat<=-2625000)
    jourcyclo<-unique(humzone$day)+2
    
    if(g==2000|g==2004|g==2008|g==2012|g==2016){ 
      if (i==1){
        dur<-31
      }else if(i==2){
        dur<-29
      }else if(i==10){
        dur<-31
      }else if(i==11){
        dur<-30
      }else{
        dur<-31
      }
    }else{# année standard
      if (i==1){
        dur<-31
      }else if(i==2){
        dur<-28
      }else if(i==10){
        dur<-31
      }else if(i==11){
        dur<-30
      }else{
        dur<-31
      } 
    }
    
    for (d in 1:dur){
      print(d)
      start <- rep(1,3)	
      if(g==2000|g==2004|g==2008|g==2012|g==2016){ 
        if (i==1){
          x<-d*4-3  
        }else if(i==2){
          x<-31*4+d*4-3 
        }else if(i==10){
          x<-274*4+d*4-3
        }else if(i==11){
          x<-305*4+d*4-3
        }else{
          x<-335*4+d*4-3
        }
      }else{
        if (i==1){
          x<-d*4-3 
        }else if(i==2){
          x<-31*4+d*4-3
        }else if(i==10){
          x<-273*4+d*4-3
        }else if(i==11){
          x<-304*4+d*4-3
        }else{
          x<-334*4+d*4-3
        } 
      }
      start <-c(1,1,x)
      count<- c(192,94,4)	
      #u
      temp <-ncvar_get(Temp,"uwnd", start=start, count=count)
      temp_long <- as.vector(temp)
      Temp_mat <- matrix(temp_long, nrow=dim(lon)*dim(lat), ncol=4)
      lonlat <- as.matrix(expand.grid(lon,lat))
      Temp_df <- data.frame(cbind(lonlat,Temp_mat))              
      names(Temp_df) <- c("lon","lat","U0","U6","U12","U18")
      #nc_close(Temp)
      #v
      tempv <-ncvar_get(Tempv,"vwnd", start=start, count=count)
      temp_longv <- as.vector(tempv)
      Temp_matv <- matrix(temp_longv, nrow=dim(lonv)*dim(latv), ncol=4)
      lonlatv <- as.matrix(expand.grid(lonv,latv))
      Temp_dfv <- data.frame(cbind(lonlatv,Temp_matv))              
      names(Temp_dfv) <- c("lon","lat","V0","V6","V12","V18")
      #nc_close(Tempv)
      # Calculation of the corresponding windspeed 
      Temp_dfws<- data.frame(matrix(NA,ncol=6,nrow=nrow(Temp_df)))
      colnames(Temp_dfws)<-c("lon","lat","WS0","WS6","WS12","WS18") 
      
      Temp_dfws$lon<-Temp_df$lon
      Temp_dfws$lat<-Temp_df$lat
      Temp_dfws$WS0<-sqrt(Temp_df$U0^2+Temp_dfv$V0^2)
      Temp_dfws$WS6<-sqrt(Temp_df$U6^2+Temp_dfv$V6^2)
      Temp_dfws$WS12<-sqrt(Temp_df$U12^2+Temp_dfv$V12^2)
      Temp_dfws$WS18<-sqrt(Temp_df$U18^2+Temp_dfv$V18^2)
      
      Temp_dfws$Windmin<-apply(Temp_dfws[,3:6],1,min)
      Temp_dfws$Windmax<-apply(Temp_dfws[,3:6],1,max)
      coordinates(Temp_dfws) <- ~lon + lat 
      bfrastmin<-rasterize(Temp_dfws,rast,Temp_dfws$Windmin)
      bfrastmax<-rasterize(Temp_dfws,rast,Temp_dfws$Windmax)
      bfrotmin<-rotate(bfrastmin)
      bfrotmax<-rotate(bfrastmax)
      proj4string(bfrotmin)<-proj4string(terre)
      proj4string(bfrotmax)<-proj4string(terre) 
      names(bfrotmin)<-paste("Windmin",d,sep="")
      names(bfrotmax)<-paste("Windmax",d,sep="")
      bfcrmin<-crop(bfrotmin,extent(c(-180,180,0,90)))
      bfcrmax<-crop(bfrotmax,extent(c(-180,180,0,90)))
      bfrotminfin<-projectRaster(bfcrmin,rastease250)
      bfrotmaxfin<-projectRaster(bfcrmax,rastease250)
      
      if(d==1){
        stmin<-bfrotminfin
        stmax<-bfrotmaxfin
      }else{
        stmin<-stack(stmin,bfrotminfin)
        stmax<-stack(stmax,bfrotmaxfin) 
      }
    }
    stmincr<-crop(stmin,extent(c(-3675000,-2825000,-3375000,-2625000)))
    stminfr<-as.data.frame(stmincr,xy=TRUE)
    stminfrclean<-stminfr[,-jourcyclo] 
    stmaxcr<-crop(stmax,extent(c(-3675000,-2825000,-3375000,-2625000)))
    stmaxfr<-as.data.frame(stmaxcr,xy=TRUE)
    stmaxfrclean<-stmaxfr[,-jourcyclo]
    
    #Average calculation
    stminfrclean$Wminpix<-apply(stminfrclean[,3:ncol(stminfrclean)],1,mean,na.rm=T)
    Tempmin<-stminfrclean[,c(1,2,ncol(stminfrclean))]
    
    stmaxfrclean$Wmaxpix<-apply(stmaxfrclean[,3:ncol(stmaxfrclean)],1,mean,na.rm=T)
    Tempmax<-stmaxfrclean[,c(1,2,ncol(stmaxfrclean))]
    
    nb$Windmin<-mean(Tempmin$Wminpix,na.rm=T)
    nb$Windmax<-mean(Tempmax$Wmaxpix,na.rm=T)
    
    
    if(i==1){
      annee[1,]<-nb
    }else if (i==2){
      annee[2,]<-nb   
    }else if (i==10){
      annee[3,]<-nb   
    }else if (i==11){
      annee[4,]<-nb  
    }else if (i==12){
      annee[5,]<-nb  
    } 
    
  }
  
  mois[((j*5)+1):((j*5)+5),]<-annee
  
}
write.table(mois,"DescriptionBaselineTerreNeuveWind.csv",sep=";",row.names=FALSE)
#Melt of the files obtained for the non-cyclonic conditions
sst<-read.csv("/Volumes/DISQUE DUR/Thèse/Tempête/cyclone/DescriptionBaselineTerreNeuveSST.csv",sep=";")
temp<-read.csv("/Volumes/DISQUE DUR/Thèse/Tempête/cyclone/DescriptionBaselineTerreNeuveTemp.csv",sep=";")
hum<-read.csv("/Volumes/DISQUE DUR/Thèse/Tempête/cyclone/DescriptionBaselineTerreNeuveHum.csv",sep=";")
wd<-read.csv("/Volumes/DISQUE DUR/Thèse/Tempête/cyclone/DescriptionBaselineTerreNeuveWind.csv",sep=";")

dbTN<-sst
dbTN$Tmin<-temp$Tmin
dbTN$Tmax<-temp$Tmax
dbTN$Tmean<-temp$Tmean
dbTN$Hmin<-hum$Hmin
dbTN$Hmax<-hum$Hmax
dbTN$Windmin<-wd$Windmin
dbTN$Windmax<-wd$Windmax
write.table(dbTN,"DescriptionBaselineTerreNeuve.csv",sep=";",row.names=FALSE)

cyTN<-read.csv("/Volumes/DISQUE DUR/Thèse/Tempête/cyclone/DescriptioncycloneTerreNeuve.csv",sep=";")

#Average calculation between 2000-2016 of environmental carateristics for non-cyclonic conditions 
nb<- data.frame(matrix(NA,ncol=18,nrow=5))
colnames(nb)<-c("Classe","Mois","SST","Tmin","Tmax","Tmean","Hmin","Hmax","Windmin","Windmax","sdsst","sdtmin","sdtmax","sdtmean","sdhmin","sdhmax","sdwdmin","sdwdmax") 

nb$Classe<-"Baseline"
nb$Mois<-c(1,2,10,11,12)
nb$Mois<-as.factor(nb$Mois)
for(i in c(1,2,10:12)){
  m<-subset(dbTN,dbTN$Mois==i)
  mean<-apply(m[,3:10],2,mean,na.rm=TRUE)
  sd<-apply(m[,3:10],2,sd,na.rm=TRUE)
  
  if(i==1){
    nb[1,3:10]<-mean
    nb[1,11:18]<-sd
  }else if (i==2){
    nb[2,3:10]<-mean
    nb[2,11:18]<-sd  
  }else if (i==10){
    nb[3,3:10]<-mean
    nb[3,11:18]<-sd 
  } else if (i==11){ 
    nb[4,3:10]<-mean
    nb[4,11:18]<-sd
  }else{ 
    nb[5,3:10]<-mean
    nb[5,11:18]<-sd 
  }
}
##Average calculation between 2000-2016 of environmental carateristics cyclonic conditions (for each category)
cy<- data.frame(matrix(NA,ncol=18,nrow=1))
colnames(cy)<-c("Classe","Mois","SST","Tmin","Tmax","Tmean","Hmin","Hmax","Windmin","Windmax","sdsst","sdtmin","sdtmax","sdtmean","sdhmin","sdhmax","sdwdmin","sdwdmax") 


for(i in c(1,2,10:12)){
  cy$Mois<-i
  cy$Mois<-as.factor(cy$Mois)
  m<-subset(cyTN,cyTN$Mois==i)
  for(j in c(1:4)){
    cy$Classe<-j
    cl<-subset(m,m$Classe==j)
    mean<-apply(cl[,4:9],2,mean,na.rm=TRUE)
    sd<-apply(cl[,4:9],2,sd,na.rm=TRUE)
    
    cy[,3:8]<-mean
    cy[,11:16]<-sd 
    
    cy[,17:18]<-0
    if (j==1){
      cy[,10]<-13 
    }else if(j==2){
      cy[,10]<-17
    }else if(j==3){
      cy[,10]<-32.5
    }else{
      cy[,10]<-70 
    }
    
    if(j==1){
      cltot<-cy
    }else{
      cltot<-rbind(cltot,cy)  
    }  
    
  }
  
  if(i==1){
    annee<-cltot
  }else{
    annee<-rbind(annee,cltot)  
  }
}

#Melting of the dataframes obtained
tot<-rbind(nb,annee)
write.table(tot,"CompTerreNeuve.csv",sep=";",row.names=FALSE)

#Graphical representation (Figure S4)
#SST
library(ggplot2)
library(ggthemes) 
tot$Month2 <- factor(tot$Mois, levels=unique(tot$Mois))
tot$Month2 <- factor(tot$Month2,levels(tot$Month2)[c(3,4,5,1,2)])
tot$Classe <- as.factor(tot$Classe)

p<-ggplot(tot, aes(x=Month2, y=SST,color=Classe,group=Classe))+
  scale_colour_manual(name=c("Conditions"),values=c("#41BEC4","#1D91C0","#225EA8","#081D58","#7FCDBB"),labels=c("Non cyclonic","Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("Baseline","1","2","3","4"))+ 
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=SST, ymax=SST+sdsst),colour="black", width=.1,alpha=0.4,position=position_dodge(.05))+
  xlab("Month") +
  ylab("Sea Surface Temperature (°C)\n") +
  scale_x_discrete(labels=c("Oct","Nov","Dec","Jan","Feb"))+
  theme_classic() 


ggsave(filename ="SSTTNfig1.jpeg",width = 20, height = 20, units = "cm")

#Tmin
p<-ggplot(tot, aes(x=Month2, y=Tmin,color=Classe,group=Classe))+
  scale_colour_manual(name=c("Conditions"),values=c("#41BEC4","#1D91C0","#225EA8","#081D58","#7FCDBB"),labels=c("Non cyclonic","Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("Baseline","1","2","3","4"))+  
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Tmin, ymax=Tmin+sdtmin),colour="black", width=.1,alpha=0.4,position=position_dodge(.05))+
  xlab("Month") +
  ylab("Minimum Air Temperature (°C)\n") +
  scale_x_discrete(labels=c("Oct","Nov","Dec","Jan","Feb"))+
  theme_classic() 


ggsave(filename ="TminTNfig1.jpeg",width = 20, height = 20, units = "cm")
#Tmax
p<-ggplot(tot, aes(x=Month2, y=Tmax,color=Classe,group=Classe))+
  scale_colour_manual(name=c("Conditions"),values=c("#41BEC4","#1D91C0","#225EA8","#081D58","#7FCDBB"),labels=c("Non cyclonic","Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("Baseline","1","2","3","4"))+  
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Tmax, ymax=Tmax+sdtmax),colour="black", width=.1,alpha=0.4,position=position_dodge(.05))+
  xlab("Month") +
  ylab("Maximum Air Temperature (°C)\n") +
  scale_x_discrete(labels=c("Oct","Nov","Dec","Jan","Feb"))+
  theme_classic() 


ggsave(filename ="TmaxTNfig1.jpeg",width = 20, height = 20, units = "cm")
#Hmin
p<-ggplot(tot, aes(x=Month2, y=Hmin,color=Classe,group=Classe))+
  scale_colour_manual(name=c("Conditions"),values=c("#41BEC4","#1D91C0","#225EA8","#081D58","#7FCDBB"),labels=c("Non cyclonic","Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("Baseline","1","2","3","4"))+
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Hmin, ymax=Hmin+sdhmin),colour="black", width=.1,alpha=0.4,position=position_dodge(.05))+
  xlab("Month") +
  ylab("Minimum relative humidity (%)\n") +
  scale_x_discrete(labels=c("Oct","Nov","Dec","Jan","Feb"))+
  theme_classic() 


ggsave(filename ="HminTNfig1.jpeg",width = 20, height = 20, units = "cm")
#Hmax
p<-ggplot(tot, aes(x=Month2, y=Hmax,color=Classe,group=Classe))+
  scale_colour_manual(name=c("Conditions"),values=c("#41BEC4","#1D91C0","#225EA8","#081D58","#7FCDBB"),labels=c("Non cyclonic","Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("Baseline","1","2","3","4"))+
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Hmax, ymax=Hmax+sdhmax),colour="black", width=.1,alpha=0.4,position=position_dodge(.05))+
  xlab("Month") +
  ylab("Maximum relative humidity (%)\n") +
  scale_x_discrete(labels=c("Oct","Nov","Dec","Jan","Feb"))+
  theme_classic() 


ggsave(filename ="HmaxTNfig1.jpeg",width = 20, height = 20, units = "cm")
#Windmax
p<-ggplot(tot, aes(x=Month2, y=Windmax,color=Classe,group=Classe))+
  scale_colour_manual(name=c("Conditions"),values=c("#CCFFFF","#99CCFF","#6699CC","#336699","#003366"),labels=c("Usual","Cyclone category 1", "Cyclone category 2","Cyclone category 3","Cyclone category 4"),breaks=c("Baseline","1","2","3","4"))+ 
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Windmax, ymax=Windmax+sdwdmax),colour="black", width=.1,alpha=0.4,position=position_dodge(.05))+
  xlab("Month") +
  ylab("Mean maximal wind speed values experienced between 2000 and 2016 
       off Newfounland\n") +
  scale_x_discrete(labels=c("Jan","Feb","Oct","Nov","Dec"))+
  theme_classic() 


ggsave(filename ="WSmaxTN.jpeg",width = 20, height = 20, units = "cm")

####Calculation of UD and core use areas for each seabird colony during winter####
#Example for little auk 
merg<-read.csv("/Volumes/Totmergoceanperiod.csv",sep=";",stringsAsFactors=T)#This dataframe contained all little auk locations: see KEY ressource Table
merg$colony<-as.character(merg$colony)
#Formatting for use
merg$DateT<-paste(merg$Date,merg$Time,sep=" ")
merg$DateT<-gsub("\\-", "/",merg$DateT) 
merg$DateT<-as.POSIXlt(merg$DateT,"%d/%m/%Y %H:%M",tz="GMT")
merg$BA<-paste(merg$ID,merg$Year,sep="_")

#For each individual, we counted the number of locations for each wintering month
for( c in c(unique(merg$colony))){
  print(c)
  subcol<-subset(merg, merg$colony==c)
  for(i in c(1:2,10:12)){
    print(i)
    submois<-subset(subcol,subcol$m==i)
    if(nrow(submois)!=0){
      submois$BA<-as.character(submois$BA) 
      num<-c(unique(submois$BA)) # number of individual for the colony "c" and the month "i"
      # creation of a dataframe implemented with the number of location per individual 
      numtab<- data.frame(matrix(NA,ncol=2,nrow=length(num)))
      colnames(numtab)<-c("ID","nbpoint")
      for(j in 1:length(num)){
        numtab[j,1]<-num[j]
        numtab[j,2]<-nrow(subset(submois,submois$BA==num[j]))
      }
      write.table(numtab,paste("nbpoint_",c,"_",i,".csv",sep="")) 
      #graphical representation
      t<-paste("Histo_nombrepoint_",c,"_",i,".jpeg",sep="")
      jpeg(filename =t, width = 2000, height = 2000, units = "px", res=300)
      hist(numtab$nbpoint,breaks="Sturges",main="nombre d'individu pour chaque nombre de points",col="lightblue",border="red",xlab="nombre de points",ylab="nombre d'individu",freq=TRUE,labels=T)
      dev.off() 
    }
  }
}

#Individual kernel calculation with BRB-MKDE.exe (long lat version:https://www.cefe.cnrs.fr/fr/recherche/ee/ec/216-simon-benhamou?tmpl=component&type=raw)
for( c in c(unique(merg$colony))){
  print(c)
  subcol<-subset(merg, merg$colony==c)
  for(i in c(1:2,10:12)){
    print(i)
    submois<-subset(subcol,subcol$m==i)
    if(nrow(submois)!=0){
      submois$BA<-as.character(submois$BA)               
      for(j in c(unique(submois$BA))){
        print(j)
        subid<-subset(submois,submois$BA==j)
        subidMKDE<-subid[,c(19,8,7)]
        subidMKDE<-subidMKDE[order(subidMKDE$DateT),]
        nom3<-paste("Merg",c,"_",i,"_",j,"MKDE.ASCII.txt",sep="")
        write.table(subidMKDE,nom3,sep=" ",row.names=FALSE)
        F <-paste("/Volumes/Merg",c,"_",i,"_",j,"MKDE.ASCII",sep="")                
        parametersVector <- paste("BRB_MKDE_DD_nograph",F,"0 250 0 0 0 0 0 0 2 0 1.5 0 0 -181.5 -1.5 1.5",sep=" ") #parametrization
        outputCommand <- system(as.character(parametersVector),show.output.on.console = FALSE, invisible=TRUE,wait=FALSE, intern=TRUE)
      } 
    }  
  }
} 

#Perform stability analysis for each individual (see STAR METHOD)
for( c in c(unique(merg$colony))){
  print(c)
  for(i in c(1,2,10:12)){
    print(i)
    nbtab<-read.csv(paste("/Volumes/nbpoint_",c,"_",i,".csv",sep=""),sep="")
    sel<-nbtab[sample(nrow(nbtab),nrow(nbtab)), ]
    sel$ID<-as.character(sel$ID)
    # For each individual, we created a folder with files containing more and more locations in order to calculate the overlap 
    for(j in c(unique(sel$ID))){
      if (file.exists(paste("/Volumes/UD_Merg",c,"_",i,"_",j,"MKDE.ASCII.asc",sep=""))==TRUE){
        ind<-read.delim(paste("/Volumes/Merg",c,"_",i,"_",j,"MKDE.ASCII.txt",sep=""),header=TRUE)
        #one folder per individual
        dir.create(paste("/Volumes/Merg",c,"_",i,"_",j,"",sep=""))
        file.copy("/Volumes/BRB_MKDE_DD_nograph.exe", paste("/Volumes/Merg",c,"_",i,"_",j,"",sep=""))
        file.copy("/Volumes/Overlap.exe", paste("/Volumes/Merg",c,"_",i,"_",j,"",sep=""))
        file.copy("/Volumes/BRB_MKDE_DD.txt", paste("/Volumes/Merg",c,"_",i,"_",j,"",sep=""))
        setwd(paste("/Volumes/Merg",c,"_",i,"_",j,"",sep=""))
        write.table(ind,paste("Merg",c,"_",i,"_",j,"_0_MKDE.ASCII.txt",sep=""),sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
        
        #We determine the number of whole turns of the loop that we can do
        arr<-floor(nrow(ind)/4)
        for(r in 1:arr){
          n=r*4
          imp<-as.data.frame(ind[1:n,1])
          write.table(imp,paste("Merg",c,"_",i,"_",j,"_",r,"_MKDE.ASCII.txt",sep=""),sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
        }
        
        setwd(paste("/Volumes/Merg",c,"_",i,"_",j,"",sep=""))
        for (g in 0:arr){
          print(g)
          F <-paste("/Volumes/Merg",c,"_",i,"_",j,"/Merg",c,"_",i,"_",j,"_",g,"_MKDE.ASCII",sep="")                
          parametersVector <- paste("BRB_MKDE_DD_nograph",F,"0 250 0 0 0 0 0 0 2 0 1.5 0 0 -181.5 -1.5 1.5",sep=" ")
          outputCommand <- system(as.character(parametersVector),show.output.on.console = FALSE, invisible=TRUE,wait=FALSE, intern=TRUE)
        }
        overlap<-system(as.character("Overlap"),show.output.on.console = FALSE, invisible=TRUE,wait=FALSE, intern=TRUE)#overlap calculation
        
        x<-list.files(paste("/Volumes/Merg",c,"_",i,"_",j,"/",sep=""),pattern="^UD_Merg")
        if (length(x)>1){
          over<-read.delim(paste("/Volumes/Merg",c,"_",i,"_",j,"/over.txt",sep=""),sep="\t",header=TRUE)
          st<- data.frame(matrix(NA,ncol=3,nrow=arr+1))
          colnames(st)<-c("ID","nbpoint","%overlap") 
          
          st$ID<-as.character(j)
          for(n in 1:(arr+1)){
            if(n==1){
              st$nbpoint[n]<-4
            }else{
              st[n,2]<-st[n-1,2]+4
            }
          }
          
          st[arr+1,3]<-100
          for(w in 1:arr){
            if(w<10){
              sto<-strsplit(as.vector(over[(arr-(10-w))*5+1,1]),split=" ")
              st[w,3]<- as.numeric(sto[[1]][length(sto[[1]])])
            }else{         
              sto<-strsplit(as.vector(over[(w-10)*5+1,1]),split=" ")
              st[w,3]<-as.numeric(sto[[1]][length(sto[[1]])])
            }
          }
        }else{
          st<- data.frame(matrix(NA,ncol=3,nrow=1))
          colnames(st)<-c("ID","nbpoint","%overlap")
          st$ID<-as.character(j) 
        }
      } else{
        st<- data.frame(matrix(NA,ncol=3,nrow=1))
        colnames(st)<-c("ID","nbpoint","%overlap")
        st$ID<-as.character(j) 
      } 
      
      if(j==c(unique(sel$ID))[1]){
        stab<-st
      }else{
        stab<-rbind(stab,st)
      }
    }
    write.table(stab,paste("Stab_Merg",c,"_",i,".csv",sep=""),sep=";") 
  }
}

# graphical analysis of the overlaps previously created: example for little auk from Kap Hoegh in december
st<-read.table("/Volumes/Stab_MergKapHoegh_12.csv",sep=";")
st$nbpoint2<-factor(st$nbpoint,levels=unique(st$nbpoint))
st$ID<-as.factor(st$ID)
p<-ggplot(st, aes(x=nbpoint2, y=X.overlap, color=ID, group=ID))+
  geom_line(size=1)
p
#It helps us to dtermine from how many locations taken into account in the UD calculation, the % overlap remained stable 
# In the little auk case, we only considered individuals having more than 30 locations 

#For each colony and each month, we calculated a mean UD
e<-extent(-180,180,0,90)
for( c in c(unique(merg$colony))){
  print(c)
  for(i in c(1:2,10:12)){
    print(i)
    nbtab<-read.csv(paste("/Volumes/nbpoint_",c,"_",i,".csv",sep=""),sep="")
    nbtabred<-subset(nbtab,30<=nbtab$nbpoint)# selection of individual having more than 30 locations (see above) 
    sel<-nbtabred[sample(nrow(nbtabred),nrow(nbtabred)), ]
    sel$ID<-as.character(sel$ID)
    
    for(j in c(unique(sel$ID))){
      print(j)
      #if the UD file exists, open it, if not, create an empty raster 
      if (file.exists(paste("/Volumes/UD_Merg",c,"_",i,"_",j,"MKDE.ASCII.asc",sep=""))==TRUE){
        rastmean<-raster(paste("/Volumes/UD_Merg",c,"_",i,"_",j,"MKDE.ASCII.asc",sep=""))
        proj4string(rastmean)<-proj4string(terre)
        rastmeane<-extend(rastmean,e,0)
        rastmeane<-crop(rastmeane,e)
      }else{
        print(NA)
        rastmeane<-raster(nrows=60, ncols=240,e, 1.5, vals=NA)
        proj4string(rastmeane)<-proj4string(terre)
      }
      
      if(j==c(unique(sel$ID))[1]){
        rastmeanfinal<-rastmeane
      }else{
        rastmeanfinal<-stack(rastmeanfinal,rastmeane)
      }
    }
    rastmeanfinalval<-mean(rastmeanfinal,na.rm=TRUE)# average calculation
    fin<-paste("UDmean_",c,"_",i,".grd", sep="")
    writeRaster(rastmeanfinalval,fin)
  }  
}

#Core use areas calculation (25% kernel) for each colony and each wintering month
rast<- raster(ext=extent(c(-180,180,0,90)), resolution=1.5)

for (i in c(1,2,10:12)) {
  for (co in c(unique(merg$colony))){
    ud<-raster(paste("/Volumes/UDmean_",co,"_",i,".gri",sep=""))
    udt<-as.data.frame(ud,xy=T)
    #We checked that the sum of UD values is equal to 1
    total<-sum(udt$layer)
    #We deleted 0 before we sorted the remaining UD values 
    udt<-subset(udt,udt$layer>0)
    udtri<- udt[order(udt$layer,decreasing=T),] #Sorting 
    #We kept the UD values of the first 25%
    udt25<-udtri[c(1:(nrow(udtri)*25/100)),]
    min(udt25$layer)
    udt25<-subset(udt25,udt25$layer>min(udt25$layer))
    coordinates(udt25)<- ~x+y
    udt25r<-rasterize(udt25,rast,udt25$layer,background=0)
    r25<-rasterToPolygons(clump(udt25r>0), na.rm=TRUE, dissolve=TRUE) # creation and saving of the core use areas delimitation
    nom<-paste("/Volumes/kernel25",co,"_",i,".shp",sep="")
    writeOGR(r25,nom,"r25",driver="ESRI Shapefile")
  }
}


####Graphical representation of overlaps between cyclones the studied area (100°W-100°E, 0°N-90°N) and seabird core use areas (Figure 1, Data S1A-DataS1F)####
#Example with little auk
#For each month and cyclone category
#Mapping of the overlap between cyclones and core use areas
maxv <- 3
minv <- 0
#brewer.pal(n = 9, name = "Greys")
cols <-c("#F0F0F0","#BDBDBD","#737373","#525252","#000000")
pal <- colorRampPalette(cols)
brks <- seq(minv,maxv,by=0.3)
nbrks <- length(brks)-1

cote<-readOGR("/Volumes/reg10coastlinepm.shp")# coast lines shapefile
cotep<-spTransform(cote, CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))#projection in North Pole Lambert Equal Area

grat<- readOGR("/Volumes/ne_50m_graticules_15.shp")# 15° graticules
gratp<- spTransform(grat, CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))#projection in North Pole Lambert Equal Area

for (i in c(1,2,10:12)){ #For each month
  # Opening of the corresponding core use areas
  KH25<-readOGR(paste("/Volumes/kernel25KapHoegh_",i,".shp",sep=""))
  KH25p<-spTransform(KH25, CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))#projection in North Pole Lambert Equal Area
  
  P25<-readOGR(paste("/Volumes/kernel25Qoororsuaq_",i,".shp",sep=""))
  P25p<-spTransform(P25, CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))
  
  for (c in c(1:4)){ #For each cyclone category)
    #Opening of the corresponding raster (Average number of cyclone of the class "c" during the month "i")
    densmoy<-raster(paste("/Volumes/NBmoycyclone_classe",c,"_mois",i,".tif"))
    
    #Graphical representation
    t<-paste("Overlap_mergule_mois_",i,"classe_",c,".jpeg",sep="")
    jpeg(filename =t, width = 3200, height = 2800, units = "px", res=300)
    plot(densmoy,breaks=brks,col=pal(nbrks),zlim=c(minv,maxv),legend=T,useRaster=FALSE,box=FALSE,axes=F,main=paste("month_",i,"classe_",c,sep=""))
    plot(cotep,col="#333333",border=NA,cex=0.1,add=T)#coastline delimitation
    plot(gratp,col="#FFFFFF",border="#FFFFFF",cex=0.1,add=T)#graticules
    plot(KH25p,border="#0099FF",cex=1,add=T)
    plot(P25p,border="#66FF66",cex=1,add=T)
    
    
    legend("right",pch="-",cex=1,pt.bg=c("#33CCFF","#66FF66"),legend=c("Kap Hoegh","Qoororsuaq"),col=c("#33CCFF","#66FF66"),box.lty=0,bg="white", ncol=1)
    dev.off()  
  }
}

#We then counted for each wintering month, how many cyclones of each category happened in the core use areas of each colony
hum<-read.csv("/Volumes/CycloneNAClasse.csv",sep=",")

nb<- data.frame(matrix(NA,ncol=5,nrow=1))
colnames(nb)<-c("Colonies","Classe","Mois","Nombremoyen","SD") 
moy<- data.frame(matrix(NA,ncol=1,nrow=1))
colnames(moy)<-c("Nombremoyen") 

for (c in c("KapHoegh","Qoororsuaq")){
  nb$Colonies<-c
  for (m in c(1,2,10:12)){
    mois<-subset(hum,hum$month==m)
    nb$Mois<-m
    #For each cyclone category
    for (i in c(1:4)){
      nb$Classe<-i
      col<-readOGR(paste("/Volumes/kernel25",c,"_",m,".shp",sep=""))
      colp<-spTransform(col, CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))
      #selection of the cyclone concerned
      cl<-subset(mois,mois$classe==i)
      #Identification of the cyclones falling within the core use ares considered 
      clsp<-cl
      coordinates(clsp) <- ~longbis+lat
      proj4string(clsp)<-proj4string(terre)
      clsp<-spTransform(clsp, CRS=CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"))#projection in North Pole Lambert Equal Area
      overlap<-over(colp,clsp,returnList = TRUE)
      
      #We selected them and calculated their mean number between 2000 and 2016 
      over<-overlap$`0`
      for (y in 0:16){
        yr<-subset(over,over$year==y)
        un<-unique(yr[c("sys_num")])
        moy$Nombremoyen<-nrow(un)
        if(y==0){
          moyannee<-moy  
        }else{
          moyannee<-rbind(moyannee,moy)  
        }
      }
      nb$Nombremoyen<-mean(moyannee$Nombremoyen)
      nb$SD<-sd(moyannee$Nombremoyen)
      
      if(i==1){
        clfin<-nb  
      }else{
        clfin<-rbind(clfin,nb)  
      }
    }
    
    if(m==1){
      moisfin<-clfin  
    }else{
      moisfin<-rbind(moisfin,clfin)  
    }
  }
  
  if(c=="KapHoegh"){
    nbfinal<-moisfin  
  }else{
    nbfinal<-rbind(nbfinal,moisfin)  
  }
}
# Saving
write.table(nbfinal,"NombremoyencycloneMergulecol.csv",sep=";",row.names=FALSE)

#Grpahical representation
#for each month
for (i in c(1,2,10:12)){
  moisbis<-subset(nbfinal,nbfinal$Mois==i)
  moisbis$Colonies <- factor(moisbis$Colonies)
  moisbis$Classe <- factor(moisbis$Classe)
  
  p<-ggplot(moisbis, aes(x=Colonies, y=Nombremoyen,fill=Classe))+
    geom_bar(stat = "identity",position = position_dodge())+
    geom_errorbar(aes(ymin=Nombremoyen-SD, ymax=Nombremoyen+SD),colour="black", width=.1,alpha=0.4,
                  position=position_dodge(.9))+
    scale_fill_manual(name=c("Cyclones intensities"),values=c("#99CCFF","#3399CC","#0066CC","#003366"),labels=c("Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("1","2","3","4"))+
    xlab("Colonies") +
    ylab("Mean number of cyclones in the core areas 
       between 2000 and 2016\n") +
    theme_classic() 
  
  ggsave(filename =paste("NbmoyencycloneMergulebarplot",i,".jpeg",sep=""),width = 20, height = 20, units = "cm")  
}

####Calculation of the cyclone exposure index####
#The two previous step were made for each colony of each species

for (i in c(1,2,10:12)){
  r<-raster(paste("/Volumes/HIVNBmoycyclone_mois ",i," .tif",sep=""))#Average number of cyclones within the studied area (100°W-100°E, 0°N-90°N) for the month "i" 
  #Opening of all the average UD for the correspondinh month, reprojection and resampling 
  #Little auk
  kh<-raster(paste("/Volumes/UDmean_KapHoegh_",i,".grd",sep=""))
  khea<-projectRaster(kh,r)
  pa<-raster(paste("/Volumes/UDmean_Qoororsuaq_",i,".grd",sep=""))
  paea<-projectRaster(pa,r)
  #Atlantic puffin
  fl<-raster(paste("/Volumes/UDmean_Flatey_",i,".grd",sep=""))
  flea<-projectRaster(fl,r)
  st<-raster(paste("/Volumes/UDmean_Storholdi_",i,".grd",sep=""))
  stea<-projectRaster(st,r)
  si<-raster(paste("/Volumes/UDmean_SkomerIsland_",i,".grd",sep=""))
  siea<-projectRaster(si,r)
  sm<-raster(paste("/Volumes/UDmean_SkelligMichael_",i,".grd",sep=""))
  smea<-projectRaster(sm,r)
  gu<-raster(paste("/Volumes/UDmean_GullIsland_",i,".grd",sep=""))
  guea<-projectRaster(gu,r)
  ms<-raster(paste("/Volumes/UDmean_MachiasSealIsland_",i,".grd",sep=""))
  msea<-projectRaster(ms,r)
  #Black-legged kittiwake
  kg<-raster(paste("/Volumes/BLK/UDmean_Kongsfjorden_",i,".grd",sep=""))
  kgea<-projectRaster(kg,r)
  ki<-raster(paste("/Volumes/BLK/UDmean_Kippaku_",i,".grd",sep=""))
  kiea<-projectRaster(ki,r)
  sk<-raster(paste("/Volumes/BLK/UDmean_Sklinna_",i,".grd",sep=""))
  skea<-projectRaster(sk,r)
  raa<-raster(paste("/Volumes/BLK/UDmean_RundeandAlesund_",i,".grd",sep=""))
  raaea<-projectRaster(raa,r)
  ro<-raster(paste("/Volumes/BLK/UDmean_Rost_",i,".grd",sep=""))
  roea<-projectRaster(ro,r)
  las<-raster(paste("/Volumes/BLK/UDmean_LanganesandSkjalfandi_",i,".grd",sep=""))
  lasea<-projectRaster(las,r)
  iom<-raster(paste("/Volumes/BLK/UDmean_IsleofMay_",i,".grd",sep=""))
  iomea<-projectRaster(iom,r)
  is<-raster(paste("/Volumes/BLK/UDmean_Isfjorden_",i,".grd",sep=""))
  isea<-projectRaster(is,r)
  ho<-raster(paste("/Volumes/BLK/UDmean_Hornoya_",i,".grd",sep=""))
  hoea<-projectRaster(ho,r)
  fjl<-raster(paste("/Volumes/BLK/UDmean_FranzJosefLand_",i,".grd",sep=""))
  fjlea<-projectRaster(fjl,r)
  fi<-raster(paste("/Volumes/BLK/UDmean_FaroeIslands_",i,".grd",sep=""))
  fiea<-projectRaster(fi,r)
  ka<-raster(paste("/Volumes/BLK/UDmean_FaroeIslands_",i,".grd",sep=""))
  kaea<-projectRaster(ka,r)
  ck<-raster(paste("/Volumes/BLK/UDmean_CapeKrutik_",i,".grd",sep=""))
  ckea<-projectRaster(ck,r)
  bj<-raster(paste("/Volumes/BLK/UDmean_Bjornoya_",i,".grd",sep=""))
  bjea<-projectRaster(bj,r)
  an<-raster(paste("/Volumes/BLK/UDmean_Anda_",i,".grd",sep=""))
  anea<-projectRaster(an,r)
  ak<-raster(paste("/Volumes/BLK/UDmean_Alkefjellet_",i,".grd",sep=""))
  akea<-projectRaster(ak,r)
  #Common Guillemot
  bjgt<-raster(paste("/Volumes/CG/UDmean_Bjornoya_",i,".grd",sep=""))
  bjgtea<-projectRaster(bjgt,r)
  figt<-raster(paste("/Volumes/CG/UDmean_FaroeIslands_",i,".grd",sep=""))
  figtea<-projectRaster(figt,r)
  cg<-raster(paste("/Volumes/CG/UDmean_CapeGorodetskiy_",i,".grd",sep=""))
  cgea<-projectRaster(cg,r)
  gr<-raster(paste("/Volumes/CG/UDmean_Grimsey_",i,".grd",sep=""))
  grea<-projectRaster(gr,r)
  hj<-raster(paste("/Volumes/CG/UDmean_Hjelmsoya_",i,".grd",sep=""))
  hjea<-projectRaster(hj,r)
  skgt<-raster(paste("/Volumes/CG/UDmean_Sklinna_",i,".grd",sep=""))
  skgtea<-projectRaster(skgt,r)
  hogt<-raster(paste("/Volumes/CG/UDmean_Hornoya_",i,".grd",sep=""))
  hogtea<-projectRaster(hogt,r)
  la<-raster(paste("/Volumes/CG/UDmean_Latrabjarg_",i,".grd",sep=""))
  laea<-projectRaster(la,r)
  lasgt<-raster(paste("/Volumes/CG/UDmean_LanganesandSkjalfandi_",i,".grd",sep=""))
  lasgtea<-projectRaster(lasgt,r)
  jm<-raster(paste("/Volumes/CG/UDmean_JanMayen_",i,".grd",sep=""))
  jmea<-projectRaster(jm,r)
  #Brunnich's guillemot
  ##Kernels for Alkefjellet and Franz Josef Land  were not available for all wintering months.
  sa<-raster(paste("/Volumes/BG/UDmean_Saunders_",i,".grd",sep=""))
  saea<-projectRaster(sa,r)
  ri<-raster(paste("/Volumes/BG/UDmean_Ritenbenk_",i,".grd",sep=""))
  riea<-projectRaster(ri,r)
  psb<-raster(paste("/Volumes/BG/UDmean_PSB_",i,".grd",sep=""))
  psbea<-projectRaster(psb,r)
  yk<-raster(paste("/Volumes/BG/UDmean_YK_",i,".grd",sep=""))
  ykea<-projectRaster(yk,r)
  kigb<-raster(paste("/Volumes/BG/UDmean_Kippaku_",i,".grd",sep=""))
  kigbea<-projectRaster(kigb,r)
  gi<-raster(paste("/Volumes/BG/UDmean_GannetIslands_",i,".grd",sep=""))
  giea<-projectRaster(gi,r)
  pli<-raster(paste("/Volumes/BG/UDmean_PrinceLeopoldIsland_",i,".grd",sep=""))
  pliea<-projectRaster(pli,r)
  co<-raster(paste("/Volumes/BG/UDmean_CoatsIslands_",i,".grd",sep=""))
  coea<-projectRaster(co,r)
  mi<-raster(paste("/Volumes/BG/UDmean_Minarets_",i,".grd",sep=""))
  miea<-projectRaster(mi,r)
  if(i!=12){ 
    akgb<-raster(paste("/Volumes/BG/UDmean_Alkefjellet_",i,".grd",sep=""))
    akgbea<-projectRaster(akgb,r)
  }
  bjgb<-raster(paste("/Volumes/BG/UDmean_Bjornoya_",i,".grd",sep=""))
  bjgbea<-projectRaster(bjgb,r)
  cggb<-raster(paste("/Volumes/BG/UDmean_CapeGorodetskiy_",i,".grd",sep=""))
  cggbea<-projectRaster(cggb,r)
  if(i!=12){
    fjlgb<-raster(paste("/Volumes/BG/UDmean_FranzJosefLand_",i,".grd",sep=""))
    fjlgbea<-projectRaster(fjlgb,r)
  }
  kagb<-raster(paste("/Volumes/BG/UDmean_KaraGate_",i,".grd",sep=""))
  kagbea<-projectRaster(kagb,r)
  di<-raster(paste("/Volumes/BG/UDmean_DiggesIslands_",i,".grd",sep=""))
  diea<-projectRaster(di,r)
  oi<-raster(paste("/Volumes/BG/UDmean_OranskieIslands_",i,".grd",sep=""))
  oiea<-projectRaster(oi,r)
  lasgb<-raster(paste("/Volumes/BG/UDmean_LanganesandSkjalfandi_",i,".grd",sep=""))
  lasgbea<-projectRaster(lasgb,r)
  lagb<-raster(paste("/Volumes/BG/UDmean_Latrabjarg_",i,".grd",sep=""))
  lagbea<-projectRaster(lagb,r)
  jmgb<-raster(paste("/Volumes/BG/UDmean_JanMayen_",i,".grd",sep=""))
  jmgbea<-projectRaster(jmgb,r)
  isgb<-raster(paste("/Volumes/BG/UDmean_Isfjorden_",i,".grd",sep=""))
  isgbea<-projectRaster(isgb,r)
  hogb<-raster(paste("/Volumes/BG/UDmean_Hornoya_",i,".grd",sep=""))
  hogbea<-projectRaster(hogb,r)
  grgb<-raster(paste("/Volumes/BG/UDmean_Grimsey_",i,".grd",sep=""))
  grgbea<-projectRaster(grgb,r)
  #Cyclone exposure index caluclation
  if(i==1){
    ind<-r*(khea+paea+flea+stea+siea+smea+guea+msea+kgea+kiea+skea+raaea+roea+lasea+iomea+isea+hoea+fjlea+fiea+kaea+ckea+bjea+anea+akea+
              bjgtea+figtea+cgea+grea+hjea+skgtea+hogtea+laea+lasgtea+jmea+
              saea+riea+psbea+ykea+kigbea+giea+pliea+coea+miea+akgbea+bjgbea+cggbea+fjlgbea+
              kagbea+diea+oiea+lasgbea+lagbea+jmgbea+isgbea+hogbea+grgbea)
    
  }else if(i==10){
    ind<-r*(khea+paea+flea+stea+siea+smea+guea+msea+kgea+kiea+skea+raaea+roea+lasea+iomea+isea+hoea+fjlea+fiea+kaea+ckea+bjea+anea+akea+
              bjgtea+figtea+cgea+grea+skgtea+hogtea+laea+lasgtea+jmea+
              saea+riea+ykea+kigbea+giea+pliea+coea+miea+bjgbea+cggbea+
              kagbea+diea+lasgbea+lagbea+jmgbea+isgbea+hogbea+grgbea)
  }else{
    ind<-r*(khea+paea+flea+stea+siea+smea+guea+msea+kgea+kiea+skea+raaea+roea+lasea+iomea+isea+hoea+fjlea+fiea+kaea+ckea+bjea+anea+akea+
              bjgtea+figtea+cgea+grea+skgtea+hjea+hogtea+laea+lasgtea+jmea+
              saea+riea+psbea+ykea+kigbea+giea+pliea+coea+miea+bjgbea+cggbea+
              kagbea+diea+lasgbea+lagbea+jmgbea+isgbea+hogbea+grgbea)
    
  }
  
  writeRaster(ind,paste("INDEX_mois",i,".tif"),overwrite=TRUE)

  
  if(i==1){
    index<-ind
  }else{
    index<-stack(index,ind)
  }
}
moyind<-mean(index)
writeRaster(moyind,"MOYINDEX.tif",overwrite=TRUE)

#Mapping
jpeg(filename ="MOYINDEX.jpeg", width = 3200, height = 2800, units = "px", res=300)
plot(moyind,col=plasma(256),legend=T,useRaster=FALSE,box=FALSE,axes=F)
plot(redcp,col="lightgrey",border=c("#999999"),add=T)
plot(gratp,col="#FFFFFF",border=c("#FFFFFF"),cex=0.2,add=T)
dev.off() 

####Graphical representation of Niche MapperTM results and statistical tests####
#NicheMapperTM calculation were made on the dedicated software (see Key Ressource table). The following script only corresponds to the graphical representation of the results and the statistical analysis to test the difference in energy requirements between cyclonic and non-cyclonic conditions.
#Example for little auk
re<-read.csv("/Volumes/Result.csv",sep=";") #outputs from NicheMapperTM
rem<-subset(re,re$Espece=="Mergule")# we select the output obtain for little auk
rem$Espece<-factor(rem$Espece)
# calculation of mean and standard deviation for each category and month
cy<- data.frame(matrix(NA,ncol=4,nrow=1))
colnames(cy)<-c("Mois","Categorie","MeanDep","SD") #Month, category, average, standard deviation

for(i in c("Janvier","Fevrier","Octobre","Novembre","Decembre")){
  cy$Mois<-i
  cy$Mois<-as.factor(cy$Mois)
  m<-subset(rem,rem$Mois==i)
  for(j in c("Baseline","1","2","3","4","BaselineINAC")){
    cy$Categorie<-j
    cl<-subset(m,m$Categorie==j)
    cy$MeanDep<-mean(cl$Dep,na.rm=TRUE)
    cy$SD<-sd(cl$Dep,na.rm=TRUE)
    
    if(j=="Baseline"){
      cltot<-cy
    }else{
      cltot<-rbind(cltot,cy)  
    }  
    
  }
  
  if(i=="Janvier"){
    annee<-cltot
  }else{
    annee<-rbind(annee,cltot)  
  }
}
write.table(annee,"MeanNMmergule.csv",sep=";",row.names=FALSE)
#graphical representation: Figure 3 and Data S1G-S1K
tot<-annee

tot$Month2 <- factor(tot$Mois, levels=unique(tot$Mois))
tot$Month2 <- factor(tot$Month2,levels(tot$Month2)[c(3,4,5,1,2)])
tot$Categorie <- as.factor(tot$Categorie)

p<-ggplot(tot, aes(x=Month2, y=MeanDep,color=Categorie,group=Categorie))+
  scale_colour_manual(name=c("Conditions"),values=c("#41BEC4","#1D91C0","#225EA8","#081D58","#7FCDBB","#C7E9B4"),labels=c("Non cyclonic without diving/flying activity","Non cyclonic","Class 1 cyclone", "Class 2 cyclone","Class 3 cyclone","Class 4 cyclone"),breaks=c("BaselineINAC","Baseline","1","2","3","4"))+ 
  geom_line(size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=MeanDep, ymax=MeanDep+SD),colour="black", width=.1,alpha=0.4,position=position_dodge(.05))+
  xlab("Month") +
  ylab("Energy requirements in kJ.day-1\n") +
  scale_x_discrete(labels=c("Oct","Nov","Dec","Jan","Feb"))+
  theme_classic() 


ggsave(filename ="MerguleNMfig3.jpeg",width = 20, height = 20, units = "cm")

#Test of the difference between energy requirements under cyclonic and non-cyclonic conditions
#Example for January 
j<-subset(rem,rem$Mois=="Janvier")
kruskal.test(j$Dep, j$Categorie) 

PT = dunnTest(Dep ~ Categorie,data=j,method="bh")  
RT<-PT$res
t<-cldList(P.adj ~ Comparison,data = RT,threshold = 0.05)
boxplot(Dep ~ Categorie, j
        ,xlab = "Class", ylab = "Value"
        , ylim = c(0,max(j$Dep,na.rm=T)*1.5)
        , notch = F, pch = ".")
mtext(side=3,text=t$Letter,at=1:length(t$Letter),cex=0.8)                          


