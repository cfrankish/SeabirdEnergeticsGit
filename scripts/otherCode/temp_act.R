# make map 1
flight1<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeFlight)) +
   #geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent flight (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub)) 

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_flightMean.pdf"))
plot(flight1)
dev.off()

# SD #

flight2<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeFlight_sd)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent flight (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_flightSD.pdf"))
plot(flight2)
dev.off()

# make map 1
RestWater1<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeRestWater)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
  scale_fill_gradientn('Time spent rest (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #scale_fill_viridis_c() +
  #scale_fill_gradientn('Time spent RestWater (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_RestWaterMean.pdf"))
plot(RestWater1)
dev.off()

# SD #

RestWater2<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeRestWater_sd)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
   geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow",cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
   scale_fill_gradientn('Time spent rest (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_RestWaterSD.pdf"))
plot(RestWater2)
dev.off()

# make map 1
Land1<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeLand)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
   geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent Land (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
   coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_LandMean.pdf"))
plot(Land1)
dev.off()

# Mean behaviour #

Land2<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeLand_sd)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
   geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow",cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent Land (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
   coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_LandSD.pdf"))
plot(Land2)
dev.off()

if (energySub$species[1] %in% c("Black-legged kittiwake", "Northern fulmar")) {

# make map 1
Forage1<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeForage)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
   geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow",cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent foraging (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  labs(colour="") +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_ForageMean.pdf"))
plot(Forage1)
dev.off()

# Mean behaviour #

Forage2<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeForage_sd)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
   geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow",cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent foraging (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  labs(colour="") +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_ForageSD.pdf"))
plot(Forage2)
dev.off()

} else {

# make map 1
Active1<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeActive)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
   geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow",cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent Active (hours)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  #coord_fixed() +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  xlab("") +
  ylab("") +
  labs(colour="") +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_ActiveMean.pdf"))
plot(Active1)
dev.off()

# Mean behaviour #

Active2<-ggplot() +
  geom_tile(data=actSub, aes(x=x, y=y, fill=timeActive_sd)) +
   geom_text(data=actSub, aes(x=x, y=y, label=NoBirds), cex=1) +
   geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow",cex=1) +
  # scale_fill_gradient2('Non-breeding SST (°C)', low = "#364B9A", mid ="white" , high = "#c1121f", midpoint = 10, limits=c(-2.5, 26)) +
  scale_fill_gradientn('Time spent Active (SD)', colors=c("#364B9A", "#4393C3","#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#A50026")) +
  #geom_sf(data=allLocations_polygon, aes(color=species), fill=NA) +
  geom_sf(data=world, color = "#E5E5E5", fill = "#E5E5E5") +
  geom_sf(data=coast) +
  geom_point(data=colonySub, aes(x=coords.x1, y=coords.x2), color="yellow", cex=1) +
  coord_sf(crs=projection_NA, xlim=c(-3139249 - 2800000, 1492751 + 100000), ylim=c(-577185.7 - 2950000, 1962814.3 + 1050000)) +
  #coord_fixed() +
  xlab("") +
  ylab("") +
  labs(colour="") +
  facet_wrap(~month, nrow=4) +
  ggtitle(paste0(energySub$species[1], ":", modelColSub))

pdf(paste0("results/figures/popMaps/activity/", energySub$species[1], "_", modelColSub, "_ActiveSD.pdf"))
plot(Active2)
dev.off()

}