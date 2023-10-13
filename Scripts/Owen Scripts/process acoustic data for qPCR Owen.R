#### Process Acoustic Data.
# Libraries
library(tidyverse)
library(marmap)
library(ggplot2)
library(rstan)
library(lubridate)
library(reshape2)
library(geosphere)
library(raster)
library(rgdal)
library(sp)
library(brms)
library(loo)

# run post-process qPCR STAN Output.R to get base maps and some other items of interest.

# Working directories
base.dir <- here()
data.dir <- here('Data','acoustics 2019')
script.dir <- here('Scripts')
plot.dir <- here("Plots and figures")

setwd(data.dir)
dat.acoustic.raw <- read.csv("EchoPro_un-kriged_output-05-Dec-2019_0.csv")

### WORK WITH THE ACOUSTIC DATA.
dat.acoustic <- dat.acoustic.raw %>% 
                      rename(transect=Transect,
                             lat=Lat,
                             lon=Lon,
                             mean.depth.m=Layer.Depth,
                             layer.thickness.m=Layer.height,
                             abundance=ntk_total,
                             biomass_kg=wgt_total,
                             abund_dens_n_nm2 = nntk_total,
                             weight_dens_kg_nm2 = nwgt_total
                             ) %>%
                      dplyr::select(transect,
                                    VL_start,
                                    VL_end,
                                    lat,
                                    lon,
                                    mean.depth.m,
                                    layer.thickness.m,
                                    abundance,
                                    biomass_kg,
                                    abund_dens_n_nm2 ,
                                    weight_dens_kg_nm2,
                                    sig_b,
                                    NASC) %>% 
                      arrange(transect,VL_start)

n_samps_acoustic <- dat.acoustic %>% group_by(transect) %>% summarise(N=length(transect))
 
#Make ID for each observation.
samp_id <- NULL
for(i in 1:nrow(n_samps_acoustic)){
  samp_id <- c(samp_id,1:n_samps_acoustic$N[i])
}
dat.acoustic$ID <- samp_id
dat.acoustic$trans_id <- paste0(dat.acoustic$transect,"_",dat.acoustic$ID)
####

max.lon.by.trans <- dat.acoustic %>% group_by(transect) %>% summarise(max.lon=max(lon),max.lat=max(lat)) %>% mutate(near.coast=1)


dat.acoustic <- dat.acoustic %>% left_join(.,max.lon.by.trans) %>% mutate(near.coast=ifelse(lon==max.lon,1,0))
# There are two transects which have two locations that correspond to max.lon. (121,127) both are up in Canada,
# So it doesn't really matter for our purposes.  But this is something to keep an eye on.
# I drop 121 and 127 here to make it obvious something is amiss with them. 
dat.acoustic <- dat.acoustic %>% filter(!transect %in% c(121,127))

TRANS <- unique(dat.acoustic$transect)

### Pull in Blake's 5km grid and projection, convert lat-lon to that coordinate system.

str_name <- here("Data/raster_grid_blake/fivekm_grid.tif")
dat_raster=raster(str_name)
dat_raster_extracted <- rasterToPoints(dat_raster)

# Get depth information.
raster_depth <- read.csv(here("Data/raster_grid_blake/weighted_mean_NGDC_depths_for_5km_gridcells.csv"))
raster_depth$depth_m <-  - raster_depth$WM_depth_m

####### 
### Get observations from the station IDs and convert to the new spatial coordinate systyem that Blake uses.
#######
# Use a projection derived by Blake.
## "+proj=laea +lat_0=30.5 +lon_0=-122.6 +x_0=1000000 +y_0=0 +datum=WGS84 +units=m +no_defs"
PROJ.txt <- crs(dat_raster)
proj      <- SpatialPointsDataFrame(coords = dat.acoustic %>% ungroup() %>% dplyr::select(lon,lat),
                                    data=dat.acoustic,
                                    proj4string = CRS("+proj=longlat"))
proj.utm <- spTransform(proj, CRSobj = CRS(as.character(PROJ.txt)))

dat.utm <- (proj.utm@coords / 1000) %>% as.data.frame() %>% rename(utm.lon=lon,utm.lat=lat)

dat.acoustic <- cbind(dat.acoustic %>% ungroup(),dat.utm)

temp.all <- NULL
for(i in 1:length(TRANS)){
  #print(i)
  temp <- dat.acoustic %>% filter(transect == TRANS[i])
  ref  <- temp %>% filter(near.coast==1) %>% dplyr::select(max.lon,max.lat)
  D <- distGeo(p1=as.matrix(data.frame(lon=temp$lon,lat=temp$lat)),p2=ref)
  D <- D / 1000
  temp <- temp %>% mutate(dist.km = D) %>% arrange(dist.km) %>% mutate(id.numb = 1:nrow(temp))
  temp.all <- bind_rows(temp.all,temp)
}

dat.acoustic <- left_join(dat.acoustic,temp.all)
dat.acoustic$biomass_mt <- dat.acoustic$biomass/1000
dat.acoustic$log10_biomass_mt <- log10(dat.acoustic$biomass_mt)
 

# Combine the data in a way to make plotting more friendly.
dat.acoustic.binned <- NULL

BY <- 4
for(i in 1:length(TRANS)){
  temp <- dat.acoustic %>% filter(transect==TRANS[i]) %>% mutate(bin.id=as.integer(0))
  MAX <- temp %>% dplyr::select(id.numb) %>% max(.)
  SEQ <- seq(BY,MAX+BY,by=BY)
  SEQ.id <- 1:length(SEQ)
  
  for(j in 1:length(SEQ.id)){
    temp <- temp %>% mutate(bin.id = case_when(id.numb >= SEQ[j]-(BY-1) & id.numb <= SEQ[j] ~ SEQ.id[j],
                              TRUE ~ bin.id ))
  }
  
  temp <- temp %>% group_by(transect,bin.id) %>% 
                          summarise(lat = mean(lat),
                                    lon = mean(lon),
                                    depth.mean = mean(mean.depth.m),
                                    thickness.mean = mean(layer.thickness.m),
                                    dist.km = mean(dist.km),
                                    abund= sum(abundance),
                                    biomass_kg = sum(biomass_kg),
                                    biomass_mt = sum(biomass_mt),
                                    mean.id.numb= mean(id.numb),
                                    near.coast=ifelse(sum(near.coast,na.rm=T)==1,1,0))
  
  dat.acoustic.binned <- bind_rows(dat.acoustic.binned,temp)
}

########### MAKE SOME PLOTS OF THIS DATA.
# Latitude
lat.lims <- c(34.39501,48.56390)
lon.lims <- c(-126.5501, -120.6928)

lat.lims.trim <- c(37.5,48.56390)
lon.lims.trim <- c(-126.5501,-122.5)

# Call Base_map.R
setwd(script.dir)
source("Base_map.R",local=T)
  
# Make some basic plots of lat-longs of biomass
lower.lim = 1
BREAK = 10^(1:5)
LAB= 10^(1:5)
hake.acoustic.p1 <- base_map_trim + 
  geom_point(data=dat.acoustic,aes(x=lon,y=lat,size=biomass_mt),shape=21,color="blue") +
  scale_size("log10(Biomass)",labels=LAB,breaks=BREAK,range=c(0.01,10),limits=c(lower.lim,NA))+
  geom_point(data=dat.acoustic %>% filter(is.infinite(log10_biomass_mt)),aes(x=lon,y=lat),shape="x",size=0.5,color="black") +
  scale_shape("log10(Biomass)",solid=FALSE)
  
print(hake.acoustic.p1)


# Binned plot
lower.lim = 1
BREAK = 10^(1:5)
LAB= 10^(1:5)
hake.acoustic.binned.p1 <- base_map_trim + 
  geom_point(data=dat.acoustic.binned,aes(x=lon,y=lat,size=biomass_mt),shape=21,color="blue") +
  scale_size("Biomass (mt)",labels=LAB,breaks=BREAK,range=c(0.01,10),limits=c(lower.lim,NA))+
  geom_point(data=dat.acoustic %>% filter(is.infinite(log10_biomass_mt)),aes(x=lon,y=lat),shape="x",size=0.5,color="black") +
  scale_shape("Biomass (mt)",solid=FALSE)

print(hake.acoustic.binned.p1)


##################################################3
#### ---------------------------------------------
#### Process the projection matrix to include areas that are surveyed by the acoustic surveys.
#### ---------------------------------------------
##################################################3
# ---- SPTAIAL DATA IS PULLED IN ABOVE

dat_raster_extracted %>% head()


# Make a set of inshore and offshore boundaries based on dat.acoustic


acoustic.lims <- dat.acoustic %>% group_by(transect) %>% 
                summarise(min.lon.utm = min(utm.lon), 
                           max.lon.utm = max(utm.lon),
                           mean.lat.utm = mean(utm.lat)) %>% 
                  arrange(transect)

ggplot(acoustic.lims) + 
    geom_point(aes(x=min.lon.utm,y=mean.lat.utm)) +
    geom_point(aes(x=max.lon.utm,y=mean.lat.utm),col=2)



dat_raster_trim <- NULL
dat_raster_extracted <- dat_raster_extracted %>% as.data.frame()
dat_raster_extracted$y.km <- dat_raster_extracted$y/1000
dat_raster_extracted$x.km <- dat_raster_extracted$x/1000
for(i in 2:nrow(acoustic.lims)){
  filt.temp <-  acoustic.lims[(i-1):i,] %>% arrange(mean.lat.utm)
  temp <- dat_raster_extracted %>% filter(y.km >= filt.temp$mean.lat.utm[1], 
                                          y.km < filt.temp$mean.lat.utm[2])
  uni.lat <- unique(temp$y.km) 
  
  SLOPE.min <-    (filt.temp$min.lon.utm[2] - filt.temp$min.lon.utm[1]) /
                  (filt.temp$mean.lat.utm[2] - filt.temp$mean.lat.utm[1])
  SLOPE.max <-    (filt.temp$max.lon.utm[2] - filt.temp$max.lon.utm[1]) /
                  (filt.temp$mean.lat.utm[2] - filt.temp$mean.lat.utm[1])
  
  uni.lat <- uni.lat %>% as.data.frame() %>% rename(y.km = ".") %>% 
              mutate(lon.min = SLOPE.min * (y.km -  filt.temp$mean.lat.utm[1]) + filt.temp$min.lon.utm[1]) %>%
              mutate(lon.max = SLOPE.max * (y.km -  filt.temp$mean.lat.utm[1]) + filt.temp$max.lon.utm[1])
  
  temp <- left_join(temp,uni.lat) %>% filter(x.km > lon.min,x.km < lon.max)
  dat_raster_trim <- bind_rows(dat_raster_trim,temp)
}

ggplot(acoustic.lims) + 
  geom_point(aes(x=min.lon.utm,y=mean.lat.utm)) +
  geom_point(aes(x=max.lon.utm,y=mean.lat.utm),col=2) +
  geom_point(data=dat_raster_trim, aes(x=x.km,y=y.km),col="blue") 
 
# Trim everything south of SF Bay, add in the depth information
dat_raster_fin <- dat_raster_trim %>% dplyr::select(fivekm_grid,y,x) %>% 
                      rename(Gridcell_ID=fivekm_grid)
dat_raster_fin <- left_join(dat_raster_fin,raster_depth) %>% filter(is.na(depth_m)==F,depth_m>30)

# Add lat and lon to dat_raster_fin for help with plotting later
proj      <- SpatialPointsDataFrame(coords = dat_raster_fin %>% ungroup() %>% dplyr::select(x,y),
                                    data = dat_raster_fin,
                                    proj4string = crs(PROJ.txt))
proj.latlon <- spTransform(proj, CRSobj = CRS("+proj=longlat"))

dat_raster_fin <- left_join(dat_raster_fin,
                            data.frame(Gridcell_ID= proj.latlon$Gridcell_ID,
                                       data.frame(proj.latlon@coords)%>% rename(lon=x,lat=y)))
dat_raster_fin <- dat_raster_fin %>% mutate(x= x/1000, y=y/1000)


###### TRIM THE PROJECTION POINTS TO CONTRAIN THE DISTRIBUTION EXAMINED LATER.

# This trims the points to include only points:
# North of transect 25 
# South of transect 85
# East of longitude -126.5

min.lat <- 38.3
max.lat <- 48.5
min.lon <- -126.5

dat_raster_fin <- dat_raster_fin %>% filter(lat > min.lat,lat < max.lat,lon > min.lon)


ggplot(dat_raster_fin) +
  geom_tile(aes(x=x,y=y,fill=depth_m),alpha=0.8)  +
  theme_bw()
  #geom_point(data=dat.station.id.trim,aes(x=utm.lon,y=utm.lat),col="red") 

# Save this data frame to file so the PCR and acoustics script can use it.                  
saveRDS(dat_raster_fin,file="../Data/_projection_rds/dat_raster_fin.rds")   

# Make sets of lines to divide up the coast into approximately equal sized chunks
lats.equal <- seq(38.25,max.lat,length.out = 11)
lats.rounded.0.5 <- c(38.25,seq(39,max.lat,by=0.5))
lats.rounded.1.0 <- c(38.25,seq(39,max.lat,by=1.0))


# Sort from south to north
max.lon.for.groups <- rep(0,length(lats.equal))
max.lats <- max.lon.for.groups
for(i in 1: length(max.lon.for.groups)){
  if(i==1){
    max.lon.for.groups[i] <- dat_raster_fin$lon %>% max()
  }else{
    max.lats[i-1] = lats.equal[i]
    max.lon.for.groups[i]  <- dat_raster_fin %>% filter(lat< lats.equal[i],
                                                        lat> (lats.equal[i]-0.1)) %>% 
                              dplyr::select(lon) %>% max()
  }
  max.lats[length(max.lon.for.groups)] = lats.equal[length(max.lon.for.groups)] + 1
}
lats.equal <- bind_cols(lat=c(lats.equal),lat.max = max.lats,lon.max=max.lon.for.groups,) %>% 
              mutate(lon.min=lon.max-1.75,ID = 1:length(lon.min))
lats.equal <- lats.equal %>% filter(lat.max<= max.lat)

# Sort from south to north
max.lon.for.groups <- rep(0,length(lats.rounded.0.5))
max.lats <- max.lon.for.groups
for(i in 1: length(max.lon.for.groups)){
  if(i==1){
    max.lon.for.groups[i] <- dat_raster_fin$lon %>% max()
  }else{
    max.lats[i-1] = lats.rounded.0.5[i]
    max.lon.for.groups[i]  <- dat_raster_fin %>% filter(lat< lats.rounded.0.5[i],
                                                        lat> (lats.rounded.0.5[i]-0.1)) %>% 
      dplyr::select(lon) %>% max()
  }
  max.lats[length(max.lon.for.groups)] = lats.rounded.0.5[length(max.lon.for.groups)] + 0.5
}
lats.rounded.0.5 <- bind_cols(lat=c(lats.rounded.0.5),lat.max = max.lats,lon.max=max.lon.for.groups) %>% 
                        mutate(lon.min=lon.max-1.75,ID = 1:length(lon.min))
lats.rounded.0.5 <- lats.rounded.0.5 %>% filter(lat.max <= max.lat)

# Sort from south to north
max.lon.for.groups <- rep(0,length(lats.rounded.1.0))
max.lats <- max.lon.for.groups
for(i in 1: length(max.lon.for.groups)){
  if(i==1){
    max.lon.for.groups[i] <- dat_raster_fin$lon %>% max()
  }else{
    max.lats[i-1] = lats.rounded.1.0[i]
    max.lon.for.groups[i]  <- dat_raster_fin %>% filter(lat< lats.rounded.1.0[i],
                                                        lat> (lats.rounded.1.0[i]-0.1)) %>% 
      dplyr::select(lon) %>% max()
  }
  max.lats[length(max.lon.for.groups)] = lats.rounded.1.0[length(max.lon.for.groups)] + 1
  
}
lats.rounded.1.0 <- bind_cols(lat=c(lats.rounded.1.0),lat.max = max.lats,lon.max=max.lon.for.groups) %>% 
                        mutate(lon.min=lon.max-1.75,
                               ID = 1:length(lon.min))


# latitudinal breaks for analysis
lat.breaks <- list(
              lats.equal = lats.equal,
              lats.rounded.0.5 = lats.rounded.0.5,
              lats.rounded.1.0 = lats.rounded.1.0)

save(lat.breaks,file="../Data/lat_breaks_for_projections.RData")

