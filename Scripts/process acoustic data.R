#### Process Acoustic Data.
# Libraries
library(tidyverse)
library(marmap)
library(ggplot2)
library(rstan)
library(lubridate)
library(reshape2)
library(geosphere)

# run post-process qPCR STAN Output.R to get base maps and some other items of interest.

# Working directories
base.dir <- "/Users/ole.shelton/Github/eDNA-Hake/"
data.dir <- paste0(base.dir,"Data/Biomass")
script.dir <- paste0(base.dir,"/Scripts")
plot.dir <- paste0(base.dir,"Plots and figures")

setwd(data.dir)
dat.acoustic <- read.csv("Biomass_2019 Hake Survey.csv")


### WORK WITH THE ACOUSTIC DATA.
dat.acoustic <- dat.acoustic %>% 
                      rename(transect=Transect,
                             lat=Latitude..deg.,
                             lon=Longitude..deg.,
                             mean.depth=Layer.Mean.Depth..m.,
                             layer.thickness=Layer.Thickness..m.,
                             abundance=Abundance,
                             biomass=Biomass..kg.)


max.lon.by.trans <- dat.acoustic %>% group_by(transect) %>% summarise(max.lon=max(lon),max.lat=max(lat)) %>% mutate(near.coast=1)


dat.acoustic <- dat.acoustic %>% left_join(.,max.lon.by.trans) %>% mutate(near.coast=ifelse(lon==max.lon,1,0))
# There are two transects which have two locations that correspond to max.lon. (121,127) both are up in Canada,
# So it doesn't really matter for our purposes.  But this is something to keep an eye on.
# I drop 121 and 127 here to make it obvious something is amiss with them. 
dat.acoustic <- dat.acoustic %>% filter(!transect %in% c(121,127))

TRANS <- unique(dat.acoustic$transect)

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
                                    depth.mean = mean(mean.depth),
                                    thickness.mean = mean(layer.thickness),
                                    dist.km = mean(dist.km),
                                    abund= sum(abundance),
                                    biomass = sum(biomass),
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







