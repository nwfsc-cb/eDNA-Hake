## This script combines maps for 2019, 2021 and 2022


### Libraries
library(tidyverse)
library(knitr)
library(reshape2)
library(viridis)
library(ggmap)
library(maps)
library(mapdata)
library(sp)
library(gstat)
library(ggplot2)
library(sf)
library(brms)
library(ggsci)
library(gridExtra)
library(cowplot)

base.dir <- "/Users/ole.shelton/Github/eDNA-Hake/"
script.dir <- "/Users/ole.shelton/Github/eDNA-Hake/Scripts"
results.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Stan Model Fits"
plot.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Plots and figures"

# This calls both the qPCR results from 2019 and the acoustic data.
# define species

# load and run in the acoustic data.
setwd(script.dir)
# source("process acoustic data for qPCR.R")
# Read in agreed upon dat_raster_fin that has been trimmed to appropriate boundaries.
dat_raster_fin <- readRDS(file="../Data/_projection_rds/dat_raster_fin.rds") 
# Read in a few breaks based on latitude for summarizing the projections on different spatial scales.
load(file="../Data/lat_breaks_for_projections.RData")

# get bathymetry data.
# source("pull NOAA bathy for acoustic data.R")

# Read in the output from the Stan model
setwd(results.dir)

#load(paste("qPCR 2019",SPECIES, MOD, "7_12 Fitted.RData"))
load("qPCR 2019 hake lat.long.smooth 4_10_fix_nu_T-FIN Base_Var Fitted NO.SURFACE= FALSE .RData")

setwd(base.dir)
# Read in the 2021 and 2022 data and Deep Sea Data
dat.2022 <- read.csv("./Data/2021-22/2022_all_ctd_meta.csv")
dat.2021.ctd  <-   read.csv("./Data/2021-22/sh202106_all_ctd_meta.csv")
dat.2021.samp <-   read.csv("./Data/2021-22/2021 Hake eDNA samples.csv")
dat.deep <- read.csv("./Data/2021-22/Deep-sea_eDNA.csv", fileEncoding="latin1")

# Change name of the 2021 ctd data.
dat.2021.ctd <- dat.2021.ctd %>% rename(CTD.cast = Station.ID)

#Change the CTD names from having a X to having a x
dat.2021.ctd$CTD.cast <- sub('X','x',dat.2021.ctd$CTD.cast)
dat.2021 <- dat.2021.samp %>% distinct(CTD.cast) %>%
                left_join(.,dat.2021.ctd) %>% mutate(Year = 2021)

dat.2022 <- dat.2022  %>% mutate(Year = 2022) %>% rename(CTD.cast = Station.ID)

# Trim data to include only id, Year, and Lat-Lon
dat.2019 <- Output.qpcr$dat.id %>% filter(!is.na(station)) %>%
                distinct(station,lat,lon) %>% mutate(Year=2019) %>%
                dplyr::select(Year,CTD.cast=station,Latitude=lat,Longitude=lon)
dat.2021 <- dat.2021 %>% dplyr::select(Year,CTD.cast,Latitude,Longitude)
dat.2022 <- dat.2022 %>% dplyr::select(Year,CTD.cast,Latitude,Longitude)
dat.deep.ctd <- dat.deep %>% filter(ROV_CTD_Rosette == "CTD Rosette") %>%
                    mutate(type="CTD") %>% 
                    distinct(Year,Latitude,Longitude,type) %>% filter(!Latitude == "") %>%
                    mutate(Latitude = as.numeric(Latitude)) %>%
                    mutate(Longitude = as.numeric(Longitude))

dat.CTD <- dat.2019 %>% bind_rows(.,dat.2021) %>% bind_rows(.,dat.2022) %>%
            mutate(type="CTD",source="Hake") %>% 
            bind_rows(.,dat.deep.ctd %>% mutate(source="Deep Sea"))

dat.ROV <- dat.deep %>% filter(ROV_CTD_Rosette == "ROV") %>%
                mutate(type="ROV",source="Deep Sea") %>% 
                distinct(Year,Latitude,Longitude,type,source) %>% 
                filter(!Latitude == "") %>%
                mutate(Latitude = as.numeric(Latitude)) %>%
                mutate(Longitude = as.numeric(Longitude)) %>%
                filter(Latitude>30)

dat.all <- bind_rows(dat.CTD,dat.ROV)

lat <- summary(dat.all$Latitude)
lon <- summary(dat.all$Longitude)

setwd(script.dir)
source("Base_map.R")


eDNA.map.A <-  base_map_trim_2022 +
  scale_shape_manual("Source",values=c(17,16)) +
  geom_point(data=dat.all,aes(x=Longitude,y=Latitude,color=as.factor(Year),shape=source),alpha=0.75) +
  scale_color_viridis_d("Year",option = "plasma", begin=0,end=0.8)

eDNA.map.B <-  base_map_trim_2022 +
  scale_shape_manual("Source",values=c(16,17)) +
  geom_point(data=dat.all,aes(x=Longitude,y=Latitude,color=as.factor(Year),shape=type),alpha=0.75) +
  scale_color_viridis_d("Year",option = "plasma", begin=0,end=0.8)


setwd(plot.dir)
quartz(file="Sample_map_v1.pdf",type="pdf",dpi=600,height=9,width=5.5)
  print(eDNA.map.A)
dev.off()

setwd(plot.dir)
  quartz(file="Sample_map_v2.pdf",type="pdf",dpi=600,height=9,width=5.5)
  print(eDNA.map.B)
dev.off()


eDNA.map.C <-  base_map_trim_2022 +
  scale_shape_manual("Source",values=c(17,16)) +
  geom_point(data=dat.all,aes(x=Longitude,y=Latitude,color=as.factor(Year),shape=source),alpha=0.75) +
  scale_color_viridis_d("Year",option = "plasma", begin=0,end=0.8) +
  facet_wrap(~type)

setwd(plot.dir)
quartz(file="Sample_map_v3.pdf",type="pdf",dpi=600,height=9,width=8)
print(eDNA.map.C)
dev.off()

