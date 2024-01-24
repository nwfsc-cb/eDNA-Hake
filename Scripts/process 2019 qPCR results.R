# This script is for reading in qPCR data, merging it with locations, and starting a set of plots.  
# It will likely be called in a rMarkdown file later.

# Libraries
library(tidyverse)
library(marmap)
library(ggplot2)
library(rstan)
library(lubridate)
library(reshape2)
library(gridExtra)
library(raster)
#library(rgdal)
library(sp)
library(brms)
library(loo)

###########################################################################
## DECLARED SPECIES OF INTEREST
#SP <- "hake" # options: hake, lamprey, eulachon
###########################################################################

# Pull in qPCR data, qPCR standards, sample id information

dat.all <- read_csv(here('Data','qPCR','Hake eDNA 2019 qPCR results 2021-01-04 results.csv'))
dat.stand <- read_csv(here('Data','qPCR','Hake eDNA 2019 qPCR results 2020-01-04 standards.csv'))
dat.sample.id <- read_csv(here('Data','qPCR','Hake eDNA 2019 qPCR results 2023-02-10 sample details.csv'))

dat.station.id <- read_csv(here('Data','CTD_hake_data_10-2019.csv'))

# Pull in posterior for wash_offset derived from hake
# setwd(paste0(base.dir,"Stan Model Fits/"))
# wash_offset_hake <- read.csv("wash_offset_hake.csv") 
######################################################
# Drop 25m samples? 
  TRIM.25 <- TRUE

# PROCESS CTD LOCATION DATA FIRST
# Modify and clean CTD data 
  # fix lat-lon data
  dat.station.id <- dat.station.id %>% 
                      rename(station='Station #',lat='SciGPS-Lat',lon='SciGPS-Lon') %>%                
                      filter(Button == "CTD at Depth" | station == "49-9") %>% 
                      mutate(Lat=substr(lat,1,nchar(lat)-1),Lon=substr(lon,1,nchar(lon)-1)) %>%
                      mutate(lat.deg=as.numeric(substr(Lat,1,2)),lon.deg = as.numeric(substr(Lon,1,3)),
                             lat.dec=as.numeric(substr(Lat,3,nchar(Lat)))/60,lon.dec=as.numeric(substr(Lon,4,nchar(Lon)))/60) %>%
                      mutate(lat = lat.deg + lat.dec, lon= -(lon.deg + lon.dec)) %>% 
                      dplyr::select(-Lat,-Lon,-lat.deg,-lat.dec,-lon.deg,-lon.dec)
  
  # Fix CTD station data.
  dat.station.id <- dat.station.id %>% 
                            mutate(station = case_when(station=="77-0" ~ "77-MT505",
                                                       station=="78-0" ~ "78-MT506",
                                                       grepl("CTD ",station) ~ substr(station,5,nchar(station)),
                                                       grepl("CTD",station) ~ substr(station,4,nchar(station)),
                                                       TRUE ~ station) ) %>%
                            mutate(transect=substr(station,1,2),
                                   transect=case_when(grepl("-",transect) ~ substr(transect,1,1),
                                                      TRUE ~ transect))

  dat.station.id <- dat.station.id %>% mutate(date= as.Date(Date,"%m/%d/%y"),year=year(date),month=month(date),day=day(date)) %>% 
                              rename(water.depth='EK-38kHz-DepthBelowSurface-M-VALUE') 
  
  ### HERE ADD SCRIPT FOR PULLING DEPTH FROM bathy
  setwd(script.dir)
  source(here("Scripts","pull NOAA bathy for acoustic data.R"),local=T)
 
  temp <- marmap::get.depth(b,
                    x= dat.station.id %>% dplyr::select(lon,lat),
                    locator = FALSE)  
  dat.station.id$bathy.bottom.depth <-  -temp$depth
  
        d.temp <- dat.station.id %>% group_by(date,year,month,day,station, transect) %>%
                    summarise(m.lat=mean(lat),m.lon=mean(lon),
                              m.water.depth=mean(water.depth),bathy.depth = mean(bathy.bottom.depth)) %>%
                    rename(lat=m.lat,lon=m.lon, water.depth=m.water.depth, bathy.bottom.depth=bathy.depth)
  # come up with consensus bottom depth
        d.temp$bottom.depth.consensus <- d.temp$water.depth
        d.temp <- d.temp %>% mutate(
                            bottom.depth.consensus = ifelse(water.depth<100 & bathy.bottom.depth>1000,bathy.bottom.depth,bottom.depth.consensus),
                            bottom.depth.consensus = ifelse(water.depth<600 & water.depth>400 & bathy.bottom.depth>1000,bathy.bottom.depth,bottom.depth.consensus),
                            bottom.depth.consensus = ifelse(water.depth<1100 & water.depth>1000 & bathy.bottom.depth>2000,bathy.bottom.depth,bottom.depth.consensus),
                            bottom.depth.consensus = ifelse(water.depth<0,bathy.bottom.depth,bottom.depth.consensus),
                            bottom.depth.consensus = ifelse(water.depth>1000 & bathy.bottom.depth<200,bathy.bottom.depth,bottom.depth.consensus))

  dat.station.id.trim <- dat.station.id %>% dplyr::select(date,year,month,day, station, transect) %>%
                            # this combines the latitude and longitude to make a single concensus value for each Station.
                            left_join(d.temp,.)
  # get rid of a small number of duplicate stations in this data.frame.
  A <- dat.station.id.trim %>% group_by(station) %>% summarise(N=length(date)) %>% filter(N>1)
  THESE  <- A$station[A$station!=""]
  temp   <- dat.station.id.trim %>% filter(station %in% THESE) %>% arrange(station)
  #pull out every other row to drop duplicates.
  temp <- temp[seq(1,nrow(temp),by=2),]
  
  dat.station.id.trim <- dat.station.id.trim %>% filter(!station %in% THESE) %>% bind_rows(.,temp)
    
  ##
  dat.sample.id <- dat.sample.id %>% rename(Date.UTC='Date UTC',sample='Tube #')
  
  dat.sample.control.id <- dat.sample.id %>% mutate(date= as.Date(Date.UTC,"%m/%d/%y"),year=year(date),month=month(date),day=day(date)) %>%
                                  dplyr::select(sample,date,year,month,day,drop.sample,field.negative.type,volume = water.filtered.L,)
  dat.sample.id  <- dat.sample.id %>% dplyr::select(sample,
                                                  station='CTD cast',
                                                  Niskin,
                                                  depth,
                                                  drop.sample,
                                                  field.negative.type,
                                                  volume = water.filtered.L,
                                                  Fluor,
                                                  Zymo='Zymo columns')

  
  ##################################################3
  ##################################################3
  ##################################################3
  ##################################################3
  # ---- PULL IN SPATIAL DATA TO MAKE GRID ESTIMATE ON AND 
  # ---- TO PROJECT TO 
  
  str_name <- paste0(base.dir,"Data/raster_grid_blake/fivekm_grid.tif")
  dat_raster=raster(str_name)
  dat_raster_extracted <- rasterToPoints(dat_raster)
  
  # Get depth information.
  raster_depth <- read.csv(paste0(base.dir,"Data/raster_grid_blake/weighted_mean_NGDC_depths_for_5km_gridcells.csv"))
  raster_depth$depth_m <-  - raster_depth$WM_depth_m
  
  ####### 
  ### Get observations from the station IDs and convert to the new spatial coordinate systyem that Blake uses.
  #######
  # Use a projection derived by Blake.
    ## "+proj=laea +lat_0=30.5 +lon_0=-122.6 +x_0=1000000 +y_0=0 +datum=WGS84 +units=m +no_defs"
  PROJ.txt <- dat_raster@srs %>% as.character()
  proj      <- SpatialPointsDataFrame(coords = dat.station.id.trim %>% ungroup() %>% dplyr::select(lon,lat),
                                      data=dat.station.id.trim,
                                      proj4string = CRS("+proj=longlat"))
  proj.utm <- spTransform(proj, CRSobj = CRS(PROJ.txt))
  
  dat.utm <- (proj.utm@coords / 1000) %>% as.data.frame() %>% rename(utm.lon=coords.x1,utm.lat=coords.x2)
  
  dat.station.id.trim <- cbind(dat.station.id.trim %>% ungroup(),dat.utm)
  
  # calculate offshore distance from the acoustic survey start point reference 
  TRANS <- dat.station.id %>% distinct(year,transect)
  THESE <- c(dat.acoustic %>% filter(year==2019) %>% distinct(transect) %>% pull(transect))
  TRANS <- TRANS %>% filter(transect %in% THESE) 
  
  temp.all <- NULL
  for(i in 1:length(TRANS$transect)){
    #print(i)
    temp <- dat.station.id.trim %>% filter(transect == TRANS$transect[i])
    ref  <- dat.acoustic %>% filter(year==TRANS$year[i],transect == TRANS$transect[i],near.coast==1) %>% dplyr::select(max.lon,max.lat)
    D <- distGeo(p1=as.matrix(data.frame(lon=temp$lon,lat=temp$lat)),p2=ref)
    D <- D / 1000
    temp <- temp %>% mutate(dist.km = D) %>% arrange(dist.km) %>% mutate(id.numb = 1:nrow(temp))
    temp.all <- bind_rows(temp.all,temp)
  }
  
  # merge in the distances from the start of the transect.
  dat.station.id.trim <- left_join(dat.station.id.trim,temp.all %>% dplyr::select(station,transect,dist.km)) %>%
                                      rename(transect.dist.km=dist.km)
  
  ### Go get the projected points for the coast from Blake's data
  dat_raster_fin <- readRDS(file="../Data/_projection_rds/dat_raster_fin.rds")  
  

  ##################################################3
  ##################################################3  
    dat.id <- full_join(dat.sample.id,dat.station.id.trim %>% filter(!is.na(station))) %>% 
              mutate(station=ifelse(sample=="412","-",station)) %>%
              mutate(control = case_when(station=="" ~ "extraction",
                                        station=="-" ~ "field",
                                        TRUE ~ "no")) %>%
              mutate(depth=as.character(depth)) %>%
              mutate(depth = case_when(depth=="-" ~ "-99",
                                       depth=="sfc" ~ "0",
                                       depth=="300/150" ~ "300",
                                       depth=="" ~ "0",
                                       TRUE ~depth))

dat.id.control <- dat.id %>% filter(!control == "no") %>% dplyr::select(sample,volume,control) %>% left_join(.,dat.sample.control.id)
dat.id.samp    <- dat.id %>% filter(control == "no") %>% mutate(depth=as.numeric(as.character(depth)))

#write.csv(dat.id.samp, file = "Hake_2019_samples_w_CTD.csv")

##################################################3
##################################################3
# STANDARDS 
THESE <- c("qPCR","well","sample","type","IPC_Ct","inhibition_rate")
##################################################3##################################################
cut.plates <- c("H8")

dat.stand <- dat.stand %>% dplyr::select(all_of(THESE),grep(SP,colnames(dat.stand))) %>% 
                      rename(task=paste0(SP,"_task"),Ct=paste0(SP,"_Ct"),copies_ul = paste0(SP,"_copies_ul")) %>%
                      filter(task=="STANDARD",!qPCR%in%cut.plates) %>%
                      mutate(Ct=as.character(Ct),
                             Ct=case_when(Ct=="Undetermined"~"-99",
                                          TRUE~Ct)) %>%
                      mutate(Ct=as.numeric(Ct),Ct_bin=ifelse(Ct>0,1,0),log10_copies=log10(copies_ul))

    # Make sure there are sufficient samples for each qPCR
    stand.summary <- dat.stand %>% group_by(qPCR,copies_ul) %>% summarise(N= length(copies_ul)) %>% as.data.frame()
    PCR <- data.frame(qPCR=unique(dat.stand$qPCR),plate_idx= 1:length(unique(dat.stand$qPCR)))
    dat.stand <- left_join(dat.stand,PCR)
    
# SAMPLES, ntc and control
dat.samp <- dat.all %>% dplyr::select(all_of(THESE),useful,Zymo,dilution,grep(SP,colnames(dat.all))) %>%
                  filter(!qPCR%in%cut.plates) %>%
                  mutate(sample=as.character(sample),
                         sample=case_when(grepl("a",sample)~substr(sample,1,4),
                                                      TRUE ~ sample)) %>%
                  rename(Ct=paste0(SP,"_Ct"),copies_ul = paste0(SP,"_copies_ul")) %>%
                  dplyr::select(-grep("copies_l",colnames(.))) %>%
                  filter(!qPCR %in% cut.plates) %>%
                  mutate(Ct = as.character(Ct),
                         Ct=ifelse(Ct=="",-99,Ct),
                         Ct=ifelse(Ct=="Undetermined",-99,Ct),
                         Ct=ifelse(is.na(Ct)==T,-99,Ct), 
                         Ct=as.numeric(Ct)) %>%
                  mutate(IPC_Ct = as.character(IPC_Ct),
                        IPC_Ct=ifelse(IPC_Ct=="",-99,IPC_Ct),
                        IPC_Ct=ifelse(IPC_Ct=="Undetermined",-99,IPC_Ct),
                        IPC_Ct=ifelse(is.na(IPC_Ct)==T,-99,IPC_Ct), 
                        IPC_Ct=as.numeric(IPC_Ct)) %>%
                  filter(type=="unknowns")

dat.samp <- left_join(dat.samp,dat.id.samp,by="sample")

# Classify each depth into one of a few categories.
dat.samp <- dat.samp %>% mutate(depth_cat=case_when(depth < 25 ~ 0,
                                                          depth ==25 ~ 25,  
                                                          depth > 25  & depth <= 60  ~ 50,
                                                          depth > 60  & depth <= 100 ~ 100,
                                                          depth > 119 & depth <= 150 ~ 150,
                                                          depth > 151 & depth <= 200 ~ 200,
                                                          depth > 240 & depth <= 350 ~ 300,
                                                          depth > 400 & depth <= 500 ~ 500))

#### THIS IS A SWITCH.  IF NO.SURFACE == "TRUE", drop all surface samples.
# if(NO.SURFACE=="TRUE"){
#   dat.samp <- dat.samp %>% filter(depth_cat!=0)
# }



# Do some modifications to all of the samples first.
# Calculate the volume sampled for each sample
dat.samp <- dat.samp %>% mutate(vol.standard = volume/2.5)

# drop samples that hit the bench top
dat.samp <- dat.samp %>% filter(!drop.sample %in% c("Y","Y?"))#,!useful=="NO") 

dat.samp.begin <- dat.samp %>% group_by(sample) %>% summarise(N=length(sample)) %>% as.data.frame()

### CHECK THE SAMPLES THAT HAD TROUBLE WITH WASHING (drop.sample == "30EtOH" or "30EtOHpaired")
dat.wash <- dat.samp %>% filter(drop.sample %in% c("30EtOH","30EtOHpaired"))

  # find unique depth-station combinations among these stations.
  dat.wash$station <- as.factor(dat.wash$station)
  uni.wash <- dat.wash %>% group_by(station,depth,Niskin,drop.sample) %>% 
                summarise(N=length(station)) %>% mutate(status="washed") %>% ungroup()

  pairs.wash <- dat.samp %>% filter(!drop.sample %in% c("30EtOH","30EtOHpaired"))
  pairs.wash <- uni.wash %>% dplyr::select(station,depth) %>% left_join(.,pairs.wash) %>% 
                filter(is.na(Niskin)==F) %>% mutate(status="unwashed")

  dat.wash <- dat.wash %>% mutate(status="washed")

  dat.wash.all <- bind_rows(dat.wash,pairs.wash) %>% arrange(station,depth)

  dat.wash.summary <- dat.wash.all %>% group_by(station,depth,Niskin,status) %>% summarise(N=length(status)) %>% arrange(station,depth,status)
  dat.wash.summary <- dat.wash.summary %>% as.data.frame() 

  # There are 24 paired samples with which to estimate the effect of the 30% EtOH treatment
  print(nrow(dat.wash.summary[dat.wash.summary$status=="unwashed",]))

  # add indicator for membership in 30EtOH club and associated pairs
  # 0 = normal sample. 1 = washed with 30% etoh. 2= pair of washed with 30% EtOH sample
    dat.wash.all <- dat.wash.all %>% mutate(wash.indicator = ifelse(status == "washed",1,2)) %>%
                    dplyr::select(-status)
    # This indicator variable gets used in the STAN code.
    dat.wash.all <- dat.wash.all %>% mutate(wash_idx = ifelse(wash.indicator==1,1,0))

  # find samples that were washed with 30% EtOH, exclude them from dat.samp, 
  # then add them back in with needed indicator variables
  exclude  <- unique(dat.wash.all$sample)
  dat.samp <- dat.samp %>% mutate(wash.indicator=0,wash_idx=0) %>% filter(!sample %in% exclude) %>% 
              bind_rows(.,dat.wash.all)
  
# Separate out the non-template controls and the field controls.
dat.ntc <- dat.all %>% filter(type=="ntc",!qPCR %in% cut.plates) %>% mutate(IPC_Ct = as.numeric(as.character(IPC_Ct))) %>%
                      group_by(qPCR) %>% 
                      dplyr::summarise(mean.ntc = mean(IPC_Ct),sd.ntc=sd(IPC_Ct))

# get rid of one duplicate.
dat.sample.control.id <- dat.sample.control.id %>% 
                          distinct(sample,date,year,month,day,
                                   drop.sample, field.negative.type,volume)

dat.control <- dat.all %>% filter(grepl("neg",type)|type=="extc") %>%
                    dplyr::select(all_of(THESE),useful,Zymo,dilution,grep(SP,colnames(dat.all))) %>%
                    rename(Ct=paste0(SP,"_Ct"),copies_ul = paste0(SP,"_copies_ul")) %>%
                     mutate(Ct = as.character(Ct),
                              Ct=ifelse(Ct=="",-99,Ct),
                              Ct=ifelse(Ct=="Undetermined",-99,Ct),
                              Ct=ifelse(is.na(Ct)==T,-99,Ct), 
                              Ct=as.numeric(Ct)) %>%
                    mutate(IPC_Ct = as.character(IPC_Ct),
                              IPC_Ct=ifelse(IPC_Ct=="",-99,IPC_Ct),
                              IPC_Ct=ifelse(IPC_Ct=="Undetermined",-99,IPC_Ct),
                              IPC_Ct=ifelse(is.na(IPC_Ct)==T,-99,IPC_Ct), 
                              IPC_Ct=as.numeric(IPC_Ct))
                  
dat.control <- left_join(dat.control,dat.sample.control.id,by="sample") %>% 
                  filter(!grepl("Y",drop.sample)) # drop samples that hit the bench top

dat.control <- dat.control %>% mutate(vol.standard = volume/2.5)
dat.control.field.neg <-   dat.control %>% filter(type=="field_neg")

# Get rid of samples with dilution == 1 if a dilution series was run on a sample and those that were inhibited
dat.samp <- left_join(dat.samp,dat.ntc) %>% mutate(mean.ntc = as.numeric(as.character(mean.ntc))) %>%
                  mutate(inhibit.val = IPC_Ct-mean.ntc,
                         #inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT & inhibit.val > -INHIBIT.LIMIT,0,1),
                         inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT ,0,1),
                         Ct_bin = ifelse(Ct>0,1,0)) %>%
                  left_join(.,PCR) %>%
                  mutate(station.depth = paste0(station,".",depth))

these.samps <- dat.samp %>% dplyr::select(sample,dilution) %>% filter(dilution<1) %>% 
                dplyr::select(sample) %>%
                unique() %>% c(t(.))

dat.samp.dil <- dat.samp %>% filter(sample %in% these.samps,dilution<1)
dat.samp <- dat.samp %>% filter(!sample %in% these.samps) %>% bind_rows(.,dat.samp.dil)
    
### MANUALLY CUT A COUPLE OF SAMPLES THAT ARE KILLING THE MODEL FIT and after inspection are clearly causing troubles.

if(SP == "hake"){
  dat.samp <- dat.samp %>% filter(!sample %in% c(1298, # This one must have dropped on the tabletop
                                               1622,   # This one must have dropped on the tabletop 30EtOH too.
                                               #1552,   # Two orders of magnitude large than pair.
                                               #1419,   # Two orders of magnitude large than pair.
                                               1326,   # This is a 25m deep spot so gets dropped anyway.
                                               536,535    # Removed both of this pair() 
                                               #55,    # This is a 25m deep spot so gets dropped anyway.
                                               #1807 # THIS IS outlier among washed samples (dropped on tabletop)
  )) 
}
######

dat.control.field.neg <- left_join(dat.control.field.neg,dat.ntc) %>% mutate(mean.ntc = as.numeric(as.character(mean.ntc))) %>%
                          mutate(inhibit.val = IPC_Ct-mean.ntc,
                                  #inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT & inhibit.val > -INHIBIT.LIMIT,0,1),
                                 inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT ,0,1),
                                  Ct_bin = ifelse(Ct>0,1,0)) %>%
                          left_join(.,PCR)
# add wash_idx to field negatives
dat.control.field.neg <- dat.control.field.neg %>% 
                          mutate(wash_idx = ifelse(drop.sample=="30etOH",1,0))

# cull dat.samp of inhibited samples
dat.inhibit <- dat.samp %>% filter(useful=="NO" | inhibit.bin==1) 
dat.samp    <- dat.samp %>% filter(inhibit.bin==0) 
dat.control.field.neg    <- dat.control.field.neg %>% filter(inhibit.bin==0) 

####### THIS SECTION IS FOR EXPLORING SOME OF THE STRNAGENESS WITH
####### UNDERSTANDING THE EFFECT OF DILUTIONS.
dat.few <- dat.samp 
dat.few$copies_ul[dat.few$copies_ul == ""] <- NA
dat.few$copies_ul <- as.numeric(as.character(dat.few$copies_ul))
B <- dat.few %>% group_by(sample,dilution) %>% 
      summarise(count=mean(copies_ul,na.rm=T),N=sum(Ct_bin))
  
C <- B %>% dplyr::select(-N) %>% ungroup() %>% pivot_wider(.,names_from = c("dilution"),values_from = "count") %>% as.data.frame()
colnames(C)[2:5] <- c("x1","x.1","x.2","x.5")

# p.1 <- ggplot(C) +
#   geom_point(aes(x=x.1,y=x.2),alpha=0.5) +
#   geom_abline(intercept=c(0,0),slope = c(2,1),color="red",linetype=c("solid","dashed")) +
#   xlab("DNA copies/ul at 1:10 dilution") +
#   ylab("DNA copies/ul at 1:5 dilution")
# 
# p.2 <- ggplot(C) +
#   geom_point(aes(x=x.1,y=x.5),alpha=0.5) +
#   geom_abline(intercept=c(0,0),slope = c(5,1),color="red",linetype=c("solid","dashed")) +
#   xlab("DNA copies/ul at 1:10 dilution") +
#   ylab("DNA copies/ul at 1:2 dilution") #+
#   # xlim(c(0,120)) +
#   # ylim(c(0,120))
# 
# p.3 <- ggplot(C) +
#   geom_point(aes(x=x.2,y=x.5),alpha=0.5) +
#   geom_abline(intercept=c(0,0),slope = c(2.5,1),color="red",linetype=c("solid","dashed")) +
#   xlab("DNA copies/ul at 1:5 dilution") +
#   ylab("DNA copies/ul at 1:2 dilution") # +
#   # xlim(c(0,160)) +
#   # ylim(c(0,160))
# 
# dilution_plot <- grid.arrange(p.1,p.2,p.3,nrow=2)

# This strongly suggests that the 0.5 dilution did not work all that well (lower copies than expected)
# So if we exclude all of the 0.5 and inspect the results
A <- dat.samp %>% group_by(qPCR,dilution) %>% summarise(N=length(dilution))

if(SP == "lamprey"){ # ONLY KEEP 0.1 for lamprey
  dat.samp <- dat.samp %>% filter(!dilution %in% c(0.5,0.2))
}
if(SP == "eulachon"){ # ONLY KEEP 0.2 for Eulachon.
  dat.samp <- dat.samp %>% filter(!dilution %in% c(0.5,0.1))
}
if(SP == "hake"){ # ONLY KEEP 0.2 for hake
  dat.samp <- dat.samp %>% filter(!dilution %in% c(0.5,0.1))
}

B <- dat.samp %>% group_by(qPCR,dilution) %>% summarise(N=length(dilution))  %>% as.data.frame()

########## CHECK ON HOW MANY SAMPLES WE HAVE ZERO REPLICATES FOR  
# Examine how many samples are retained out of the intial number done.
dat.samp.ok <- dat.samp %>% group_by(sample,inhibit.bin) %>% summarise(N=length(sample)) %>% as.data.frame()

dat.samp.begin %>% nrow() #1844
dat.samp.ok %>% nrow() #1833 .... lost 3 samples through filtering, an extra 8 were excluded as outliers.

#dat.all %>% filter(sample %in% c(dat.samp.begin %>% filter(!sample %in% dat.samp.ok$sample) %>%pull(sample))) %>% as.data.frame()


### GET RID OF 25 m deep samples and 200m samples if doing smoothes by depth category.
if(TRIM.25 ==TRUE) {
  dat.samp <- dat.samp %>% filter(!depth_cat %in% c(25,200))
}


## make clear that these are 2019 specific objects
dat.samp <- dat.samp %>% dplyr::select(-useful,-Zymo.x,-Zymo.y,-Fluor,-water.depth)
COL <- names(dat.samp)
THESE <- c(which(COL=="year"),which(COL!="year"))
# reorder columns.
dat.samp.2019 <- dat.samp %>% dplyr::select(COL[THESE])

dat.stand.2019 <- dat.stand
PCR.2019 <- PCR
dat.control.field.neg.2019 <- dat.control.field.neg %>% dplyr::select(-useful,-Zymo)
dat.control.2019 <- dat.control %>% dplyr::select(-useful,-Zymo)

dat.inhibit <- dat.inhibit %>% dplyr::select(-useful,-Zymo.x,-Zymo.y,-Fluor,-water.depth)
COL <- names(dat.inhibit)
THESE <- c(which(COL=="year"),which(COL!="year"))
# reorder columns.
dat.inhibit.2019 <- dat.inhibit %>% dplyr::select(COL[THESE])



