
# Pull in qPCR data, qPCR standards, sample id information
dat.all <- read_csv(here('Data','qPCR','Hake eDNA 2021 qPCR results 10.13.2023.csv'))
dat.stand <- read_csv(here('Data','qPCR','Hake eDNA 2021 qPCR standards 10.13.2023.csv'))
dat.sample.id <- read_csv(here('Data','qPCR','Hake eDNA 2021 qPCR sample details 10.16.2023.csv'),skip=2)

dat.station.id <- read_csv(here('Data','CTD meta 2021.csv'))

# load and run the acoustic data. this is needed to reference the offshore-ness of 
setwd(script.dir)
source("process acoustic data for qPCR 2019 on.R",local=T)
# dat.acoustic and dat.acoustic.binned are the relevant data frames

# Pull in posterior for wash_offset derived from hake
# setwd(paste0(base.dir,"Stan Model Fits/"))
# wash_offset_hake <- read.csv("wash_offset_hake.csv") 
######################################################
if(MODEL.TYPE == "Base"){
  TRIM.25 <- FALSE
} 
if(MODEL.TYPE == "lat.long.smooth"|MODEL.TYPE == "lat.long.smooth.base"){
  TRIM.25 <- TRUE
}  

# Merge in ntc to data to flag inhibited samples 
INHIBIT.LIMIT <- 0.5

# Fix naming differences between station data.
# make all stations that start with a capital X, start with a lowercase x
dat.station.id <- dat.station.id %>% rename(station='Station ID',lat=Latitude,lon=Longitude) %>%
                      filter(grepl('X|x',station)) %>% 
                      mutate(station=sub('X','x',station)) %>%
                      # Extract transect from station number.
                      mutate(transect = gsub("\\_.*","", station) ) %>%
                      mutate(transect = gsub(".*x","", transect) )
# Fix dates.
dat.station.id <- dat.station.id %>% 
                    mutate(date= as.Date(UTCDate,"%m/%d/%y"),year=year(date),month=month(date),day=day(date)) #%>% 
                    #rename(water.depth=EK.38kHz.DepthBelowSurface.M.VALUE) 

### HERE ADD SCRIPT FOR PULLING DEPTH FROM bathy
setwd(script.dir)
source("pull NOAA bathy for acoustic data.R",local=T)

temp <- marmap::get.depth(b,
                          x= dat.station.id %>% dplyr::select(lon,lat),
                          locator = FALSE)  
dat.station.id$bathy.bottom.depth <-  -temp$depth

d.temp <- dat.station.id %>% group_by(date,year,month,day,station, transect) %>%
  summarise(m.lat=mean(lat),m.lon=mean(lon),
            #m.water.depth=mean(water.depth),
            bathy.depth = mean(bathy.bottom.depth)) %>%
  rename(lat=m.lat,lon=m.lon, 
         #water.depth=m.water.depth, 
         bathy.bottom.depth=bathy.depth)
# come up with consensus bottom depth (vestigal from 2019 analysis)
# 2021 depth comes entirely from NOAA.
d.temp$bottom.depth.consensus <- d.temp$bathy.bottom.depth
# d.temp <- d.temp %>% mutate(
#   bottom.depth.consensus = ifelse(water.depth<100 & bathy.bottom.depth>1000,bathy.bottom.depth,bottom.depth.consensus),
#   bottom.depth.consensus = ifelse(water.depth<600 & water.depth>400 & bathy.bottom.depth>1000,bathy.bottom.depth,bottom.depth.consensus),
#   bottom.depth.consensus = ifelse(water.depth<1100 & water.depth>1000 & bathy.bottom.depth>2000,bathy.bottom.depth,bottom.depth.consensus),
#   bottom.depth.consensus = ifelse(water.depth<0,bathy.bottom.depth,bottom.depth.consensus),
#   bottom.depth.consensus = ifelse(water.depth>1000 & bathy.bottom.depth<200,bathy.bottom.depth,bottom.depth.consensus))

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

## WORK WITH DAT.SAMPLE.ID
dat.sample.control.id <- dat.sample.id %>% rename(Date.UTC ='Date UTC') %>% 
                        mutate(date= as.Date(Date.UTC,"%m/%d/%y"),year=year(date),month=month(date),day=day(date)) %>%
                        dplyr::select(sample='Tube #',date,year,month,day,drop.sample,field.negative.type,volume = water.filtered.L)
dat.sample.id  <- dat.sample.id %>% mutate(year= 2021) %>%
                        dplyr::select(year,
                                                  sample='Tube #',
                                                  station='CTD cast',
                                                  Niskin,
                                                  depth,
                                                  drop.sample,
                                                  field.negative.type,
                                                  volume = water.filtered.L)
                                                  #Fluor,
                                                  #Zymo=Zymo.columns)
# replace any CTD cast with capital X name = x
dat.sample.id <- dat.sample.id %>% mutate(station=sub('X','x',station))

##################################################3
##################################################3
##################################################3
##################################################3
# ---- PULL IN SPATIAL DATA TO MAKE GRID ESTIMATE ON AND 
# ---- TO PROJECT TO 

str_name <- paste0(base.dir,"/Data/raster_grid_blake/fivekm_grid.tif")
dat_raster=raster(str_name)
dat_raster_extracted <- rasterToPoints(dat_raster)

# Get depth information.
raster_depth <- read.csv(paste0(base.dir,"/Data/raster_grid_blake/weighted_mean_NGDC_depths_for_5km_gridcells.csv"))
raster_depth$depth_m <-  - raster_depth$WM_depth_m

####### 
### Get observations from the station IDs and convert to the new spatial coordinate systyem that Blake uses.
#######
# Use a projection derived by Blake.
## "+proj=laea +lat_0=30.5 +lon_0=-122.6 +x_0=1000000 +y_0=0 +datum=WGS84 +units=m +no_defs"
PROJ.txt <- crs(dat_raster)
proj      <- SpatialPointsDataFrame(coords = dat.station.id.trim %>% ungroup() %>% dplyr::select(lon,lat),
                                    data=dat.station.id.trim,
                                    proj4string = CRS("+proj=longlat"))
proj.utm <- spTransform(proj, CRSobj = CRS(as.character(PROJ.txt)))

dat.utm <- (proj.utm@coords / 1000) %>% as.data.frame() %>% rename(utm.lon=lon,utm.lat=lat)

dat.station.id.trim <- cbind(dat.station.id.trim %>% ungroup(),dat.utm)

# calculate offshore distance from the acoustic survey start point reference 
TRANS <- dat.station.id %>% distinct(year,transect)
THESE <- c(dat.acoustic %>% filter(year==2021) %>% distinct(transect) %>% pull(transect))
TRANS <- TRANS %>% filter(transect %in% THESE) 

temp.all <- NULL
for(i in 1:length(TRANS$transect)){
  #print(i)
  temp <- dat.station.id.trim %>% filter(transect == TRANS$transect[i])
  ref  <- dat.acoustic %>% filter(year==TRANS$year[i], transect == TRANS$transect[i],near.coast==1) %>% dplyr::select(max.lon,max.lat)
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
dat.id <- full_join(dat.sample.id,dat.station.id.trim) %>% 
  #mutate(station=ifelse(sample=="412","-",station)) %>%
  mutate(control = case_when(station=="N/A" ~ "control",
                             station=="-" ~ "extraction",
                             TRUE ~ "no")) %>%
  mutate(depth=as.character(depth)) %>%
  mutate(depth = case_when(depth=="-" ~ "-99",   # extraction controls.
                           depth=="N/A" ~ "-99", # these are all the various controls and aerial samples.
                           is.na(depth) ~ "-99", # This is a single station that wasn't sampled for eDNA.
                           depth=="sfc" ~ "0",
                           TRUE ~depth))

dat.id.control <- dat.id %>% filter(!control == "no") %>% dplyr::select(sample,volume,control) %>% left_join(.,dat.sample.control.id)
dat.id.samp    <- dat.id %>% filter(control == "no") %>% mutate(depth=as.numeric(as.character(depth)))

#write.csv(dat.id.samp, file = "Hake_2019_samples_w_CTD.csv")

####################################################
###  END MAKING DAT.ID.
####################################################

##################################################3
##################################################3
# STANDARDS 
THESE <- c("qPCR","well","sample","type","IPC_Ct","inhibition_rate")
##################################################3##################################################
cut.plates <- c("H21_S01") # Cut the Sadie only qPCR

dat.stand <- dat.stand %>% dplyr::select(all_of(THESE),grep(SP,colnames(dat.stand))) %>% 
  mutate(task="STANDARD") %>%
  rename(Ct=paste0(SP,"_Ct"),copies_ul = paste0(SP,"_copies_ul")) %>%
  filter(task=="STANDARD") %>%
  filter(!qPCR%in%cut.plates) %>%
  mutate(Ct=as.character(Ct),
         Ct=case_when(Ct=="Undetermined"~"-99",
                      TRUE~Ct)) %>%
  mutate(Ct=as.numeric(Ct),Ct_bin=ifelse(Ct>0,1,0),log10_copies=log10(copies_ul))

# Make sure there are sufficient samples for each qPCR
stand.summary <- dat.stand %>% group_by(qPCR,copies_ul) %>% summarise(N= length(copies_ul)) %>% as.data.frame()
PCR <- data.frame(qPCR=unique(dat.stand$qPCR),plate_idx= 1:length(unique(dat.stand$qPCR)))
dat.stand <- left_join(dat.stand,PCR)

# SAMPLES, ntc and control
dat.samp <- dat.all %>% dplyr::select(all_of(THESE),dilution,grep(SP,colnames(dat.all))) %>%
  filter(!qPCR %in% cut.plates) %>% # Get rid of SADIE samples.
  filter(!grepl("S",sample)) %>% # Get rid of SADIE samples (some are mixed in with the basic samples.)
  # mutate(sample=as.character(sample),
  #        sample=case_when(grepl("a",sample)~substr(sample,1,4),
  #                         TRUE ~ sample)) %>%
  rename(Ct=paste0(SP,"_Ct"),copies_ul = paste0(SP,"_copies_ul")) %>%
  dplyr::select(-grep("copies_l",colnames(.))) %>%
  
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


# Filter out samples that are at expected depth 
dpths <- c(0,50,150,200,300,500)
A <- dat.samp %>% filter(!depth %in% dpths)
# manually look through these samples 

# Classify each depth into one of a few categories.
# Categories should be 0, 50, 150,200,300,500
dat.samp <- dat.samp %>% mutate(depth_cat=case_when(depth < 25 ~ 0,
                                                    depth > 25  & depth <= 50  ~ 50,
                                                    depth > 60  & depth <= 160 ~ 150,
                                                    depth > 160 & depth <= 200 ~ 200,
                                                    depth > 225 & depth <= 339 ~ 300,
                                                    depth > 340 & depth <= 500 ~ 500))

#### THIS IS A SWITCH.  IF NO.SURFACE == "TRUE", drop all surface samples.
if(NO.SURFACE=="TRUE"){
  dat.samp <- dat.samp %>% filter(depth_cat!=0)
}

# Do some modifications to all of the samples first.

# drop samples that hit the bench top
dat.samp <- dat.samp %>% filter(!grepl("Y",drop.sample))#,!useful=="NO") 

# Calculate the volume sampled for each sample
dat.samp <- dat.samp %>% mutate(vol.standard = as.numeric(volume)/2.5)

dat.samp.begin <- dat.samp %>% group_by(sample) %>% summarise(N=length(sample)) %>% as.data.frame()

### Add field to make these samples comparable to 2019. Add field for no wash error in 2021
dat.samp <- dat.samp %>% mutate(wash.indicator=0,wash_idx=0) 

# Separate out the non-template controls and the field controls.
dat.ntc <- dat.all %>% filter(type=="ntc") %>% mutate(IPC_Ct = as.numeric(as.character(IPC_Ct))) %>%
  group_by(qPCR) %>% 
  dplyr::summarise(mean.ntc = mean(IPC_Ct),sd.ntc=sd(IPC_Ct))

dat.control <- dat.all %>% filter(grepl("extc",type)|type=="extc"|grepl("neg",type)) %>%
  dplyr::select(all_of(THESE),dilution,grep(SP,colnames(dat.all))) %>%
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

dat.control <- dat.control %>% mutate(volume = ifelse(volume=="N/A",NA,volume)) %>%  
                mutate(vol.standard = as.numeric(volume)/2.5)
dat.control.field.neg <-   dat.control %>% filter(type=="field_neg")


# Get rid of samples with dilution == 1 if a dilution series was run on a sample and those that were inhibited
dat.samp <- left_join(dat.samp,dat.ntc) %>% mutate(mean.ntc = as.numeric(as.character(mean.ntc))) %>%
  mutate(inhibit.val = IPC_Ct-mean.ntc,
         #inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT & inhibit.val > -INHIBIT.LIMIT,0,1),
         inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT ,0,1),
         Ct_bin = ifelse(Ct>0,1,0)) %>%
  left_join(.,PCR) %>%
  mutate(station.depth = paste0(station,".",depth))

these.samps <- dat.samp %>% dplyr::select(sample,inhibit.val,dilution) %>% 
              filter(inhibit.val < INHIBIT.LIMIT,dilution<1) %>% 
              dplyr::select(sample,dilution) %>%
              group_by(sample,dilution) %>% summarise(N=length(sample)) %>% arrange(sample) %>%
              as.data.frame()

dat.samp.dil <- dat.samp %>% filter(sample %in% these.samps,inhibit.val < INHIBIT.LIMIT,dilution<1)
dat.samp <- dat.samp %>% filter(!sample %in% these.samps) %>% bind_rows(.,dat.samp.dil)

# dat.samp includes a few samples with multiple different dilutions (0.1 and 0.2) I have left them this way for now.

### MANUALLY CUT A COUPLE OF SAMPLES THAT ARE KILLING THE MODEL FIT and after inspection are clearly causing troubles.
if(SP == "hake"){
  dat.samp <- dat.samp 
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
dat.inhibit <- dat.samp %>% filter( inhibit.bin==1) 
dat.samp    <- dat.samp %>% filter(inhibit.bin==0) 
dat.control.field.neg    <- dat.control.field.neg %>% filter(inhibit.bin==0) 

####### THIS SECTION IS FOR EXPLORING SOME OF THE STRNAGENESS WITH
####### UNDERSTANDING THE EFFECT OF DILUTIONS.  IGNORE
# dat.few <- dat.samp 
# dat.few$copies_ul[dat.few$copies_ul == ""] <- NA
# dat.few$copies_ul <- as.numeric(as.character(dat.few$copies_ul))
# B <- dat.few %>% group_by(sample,dilution) %>% 
#   summarise(count=mean(copies_ul,na.rm=T),N=sum(Ct_bin))
# 
# C <- B %>% dplyr::select(-N) %>% ungroup() %>% pivot_wider(.,names_from = c("dilution"),values_from = "count") %>% as.data.frame()
# colnames(C)[2:5] <- c("x1","x.1","x.2","x.5")
# 
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

# if(SP == "lamprey"){ # ONLY KEEP 0.1 for lamprey
#   dat.samp <- dat.samp %>% filter(!dilution %in% c(0.5,0.2))
# }
# if(SP == "eulachon"){ # ONLY KEEP 0.2 for Eulachon.
#   dat.samp <- dat.samp %>% filter(!dilution %in% c(0.5,0.1))
# }
# if(SP == "hake"){ # ONLY KEEP 0.2 for hake
#   dat.samp <- dat.samp %>% filter(!dilution %in% c(0.5,0.1))
# }
B <- dat.samp %>% group_by(qPCR,dilution) %>% summarise(N=length(dilution))  %>% as.data.frame()

########## CHECK ON HOW MANY SAMPLES WE HAVE ZERO REPLICATES FOR  
# Examine how many samples are retained out of the intial number done.
dat.samp.ok <- dat.samp %>% group_by(sample,inhibit.bin) %>% summarise(N=length(sample)) %>% as.data.frame()

dat.samp.begin %>% nrow() #1916
dat.samp.ok %>% nrow() #1901 .... lost 15 samples through filtering. Most are 2 transect 25 and most ar e due to inhibition... weird.

dat.all %>% filter(sample %in% c(dat.samp.begin %>% filter(!sample %in% dat.samp.ok$sample) %>%pull(sample))) %>% as.data.frame()

### GET RID OF 25 m deep samples and 200m samples if doing smoothes by depth category.
if(TRIM.25 ==TRUE) {
  dat.samp <- dat.samp %>% filter(!depth_cat %in% c(25))
}

































################################### Making indicators and indices for stan.
STATION.DEPTH <- data.frame(station.depth=unique(dat.samp$station.depth),
                            station_depth_idx = 1:length(unique(dat.samp$station.depth)))
# add depth and depth category
STATION.DEPTH <- STATION.DEPTH %>% left_join(.,
                                             dat.samp %>% dplyr::select(station.depth,depth,depth_cat,transect,station) %>%
                                               group_by(station.depth,depth,depth_cat,transect,station) %>% summarise(N=length(depth)))
# Add utm.lat, utm.lon, observed water depth.
STATION.DEPTH <- STATION.DEPTH %>% left_join(.,
                                             dat.station.id.trim %>% dplyr::select(transect,station,utm.lon,utm.lat,lat,lon,
                                                                                   bottom.depth=water.depth,
                                                                                   bathy.bottom.depth,
                                                                                   bottom.depth.consensus,
                                                                                   transect.dist.km)) %>%
  dplyr::select(-N)

STATION.DEPTH$depth_cat_factor <- factor(STATION.DEPTH$depth_cat)

########################################################################
## SAMPLES
########################################################################
# Add indices associated with each sample.
SAMPLES <- dat.samp %>% dplyr::select(station.depth,sample) %>% group_by(station.depth,sample) %>% 
  summarise(N=length(sample)) %>% dplyr::select(-N) %>% left_join(.,STATION.DEPTH) %>% arrange(station_depth_idx) %>%
  ungroup()
SAMPLES <- SAMPLES %>% mutate(sample_idx = 1:nrow(SAMPLES))
SAMPLES <- SAMPLES %>% group_by(station.depth) %>% summarise(N=length(sample)) %>% 
  mutate(singleton_idx=ifelse(N==1,0,1)) %>% dplyr::select(station.depth,singleton_idx) %>%
  left_join(SAMPLES,.) %>% rename(samp_station_depth_idx =station_depth_idx)
SAMPLES <- dat.samp %>% group_by(sample) %>% 
  summarise(vol_standard=mean(vol.standard),log_vol_standard=log10(vol_standard),wash_idx=mean(wash_idx)) %>%
  left_join(SAMPLES,.)

these_samp <- dat.samp %>% group_by(station.depth) %>% summarise(Ct_bin_idx = max(Ct_bin))
SAMPLES <- left_join(SAMPLES,these_samp) 

# Make index for depth category 
if(NO.SURFACE!="TRUE"){
  SAMPLES <- SAMPLES %>% mutate(depth_idx = case_when(depth_cat == 0 ~ 1,
                                                      depth_cat == 50 ~ 2,
                                                      depth_cat == 100 ~ 3,
                                                      depth_cat == 150 ~ 4,
                                                      depth_cat == 300 ~ 5,
                                                      depth_cat == 500 ~ 6))
}else{
  SAMPLES <- SAMPLES %>% mutate(depth_idx = case_when(depth_cat == 50 ~ 1,
                                                      depth_cat == 100 ~ 2,
                                                      depth_cat == 150 ~ 3,
                                                      depth_cat == 300 ~ 4,
                                                      depth_cat == 500 ~ 5))
}

N_depth <- length(unique(SAMPLES$depth_cat))

SAMPLES.CONTROL <- data.frame(sample=unique(dat.control.field.neg$sample),
                              sample_control_idx = 1:length(unique(dat.control.field.neg$sample)))

dat.samp <- left_join(dat.samp,SAMPLES) %>% left_join(.,STATION.DEPTH) 
dat.samp$loo_idx <- 1:nrow(dat.samp)
dat.control.field.neg <- left_join(dat.control.field.neg,SAMPLES.CONTROL)

# This identifies where there is only one sample for each station-depth combination.  This is important for model identifiability.
# dat.samp   <-  dat.samp %>% group_by(station_depth_idx) %>% summarise(N=length(unique(sample))) %>% 
#               mutate(singleton_idx=ifelse(N==1,0,1)) %>% dplyr::select(-N) %>%
#               left_join(dat.samp,.) %>% mutate()

### THERE IS SOME WEIRDNESS WITH SAMPLES 1306, 1361, 1362, and 412.  I think I have resolved them.
dat.stand.bin <- dat.stand
dat.stand.pos <- dat.stand %>% filter(Ct_bin==1)
N_stand_bin <- nrow(dat.stand.bin)
N_stand_pos <- nrow(dat.stand.pos)

dat.obs.bin <- dat.samp
dat.obs.pos <- dat.samp %>% filter(Ct_bin==1)
N_obs_bin <- nrow(dat.obs.bin)
N_obs_pos <- nrow(dat.obs.pos)

bin_log_dilution_obs   <- log10(dat.obs.bin$dilution)
pos_log_dilution_obs   <- log10(dat.obs.pos$dilution)
log_vol_obs            <- SAMPLES$log_vol_standard
singleton_idx          <- SAMPLES$singleton_idx
samp_station_depth_idx <- SAMPLES$samp_station_depth_idx
wash_idx               <- SAMPLES$wash_idx


dat.control.bin <- dat.control.field.neg
dat.control.pos <- dat.control.field.neg %>% filter(Ct_bin==1)
N_control_bin <- nrow(dat.control.bin)
N_control_pos <- nrow(dat.control.pos)
bin_log_vol_control      <- log10(dat.control.bin$vol.standard)
pos_log_vol_control      <- log10(dat.control.pos$vol.standard)
bin_log_dilution_control <- log10(dat.control.bin$dilution)
pos_log_dilution_control <- log10(dat.control.pos$dilution)

N_pcr <- max(PCR$plate_idx)
N_sample  <- max(SAMPLES$sample_idx)
N_station_depth <- max(STATION.DEPTH$station_depth_idx)
N_control_sample <- max(SAMPLES.CONTROL$sample_control_idx)

if(SP=="hake"){
  wash_offset_prior <- c(-1,1)                          
}else{
  wash_offset_prior <- c(wash_offset_hake$Mean,wash_offset_hake$SD*0.01)                           
}

# Deal with singleton in control positive for lamprey.
if(nrow(dat.control.pos)==1){
  N_control_pos = 2
  single_n_control_pos <- 1
  pos_control    = c(dat.control.pos$Ct,0)
  pos_log_vol_control= c(pos_log_vol_control,0)
  pos_log_dilution_control = c(pos_log_dilution_control,0)
  pcr_control_pos_idx    <- c(dat.control.pos$plate_idx,0)
  sample_control_pos_idx <- c(dat.control.pos$sample_control_idx,0) 
  wash_control_pos_idx   <- c(dat.control.pos$wash_idx,0)
}else{
  single_n_control_pos = 0
  pos_control         = dat.control.pos$Ct
  pos_log_vol_control = pos_log_vol_control
  pos_log_dilution_control = pos_log_dilution_control
  pcr_control_pos_idx    <- dat.control.pos$plate_idx
  sample_control_pos_idx <- dat.control.pos$sample_control_idx
  wash_control_pos_idx   <- dat.control.pos$wash_idx
}
#################################################################### 
#################################################################### 
#################################################################### 
# Make the code for smoothes using the data in the STATION.DEPTH data.frame


if(MODEL.TYPE == "lat.long.smooth"){
  # LAT-LON TENSOR SMOOTHES
  utm.lon.lims <- c(min(STATION.DEPTH$utm.lon), max(STATION.DEPTH$utm.lon))
  utm.lat.lims <- c(min(STATION.DEPTH$utm.lat), max(STATION.DEPTH$utm.lat))
  knots.lon    <- seq(utm.lon.lims[1],utm.lon.lims[2],length.out=N.knots.lon)
  knots.lat    <- seq(utm.lat.lims[1],utm.lat.lims[2],length.out=N.knots.lat)
  
  # bottom.depth smooth
  #N.knots.bd <- 6
  bd.lims <- c(min(STATION.DEPTH$bottom.depth.consensus), max(STATION.DEPTH$bottom.depth.consensus))
  knots.bd <- seq(bd.lims[1],bd.lims[2],length.out=N.knots.bd)
  
  # depth smooth
  #N.knots.depth <- 4
  depth.lims <- c(min(STATION.DEPTH$depth_cat), max(STATION.DEPTH$depth_cat))
  knots.depth <- seq(depth.lims[1],depth.lims[2],length.out=N.knots.depth)
  
  STATION.DEPTH.temp <- STATION.DEPTH %>% mutate(Y=rnorm(nrow(STATION.DEPTH),0,1))
  
  brms.object <- brm(Y ~ t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)+
                       s(bottom.depth.consensus,k=N.knots.bd) +
                       depth_cat_factor, 
                     #knots=list(#utm.lon=knots.lon,
                     #utm.lat=knots.lat,
                     #      bottom.depth = knots.bd),
                     data=STATION.DEPTH.temp,
                     chains=0)
  
  smooth.dat <- standata(brms.object)
  code.dat   <- stancode(brms.object)
  # THIS IS THE MAGIC SAUCE FOR MAKING PREDICTIONS ON THE SAME SET OF BASIS FUNCTIONS.
  #    new.dat <- data.frame(Y = rep(0,3),bottom.depth = Q$bottom.depth[1:3])
  #    smooth.dat.pred <- standata(brms.object,newdata=new.dat)
  
  # C <- brm(Y ~ s(bottom.depth,k=N.knots.bd),
  #          knots=list(bottom.depth = knots.bd),
  #          data=Q,
  #          chains=0)
  
  # brms.object <- 
  # smooth.datC <- standata(C)
  # code.datC   <- stancode(C)
}

if(MODEL.TYPE == "lat.long.smooth.base"){
  
  # LAT-LON TENSOR SMOOTHES
  utm.lon.lims <- c(min(STATION.DEPTH$utm.lon), max(STATION.DEPTH$utm.lon))
  utm.lat.lims <- c(min(STATION.DEPTH$utm.lat), max(STATION.DEPTH$utm.lat))
  knots.lon    <- seq(utm.lon.lims[1],utm.lon.lims[2],length.out=N.knots.lon)
  knots.lat    <- seq(utm.lat.lims[1],utm.lat.lims[2],length.out=N.knots.lat)
  
  # bottom.depth smooth
  # N.knots.bd <- 6
  # bd.lims <- c(min(STATION.DEPTH$bathy.bottom.depth), max(STATION.DEPTH$bathy.bottom.depth))
  # knots.bd <- seq(bd.lims[1],bd.lims[2],length.out=N.knots.bd)
  
  # depth smooth
  # N.knots.depth <- 4
  # depth.lims <- c(min(STATION.DEPTH$depth_cat), max(STATION.DEPTH$depth_cat))
  # knots.depth <- seq(depth.lims[1],depth.lims[2],length.out=N.knots.depth)
  
  STATION.DEPTH.temp <- STATION.DEPTH %>% mutate(Y=rnorm(nrow(STATION.DEPTH),0,1))
  
  # print("HEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHE")
  # print(N.knots.lon)
  brms.object <- brm(Y ~ t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)+
                       depth_cat_factor,
                     #knots=list(#utm.lon=knots.lon,
                     #utm.lat=knots.lat,
                     #      bottom.depth = knots.bd),
                     data=STATION.DEPTH.temp,
                     chains=0)
  # print("HERE")
  smooth.dat <- standata(brms.object)
  code.dat   <- stancode(brms.object)
  # THIS IS THE MAGIC SAUCE FOR MAKING PREDICTIONS ON THE SAME SET OF BASIS FUNCTIONS.
  #    new.dat <- data.frame(Y = rep(0,3),bottom.depth = Q$bottom.depth[1:3])
  #    smooth.dat.pred <- standata(brms.object,newdata=new.dat)
  
  # C <- brm(Y ~ s(bottom.depth,k=N.knots.bd),
  #          knots=list(bottom.depth = knots.bd),
  #          data=Q,
  #          chains=0)
  
  # brms.object <- 
  # smooth.datC <- standata(C)
  # code.datC   <- stancode(C)
}


#################################################################### 
#################################################################### 
stan_data = list(
  "bin_stand"     = dat.stand.bin$Ct_bin,
  "pos_stand"     = dat.stand.pos$Ct,
  "D_bin_stand"   = dat.stand.bin$log10_copies,
  "D_pos_stand"   = dat.stand.pos$log10_copies,
  
  #Observations
  "bin_obs"    = dat.obs.bin$Ct_bin,
  "pos_obs"    = dat.obs.pos$Ct,
  
  #Covariates associated with individual PCRs
  "bin_log_dilution_obs"= bin_log_dilution_obs,
  "pos_log_dilution_obs"= pos_log_dilution_obs,
  
  # Covariates associated with individual samples
  "log_vol_obs"     = log_vol_obs,
  "singleton_idx"   = singleton_idx,
  "wash_idx"        = wash_idx,
  
  # Indices for sample bottles combinations and singleton indicator
  "sample_idx" = SAMPLES$sample_idx,
  "samp_station_depth_idx" = SAMPLES$samp_station_depth_idx,
  "singleton_idx" = SAMPLES$singleton_idx,
  "Ct_bin_idx" = SAMPLES$Ct_bin_idx,
  "depth_idx"  = SAMPLES$depth_idx,
  
  # Indices for station-depth combination
  "station_depth_idx" = STATION.DEPTH$station_depth_idx,
  
  # Control Values
  "bin_control"    = dat.control.bin$Ct_bin,
  "pos_control"    = pos_control,
  "bin_log_vol_control" = bin_log_vol_control,
  "pos_log_vol_control" = pos_log_vol_control,
  "bin_log_dilution_control"= bin_log_dilution_control,
  "pos_log_dilution_control"= pos_log_dilution_control,
  
  # Indices and counters
  "N_pcr"    = N_pcr,    # Number of PCR plates
  "N_depth"  = N_depth,
  "N_sample" = N_sample,
  "N_control_sample" = N_control_sample,
  "N_station_depth"  = N_station_depth,
  "N_stand_bin" = N_stand_bin,
  "N_stand_pos" = N_stand_pos,
  "N_obs_bin"   = N_obs_bin,
  "N_obs_pos"   = N_obs_pos,
  "N_control_bin"   = N_control_bin,
  "N_control_pos"   = N_control_pos,
  
  # Indices for qPCR plate
  "pcr_stand_bin_idx"   = dat.stand.bin$plate_idx,
  "pcr_stand_pos_idx" = dat.stand.pos$plate_idx,
  "pcr_obs_bin_idx"   = dat.obs.bin$plate_idx,
  "pcr_obs_pos_idx" =   dat.obs.pos$plate_idx,
  "pcr_control_bin_idx"   = dat.control.bin$plate_idx,
  "pcr_control_pos_idx" =   pcr_control_pos_idx,
  
  # Index for leave.one.out
  "loo_pos_idx" = dat.obs.pos$loo_idx,
  
  # Control samples  
  "sample_control_idx" = SAMPLES.CONTROL$sample_control_idx,
  
  # Indices for Samples
  "sample_pos_idx" = dat.obs.pos$sample_idx,
  "sample_bin_idx" = dat.obs.bin$sample_idx,
  "sample_control_pos_idx" = sample_control_pos_idx,
  "sample_control_bin_idx" = dat.control.bin$sample_control_idx,
  
  "station_depth_pos_idx" = dat.obs.pos$station_depth_idx,
  "station_depth_bin_idx" = dat.obs.bin$station_depth_idx,
  
  # "singleton_pos_idx" = dat.obs.pos$singleton_idx,
  # "singleton_bin_idx" = dat.obs.bin$singleton_idx,
  
  "wash_pos_idx" = dat.obs.pos$wash_idx,
  "wash_bin_idx" = dat.obs.bin$wash_idx,
  
  "wash_control_pos_idx" = wash_control_pos_idx,
  "wash_control_bin_idx" = dat.control.bin$wash_idx,
  
  #Offset of density for improving fitting characteristics
  "OFFSET" = 0,
  
  "single_n_control_pos" = single_n_control_pos,
  
  # Contamination Weight
  "wash_offset_prior" = wash_offset_prior,
  
  
  # SMOOTH COMPONENTS.
  # Data for linear effects
  "K" = smooth.dat$K, # number of population-level effects
  "X" = smooth.dat$X,  # population-level design matrix
  # data for splines
  "Ks" = smooth.dat$Ks, # number of linear effects
  "Xs" = smooth.dat$Xs, # design matrix for the linear effects
  # data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)0
  "nb_1" = smooth.dat$nb_1,  # number of bases
  "knots_1" = smooth.dat$knots_1, # number of knots
  # basis function matrices
  "Zs_1_1" = smooth.dat$Zs_2_1,
  "Zs_1_2" = smooth.dat$Zs_2_2,
  "Zs_1_3" = smooth.dat$Zs_2_3,
  # data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)50
  "nb_2" = smooth.dat$nb_2,  # number of bases
  "knots_2" = smooth.dat$knots_2,  # number of knots
  # basis function matrices
  "Zs_2_1" = smooth.dat$Zs_2_1,
  "Zs_2_2" = smooth.dat$Zs_2_2,
  "Zs_2_3" = smooth.dat$Zs_2_3,
  # data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)100
  "nb_3" = smooth.dat$nb_3,  # number of bases
  "knots_3" = smooth.dat$knots_3,   # number of knots
  # basis function matrices
  "Zs_3_1"= smooth.dat$Zs_3_1,
  "Zs_3_2"= smooth.dat$Zs_3_2,
  "Zs_3_3"= smooth.dat$Zs_3_3,
  
  # data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)150
  "nb_4"= smooth.dat$nb_4,  # number of bases
  "knots_4"= smooth.dat$knots_4,  # number of knots
  # basis function matrices
  "Zs_4_1"= smooth.dat$Zs_4_1,
  "Zs_4_2"= smooth.dat$Zs_4_2,
  "Zs_4_3"= smooth.dat$Zs_4_3,
  
  # data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)300
  "nb_5"= smooth.dat$nb_5, # number of bases
  "knots_5"= smooth.dat$knots_5, # number of knots
  # basis function matrices
  "Zs_5_1"= smooth.dat$Zs_5_1,
  "Zs_5_2"= smooth.dat$Zs_5_2,
  "Zs_5_3"= smooth.dat$Zs_5_3
) 
if(NO.SURFACE != "TRUE"){
  stan_data = c(stan_data,
                list(# data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)500
                  "nb_6" = smooth.dat$nb_6,  # number of bases
                  "knots_6"= smooth.dat$knots_6,  # number of knots
                  # basis function matrices
                  "Zs_6_1"= smooth.dat$Zs_6_1,
                  "Zs_6_2"= smooth.dat$Zs_6_2,
                  "Zs_6_3"= smooth.dat$Zs_6_3))
  if(MODEL.TYPE == "lat.long.smooth"){
    stan_data = c(stan_data,
                  list(# data for spline s(bottom.depth,k=N.knots.bd)
                    "nb_7"= smooth.dat$nb_7,  # number of bases
                    "knots_7"= smooth.dat$knots_7,  # number of knots
                    # basis function matrices
                    "Zs_7_1"= smooth.dat$Zs_7_1))
  }
}
if(NO.SURFACE == "TRUE"){
  if(MODEL.TYPE == "lat.long.smooth"){
    stan_data = c(stan_data,
                  list(# data for spline s(bottom.depth,k=N.knots.bd)
                    "nb_7"= smooth.dat$nb_7,  # number of bases
                    "knots_7"= smooth.dat$knots_7,  # number of knots
                    # basis function matrices
                    "Zs_7_1"= smooth.dat$Zs_7_1))
  }
}




stan_pars = c(
  "beta_0", # intercept for standards
  "beta_1", # slope for standards
  "phi_0",  # logit intercept for standards
  "phi_1",  # logit slope for standard,
  "wash_offset", # parameter for concentration offset from washing with 30% EtOH instead of 70%.
  
  "mu_contam", # log-mean of average contamination
  "sigma_contam", # lod-sd of contamination
  
  "D",     # Latent variable for Log10-count in each station-location combination
  "D_delta",
  #"D_error",
  "D_contam",
  "D_control", # Latent variable for log-count in field negative controls.
  
  
  "sigma_pcr",     # variability among samples, given individual bottle, site, and month 
  "sigma_stand_int",
  #"nu",
  "delta",  # random effect for sample #
  "tau_sample"# sd among samples given station-depth 
  # "pred_bin",
  # "pred_pos",
  # "log_lik"
)   

if(MODEL.TYPE=="lat.long.smooth" & NO.SURFACE !="TRUE"){
  stan_pars <- c(stan_pars,
                 # smooth and linear coefficients
                 "D",
                 "D_delta",
                 # "Intercept",
                 "b",
                 "bs",
                 "s_1_1","s_1_2","s_1_3",
                 "s_2_1","s_2_2","s_2_3",
                 "s_3_1","s_3_2","s_3_3",
                 "s_4_1","s_4_2","s_4_3",
                 "s_5_1","s_5_2","s_5_3",
                 "s_6_1","s_6_2","s_6_3",
                 "s_7_1"
  ) 
  N_knots_all <- list("N.knots.lon"=N.knots.lon,
                      "N.knots.lat"=N.knots.lat,
                      "N.knots.bd"=N.knots.bd)
}
if(MODEL.TYPE=="lat.long.smooth" & NO.SURFACE =="TRUE"){
  stan_pars <- c(stan_pars,
                 # smooth and linear coefficients
                 "D",
                 "D_delta",
                 # "Intercept",
                 "b",
                 "bs",
                 "s_1_1","s_1_2","s_1_3",
                 "s_2_1","s_2_2","s_2_3",
                 "s_3_1","s_3_2","s_3_3",
                 "s_4_1","s_4_2","s_4_3",
                 "s_5_1","s_5_2","s_5_3",
                 "s_6_1"
  ) 
  N_knots_all <- list("N.knots.lon"=N.knots.lon,
                      "N.knots.lat"=N.knots.lat,
                      "N.knots.bd"=N.knots.bd)
}


if(MODEL.TYPE=="lat.long.smooth.base"){
  stan_pars <- c(stan_pars,
                 # smooth and linear coefficients
                 # "Intercept",
                 "b",
                 "bs",
                 "s_1_1","s_1_2","s_1_3",
                 "s_2_1","s_2_2","s_2_3",
                 "s_3_1","s_3_2","s_3_3",
                 "s_4_1","s_4_2","s_4_3",
                 "s_5_1","s_5_2","s_5_3",
                 "s_6_1","s_6_2","s_6_3"
  )   
  N_knots_all <- list("N.knots.lon"=N.knots.lon,
                      "N.knots.lat"=N.knots.lat)
}
if(MODEL.VAR=="Base_Var"){
  stan_pars <- c(stan_pars)
}
if(MODEL.VAR=="Linear_Var"){
  stan_pars <- c(stan_pars,
                 "gamma0","gamma1",
                 "eta0","eta1")
}



### INTIAL VALUES
stan_init_f2 <- function(n.chain,N_pcr,N_station_depth,N_control_sample){ 
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      
      sigma_stand_int = runif(1,0.01,2),
      beta_0 = runif(N_pcr,30,40),
      beta_1 = rnorm(N_pcr,-4,1),
      wash_offset = rnorm(1,-2,1),
      
      phi_0  = runif(N_pcr,0,5),
      phi_1  = rnorm(N_pcr,3,0.5),
      sigma_pcr = runif(1,0.01,0.4),
      sigma_stand_int = runif(1,0.01,0.4),
      D_control  = rnorm(N_control_sample,0,2),
      gamma0 =runif(1,0.1,0.5),
      gamma1 =runif(1,0,0.05),
      eta0 =runif(1,0.1,0.5),
      eta1 =runif(1,0,0.05)
    )
  }
  return(A)
}

################################################################
################################################################
# STAN MODEL FOR THE UNKNOWN SAMPLES
################################################################
################################################################

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

N_CHAIN = 4
Warm = 1500
Iter = 9000
Treedepth = 13
Adapt_delta = 0.86

LOC <- paste0(base.dir,"Scripts/Stan Files/")
setwd(LOC)

if(MODEL.TYPE=="Base"){
  stanMod = stan(file = "qPCR_Hake.stan" ,data = stan_data, 
                 verbose = FALSE, chains = N_CHAIN, thin = 6, 
                 warmup = Warm, iter = Warm + Iter, 
                 control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                 pars = stan_pars,
                 boost_lib = NULL,
                 sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,".csv"),
                 init = stan_init_f2(n.chain=N_CHAIN,
                                     N_pcr= N_pcr,
                                     N_station_depth = N_station_depth,
                                     N_control_sample = N_control_sample)
  )
}

if(MODEL.TYPE=="lat.long.smooth"){
  if(MODEL.VAR=="Base_Var"){
    if(NO.SURFACE!="TRUE"){
      stanMod = stan(file = "qPCR_Hake_smoothes.stan" ,data = stan_data, 
                     verbose = FALSE, chains = N_CHAIN, thin = 4, 
                     warmup = Warm, iter = Warm + Iter, 
                     control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                     pars = stan_pars,
                     boost_lib = NULL,
                     sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,"_",MODEL.VAR,".csv"),
                     init = stan_init_f2(n.chain=N_CHAIN,
                                         N_pcr= N_pcr,
                                         N_station_depth = N_station_depth,
                                         N_control_sample = N_control_sample)
      )
    }
    if(NO.SURFACE=="TRUE"){
      stanMod = stan(file = "qPCR_Hake_smoothes_no_surf.stan" ,data = stan_data, 
                     verbose = FALSE, chains = N_CHAIN, thin = 4, 
                     warmup = Warm, iter = Warm + Iter, 
                     control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                     pars = stan_pars,
                     boost_lib = NULL,
                     sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,"_",MODEL.VAR,".csv"),
                     init = stan_init_f2(n.chain=N_CHAIN,
                                         N_pcr= N_pcr,
                                         N_station_depth = N_station_depth,
                                         N_control_sample = N_control_sample)
      )
    }
  }
  if(MODEL.VAR=="Linear_Var"){
    stanMod = stan(file = "qPCR_Hake_smoothes_VAR2.stan" ,data = stan_data, 
                   verbose = FALSE, chains = N_CHAIN, thin = 2, 
                   warmup = Warm, iter = Warm + Iter, 
                   control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                   pars = stan_pars,
                   boost_lib = NULL,
                   sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,"_",MODEL.VAR,".csv"),
                   init = stan_init_f2(n.chain=N_CHAIN,
                                       N_pcr= N_pcr,
                                       N_station_depth = N_station_depth,
                                       N_control_sample = N_control_sample)
    )
  }
}

if(MODEL.TYPE=="lat.long.smooth.base"){
  if(MODEL.VAR=="Base_Var"){
    stanMod = stan(file = "qPCR_Hake_smoothes_no_depth.stan" ,data = stan_data, 
                   verbose = FALSE, chains = N_CHAIN, thin = 2, 
                   warmup = Warm, iter = Warm + Iter, 
                   control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                   pars = stan_pars,
                   boost_lib = NULL,
                   sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,"_",MODEL.VAR,".csv"),
                   init = stan_init_f2(n.chain=N_CHAIN,
                                       N_pcr= N_pcr,
                                       N_station_depth = N_station_depth,
                                       N_control_sample = N_control_sample)
    )
  }
  if(MODEL.VAR=="Linear_Var"){
    stanMod = stan(file = "qPCR_Hake_smoothes_no_depth_VAR.stan" ,data = stan_data, 
                   verbose = FALSE, chains = N_CHAIN, thin = 2, 
                   warmup = Warm, iter = Warm + Iter, 
                   control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                   pars = stan_pars,
                   boost_lib = NULL,
                   sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,"_",MODEL.VAR,".csv"),
                   init = stan_init_f2(n.chain=N_CHAIN,
                                       N_pcr= N_pcr,
                                       N_station_depth = N_station_depth,
                                       N_control_sample = N_control_sample)
    )
  }
}

# get_adaptation_info(stanMod)
pars <- rstan::extract(stanMod, permuted = TRUE)
# log_lik_1 <- extract_log_lik(stanMod, merge_chains = FALSE)
# r_eff <- relative_eff(exp(log_lik_1))
# loo_val <- loo(log_lik_1, r_eff = r_eff)
# print(loo_val)
# plot(loo_val)

# bad <- loo_val$pointwise %>% as.data.frame()
# bad$ID <- 1:nrow(bad)
# bad2 <- bad %>% filter(influence_pareto_k>1)
# 
# bad.dat <- dat.samp[bad2$ID,]
# dim(bad.dat)
# 
# ggplot()+
#     geom_histogram(data=bad.dat %>% filter(Ct>0),aes(Ct),fill="red",alpha=0.5) +
#    geom_histogram(data=dat.samp %>% filter(Ct>0),aes(Ct),fill="blue",alpha=0.5) + theme_bw()
###############3333

samp_params <- get_sampler_params(stanMod)
#samp_params 
stanMod_summary <- summary(stanMod)$summary
stanMod_summary_parts <- list()
stanMod_summary_parts[[as.name("D")]] <- summary(stanMod,pars="D")$summary  
stanMod_summary_parts[[as.name("reg")]] <- summary(stanMod,pars=c("beta_0",
                                                                  "beta_1",
                                                                  "phi_0",
                                                                  "phi_1"))$summary
stanMod_summary_parts[[as.name("param")]] <- summary(stanMod,pars=c(
  "wash_offset",
  "mu_contam",
  "sigma_contam",
  "tau_sample",
  "sigma_stand_int", # variability among standard regression.
  "sigma_pcr"#,"nu"
))$summary     # variability among samples, given individual bottle, site, and month ))$summary
if(MODEL.TYPE == "lat.long.smooth"){
  if(NO.SURFACE != "TRUE"){
    stanMod_summary_parts[[as.name("smoothes")]] <- summary(stanMod, pars=c("b", "bs",
                                                                            "s_1_1","s_1_2","s_1_3",  
                                                                            "s_2_1","s_2_2","s_2_3",
                                                                            "s_3_1","s_3_2","s_3_3",
                                                                            "s_4_1","s_4_2","s_4_3",
                                                                            "s_5_1","s_5_2","s_5_3",
                                                                            "s_6_1","s_6_2","s_6_3",
                                                                            "s_7_1"))$summary
  }else{
    stanMod_summary_parts[[as.name("smoothes")]] <- summary(stanMod, pars=c("b", "bs",
                                                                            "s_1_1","s_1_2","s_1_3",  
                                                                            "s_2_1","s_2_2","s_2_3",
                                                                            "s_3_1","s_3_2","s_3_3",
                                                                            "s_4_1","s_4_2","s_4_3",
                                                                            "s_5_1","s_5_2","s_5_3",
                                                                            "s_6_1"))$summary
  }
}
if(MODEL.TYPE == "lat.long.smooth.base"){
  stanMod_summary_parts[[as.name("smoothes")]] <- summary(stanMod, pars=c("b", "bs",
                                                                          "s_1_1","s_1_2","s_1_3",  
                                                                          "s_2_1","s_2_2","s_2_3",
                                                                          "s_3_1","s_3_2","s_3_3",
                                                                          "s_4_1","s_4_2","s_4_3",
                                                                          "s_5_1","s_5_2","s_5_3",
                                                                          "s_6_1","s_6_2","s_6_3"))$summary
  
}




ID <- rownames(stanMod_summary_parts$D)
stanMod_summary_D <- stanMod_summary_parts$D %>%  as.data.frame() %>% mutate(ID=ID) %>% arrange(desc(Rhat))

####
base_params <- c(
  "beta_0",
  "beta_1",
  "phi_0",
  "phi_1",
  "wash_offset",
  
  "mu_contam",
  "sigma_contam",
  
  #"D",
  #"delta",
  "tau_sample",
  
  "sigma_stand_int", # variability among standard regression.
  "sigma_pcr"     # variability among samples, given individual bottle, site, and month 
  
) 

##### MAKE SOME DIAGNOSTIC PLOTS
TRACE <- list()
TRACE[[as.name("D")]] <- traceplot(stanMod,pars=c("lp__","D[12]","D[234]","D[456]","D[855]"),inc_warmup=FALSE)
#TRACE[[as.name("mu_smooth")]] <- traceplot(stanMod,pars=c("lp__","mu_smooth[12]","mu_smooth[234]","mu_smooth[456]","mu_smooth[878]"),inc_warmup=FALSE)

TRACE[[as.name("Phi0")]] <- traceplot(stanMod,pars=c("lp__","phi_0"),inc_warmup=FALSE)
TRACE[[as.name("Phi1")]] <- traceplot(stanMod,pars=c("lp__","phi_1"),inc_warmup=FALSE)
TRACE[[as.name("Beta0")]] <- traceplot(stanMod,pars=c("lp__","beta_0"),inc_warmup=FALSE)
TRACE[[as.name("Beta1")]] <- traceplot(stanMod,pars=c("lp__","beta_1"),inc_warmup=FALSE)
TRACE[[as.name("Contam")]] <- traceplot(stanMod,pars=c("lp__","wash_offset","mu_contam","sigma_contam"),inc_warmup=FALSE)

if(MODEL.TYPE=="lat.long.smooth"){
  TRACE[[as.name("b")]] <- traceplot(stanMod,pars=c("b"),inc_warmup=FALSE)
  TRACE[[as.name("bs")]] <- traceplot(stanMod,pars=c("bs"),inc_warmup=FALSE)
  
  TRACE[[as.name("s_1")]] <- traceplot(stanMod,pars=c("s_1_1","s_1_2","s_1_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_2")]] <- traceplot(stanMod,pars=c("s_2_1","s_2_2","s_2_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_3")]] <- traceplot(stanMod,pars=c("s_3_1","s_3_2","s_3_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_4")]] <- traceplot(stanMod,pars=c("s_4_1","s_4_2","s_4_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_5")]] <- traceplot(stanMod,pars=c("s_5_1","s_5_2","s_5_3"),inc_warmup=FALSE)
  
  if(NO.SURFACE!="TRUE"){
    TRACE[[as.name("s_6")]] <- traceplot(stanMod,pars=c("s_6_1","s_6_2","s_6_3"),inc_warmup=FALSE)
    TRACE[[as.name("s_7")]] <- traceplot(stanMod,pars=c("s_7_1"),inc_warmup=FALSE)
  }else{
    TRACE[[as.name("s_6")]] <- traceplot(stanMod,pars=c("s_6_1"),inc_warmup=FALSE)
  }
}else if(MODEL.TYPE=="lat.long.smooth.base"){
  TRACE[[as.name("b")]] <- traceplot(stanMod,pars=c("b"),inc_warmup=FALSE)
  TRACE[[as.name("bs")]] <- traceplot(stanMod,pars=c("bs"),inc_warmup=FALSE)
  
  TRACE[[as.name("s_1")]] <- traceplot(stanMod,pars=c("s_1_1","s_1_2","s_1_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_2")]] <- traceplot(stanMod,pars=c("s_2_1","s_2_2","s_2_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_3")]] <- traceplot(stanMod,pars=c("s_3_1","s_3_2","s_3_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_4")]] <- traceplot(stanMod,pars=c("s_4_1","s_4_2","s_4_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_5")]] <- traceplot(stanMod,pars=c("s_5_1","s_5_2","s_5_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_6")]] <- traceplot(stanMod,pars=c("s_6_1","s_6_2","s_6_3"),inc_warmup=FALSE)
  #TRACE[[as.name("s_7")]] <- traceplot(stanMod,pars=c("s_7_1"),inc_warmup=FALSE)
}

if(MODEL.VAR=="Base_Var"){
  TRACE[[as.name("Var")]] <- traceplot(stanMod,pars=c("lp__","tau_sample","sigma_stand_int","sigma_pcr"),inc_warmup=FALSE)
}
if(MODEL.VAR=="Linear_Var"){
  TRACE[[as.name("Var")]] <- traceplot(stanMod,pars=c("lp__","tau_sample","gamma0","gamma1","sigma_pcr"),inc_warmup=FALSE)
}
TRACE$Var
TRACE$Contam
TRACE$b
TRACE$bs

#### Get rid of all the predictions for likelihood and for individual observations
#### from pars.  
#### Save only the mean prediction for each observation and the loo summary
#### to maintain a reasonable file size.

# pred_bin <- apply(pars$pred_bin,2,quantile,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# pred_bin_mean <- colMeans(pars$pred_bin)
# pred_bin <- pred_bin %>% t() %>% as.data.frame() %>% bind_cols(pred_bin_mean,.) %>% as.data.frame()
# colnames(pred_bin) <- c("Mean",paste0("X.",c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)))
# 
# pred_pos <- apply(pars$pred_pos,2,quantile,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# pred_pos_mean <- colMeans(pars$pred_pos)
# pred_pos <- pred_pos %>% t() %>% as.data.frame() %>% bind_cols(pred_pos_mean,.) %>% as.data.frame()
# colnames(pred_pos) <- c("Mean",paste0("X.",c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)))

#pars2 <- pars
# pars$pred_bin <- NULL
# pars$pred_pos <- NULL
# pars$log_lik <- NULL

#### WRITE TO FILE
Output.qpcr <- list(
  # STAN MODEL ASSOCIATED THINGS
  stanMod = stanMod, 
  stanMod_summary = stanMod_summary,
  loo_val = loo_val, #loo summary
  log_lik_1 = log_lik_1, #log_likelihood object for summary
  r_eff = r_eff,
  Iter=Iter,
  N_knots_all = N_knots_all,
  #stanMod_summary_reg = stanMod_summary_reg,
  stanMod_summary_parts = stanMod_summary_parts,
  #stanMod_summary_D = stanMod_summary_D,
  samp = pars, samp_params=samp_params,
  #TRACE = TRACE,
  SPECIES = SP,
  MODEL.TYPE = MODEL.TYPE,
  MODEL.VAR = MODEL.VAR,
  MODEL.ID = MODEL.ID,
  # Input Data
  brms.object =brms.object, # for constructing the model.
  N_knots_all = N_knots_all,# for constructing the model.
  dat.station.id.trim=dat.station.id.trim,
  dat.sample.id=dat.sample.id,
  dat.id =dat.id,
  dat.id.control=dat.id.control,
  dat.inhibit=dat.inhibit,
  dat.control.field.neg = dat.control.field.neg,
  dat.control=dat.control,
  dat.stand.bin = dat.stand.bin,
  dat.stand.pos = dat.stand.pos,
  dat.obs.bin =dat.obs.bin,
  dat.obs.pos = dat.obs.pos,
  INHIBIT.LIMIT = INHIBIT.LIMIT,
  STATION.DEPTH=STATION.DEPTH,
  SAMPLES=SAMPLES,
  SAMPLES.CONTROL=SAMPLES.CONTROL,
  PCR=PCR,
  N_station_depth=N_station_depth,
  N_sample = N_sample,   # Number of site-month-bottle combinations.
  N_pcr    = N_pcr,    # Number of PCR plates
  dat_raster_fin = dat_raster_fin,
  NO.SURFACE = NO.SURFACE
  #dat_raster_trim = dat_raster_trim
  #Predictions
  # pred_bin = pred_bin,
  # pred_pos = pred_pos
)

base.dir <- "/Users/ole.shelton/Github/eDNA-Hake/"
setwd(base.dir)
setwd("./Stan Model Fits/")
save(Output.qpcr,file=paste("qPCR 2019",SP,MODEL.TYPE,MODEL.ID,MODEL.VAR,"Fitted NO.SURFACE=",NO.SURFACE,".RData"))

