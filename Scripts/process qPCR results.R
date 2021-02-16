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

# Working directories
base.dir <- "/Users/ole.shelton/Github/eDNA-Hake/"
data.dir <- paste0(base.dir,"Data")
plot.dir <- paste0(base.dir,"Plots and figures")

# Pull in qPCR data, qPCR standards, sample id information
setwd(data.dir)
dat.all <- read.csv("./qPCR/Hake eDNA 2019 qPCR results 2021-01-04 results.csv")
dat.stand <- read.csv("./qPCR/Hake eDNA 2019 qPCR results 2020-01-04 standards.csv")
dat.sample.id <- read.csv("./Hake eDNA 2019 qPCR results 2020-12-15 sample details.csv")
dat.station.id <- read.csv("./CTD_hake_data_10-2019.csv")

######################################################
# PROCESS CTD LOCATION DATA FIRST
# Modify and clean CTD data 
  # fix lat-lon data
  dat.station.id <- dat.station.id %>% filter(Button == "CTD at Depth" | Station.. == "49-9") %>% 
                      mutate(lat=as.character(SciGPS.Lat), lon=as.character(SciGPS.Lon)) %>%
                      mutate(Lat=substr(lat,1,nchar(lat)-1),Lon=substr(lon,1,nchar(lon)-1)) %>%
                      mutate(lat.deg=as.numeric(substr(Lat,1,2)),lon.deg = as.numeric(substr(Lon,1,3)),
                             lat.dec=as.numeric(substr(Lat,3,nchar(Lat)))/60,lon.dec=as.numeric(substr(Lon,4,nchar(Lon)))/60) %>%
                      mutate(lat = lat.deg + lat.dec, lon= -(lon.deg + lon.dec)) %>% 
                      dplyr::select(-Lat,-Lon,-lat.deg,-lat.dec,-lon.deg,-lon.dec)
  
  # Fix CTD station data.
  dat.station.id <- dat.station.id %>% mutate(station=as.character(Station..)) %>% 
                            mutate(station = case_when(station=="77-0" ~ "77-MT505",
                                                       station=="78-0" ~ "78-MT506",
                                                       grepl("CTD ",station) ~ substr(station,5,nchar(station)),
                                                       grepl("CTD",station) ~ substr(station,4,nchar(station)),
                                                       TRUE ~ station) ) %>%
                            mutate(transect=substr(station,1,2),
                                   transect=case_when(grepl("-",transect) ~ substr(transect,1,1),
                                                      TRUE ~ transect))

  dat.station.id <- dat.station.id %>% mutate(date= as.Date(Date,"%m/%d/%y"),year=year(date),month=month(date),day=day(date)) %>% 
                              rename(water.depth=EK.38kHz.DepthBelowSurface.M.VALUE) 
    
        d.temp <- dat.station.id %>% group_by(date,year,month,day,station, transect) %>%
                    summarise(m.lat=mean(lat),m.lon=mean(lon),m.water.depth=mean(water.depth)) %>%
                    rename(lat=m.lat,lon=m.lon,water.depth=m.water.depth)

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
  dat.sample.control.id <- dat.sample.id %>% mutate(date= as.Date(Date.UTC,"%m/%d/%y"),year=year(date),month=month(date),day=day(date)) %>%
                                  dplyr::select(sample=Tube..,date,year,month,day,drop.sample,field.negative.type,volume = water.filtered.L,)
  dat.sample.id  <- dat.sample.id %>% dplyr::select(sample=Tube..,
                                                  station=CTD.cast,
                                                  Niskin,
                                                  depth,
                                                  drop.sample,
                                                  field.negative.type,
                                                  volume = water.filtered.L,
                                                  Fluor,
                                                  Zymo=Zymo.columns)

    dat.id <- full_join(dat.sample.id,dat.station.id.trim) %>% 
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

###########################################################################
## DECLARED SPECIES OF INTEREST
SP <- "hake" # options: hake, lamprey, eulachon
###########################################################################

# STANDARDS 
THESE <- c("qPCR","well","sample","type","IPC_Ct","inhibition_rate")
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
  uni.wash <- dat.wash %>% group_by(station,depth,Niskin,drop.sample) %>% summarise(N=length(station)) %>% mutate(status="washed") %>% ungroup()

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
dat.ntc <- dat.all %>% filter(type=="ntc") %>% mutate(IPC_Ct = as.numeric(as.character(IPC_Ct))) %>%
                      group_by(qPCR) %>% 
                      dplyr::summarise(mean.ntc = mean(IPC_Ct),sd.ntc=sd(IPC_Ct))

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
                  filter(!drop.sample == "Y") # drop samples that hit the bench top

dat.control <- dat.control %>% mutate(vol.standard = volume/2.5)
dat.control.field.neg <-   dat.control %>% filter(type=="field_neg")


# Merge in ntc to data to flag inhibited samples 
INHIBIT.LIMIT <- 0.5

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
dat.samp <- dat.samp %>% filter(!sample %in% c(1298,
                                               1622,
                                               1552,
                                               1419,
                                               1326,
                                               536,
                                               55
                                               )) # THESE ARE CRAZY OUTLIERS ()

######

dat.control.field.neg <- left_join(dat.control.field.neg,dat.ntc) %>% mutate(mean.ntc = as.numeric(as.character(mean.ntc))) %>%
                          mutate(inhibit.val = IPC_Ct-mean.ntc,
                                  #inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT & inhibit.val > -INHIBIT.LIMIT,0,1),
                                 inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT ,0,1),
                                  Ct_bin = ifelse(Ct>0,1,0)) %>%
                          left_join(.,PCR)

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

p.1 <- ggplot(C) +
  geom_point(aes(x=x.1,y=x.2),alpha=0.5) +
  geom_abline(intercept=c(0,0),slope = c(2,1),color="red",linetype=c("solid","dashed")) +
  xlab("DNA copies/ul at 1:10 dilution") +
  ylab("DNA copies/ul at 1:5 dilution")

p.2 <- ggplot(C) +
  geom_point(aes(x=x.1,y=x.5),alpha=0.5) +
  geom_abline(intercept=c(0,0),slope = c(5,1),color="red",linetype=c("solid","dashed")) +
  xlab("DNA copies/ul at 1:10 dilution") +
  ylab("DNA copies/ul at 1:2 dilution") +
  xlim(c(0,120)) +
  ylim(c(0,120))

p.3 <- ggplot(C) +
  geom_point(aes(x=x.2,y=x.5),alpha=0.5) +
  geom_abline(intercept=c(0,0),slope = c(2.5,1),color="red",linetype=c("solid","dashed")) +
  xlab("DNA copies/ul at 1:5 dilution") +
  ylab("DNA copies/ul at 1:2 dilution") +
  xlim(c(0,160)) +
  ylim(c(0,160))

grid.arrange(p.1,p.2,p.3,nrow=2)

# This strongly suggests that the 0.5 dilution did not work all that well (lower copies than expected)
# So if we exclude all of the 0.5 and inspect the results
A <- dat.samp %>% group_by(qPCR,dilution) %>% summarise(N=length(dilution))

dat.samp <- dat.samp %>% filter(dilution !=0.5)
B <- dat.samp %>% group_by(qPCR,dilution) %>% summarise(N=length(dilution))  %>% as.data.frame()

########## CHECK ON HOW MANY SAMPLES WE HAVE ZERO REPLICATES FOR  

# Examine how many samples are retained out of the intial number done.
dat.samp.ok <- dat.samp %>% group_by(sample,inhibit.bin) %>% summarise(N=length(sample)) %>% as.data.frame()

dat.samp.begin %>% nrow() #1844
dat.samp.ok %>% nrow() #1841 .... lost 3 samples through filtering.

################################### Making indicators and indexs for stan.
STATION.DEPTH <- data.frame(station.depth=unique(dat.samp$station.depth),station_depth_idx = 1:length(unique(dat.samp$station.depth)))

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

SAMPLES.CONTROL <- data.frame(sample=unique(dat.control.field.neg$sample),
                              sample_control_idx = 1:length(unique(dat.control.field.neg$sample)))

dat.samp <- left_join(dat.samp,SAMPLES) %>% left_join(.,STATION.DEPTH) 
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

#####################################3
#####################################3
#####################################3
#####################################3
#####################################3
#####################################3

# Derive Estimates of Concentration 

OFFSET = 0 # Value to improve Fitting in STAN

stan_data = list(
  "bin_stand"     = dat.stand.bin$Ct_bin,
  "pos_stand"     = dat.stand.pos$Ct,
  "D_bin_stand"   = dat.stand.bin$log10_copies,
  "D_pos_stand" =   dat.stand.pos$log10_copies,
  
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
  
  # Indices for station-depth combination
  "station_depth_idx" = STATION.DEPTH$station_depth_idx,
  
  # Control Values
  "bin_control"    = dat.control.bin$Ct_bin,
  "pos_control"    = dat.control.pos$Ct,
  "bin_log_vol_control"= bin_log_vol_control,
  "pos_log_vol_control"= pos_log_vol_control,  
  "bin_log_dilution_control"= bin_log_dilution_control,
  "pos_log_dilution_control"= pos_log_dilution_control,  
  
  # Indices and counters
  "N_pcr"    = N_pcr,    # Number of PCR plates
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
  "pcr_control_pos_idx" =   dat.control.pos$plate_idx,

  # Control samples  
  "sample_control_idx" = SAMPLES.CONTROL$sample_control_idx,
  
  # Indices for Samples
  "sample_pos_idx" = dat.obs.pos$sample_idx,
  "sample_bin_idx" = dat.obs.bin$sample_idx,
  "sample_control_pos_idx" = dat.control.pos$sample_control_idx,
  "sample_control_bin_idx" = dat.control.bin$sample_control_idx,

  "station_depth_pos_idx" = dat.obs.pos$station_depth_idx,
  "station_depth_bin_idx" = dat.obs.bin$station_depth_idx,
  
  # "singleton_pos_idx" = dat.obs.pos$singleton_idx,
  # "singleton_bin_idx" = dat.obs.bin$singleton_idx,

  "wash_pos_idx" = dat.obs.pos$wash_idx,
  "wash_bin_idx" = dat.obs.bin$wash_idx,
    
  #Offset of density for improving fitting characteristics
  "OFFSET" = OFFSET
)

stan_pars = c(
  "beta_0", # intercept for standards
  "beta_1", # slope for standards
  "phi_0",  # logit intercept for standards
  "phi_1",  # logit slope for standard,
  "wash_offset", # parameter for concentration offset from washing with 30% EtOH instead of 70%.
  
  "mu_contam", # log-mean of average contamination
  "sigma_contam", # lod-sd of contamination

  "D",     # Latent variable for Log-count in each station-location combination
  "D_error",
  "D_contam",
  "D_control", # Latent variable for log-count in field negative controls.
  
  "sigma_stand_int", # variability among standards regression.
  "sigma_pcr",     # variability among samples, given individual bottle, site, and month 
  
  "delta",  # random effect for sample #
  "tau_sample"   # sd among samples given station-depth 

  # random effect for each bottle around the site-month mean.
  
  #"sigma_stand_slope", # variability among standards regression.
  #"sigma_stand_slope2", # variability among standards regression.
  
  # "sigma_stand_int_bar", # variability among standards regression.
  # "sigma_stand_slope_bar", # variability among standards regression.
  # "sigma_stand_int_sd", # variability among standards regression.
  # "sigma_stand_slope_sd", # variability among standards regression.

  # "geom_month_index",
  # "arith_month_index"
)   

### INTIAL VALUES
  stan_init_f1 <- function(n.chain,N_pcr,N_station_depth,N_control_sample){ 
    A <- list()
    for(i in 1:n.chain){
      A[[i]] <- list(
        
        sigma_stand_int = runif(1,0.01,2),
        #sigma_stand_slope = runif(1,-1,-0.1),
        # beta_0_bar = runif(1,20,40),
        # beta_1_bar = rnorm(1,-3,1),
        beta_0 = runif(N_pcr,30,40),
        beta_1 = rnorm(N_pcr,-4,1),
        wash_offset = rnorm(1,-2,1),
        
        # phi_0_bar = runif(1,10,25),
        # phi_1_bar = rnorm(1,5,1),
        phi_0  = runif(N_pcr,0,5),
        phi_1  = rnorm(N_pcr,3,0.5),
        sigma_pcr = runif(1,0.01,0.4),
        tau_sample = runif(1,0.01,0.2),
        D      = rnorm(N_station_depth,0,2),
        D_control  = rnorm(N_control_sample,0,2)
        # gamma  = rnorm(N_hybrid,0,2)

      )
    }
    return(A)
  }

#################################################################### 
#################################################################### 
##### STAN
#################################################################### 
#################################################################### 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
  
N_CHAIN = 3
Warm = 1000
Iter = 5000
Treedepth = 10
Adapt_delta = 0.60

LOC <- paste0(base.dir,"Scripts/Stan Files/")
setwd(LOC)

stanMod = stan(file = "qPCR_Hake.stan" ,data = stan_data, 
               verbose = FALSE, chains = N_CHAIN, thin = 2, 
               warmup = Warm, iter = Warm + Iter, 
               control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
               pars = stan_pars,
               boost_lib = NULL,
               sample_file = "./Output files/test.csv",
               init = stan_init_f1(n.chain=N_CHAIN,
                                   N_pcr= N_pcr,
                                   N_station_depth = N_station_depth,
                                   N_control_sample = N_control_sample)
)

#################################################################### 
#################################################################### 
#################################################################### 

pars <- rstan::extract(stanMod, permuted = TRUE)
# get_adaptation_info(stanMod)
samp_params <- get_sampler_params(stanMod)
#samp_params 
stanMod_summary <- summary(stanMod)$summary
stanMod_summary_reg <- summary(stanMod,pars=c("beta_0",
                               "beta_1",
                               "phi_0",
                               "phi_1"))$summary
stanMod_summary_param <- summary(stanMod,pars=c(
                                         "wash_offset",
                                         "mu_contam",
                                         "sigma_contam",
                                         "tau_sample",
                                         "sigma_stand_int", # variability among standard regression.
                                         "sigma_pcr"))$summary     # variability among samples, given individual bottle, site, and month ))$summary
stanMod_summary_D <- summary(stanMod,pars="D")$summary
ID <- rownames(stanMod_summary_D)
stanMod_summary_D <- stanMod_summary_D %>%  as.data.frame() %>% mutate(ID=ID) %>% arrange(desc(Rhat))

#round(stanMod_summary,2)

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
TRACE[[as.name("D")]] <- traceplot(stanMod,pars=c("lp__",paste0("D[",sample(1:N_sample,8),"]")),inc_warmup=FALSE)
TRACE[[as.name("Phi0")]] <- traceplot(stanMod,pars=c("lp__","phi_0"),inc_warmup=FALSE)
TRACE[[as.name("Phi1")]] <- traceplot(stanMod,pars=c("lp__","phi_1"),inc_warmup=FALSE)
TRACE[[as.name("Beta0")]] <- traceplot(stanMod,pars=c("lp__","beta_0"),inc_warmup=FALSE)
TRACE[[as.name("Beta1")]] <- traceplot(stanMod,pars=c("lp__","beta_1"),inc_warmup=FALSE)
TRACE[[as.name("Var")]] <- traceplot(stanMod,pars=c("lp__","tau_sample","sigma_stand_int","sigma_pcr"),inc_warmup=FALSE)
TRACE[[as.name("Contam")]] <- traceplot(stanMod,pars=c("lp__","wash_offset","mu_contam","sigma_contam"),inc_warmup=FALSE)

# pull 9 random D values.


TRACE$Var
TRACE$Contam

#pairs(stanMod, pars = c(base_params), log = FALSE, las = 1)

B1 <- apply(pars$beta_0,2,mean)
B2 <- apply(pars$beta_1,2,mean)

P0 <- apply(pars$phi_0,2,mean)
P1 <- apply(pars$phi_1,2,mean)

V0 <-  mean(pars$sigma_stand_int)
#V1 <- mean(pars$sigma_stand_slope)
#V2 <- mean(pars$sigma_stand_slope2)

# Plot regression against Standard
X <- seq(0,5,length.out=1000) - OFFSET 
Y <- t(B1 + B2 %*% t(X ))

STAND.REG <- data.frame(X=X,Y=Y)
STAND.REG <- melt(STAND.REG,id.vars="X",value.name="Y")

x.lim=c(min(X),max(X))
y.lim=c(20,40)
for(i in 1:ncol(Y)){
  plot(Y[,i]~X,xlim=x.lim,ylim=y.lim,type="l",col=2)
  par(new=T)
}
plot(dat.stand.pos$Ct~ dat.stand.pos$log10_copies,xlim=x.lim,ylim=y.lim)


#BREAKS <- c(-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5)
#LABS   <- BREAKS + OFFSET


lab.temp<- data.frame(variable=paste0("Y.",PCR$plate_idx), qPCR=PCR$qPCR)
STAND.REG <-left_join(STAND.REG,lab.temp)

stand.plot <- ggplot(dat.stand.pos) +
  geom_point(aes(x=log10_copies - OFFSET ,y=Ct,color=qPCR),alpha=0.75) +
  scale_shape_discrete(name="qPCR Plate") +
  theme_bw()
stand.plot <- stand.plot +
  geom_line(data=STAND.REG,aes(x=X,y=Y,color=qPCR)) +
  scale_color_discrete(name="qPCR Plate") +
  ylab("PCR cycle")  +
  xlab("log10 copies DNA") +
  facet_wrap(~qPCR)
#scale_x_continuous(labels = paste0("1e",LABS),limits = c(-2,4))
stand.plot

# Plot occurrence of standard
Y <- t(P0 + P1 %*% t(X ))
LOGIT <- data.frame(X=X,Y=plogis(Y))
LOGIT <- melt(LOGIT,id.vars="X",value.name="Y")

LOGIT <-left_join(LOGIT,lab.temp)

stand.plot.pres <- ggplot(dat.stand.bin) +
  geom_jitter(aes(x=log10_copies - OFFSET,y=Ct_bin,color=qPCR),alpha=0.75,width=0,height=0.05) +
  geom_line(data=LOGIT,aes(y=Y,x=X,color=qPCR)) +
  theme_bw() +
  scale_shape_discrete(name="qPCR Plate") +
  scale_color_discrete(name="qPCR plate") +
  xlab("log10 copies DNA") +
  ylab("Amplification success") +
  scale_y_continuous(breaks=c(0,1),labels = c("No","Yes"))
#scale_x_continuous(breaks=BREAKS,labels = paste0("1e",LABS),limits = c(-2.5,4))

stand.plot.pres

################################################################################################
################################################################################################
################################################################################################
################################################################################################
##### Extract data of interest, save to file for use elsewhere
################################################################################################
################################################################################################
################################################################################################
###############################################################################################

# THIS IS THE FACTOR FOR EXPANDING FROM COPIES PER ul to COPIES PER L 
# Based on sampling volume of 2.5 L

EXPAND <- 40
LOG.EXPAND <- log10(EXPAND)

PROBS <- c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975)
Log   <- paste0("log",PROBS)
station_depth_out <- data.frame(station_depth_idx= 1:ncol(pars$D), 
                         Mean.log=apply(pars$D,2,mean),
                         Sd.log=apply(pars$D,2,sd),
                         Log.Val=data.frame(t(apply(pars$D,2,quantile,probs=PROBS))),
                         Mean=apply(10^pars$D,2,mean),
                         Sd=apply(10^pars$D,2,sd),
                         Val=data.frame(t(apply(10^pars$D,2,quantile,probs=PROBS))))

station_depth_out_liter <- data.frame(station_depth_idx= 1:ncol(pars$D), 
                                Mean.log=apply(pars$D+LOG.EXPAND,2,mean),
                                Sd.log=apply(pars$D+LOG.EXPAND,2,sd),
                                Log.Val=data.frame(t(apply(pars$D+LOG.EXPAND,2,quantile,probs=PROBS))),
                                Mean=apply(10^(pars$D+LOG.EXPAND),2,mean),
                                Sd=apply(10^(pars$D+LOG.EXPAND),2,sd),
                                Val=data.frame(t(apply(10^(pars$D+LOG.EXPAND),2,quantile,probs=PROBS))))

sample_contam_error_out <- data.frame(sample_idx= 1:ncol(pars$D_error), 
                                Mean.log=apply(pars$D_error,2,mean),
                                Sd.log=apply(pars$D_error,2,sd),
                                Log.Val=data.frame(t(apply(pars$D_error,2,quantile,probs=PROBS))),
                                Mean=apply(10^pars$D_error,2,mean),
                                Sd=apply(10^pars$D_error,2,sd),
                                Val=data.frame(t(apply(10^pars$D_error,2,quantile,probs=PROBS))))

sample_contam_error_out_liter <- data.frame(sample_idx= 1:ncol(pars$D_error), 
                                      Mean.log=apply(pars$D_error+LOG.EXPAND,2,mean),
                                      Sd.log=apply(pars$D_error+LOG.EXPAND,2,sd),
                                      Log.Val=data.frame(t(apply(pars$D_error+LOG.EXPAND,2,quantile,probs=PROBS))),
                                      Mean=apply(10^(pars$D_error+LOG.EXPAND),2,mean),
                                      Sd=apply(10^(pars$D_error+LOG.EXPAND),2,sd),
                                      Val=data.frame(t(apply(10^(pars$D_error+LOG.EXPAND),2,quantile,probs=PROBS))))

sample_contam_total_out <- data.frame(sample_idx= 1:ncol(pars$D_contam), 
                                      Mean.log=apply(pars$D_contam,2,mean),
                                      Sd.log=apply(pars$D_contam,2,sd),
                                      Log.Val=data.frame(t(apply(pars$D_contam,2,quantile,probs=PROBS))),
                                      Mean=apply(10^pars$D_contam,2,mean),
                                      Sd=apply(10^pars$D_contam,2,sd),
                                      Val=data.frame(t(apply(10^pars$D_contam,2,quantile,probs=PROBS))))

sample_contam_total_out_liter <- data.frame(sample_idx= 1:ncol(pars$D_contam), 
                                            Mean.log=apply(pars$D_contam+LOG.EXPAND,2,mean),
                                            Sd.log=apply(pars$D_contam+LOG.EXPAND,2,sd),
                                            Log.Val=data.frame(t(apply(pars$D_contam+LOG.EXPAND,2,quantile,probs=PROBS))),
                                            Mean=apply(10^(pars$D_contam+LOG.EXPAND),2,mean),
                                            Sd=apply(10^(pars$D_contam+LOG.EXPAND),2,sd),
                                            Val=data.frame(t(apply(10^(pars$D_contam+LOG.EXPAND),2,quantile,probs=PROBS))))

field_neg_out <- data.frame(sample_control_idx= 1:ncol(pars$D_control), 
                                Mean.log=apply(pars$D_control,2,mean),
                                Sd.log=apply(pars$D_control,2,sd),
                                Log.Val=data.frame(t(apply(pars$D_control,2,quantile,probs=PROBS))),
                                Mean=apply(10^pars$D_control,2,mean),
                                Sd=apply(10^pars$D_control,2,sd),
                                Val=data.frame(t(apply(10^pars$D_control,2,quantile,probs=PROBS))))

field_neg_out_liter <- data.frame(sample_control_idx= 1:ncol(pars$D_control), 
                                  Mean.log=apply(pars$D_control+LOG.EXPAND,2,mean),
                                  Sd.log=apply(pars$D_control+LOG.EXPAND,2,sd),
                                  Log.Val=data.frame(t(apply(pars$D_control+LOG.EXPAND,2,quantile,probs=PROBS))),
                                  Mean=apply(10^(pars$D_control+LOG.EXPAND),2,mean),
                                  Sd=apply(10^(pars$D_control+LOG.EXPAND),2,sd),
                                  Val=data.frame(t(apply(10^(pars$D_control+LOG.EXPAND),2,quantile,probs=PROBS))))

delta_out <- data.frame(sample_idx= 1:ncol(pars$delta), 
                                Mean.log=apply(pars$delta,2,mean) ,
                                Sd.log=apply(pars$delta,2,sd),
                                Log.Val=data.frame(t(apply(pars$delta,2,quantile,probs=PROBS))),
                                Mean=apply(10^pars$delta,2,mean),
                                Sd=apply(10^pars$delta,2,sd),
                                Val=data.frame(t(apply(10^pars$delta,2,quantile,probs=PROBS))))


Output.qpcr <- list(stanMod = stanMod, stanMod_summary = stanMod_summary,samp = pars, samp_params=samp_params,
                    TRACE = TRACE,
                    SPECIES = SP,
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
                    OFFSET = OFFSET,
                    base_params =base_params,
                    STATION.DEPTH=STATION.DEPTH,
                    SAMPLES=SAMPLES,
                    SAMPLES.CONTROL=SAMPLES.CONTROL,
                    PCR=PCR,
                    N_station_depth=N_station_depth,
                    N_sample = N_sample,   # Number of site-month-bottle combinations.
                    N_pcr    = N_pcr,    # Number of PCR plates
                    stand.plot = stand.plot,
                    stand.plot.pres = stand.plot.pres,
                    station_depth_out=station_depth_out,
                    station_depth_out_liter=station_depth_out_liter,
                    field_neg_out=field_neg_out,
                    field_neg_out_liter=field_neg_out_liter,
                    delta_out=delta_out)

setwd(base.dir)
setwd("./Stan Model Fits/")
save(Output.qpcr,file=paste("qPCR 2019",SP,"Fitted.RData"))
#################################################################
##################################
