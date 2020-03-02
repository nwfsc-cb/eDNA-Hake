# This script is for reading in qPCR data, merging it with locations, and starting a set of plots.  
# It will likely be called in a rMarkdown file later.

# Libraries
library(tidyverse)
library(marmap)
library(ggplot2)
library(rstan)
library(lubridate)
library(reshape2)

# Working directories
base.dir <- "/Users/ole.shelton/Github/eDNA-Hake/"
data.dir <- paste0(base.dir,"Data")
plot.dir <- paste0(base.dir,"Plots and figures")

# Pull in qPCR data, qPCR standards, sample id information
setwd(data.dir)
dat.all <- read.csv("./qPCR/Hake_eDNA_2019_qPCR_results.csv")
dat.stand <- read.csv("./qPCR/Hake_eDNA_2019_qPCR_standards.csv")
dat.station.id <- read.csv("./CTD_hake_data_10-2019.csv")
dat.sample.id  <- read.csv("./2019 Hake- Shimada cruise - eDNA water sampling.csv")

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

  dat.station.id <- dat.station.id %>% mutate(date= as.Date(Date,"%m/%d/%y"),year=year(date),month=month(date),day=day(date))
  
  dat.station.id.trim <- dat.station.id %>% dplyr::select(date,year,month,day, lat,lon,station, transect)
                            # mutate(transect = as.numeric(substr(Station,1,2)))\

  
  dat.sample.control.id <- dat.sample.id %>% mutate(date= as.Date(Date.UTC,"%m/%d/%y"),year=year(date),month=month(date),day=day(date)) %>%
                                  dplyr::select(sample=Tube..,date,year,month,day)
  dat.sample.id  <- dat.sample.id %>% dplyr::select(sample=Tube..,
                                                  station=CTD.cast,
                                                  Niskin,
                                                  depth,
                                                  volume = water.filtered..L.,
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

dat.stand <- dat.stand %>% dplyr::select(THESE,grep(SP,colnames(dat.stand))) %>% 
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
dat.samp <- dat.all %>% dplyr::select(THESE,useful,Zymo,grep(SP,colnames(dat.all))) %>%
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


dat.ntc <- dat.all %>% filter(type=="ntc") %>% mutate(IPC_Ct = as.numeric(as.character(IPC_Ct))) %>%
                      group_by(qPCR) %>% 
                      dplyr::summarise(mean.ntc = mean(IPC_Ct),sd.ntc=sd(IPC_Ct))

dat.control <- dat.all %>% filter(grepl("neg",type)|type=="extc") %>%
                    dplyr::select(THESE,useful,Zymo,grep(SP,colnames(dat.all))) %>%
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
                  
dat.control <- left_join(dat.control,dat.sample.control.id,by="sample")
dat.control.field.neg <-   dat.control %>% filter(type=="field_neg")


# Merge in ntc to data to flag inhibited samples 
INHIBIT.LIMIT <- 1

dat.samp <- left_join(dat.samp,dat.ntc) %>% mutate(mean.ntc = as.numeric(as.character(mean.ntc))) %>%
                  mutate(inhibit.val = IPC_Ct-mean.ntc,
                         inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT & inhibit.val > -INHIBIT.LIMIT,0,1),
                         Ct_bin = ifelse(Ct>0,1,0)) %>%
                  left_join(.,PCR) %>%
                  mutate(station.depth = paste0(station,".",depth))

dat.control.field.neg <- left_join(dat.control.field.neg,dat.ntc) %>% mutate(mean.ntc = as.numeric(as.character(mean.ntc))) %>%
                          mutate(inhibit.val = IPC_Ct-mean.ntc,
                                  inhibit.bin=ifelse(inhibit.val < INHIBIT.LIMIT & inhibit.val > -INHIBIT.LIMIT,0,1),
                                  Ct_bin = ifelse(Ct>0,1,0)) %>%
                          left_join(.,PCR)
                          

# cull dat.samp of inhibited samples
dat.inhibit <- dat.samp %>% filter(useful=="NO" | inhibit.bin==1) 
dat.samp    <- dat.samp %>% filter(!useful=="NO", inhibit.bin==0) 
dat.control.field.neg    <- dat.control.field.neg %>% filter(inhibit.bin==0) 

SAMPLES <- data.frame(sample=unique(dat.samp$sample),sample_idx = 1:length(unique(dat.samp$sample)))
STATION.DEPTH <- data.frame(station.depth=unique(dat.samp$station.depth),station_depth_idx = 1:length(unique(dat.samp$station.depth)))

SAMPLES.CONTROL <- data.frame(sample=unique(dat.control.field.neg$sample),
                              sample_control_idx = 1:length(unique(dat.control.field.neg$sample)))

dat.samp <- left_join(dat.samp,SAMPLES) %>% left_join(.,STATION.DEPTH) 
dat.control.field.neg <- left_join(dat.control.field.neg,SAMPLES.CONTROL)

# This identifies where there is only one sample for each station-depth combination.  This is important for model identifiability.
dat.samp   <-  dat.samp %>% group_by(station_depth_idx) %>% summarise(N=length(unique(sample))) %>% 
              mutate(singleton_idx=ifelse(N==1,0,1)) %>% dplyr::select(-N) %>%
              left_join(dat.samp,.)

### THERE IS SOME WEIRDNESS WITH SAMPLES 1306, 1361, 1362, and 412.  I think I have resolved them.
dat.stand.bin <- dat.stand
dat.stand.pos <- dat.stand %>% filter(Ct_bin==1)
N_stand_bin <- nrow(dat.stand.bin)
N_stand_pos <- nrow(dat.stand.pos)

dat.obs.bin <- dat.samp
dat.obs.pos <- dat.samp %>% filter(Ct_bin==1)
N_obs_bin <- nrow(dat.obs.bin)
N_obs_pos <- nrow(dat.obs.pos)

dat.control.bin <- dat.control.field.neg
dat.control.pos <- dat.control.field.neg %>% filter(Ct_bin==1)
N_control_bin <- nrow(dat.control.bin)
N_control_pos <- nrow(dat.control.pos)

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

OFFSET = 0 # Value to imrpove Fitting in STAN


stan_data = list(
  "bin_stand"     = dat.stand.bin$Ct_bin,
  "pos_stand"     = dat.stand.pos$Ct,
  "D_bin_stand"   = dat.stand.bin$log10_copies,
  "D_pos_stand" =   dat.stand.pos$log10_copies,
  
  "bin_obs"    = dat.obs.bin$Ct_bin,
  "pos_obs"    = dat.obs.pos$Ct,
  
  "bin_control"    = dat.control.bin$Ct_bin,
  "pos_control"    = dat.control.pos$Ct,
  
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
  
    
  # Indices for sample bottles combination
  "sample_idx" = SAMPLES$sample_idx,
  "sample_control_idx" = SAMPLES.CONTROL$sample_control_idx,
  
  # Indices for station-depth combination
  "station_depth_idx" = STATION.DEPTH$station_depth_idx,
  
  # Indices for Samples
  "sample_pos_idx" = dat.obs.pos$sample_idx,
  "sample_bin_idx" = dat.obs.bin$sample_idx,
  "sample_control_pos_idx" = dat.control.pos$sample_control_idx,
  "sample_control_bin_idx" = dat.control.bin$sample_control_idx,

  "station_depth_pos_idx" = dat.obs.pos$station_depth_idx,
  "station_depth_bin_idx" = dat.obs.bin$station_depth_idx,
  
  "singleton_pos_idx" = dat.obs.pos$singleton_idx,
  "singleton_bin_idx" = dat.obs.bin$singleton_idx,
  
  #Offset of density for improving fitting characteristics
  "OFFSET" = OFFSET
)


stan_pars = c(
  "beta_0", # intercept for standards
  "beta_1", # slope for standards
  "phi_0",  # logit intercept for standards
  "phi_1",  # logit slope for standard,
  
  "D",     # Latent variable for Log-count in each station-location combination
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
        
        # phi_0_bar = runif(1,10,25),
        # phi_1_bar = rnorm(1,5,1),
        phi_0  = runif(N_pcr,0,5),
        phi_1  = rnorm(N_pcr,3,1),
        D      = rnorm(N_station_depth,0,2),
        D_control  = rnorm(N_control_sample,0,2),
        # gamma  = rnorm(N_hybrid,0,2)
        sigma_pcr = runif(1,0.01,0.4),
        tau_sample = runif(1,0.01,0.2)
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
  
N_CHAIN = 5
Warm = 2000
Iter = 5000
Treedepth = 11
Adapt_delta = 0.70

LOC <- paste0(base.dir,"/Scripts/Stan Files/")
setwd(LOC)

stanMod = stan(file = "qPCR_Hake.stan" ,data = stan_data, 
               verbose = FALSE, chains = N_CHAIN, thin = 5, 
               warmup = Warm, iter = Warm + Iter, 
               control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
               pars = stan_pars,
               boost_lib = NULL,
               # sample_file = "./STAN_models/Output files/test.csv",
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
round(stanMod_summary,2)

base_params <- c(
  "beta_0",
  "beta_1",
  "phi_0",
  "phi_1",
  
  #"D",
  #"delta",
  "tau_sample",

  "sigma_stand_int", # variability among standard regression.
  "sigma_pcr"     # variability among samples, given individual bottle, site, and month 
  
) 

##### MAKE SOME DIAGNOSTIC PLOTS

print(traceplot(stanMod,pars=c("lp__",base_params),inc_warmup=FALSE))

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
plot(dat.stand.pos$Ct~ dat.stand.pos$log10.density,xlim=x.lim,ylim=y.lim)


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
  xlab("log10 copies DNA") 
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
################################################################################################
PROBS <- c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975)
Log   <- paste0("log",PROBS)
station_depth_out <- data.frame(station_depth_idx= 1:ncol(pars$D), 
                         Mean.log=apply(pars$D,2,mean),
                         Sd.log=apply(pars$D,2,sd),
                         Log.Val=data.frame(t(apply(pars$D,2,quantile,probs=PROBS))),
                         Mean=apply(10^pars$D,2,mean),
                         Sd=apply(10^pars$D,2,sd),
                         Val=data.frame(t(apply(10^pars$D,2,quantile,probs=PROBS))))

field_neg_out <- data.frame(sample_control_idx= 1:ncol(pars$D_control), 
                                Mean.log=apply(pars$D_control,2,mean),
                                Sd.log=apply(pars$D_control,2,sd),
                                Log.Val=data.frame(t(apply(pars$D_control,2,quantile,probs=PROBS))),
                                Mean=apply(10^pars$D_control,2,mean),
                                Sd=apply(10^pars$D_control,2,sd),
                                Val=data.frame(t(apply(10^pars$D_control,2,quantile,probs=PROBS))))


delta_out <- data.frame(sample_idx= 1:ncol(pars$delta), 
                                Mean.log=apply(pars$delta,2,mean) ,
                                Sd.log=apply(pars$delta,2,sd),
                                Log.Val=data.frame(t(apply(pars$delta,2,quantile,probs=PROBS))),
                                Mean=apply(10^pars$delta,2,mean),
                                Sd=apply(10^pars$delta,2,sd),
                                Val=data.frame(t(apply(10^pars$delta,2,quantile,probs=PROBS))))


Output.qpcr <- list(stanMod = stanMod, stanMod_summary = stanMod_summary,samp = pars, samp_params=samp_params,
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
                    field_neg_out=field_neg_out,
                    delta_out=delta_out
)

setwd(base.dir)
setwd("./Stan Model Fits/")
save(Output.qpcr,file=paste("qPCR Hake 2019",SP,"Fitted.RData"))
####
# you need to run the above code several times.  once for each species-standard combination.

#################################################################
##################################
