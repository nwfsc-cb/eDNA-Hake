###
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

script.dir <- "/Users/ole.shelton/Github/eDNA-Hake/Scripts"
results.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Stan Model Fits"
plot.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Plots and figures"

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

# CHANGE THIS FOR switching between species.
SPECIES <- "hake" # eulachon, hake 
# CHANGE THIS FOR switching among model output.
MOD <-  "lat.long.smooth"

#load(paste("qPCR 2019",SPECIES, MOD, "7_12 Fitted.RData"))
load("qPCR 2019 hake lat.long.smooth 4_10_fix_nu_T-FIN Base_Var Fitted NO.SURFACE= FALSE .RData")
#save(Output.qpcr,file=paste("qPCR 2019",SPECIES, MOD, "Fitted.RData"))

### Cacluate posterior summaries and predictive surfaces
setwd(script.dir)
source("summarize_stan_output_PCR.R")

#########
### MAKE A BUNCH OF DATA FRAMES FOR LATER USE.
#########
### Get some of the output on sample and station IDs from file.
dat.samp <- Output.qpcr$dat.obs.bin
dat.control.field.neg <- Output.qpcr$dat.control.field.neg %>% left_join(.,field_neg_out_liter)
STATION.DEPTH <- left_join(STATION.DEPTH,station_depth_out_liter)

##### Calculate some info about the field negative controls
N.station.run       <- length(unique(c(STATION.DEPTH$station)))
N.station.depth.run <- length(unique(c(STATION.DEPTH$station.depth)))
N.sample.run        <- length(unique(c(dat.samp$sample)))

# Make a map of where we have samples.
dat.id <- Output.qpcr$dat.id
dat.station.id.trim <- Output.qpcr$dat.station.id.trim
dat.sample.id <- Output.qpcr$dat.sample.id

# combine ids into stations, samples outstanding, and processed samples
dat.A <- dat.station.id.trim %>% dplyr::select(station,lat,lon,water.depth) %>% mutate(Level="CTD Stations")
dat.B <- left_join(dat.sample.id,dat.id) %>% dplyr::select(station,lat,lon,water.depth) %>%
  group_by(station,lat,lon,water.depth) %>% summarise(N=length(station)) %>% dplyr::select(-N) %>% filter(is.na(lat)==F) %>%
  mutate(Level="Sampled Stations")
dat.C <- dat.samp %>% dplyr::select(station,lat,lon) %>% group_by(station,lat,lon) %>% summarise(N=length(station)) %>%
  dplyr::select(-N) %>% mutate(Level="Processed Stations") %>% 
  left_join(.,dat.station.id.trim %>% dplyr::select(station,water.depth) )  

### THIS IS ALL OF THE CTD STATIONS.  IF YOU NEED A LAT LONG FOR EACH STATION, MERGE IN FROM HERE.
dat.process <- bind_rows(dat.A,dat.B,dat.C) %>% filter(!station=="")
dat.process$Level <- factor(dat.process$Level, levels=c("CTD Stations","Sampled Stations","Processed Stations"))

# Calculate a Category for  water depth associated with each station.
STATION.DEPTH <- STATION.DEPTH %>% 
  mutate(water.depth.cat="",water.depth.cat = 
           case_when(bottom.depth.consensus >=2000 ~ "w_2000_plus",
                     bottom.depth.consensus < 2000 & bottom.depth.consensus >=1200 ~ "w_1200_2000",
                     bottom.depth.consensus < 1250 & bottom.depth.consensus >= 750 ~ "w_750_1200",
                     bottom.depth.consensus < 750 & bottom.depth.consensus >= 400 ~ "w_400_750",
                     bottom.depth.consensus < 400 & bottom.depth.consensus >= 210 ~ "w_210_400",
                     bottom.depth.consensus < 210  & bottom.depth.consensus >= 124 ~ "w_125_210",
                     bottom.depth.consensus < 124 & bottom.depth.consensus >= 75 ~ "w_75_125",
                     bottom.depth.consensus < 75  ~ "w_0_75", 
                     TRUE ~ water.depth.cat )) %>%
  mutate(water.depth.cat.2= water.depth.cat, 
         water.depth.cat.2 = case_when(water.depth.cat %in% c("w_2000_plus","w_1200_2000","w_750_1200") ~ "1000+",
                                       #water.depth.cat %in% c("w_750_1250") ~ "1000",
                                       water.depth.cat %in% c("w_400_750") ~ "500",
                                       water.depth.cat %in% c("w_210_400") ~ "300",
                                       water.depth.cat %in% c("w_125_210") ~ "150",
                                       water.depth.cat %in% c("w_75_125") ~ "100",
                                       water.depth.cat %in% c("w_0_75") ~ "50"))

# Repeat for the Niskin level samples
SAMPLES <- left_join(SAMPLES,
                     dat.samp %>% group_by(sample) %>% summarise(dilution = min(dilution),wash = mean(wash_idx)) %>% 
                       mutate(inhibit_cat = ifelse(dilution < 1,"Yes","No"),
                              wash_cat = ifelse(wash==1,"Yes","No")))

SAMPLES <- left_join(SAMPLES,
                     dat.id %>% dplyr::select(sample,station,lat,lon))

SAMPLES <- left_join(SAMPLES,D_delta_out_liter)

SAMPLES <- SAMPLES %>% 
  mutate(water.depth.cat="",water.depth.cat = 
           case_when(bottom.depth.consensus >=2000 ~ "w_2000_plus",
                     bottom.depth.consensus < 2000 & bottom.depth.consensus >=1200 ~ "w_1200_2000",
                     bottom.depth.consensus < 1250 & bottom.depth.consensus >= 750 ~ "w_750_1200",
                     bottom.depth.consensus < 750 & bottom.depth.consensus >= 400 ~ "w_400_750",
                     bottom.depth.consensus < 400 & bottom.depth.consensus >= 210 ~ "w_210_400",
                     bottom.depth.consensus < 210  & bottom.depth.consensus >= 124 ~ "w_125_210",
                     bottom.depth.consensus < 124 & bottom.depth.consensus >= 75 ~ "w_75_125",
                     bottom.depth.consensus < 75  ~ "w_0_75", 
                     TRUE ~ water.depth.cat )) %>%
  mutate(water.depth.cat.2= water.depth.cat, 
         water.depth.cat.2 = case_when(water.depth.cat %in% c("w_2000_plus","w_1200_2000","w_750_1200") ~ "1000+",
                                       #water.depth.cat %in% c("w_750_1250") ~ "1000",
                                       water.depth.cat %in% c("w_400_750") ~ "500",
                                       water.depth.cat %in% c("w_210_400") ~ "300",
                                       water.depth.cat %in% c("w_125_210") ~ "150",
                                       water.depth.cat %in% c("w_75_125") ~ "100",
                                       water.depth.cat %in% c("w_0_75") ~ "50"))

##### MAP MAKING COMMANDS
# Plot Limits for all samples (see Base_Map.R)
# lat.lims <- c(min(dat.id$lat,na.rm=T),max(dat.id$lat,na.rm=T))
# lon.lims <- c(min(dat.id$lon,na.rm=T),max(dat.id$lon,na.rm=T))

#lat.lims.trim <- c(37.5,max(dat.id$lat,na.rm=T))
# Call Base_map.R
setwd(script.dir)
source("Base_map.R")

sample.process <- base_map + 
  scale_color_manual(values=c(grey(0.6),"blue","red"))+
  geom_point(data=dat.process,aes(x=lon,y=lat,color=Level),size=0.8) 

sample.process.facet <- sample.process +
  facet_wrap(~Level)

####################################################################
# STANDARDS, CONTAMINATION, OTHER NUTS and BOLTS PROCESSES
#####################################################################

# MAPS OF INHIBITION AND WASHED SAMPLES
inhibit.plot <- base_map_trim +
  geom_point(data=SAMPLES,aes(x=lon,y=lat,color=inhibit_cat)) +
  scale_color_viridis_d(begin=0,end=0.8)+
  facet_wrap(~inhibit_cat)
inhibit.plot.by.depth <- base_map_trim +
  geom_point(data=SAMPLES %>% filter(is.na(depth_cat)==F) %>% mutate(depth_cat = ifelse(depth_cat==0,3,depth_cat)),aes(x=lon,y=lat,color=inhibit_cat),alpha=0.5) +
  facet_wrap(~depth_cat,nrow = 2) +
  scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) 
inhibit.plot.by.depth.few <- base_map_trim +
  geom_point(data=SAMPLES %>% filter(depth_cat %in% c(0,50,100)) ,aes(x=lon,y=lat,color=inhibit_cat),alpha=0.5) +
  facet_wrap(~depth_cat,nrow = 1) +
  scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) 

wash.plot <- base_map_trim +
  geom_point(data=SAMPLES,aes(x=lon,y=lat,color=wash_cat)) +
  scale_color_viridis_d(begin=0,end=0.8)+
  facet_wrap(~wash_cat)
wash.plot.by.depth <- base_map_trim +
  geom_point(data=SAMPLES %>% filter(is.na(depth_cat)==F),aes(x=lon,y=lat,color=wash_cat),alpha=0.5) +
  facet_wrap(~depth_cat,nrow = 2) +
  scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) 

## Calculate some summaries of how many samples were inhibited 
inhibit.summ.dat <- STATION.DEPTH %>% group_by(depth_cat,inhibit) %>% summarise(N=length(inhibit)) %>% 
  pivot_wider(values_from="N",names_from=c("inhibit"))
inhibit.summ.dat[is.na(inhibit.summ.dat)] <- 0
colnames(inhibit.summ.dat)[1:3] <- c("Depth_m","No","Yes")


###----- Plots of STANDARDS (from summarize_stan_output.R)
# stand.plot
# stand.plot.facet <- stand.plot +facet_wrap(~qPCR,ncol=6)
# stand.plot.pres

#######################################################
####### ---- Calculate CONTAMINATION
#######################################################

# Calculations for Field Negative Controls.
SAMPLES.CONTROL <- Output.qpcr$SAMPLES.CONTROL
dat.id.control <-  Output.qpcr$dat.id.control
field.neg.summ <- left_join(SAMPLES.CONTROL,field_neg_out_liter) %>% left_join(.,dat.id.control)

mu_contam_l  <- mu_contam_out$Mean[mu_contam_out$type == "liter"] 
sigma_contam_l <- sigma_contam_out$Mean

BREAKS <- seq(0,max(field.neg.summ$Mean)*1.1,length.out=30)
XLIM <- c(0,max(field.neg.summ$Mean)*1.1)
X <- seq(0.25,max(XLIM),by=0.5)
DENS <- dlnorm(X,mu_contam_l,sigma_contam_l)
DENS_mod <- (10/max(DENS))*DENS
Dens <- data.frame(X=X,DENS_mod=DENS_mod)
#REAL <- rlnorm(1e6,mu_contam,sigma_contam)  

PROB.DENS <- cbind(X,plnorm(X,mu_contam_l,sigma_contam_l))

control.hist <- ggplot() +
  #geom_histogram(data=field.neg.summ,aes(Mean,fill=field.negative.type),breaks=BREAKS) +
  geom_histogram(data=field.neg.summ,aes(Mean),breaks=BREAKS) +
  xlab("Estimated copies per L") +
  scale_fill_discrete("Type")+
  ylab("Number of negative controls") +
  xlim(XLIM) +
  geom_line(data=Dens,aes(x=X,y=DENS_mod),col="blue")+
  theme_bw()

control.hist1 <- control.hist + geom_vline(xintercept = c(20,40),color="blue",linetype="dashed")

ALPHA = 0.9 
control.by.time <- ggplot(field.neg.summ) +
  geom_point(aes(x=date,y=Mean,fill=field.negative.type,color=field.negative.type),alpha=ALPHA,shape=21) +
  geom_errorbar(aes(x=date,ymin=Val.X5.,ymax=Val.X95.,color=field.negative.type),alpha=ALPHA,width=0) +
  xlab("Date") +
  ylab("Copies per L") +
  scale_fill_discrete("Type")+
  scale_color_discrete("Type")+
  geom_hline(aes(yintercept = 20),linetype="dashed")+
  theme_bw()

## Overplot the estimated hierarchical distribution:

#control.by.time.trim <- control.by.time + ylim(c(0,1000))
detect.lev = 20
N_control_positive <- field.neg.summ %>% filter(Mean > detect.lev ) %>% nrow(.)
N_control_run <- nrow(field.neg.summ)
N_control_total <- dat.id.control %>% filter(control=="field") %>% nrow(.)

#######################################################
####### ---- Compare distribution of field samples against negative controls.
#######################################################
dat.sum.v.control <- bind_rows(SAMPLES %>% dplyr::select(sample,
                                                             Mean.log,Log.Val.X5.,Log.Val.X95.,
                                                             Mean,Val.X5.,Val.X95.) %>% mutate(cat="field") %>% as.data.frame(),
                               field.neg.summ %>% dplyr::select(sample,
                                                                Mean.log,Log.Val.X5.,Log.Val.X95.,
                                                                Mean,Val.X5.,Val.X95.) %>% mutate(cat="control") %>% as.data.frame()
                    )

BREAKS <- 10^(c(-2,seq(log10(20),5,by=0.25)))
XLIM <- c(0,max(dat.sum.v.control$Mean)*1.01)
compare.hist2 <- ggplot(dat.sum.v.control) +
  geom_histogram(aes(Mean,fill=cat),breaks=BREAKS,alpha=0.5,position="identity",) +
  xlab("Estimated Copies/L") +
  ylab("Number of Samples") +
  xlim(XLIM) +
  scale_fill_viridis_d("Category",begin=0,end=0.5) +
  scale_x_log10() +
  theme_bw()
compare.hist2

x.int <- data.frame(INT=c(20))
depth.hist.facet <-  ggplot(SAMPLES %>% mutate(depth_cat=ifelse(depth_cat==0,3,depth_cat))) +
  geom_histogram(aes(x=Mean),alpha=0.5,fill=viridis(1)) +
  xlab("Estimated Copies/L") +
  #ylab("Number of Samples") +
  scale_x_log10() +          
  #scale_fill_viridis_d("Inhibition\nCategory",begin=0.0,end=0.5) +
  facet_wrap(~depth_cat,nrow = 2) +
  geom_vline(data=x.int,aes(xintercept =INT),color="black" ,linetype=c("dashed")) +
  theme_bw()

if(SPECIES == "hake"){
  X <- seq(0.25,max(XLIM),by=0.5)
  DENS <- dlnorm(X,mu_contam_l,sigma_contam_l)
  DENS_mod <- (25/max(DENS))*DENS
  Dens <- data.frame(X=X,DENS_mod=DENS_mod)
  
  depth.hist.facet <- depth.hist.facet +
    geom_line(data=Dens,aes(x=X,y=DENS_mod),col="blue")
  depth.hist.facet
}

##########################################################
### Plot marignal distributions of wash effects
##########################################################

wash.posterior <- data.frame(wash.offset =Output.qpcr$samp$wash_offset) 


wash.param.hist <- ggplot(wash.posterior) +
              geom_density(aes(wash.offset),alpha=0.5,fill=viridis(1))+
              geom_vline(xintercept = mean(wash.posterior$wash.offset),linetype="dashed") +
              xlab(expression(omega*" (EtOH wash offset)")) +
              theme_bw()

##########################################################
## Various Variance parameters.
##########################################################

sigma.stand.param <- ggplot() +
  geom_density(aes(Output.qpcr$samp$sigma_stand_int),alpha=0.5,fill=viridis(1))+
  geom_vline(xintercept = mean(Output.qpcr$samp$sigma_stand_int),linetype="dashed") +
  xlab(expression(sigma*" (Standards)")) +
  theme_bw()

sigma.pcr.param <- ggplot() +
  geom_density(aes(Output.qpcr$samp$sigma_pcr),alpha=0.5,fill=viridis(1))+
  geom_vline(xintercept = mean(Output.qpcr$samp$sigma_pcr),linetype="dashed") +
  xlab(expression(eta*" (PCR)")) +
  theme_bw()

tau_long <- Output.qpcr$samp$tau_sample %>% as.data.frame()
colnames(tau_long) <- levels(dat.obs.bin$depth_cat_factor)
tau_long <- tau_long %>% mutate(N = 1:n()) %>% pivot_longer(!N,names_to="depth_factor",values_to="tau")

tau_long_summary <- tau_long %>% group_by(depth_factor) %>% 
                            dplyr::summarise(Mean=mean(tau),SD=sd(tau),
                                             Median=median(tau),
                                             Q.0.01 = quantile(tau,probs=c(0.01)),
                                             Q.0.025 = quantile(tau,probs=c(0.025)),
                                             Q.0.05 = quantile(tau,probs=c(0.05)),
                                             Q.0.25 = quantile(tau,probs=c(0.25)),
                                             Q.0.75 = quantile(tau,probs=c(0.75)),
                                             Q.0.95 = quantile(tau,probs=c(0.95)),
                                             Q.0.975 = quantile(tau,probs=c(0.975)),
                                             Q.0.99 = quantile(tau,probs=c(0.99)))
 
tau_long_summary$depth_factor_number <- tau_long_summary$depth_factor %>% as.numeric()

# Tau histograms
tau.sample.param <- ggplot(tau_long) +
  geom_density(aes(tau),alpha=0.5,fill=viridis(1))+
  #geom_vline(xintercept = mean(tau_long),linetype="dashed") +
  xlab(expression(tau*" (sample)")) +
  facet_wrap(~depth_factor) +
  theme_bw()

tau.sample.by.depth <- ggplot(tau_long_summary) +
  geom_point(aes(x=depth_factor_number,y=Mean)) +
  geom_errorbar(aes(x=depth_factor_number,ymin=Q.0.05,ymax=Q.0.95),width=0) +
  ylab(expression(tau))+
  xlab("Depth(m)")+
  theme_bw()



# nu.param <- ggplot() +
#   geom_density(aes(Output.qpcr$samp$nu),alpha=0.5,fill=viridis(1))+
#   geom_vline(xintercept = mean(Output.qpcr$samp$nu),linetype="dashed") +
#   xlab(expression(nu*" (sample)")) +
#   theme_bw()

sd.PCR <- Output.qpcr$samp$sigma_pcr * sqrt(3/ (3 -2))

sd.PCR.param <- ggplot() +
  geom_density(aes(sd.PCR),alpha=0.5,fill=viridis(1))+
  geom_vline(xintercept = c(median(sd.PCR),mean(sd.PCR)),linetype=c("dotdash","dashed"),color=c("blue","black")) +
  xlab(expression("SD for PCR")) +
  theme_bw()

sd_param_plots <- grid.arrange(sigma.stand.param,tau.sample.param,
             sigma.pcr.param,sd.PCR.param,
             layout_matrix= rbind(c(1,2),
                                  c(3,4)))

#######################################################
#######################################################
####### ------- Make Maps
#######################################################
#######################################################
# SAMPLES has the data on an individual niskin level
# STATION.DEPTH has data on a station-depth level.
# smooth.projections %>% names()

# NISKIN LEVEL FIRST

LEV <- c("0","50","100","150","300","500","1000+") #,"1500+")
SAMPLES$water.depth.cat.2 <- factor(SAMPLES$water.depth.cat.2,levels=LEV) 
STATION.DEPTH$water.depth.cat.2 <- factor(STATION.DEPTH$water.depth.cat.2,levels=LEV) 
# drop sample sampled at 200m 
# dat.SP <- dat.SP %>% filter(!depth==200)

# BREAKS = c(-100,-2,0,1,2,3,4,5)
# SP.p1 <- base_map_trim + 
#   geom_point(data=dat.SP,aes(x=lon,y=lat,size=Mean.log,color=inhibit_cat),shape=21) +
#   scale_size("Log Copies per L",breaks=BREAKS)+
#   scale_shape("Log Copies",solid=FALSE)+
#   facet_wrap(~depth_cat)
# SP.p1


summary(SAMPLES$Mean)
summary(STATION.DEPTH$Mean)

# NISKIN level.
lower.lim <- 20
lower.lim.inhibit <- 100
LAB<-c(lower.lim,100,1000,10000,20000,40000,100000,140000)
p2 <- base_map_trim + 
  geom_point(data=SAMPLES,aes(x=lon,y=lat,size=Mean,color=inhibit_cat),shape=21) +
  scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA))+
  geom_point(data=SAMPLES%>% filter(Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
  scale_shape("Copies / L",solid=FALSE)+
  scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) +
  facet_wrap(~depth_cat,nrow=2)
p2

p3 <- base_map_trim + 
  geom_point(data=SAMPLES,aes(x=lon,y=lat,size=Mean,color=wash_cat),shape=21) +
  scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA))+
  geom_point(data=SAMPLES%>% filter(Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
  scale_shape("Copies / L",solid=FALSE)+
  scale_color_viridis_d("Wash\nCategory",begin=0,end=0.8) +
  facet_wrap(~depth_cat,nrow=2)
p3


SAMPLES.mod <- SAMPLES %>% mutate(depth_cat =ifelse(depth_cat==0,3,depth_cat))
p4 <- base_map_trim + 
  geom_point(data=SAMPLES.mod,
             aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
  geom_point(data=SAMPLES.mod %>% filter(Mean<=lower.lim),
             aes(x=lon,y=lat),shape="x",size=1,color="black") +
  scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA)) +
  scale_shape("Copies / L",solid=FALSE) +
  #scale_color_viridis_d("Wash\nCategory",begin=0,end=0.8) +
  facet_wrap(~depth_cat,nrow=2)
p4

# Station-Depth level.
q2 <- base_map_trim + 
  geom_point(data=STATION.DEPTH,aes(x=lon,y=lat,size=Mean),shape=21,color = "red") +
  scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA))+
  geom_point(data=STATION.DEPTH %>% filter(Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
  scale_shape("Copies / L",solid=FALSE)+
  scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) +
  facet_wrap(~depth_cat,nrow=2)
q2

q3 <- base_map_trim + 
  geom_point(data=STATION.DEPTH,aes(x=lon,y=lat,size=Mean),shape=21,color = "red") +
  scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA))+
  geom_point(data=STATION.DEPTH %>% filter(Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
  scale_shape("Copies / L",solid=FALSE)+
  scale_color_viridis_d("Wash\nCategory",begin=0,end=0.8) +
  facet_wrap(~depth_cat,nrow=2)
q3

# 
# p4 <- base_map_trim + 
#   geom_point(data=SAMPLES %>% filter(inhibit==0),aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
#   geom_point(data=SAMPLES %>% filter(inhibit==1,Mean > lower.lim.inhibit),aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
#   scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA))+
#   geom_point(data=SAMPLES %>% filter(inhibit==0,Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
#   geom_point(data=SAMPLES %>% filter(inhibit==1,Mean<=lower.lim.inhibit),aes(x=lon,y=lat),shape="x",size=1,color="black") +
#   scale_shape("Copies / L",solid=FALSE)+
#   scale_color_viridis_d("Wash\nCategory",begin=0,end=0.8) +
#   facet_wrap(~depth_cat,nrow=2)
# p4

# Plot that trims out only those spots shallower than 150m.
# DEPTHS <- c(0,50)
# 
# p5 <- base_map_trim + 
#   geom_point(data=SAMPLES %>% filter(depth_cat %in% c(0,50,100,150)),
#              aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
#   geom_point(data=SAMPLES s%>% filter(depth_cat %in% c(0,50,100,150))
#                ,aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
#   scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA))+
#   geom_point(data=SAMPLES %>% filter(depth %in% c(0,50,100,150)) %>%
#                filter(inhibit==0,Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
#                filter(inhibit==1,Mean<=lower.lim.inhibit),aes(x=lon,y=lat),shape="x",size=1,color="black") +
#   scale_shape("Copies / L",solid=FALSE)+
#   scale_color_viridis_d("Wash\nCategory",begin=0,end=0.8) +
#   facet_wrap(~depth_cat,nrow=1)
# p5


######################################################3
# Smooth projections.
######################################################3
SIZE = 1.3
STROKE = 0 

z.lim = c(20,2250) 
z.breaks <-  c(20,100,250,500,1000,1500,2250)
z.lim.labs <- z.breaks

DEPTH <- levels(STATION.DEPTH$depth_cat_factor)
p_D <- list()

D_pred_background <- D_pred_log_combined %>% filter(depth_cat_factor == 0)

# Color scales for viridis
OPT = "plasma" # options are "viridis"(default), "magma", "plasma", "inferno"


for(i in 1:length(DEPTH)){
  
  LAB <- paste0(DEPTH[i],"m")
  if(DEPTH[i]==0){LAB <- "3m"}
  p_D[[as.name(paste0("x_",DEPTH[i]))]] <-
    base_map_trim_proj +
    geom_point(data= D_pred_combined %>% filter(depth_cat_factor == 0),
             aes(x=lon,y=lat),size=SIZE,fill=grey(0.6),color=grey(0.6),alpha=1,stroke=STROKE,shape=15) +
    geom_point(data= D_pred_combined %>% filter(depth_cat_factor == DEPTH[i]),
            aes(x=lon,y=lat,color=Mean,fill=Mean),alpha=1,size=SIZE,stroke=STROKE,shape=15) +
    geom_point(data= D_pred_combined %>% filter(depth_cat_factor == DEPTH[i],Mean < z.lim[1]),
             aes(x=lon,y=lat),alpha=1,size=SIZE,stroke=STROKE,shape=15,color=viridis(1,begin=0,end=0.001)) +
    geom_point(data= D_pred_combined %>% filter(depth_cat_factor == DEPTH[i],Mean > z.lim[2]),
             aes(x=lon,y=lat),alpha=1,size=SIZE,stroke=STROKE,shape=15,
             color=viridis(1,begin=0.999,end=1),fill=viridis(1,begin=0.999,end=1)) +
    scale_color_viridis_c(option=OPT,trans="sqrt",limits=z.lim,breaks=z.breaks,labels=z.lim.labs,name=expression("DNA Copies L"^-1)) +
    scale_fill_viridis_c(option=OPT,trans="sqrt",limits=z.lim,breaks=z.breaks,labels=z.lim.labs,name=expression("DNA Copies L"^-1)) +
    ggtitle(LAB) +
    theme_bw() 
  
  p_D[[as.name(paste0("x_",DEPTH[i],"_dots"))]] <- 
        p_D[[as.name(paste0("x_",DEPTH[i]))]] +
        geom_point(data=STATION.DEPTH %>% filter(depth_cat_factor == DEPTH[i]),aes(x=lon,y=lat),col="red")
}  

  SIZE = 1.3
  p_D_facet <- base_map_trim_proj +
    # geom_point(data= D_pred_background,
    #            aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
    geom_point(data= D_pred_combined ,
               aes(x=lon,y=lat,color=Mean,fill=Mean),alpha=0.75,size=SIZE,stroke=STROKE,shape=15) +
    geom_point(data= D_pred_combined %>% filter(Mean < z.lim[1]),
               aes(x=lon,y=lat),alpha=0.75,size=SIZE,stroke=STROKE,shape=15,
               color=viridis(1,begin=0,end=0.001),fill=viridis(1,begin=0,end=0.001)) +
    geom_point(data= D_pred_combined %>% filter(Mean > z.lim[2]),
               aes(x=lon,y=lat),alpha=0.75,size=SIZE,stroke=STROKE,shape=15,
               color=viridis(1,begin=0.999,end=1),fill=viridis(1,begin=0.999,end=1)) +
    scale_color_viridis_c(option=OPT,trans="sqrt",limits=z.lim,breaks=z.breaks,labels=z.lim.labs,name=expression("DNA Copies L"^-1)) +
    scale_fill_viridis_c(option=OPT,trans="sqrt",limits=z.lim,breaks=z.breaks,labels=z.lim.labs,name=expression("DNA Copies L"^-1)) +
        facet_wrap(~depth_cat_factor) +
    theme_bw() 

  p_D_legend_only <- cowplot::get_legend(p_D$x_0)
  ggdraw(p_D_legend_only)
    
  ###################
  ## Plot of sum to surface result.
  ###################
  SIZE = 1.3
  STROKE = 0 
  
  z.lim = c(20,8500) 
  z.breaks <- c(20,200,500,1000,2000,4000,6000,8000)
  z.lim.labs <- z.breaks
  
  lon.lab <- -126.1
  lat.breaks$lats.rounded.1.0 <- lat.breaks$lats.rounded.1.0 %>% 
                                      mutate(mid.lat = (lat+lat.max)/2,
                                             lon.lab = lon.lab,
                                             mid.lat = ifelse(mid.lat==max(mid.lat),48.1,mid.lat),
                                             lon.lab = ifelse(mid.lat==max(mid.lat),-126.4,lon.lab))
  
  lat.breaks$lats.equal <- lat.breaks$lats.equal %>% 
                            mutate(mid.lat = (lat+lat.max)/2,
                                   lon.lab = lon.lab,
                                   mid.lat = ifelse(mid.lat==max(mid.lat),48.1,mid.lat),
                                   lon.lab = ifelse(mid.lat==max(mid.lat),-126.4,lon.lab))
  
  lat.breaks$lats.rounded.0.5 <- lat.breaks$lats.rounded.0.5 %>% 
                                  mutate(mid.lat = (lat+lat.max)/2,
                                  lon.lab = lon.lab,
                                  mid.lat = ifelse(mid.lat==max(mid.lat),48.1,mid.lat),
                                  lon.lab = ifelse(mid.lat==max(mid.lat),-126.4,lon.lab))
  
  p_log_D_final_proj <- base_map_trim_proj +
    # geom_point(data= D_pred_background,
    #            aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
    geom_point(data= D_final_projected ,
               aes(x=lon,y=lat,color=Mean,fill=Mean),alpha=1,size=SIZE,stroke=STROKE,shape=15) +
    # geom_point(data= D_final_projected %>% filter(Mean>z.lim[2]),
    #            aes(x=lon,y=lat),alpha=1,color=viridis(1,start=0.999, end=1),size=SIZE,stroke=STROKE,shape=15) +
    scale_color_viridis_c(option=OPT,trans="sqrt",  
                          limits=z.lim,breaks=z.breaks,labels=z.lim.labs,
                          name=expression("DNA Index")) +
    scale_fill_viridis_c(option=OPT,trans="sqrt",  
                          limits=z.lim,breaks=z.breaks,labels=z.lim.labs,
                          name=expression("DNA Index")) +
    theme_bw() 
  
  p_log_D_final_proj
  
  
  p_DNA_lat_1.0 <- p_log_D_final_proj +
          geom_segment(data=lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                       linetype="dashed")+
          geom_text(data=lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2) +
          theme( legend.position = c(.95, .5),
                 legend.justification = c("right", "top"),
                 legend.box.just = "right")
  
  p_DNA_lat_0.5 <- p_log_D_final_proj +
          geom_segment(data=lat.breaks$lats.rounded.0.5,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                       linetype="dashed")+
          geom_text(data=lat.breaks$lats.rounded.0.5,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2)

  p_DNA_lat_equal <- p_log_D_final_proj +
          geom_segment(data=lat.breaks$lats.equal,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                       linetype="dashed")+
          geom_text(data=lat.breaks$lats.equal,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2)
  
  ######################################################3
  ## SD in Space
  ######################################################3
  SIZE = 1.3
  STROKE = 0 
  
  z.lim = c(0,1.25) 
  z.breaks <- seq(z.lim[1],z.lim[2],length.out=5)

  DEPTH <- levels(STATION.DEPTH$depth_cat_factor)
  p_log_D_SD <- list()
  
  D_pred_background <- D_pred_log_combined %>% filter(depth_cat_factor == 0)
  
  for(i in 1:length(DEPTH)){
    LAB <- paste0(DEPTH[i],"m")
    if(DEPTH[i]==0){LAB <- "3m"}
    p_log_D_SD[[as.name(paste0("x_",DEPTH[i]))]] <-
      base_map_trim_proj +
      geom_point(data= D_pred_background,
                 aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
      geom_point(data= D_pred_log_combined %>% filter(depth_cat_factor == DEPTH[i]),
                 aes(x=lon,y=lat,color=SD),alpha=1,size=SIZE,stroke=STROKE) +
      geom_point(data= D_pred_log_combined %>% filter(SD > z.lim[2]),
                 aes(x=lon,y=lat),alpha=1,size=SIZE,stroke=STROKE,color=viridis(1,begin=0.999,end=1)) +
      scale_color_viridis_c(option=OPT,limits=z.lim,breaks=z.breaks,labels=round(z.breaks,2),name=expression("SD DNA Copies L"^-1)) +
      ggtitle(LAB) +
      theme_bw() 
    
    p_log_D_SD[[as.name(paste0("x_",DEPTH[i],"_dots"))]] <- 
      p_log_D_SD[[as.name(paste0("x_",DEPTH[i]))]] +
      geom_point(data=STATION.DEPTH %>% filter(depth_cat_factor == DEPTH[i]),aes(x=lon,y=lat),col="red")
  }  
  
  SIZE = 1.3
  p_log_D_SD_facet <- base_map_trim_proj +
    # geom_point(data= D_pred_background,
    #            aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
    geom_point(data= D_pred_log_combined ,
               aes(x=lon,y=lat,color=SD),alpha=1,size=SIZE,stroke=STROKE) +
    geom_point(data= D_pred_log_combined %>% filter(SD > z.lim[2]),
               aes(x=lon,y=lat),alpha=1,size=SIZE,stroke=STROKE,color=viridis(1,begin=0.999,end=1)) +
    scale_color_viridis_c(option=OPT,limits=z.lim,breaks=z.breaks,labels=round(z.breaks,2),name=expression("SD DNA Copies L"^-1)) +
    facet_wrap(~depth_cat_factor) +
    theme_bw() 

  ######################################################3
  ## CV in Space
  ######################################################3
  SIZE = 1.3
  STROKE = 0 
  
  z.lim = c(0.15,3.5) 
  z.breaks <- seq(z.lim[1],z.lim[2],length.out=6)
  
  DEPTH <- levels(STATION.DEPTH$depth_cat_factor)
  p_D_CV <- list()
  
  D_pred_combined <- D_pred_combined %>% mutate(CV=SD/Mean)
  
  D_pred_background <- D_pred_combined %>% filter(depth_cat_factor == 0)
  
  for(i in 1:length(DEPTH)){
    LAB <- paste0(DEPTH[i],"m")
    if(DEPTH[i]==0){LAB <- "3m"}
    p_D_CV[[as.name(paste0("x_",DEPTH[i]))]] <-
      base_map_trim_proj +
      geom_point(data= D_pred_background,
                 aes(x=lon,y=lat),size=SIZE,color=grey(0.6),fill=grey(0.6),alpha=1,stroke=STROKE,shape=15) +
      geom_point(data= D_pred_combined %>% filter(depth_cat_factor == DEPTH[i]),
                 aes(x=lon,y=lat,color=CV,fill=CV),alpha=1,size=SIZE,stroke=STROKE,shape=15) +
      geom_point(data= D_pred_combined %>% filter(CV > z.lim[2],depth_cat_factor == DEPTH[i]),
                  aes(x=lon,y=lat),alpha=1,size=SIZE,stroke=STROKE,shape=15,
                  color=viridis(1,begin=0.999,end=1),fill=viridis(1,begin=0.99999,end=1)) +
      scale_color_viridis_c(option=OPT,limits=z.lim,trans="sqrt",
                            breaks=z.breaks,labels=round(z.breaks,2),name=expression("CV")) +
      scale_fill_viridis_c(option=OPT,limits=z.lim,trans="sqrt",
                            breaks=z.breaks,labels=round(z.breaks,2),name=expression("CV")) +
      ggtitle(LAB) +
      theme_bw() 
    
    p_D_CV[[as.name(paste0("x_",DEPTH[i],"_dots"))]] <- 
      p_D_CV[[as.name(paste0("x_",DEPTH[i]))]] +
      geom_point(data=STATION.DEPTH %>% filter(depth_cat_factor == DEPTH[i]),aes(x=lon,y=lat),col="red")
  }  
  
  p_D_CV_legend_only <- cowplot::get_legend(p_D_CV$x_0)
  ggdraw(p_D_CV_legend_only)  
  
  SIZE = 1.3
  p_D_CV_facet <- base_map_trim_proj +
    # geom_point(data= D_pred_background,
    #            aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
    geom_point(data= D_pred_combined %>% 
                 mutate(depth_cat_factor=ifelse(depth_cat_factor==0,3,depth_cat_factor)),
               aes(x=lon,y=lat,color=CV),alpha=1,size=SIZE,stroke=STROKE) +
    geom_point(data= D_pred_combined %>% 
                      mutate(depth_cat_factor=ifelse(depth_cat_factor==0,3,depth_cat_factor)) %>% filter(CV > z.lim[2]),
               aes(x=lon,y=lat),alpha=1,size=SIZE,stroke=STROKE,color=viridis(1,begin=0.999999,end=1)) +
    scale_color_viridis_c(option=OPT,limits=z.lim,breaks=z.breaks,labels=round(z.breaks,2),name=expression("CV")) +
    facet_wrap(~depth_cat_factor) +
    theme_bw() 
  
  
### Residuals on a map  
  z.lim = c(log10(20),log10(5000)) 
  z.breaks <- log10(c(20,100,250,1000,2500,5000))
  z.lim.labs <- 10^z.breaks
  
  
  p_log_resid_facet <- base_map_trim_proj +
    geom_point(data= Resid ,
              aes(x=lon,y=lat,color=Resid.log,fill=Resid.log),alpha=0.4,stroke=STROKE,size=4) +
    scale_color_gradient2(low="red",mid=grey(0.8),high="blue") + 
    scale_fill_gradient2(low="red",mid=grey(0.8),high="blue") +
    facet_wrap(~depth_cat_factor) +
    theme_bw() 
  p_log_resid_facet
  
  p_resid_facet <- base_map_trim_proj +
    geom_point(data= Resid ,
               aes(x=lon,y=lat,color=Resid,fill=Resid),alpha=0.4,stroke=STROKE,size=4) +
    scale_color_gradient2(low="red",mid=grey(0.8),high="blue") + 
    scale_fill_gradient2(low="red",mid=grey(0.8),high="blue") +
    facet_wrap(~depth_cat_factor) +
    theme_bw() 
  p_resid_facet
  
  
  p_marginal_smooth <- list()
  if(MODEL.TYPE=="lat.long.smooth"){
      p_marginal_smooth[[as.name("bottom_depth")]] <- ggplot(D_pred_bottom_depth_combined) +
        geom_line(aes(x=bottom.depth.consensus,y=Mean),alpha=0.8) +
        geom_ribbon(aes(x=bottom.depth.consensus,ymin=Val.95.,ymax=Val.5.),alpha=0.4) +
        theme_bw()
    
  }
  
  #################### 
  ## Smoothes Only
  ###################
   D_pred_smooth_combined <- smooth.projections$D_pred_smooth_combined
  
   z.lim = c(min(D_pred_smooth_combined),max(D_pred_smooth_combined))
   z.breaks <- log10(c(20,100,250,500,1000,2500))
   #z.lim.labs <- 10^z.breaks
   
   
   SIZE = 1.3
   p_log_D_smoothes_only_facet <- base_map_trim_proj +
     # geom_point(data= D_pred_background,
     #            aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
     geom_point(data= D_pred_smooth_combined ,
                aes(x=lon,y=lat,color=Mean),alpha=0.75,size=SIZE,stroke=STROKE) +
     geom_point(data= D_pred_smooth_combined %>% filter(Mean < z.lim[1]),
                aes(x=lon,y=lat),alpha=0.75,size=SIZE,stroke=STROKE,color=viridis(1,begin=0,end=0.001)) +
     geom_point(data= D_pred_smooth_combined %>% filter(Mean > z.lim[2]),
                aes(x=lon,y=lat),alpha=0.75,size=SIZE,stroke=STROKE,color=viridis(1,begin=0.999,end=1)) +
     scale_color_viridis_c(option=OPT,
                           # limits=z.lim,
                           # breaks=z.breaks,
                           # labels=z.lim.labs,
                           name=expression("log"[10]*"DNA Copies L"^-1)) +
     facet_wrap(~depth_cat_factor) +
     theme_bw() 


######################################################3
######################################################3
##### --------  Make Transect Cross-section plots.
######################################################3
######################################################3

## OK.  I'm hoping this is a way of plotting he overlap between the eDNA 
## and the acoustic data.

# coast.ref.point <- dat.acoustic %>% filter(near.coast==1) %>% dplyr::select(transect, lat, lon,near.coast)
# coast.ref.point$transect <- as.character(coast.ref.point$transect)
# 
# temp.all <- NULL
# TRANS <-  sort(unique(STATION.DEPTH$transect))
# for(i in 1:length(TRANS)){
#   temp <- STATION.DEPTH %>% filter(transect == TRANS[i])
#   ref  <- coast.ref.point %>%filter(transect == TRANS[i]) %>% dplyr::select(lon,lat)
#   Disty <- distGeo(p1=as.matrix(data.frame(lon=temp$lon,lat=temp$lat)),p2=ref)
#   Disty <- Disty / 1000
#   temp <- temp %>% as.data.frame() %>% mutate(dist.km = Disty) %>% arrange(dist.km) %>% mutate(id.numb = 1:nrow(temp))
#   temp.all <- bind_rows(temp.all,temp)
# }
# 
# dat.SP      <- left_join(STATION.DEPTH,temp.all)
# dat.SP.raw  <- left_join(SAMPLES,
#                           temp.all %>% dplyr::select(station.depth,dist.km))
#                            
# 
# dat.acoustic  <- dat.acoustic %>% mutate(mean.depth = ifelse(mean.depth.m==0,NA,mean.depth.m) )
# dat.acoustic$transect <- as.character(dat.acoustic$transect)
# 
# dat.merge <- full_join(dat.SP %>% dplyr::select(transect,lat,lon,dist.km,depth,Mean) ,
#                        dat.acoustic %>% dplyr::select(transect,lat,lon,dist.km,mean.depth,biomass_mt,layer.thickness))
# dat.merge.raw <- full_join(dat.SP.raw %>% dplyr::select(transect,lat,lon,dist.km,depth,Mean) ,
#                        dat.acoustic %>% dplyr::select(transect,lat,lon,dist.km,mean.depth,biomass_mt,layer.thickness))
# 
# SP.transect.plots <- list()
# SP.transect.plots.raw <- list()
# 
# #################
# ## OK. Make some plots
# #################
# if(SPECIES=="hake"){ 
#   
#   ALP.red =0.5
#   ALP.blue = 0.3
#   ALP.blue.line = 0.8
#   
#   lower.lim.copies = 20
#   BREAK.copies = c(20,100,500,1000,2500,5000,7500,10000,20000,40000)
#   BREAK.copies = c(20,100,500,1000,2500,5000)
#   BREAK.copies.raw = c(20,100,500,1000,2500,5000,7500,10000,20000,40000)
#   LAB.copies= BREAK.copies
#   LAB.copies.raw= BREAK.copies.raw
#   
#   for(i in 1:length(TRANS)){
#     nom <- as.name(paste0("Transect.",TRANS[i]))
#     SP.transect.plots[[nom]] <- 
#       ggplot(dat.merge%>% filter( transect== TRANS[i])) +
#       geom_point(aes(x= dist.km, y=mean.depth,size=biomass_mt,shape="Biomass",color="Biomass"),alpha=ALP.blue.line) +  
#       geom_point(aes(x= dist.km, y=depth,size=Mean,shape="Copies / L",color="Copies / L"),alpha=ALP.red) +
#       geom_point(data=dat.merge %>% filter( transect== TRANS[i],Mean<lower.lim.copies),
#                  aes(x= dist.km, y=depth),shape="x") +
#       scale_size_continuous(name="Copies / L\nand\nBiomass (mt)",labels=LAB.copies,breaks=BREAK.copies,range=c(0.01,10),limits=c(lower.lim.copies,max(BREAK.copies))) +  
#       # plot the bathymetry
#       geom_polygon(data=bathy.transects %>% filter(transect==TRANS[i]),
#                    aes(x=dist.km,y=-depth),fill="tan",alpha=ALP.red) +
#       geom_ribbon(data=dat.acoustic %>% filter( transect== TRANS[i]),
#                   aes(x=dist.km,
#                       ymin=mean.depth - 0.5*layer.thickness,
#                       ymax=mean.depth + 0.5*layer.thickness),fill="blue",alpha=ALP.blue)+
#       scale_shape_manual("Type",values=c(21,16),labels=c("Biomass","Copies / L")) +
#       scale_color_manual("Type",values=c("blue","red"),labels=c("Biomass","Copies / L")) +
#       xlab("Distance from start of transect (km)") +
#       ylab("Depth(m)") +
#       ggtitle(paste0("Transect ",TRANS[i],"; Latitude ", 
#                      dat.bathy %>% filter(transect==TRANS[i]) %>% dplyr::select(lat_Mean) %>% round(.,1) %>% .$lat_Mean)) +
#       scale_y_reverse() +
#       coord_cartesian(ylim=c(500,0))+
#       scale_x_reverse() +
#       theme_bw()
#   }
#   for(i in 1:length(TRANS)){
#     nom <- as.name(paste0("Transect.",TRANS[i]))
#     SP.transect.plots.raw[[nom]] <- 
#       ggplot(dat.merge.raw %>% filter( transect== TRANS[i])) +
#       geom_point(aes(x= dist.km, y=mean.depth,size=biomass_mt,shape="Biomass",color="Biomass"),alpha=ALP.blue.line) +  
#       geom_point(aes(x= dist.km, y=depth,size=Mean,shape="Copies / L",color="Copies / L"),alpha=ALP.red) +
#       geom_point(data=dat.merge %>% filter( transect== TRANS[i],Mean<lower.lim.copies),
#                  aes(x= dist.km, y=depth),shape="x") +
#       scale_size_continuous(name="Copies / L\nand\nBiomass (mt)",labels=LAB.copies.raw,
#                             breaks=BREAK.copies.raw,range=c(0.01,10),
#                             limits=c(lower.lim.copies,max(BREAK.copies.raw))) +  
#       # plot the bathymetry
#       geom_polygon(data=bathy.transects %>% filter(transect==TRANS[i]),
#                    aes(x=dist.km,y=-depth),fill="tan",alpha=ALP.red) +
#       geom_ribbon(data=dat.acoustic %>% filter( transect== TRANS[i]),
#                   aes(x=dist.km,
#                       ymin=mean.depth - 0.5*layer.thickness,
#                       ymax=mean.depth + 0.5*layer.thickness),fill="blue",alpha=ALP.blue)+
#       scale_shape_manual("Type",values=c(21,16),labels=c("Biomass","Copies / L")) +
#       scale_color_manual("Type",values=c("blue","red"),labels=c("Biomass","Copies / L")) +
#       xlab("Distance from start of transect (km)") +
#       ylab("Depth(m)") +
#       ggtitle(paste0("Transect ",TRANS[i],"; Latitude ", 
#                      dat.bathy %>% filter(transect==TRANS[i]) %>% dplyr::select(lat_Mean) %>% round(.,1) %>% .$lat_Mean)) +
#       scale_y_reverse() +
#       coord_cartesian(ylim=c(500,0))+
#       scale_x_reverse() +
#       theme_bw()
#   }
#   
#   setwd("../Plots and figures")
#     pdf(file=paste("Hake transect profiles smoothed_",MODEL.TYPE,".pdf"),onefile = T,height=6,width=7)
#     print(SP.transect.plots)
#   dev.off()
#   
#   setwd("../Plots and figures")
#     pdf(file=paste("Hake transect profiles raw_",MODEL.TYPE,".pdf"),onefile = T,height=6,width=7)
#     print(SP.transect.plots.raw)
#   dev.off()
#   
# }else{ #### THIS IS A SECTION FOR non-hake plots
#   ALP.red =0.5
#   ALP.blue = 0.3
#   ALP.blue.line = 0.8
#   
#   lower.lim.copies = 20
#   BREAK.copies = c(20,100,500,1000,2500,5000,7500,10000,20000,40000,100000,140000)
#   LAB.copies= BREAK.copies
#   
#   for(i in 1:length(TRANS)){
#     nom <- as.name(paste0("Transect.",TRANS[i]))
#     SP.transect.plots[[nom]] <- 
#       ggplot(dat.SP%>% filter( transect== TRANS[i])) +
#       geom_point(aes(x= dist.km, y=depth,size=Mean),alpha=ALP.red,color="red",shape=16) +
#       geom_point(data=dat.SP %>% filter( transect== TRANS[i],Mean<lower.lim.copies),
#                  aes(x= dist.km, y=depth),shape="x") +
#       scale_size_continuous(name="Copies / L",labels=LAB.copies,breaks=BREAK.copies,range=c(0.01,10),limits=c(lower.lim.copies,max(BREAK.copies))) +  
#       # plot the bathymetry
#       geom_polygon(data=bathy.transects %>% filter(transect==TRANS[i]),
#                    aes(x=dist.km,y=-depth),fill="tan",alpha=ALP.red) +
#       # scale_shape_manual("Copies / L",values=c(16),labels=c("Copies / L")) +
#       # scale_color_manual("Copies / L",values=c("red"),labels=c("Copies / L")) +
#       xlab("Distance from start of transect (km)") +
#       ylab("Depth(m)") +
#       ggtitle(paste0("Transect ",TRANS[i],"; Latitude ", 
#                      dat.bathy %>% filter(transect==TRANS[i]) %>% dplyr::select(lat_Mean) %>% round(.,1) %>% .$lat_Mean)) +
#       scale_y_reverse() +
#       coord_cartesian(ylim=c(500,0))+
#       scale_x_reverse() +
#       theme_bw()
#   }
#   
#   # setwd("../Plots and figures")
#   # pdf(file=paste(SPECIES,"transect profiles.pdf"),onefile = T,height=6,width=7)
#   # print(SP.transect.plots)
#   # dev.off()
#}

######################################################3
######################################################3
######################################################3
######################################################3
######################################################3
##### --------  Make Some Marginal Plots
######################################################3
######################################################3
######################################################3
######################################################3
######################################################3
#### NEED TO DIFFERENTIATE AMONG CTDs in very deep water and ~500m.


# Model Estimates of average abundance at each depth
depth_cat_plot <- ggplot(depth_fact_summary) +
                  geom_point(aes(x=depth_id,y=Mean)) +
                  geom_errorbar(aes(x=depth_id,ymin=Val.X5.,ymax=Val.X95.),width=0)+
                  labs(x="Sample Depth(m)",y=expression("DNA Copies L"^-1)) +
                  expand_limits(y=0) +
                  theme_bw()
depth_cat_plot


water_depth_marg <-  STATION.DEPTH %>% group_by(water.depth.cat.2,depth_cat) %>% summarise(grand.mean = mean(Mean),
                                                                      grand.median = median(Mean),
                                                                      x.05 = quantile(Mean,probs=0.05),
                                                                      x.25 = quantile(Mean,probs=0.25),
                                                                      x.75 = quantile(Mean,probs=0.75),
                                                                      x.95 = quantile(Mean,probs=0.95)
                                                                      ) %>% rename(water_depth = water.depth.cat.2)


water_depth_marg <- water_depth_marg %>% 
                          mutate(depth_cat_jitt = case_when(water_depth=="50" ~ depth_cat-10,
                                                            water_depth=="100" ~ depth_cat-6,
                                                            water_depth=="150" ~ depth_cat-2,
                                                            water_depth=="300" ~ depth_cat+2,
                                                            water_depth=="500" ~ depth_cat+6,
                                                            water_depth=="1000+" ~ depth_cat+10
                                                                           ))

water_depth_marg <- water_depth_marg %>% 
  mutate(depth_range = case_when(water_depth=="50" ~ "< 75",
                                    water_depth=="100" ~ "75 - 125" ,
                                    water_depth=="150" ~ "125 - 210",
                                    water_depth=="300" ~ "210 - 400",
                                    water_depth=="500" ~ "400 - 750",
                                    water_depth=="1000+" ~ "> 750"
  ))

water_depth_marg$depth_range <- factor(water_depth_marg$depth_range,
                                       levels = c("< 75","75 - 125" ,"125 - 210","210 - 400","400 - 750","> 750"))


YLIM=c(0,max(water_depth_marg$x.95))

marginal_est_D_by_depth <- ggplot(water_depth_marg) +
    geom_point(aes(y=grand.mean,x=depth_cat_jitt,color=depth_range),position="jitter") +
    #geom_errorbar(aes(ymin=x.05,ymax=x.95,x=depth_cat_jitt,color=water_depth),width=0,alpha=0.5)+
    geom_errorbar(aes(ymin=x.25,ymax=x.75,x=depth_cat_jitt,color=depth_range),width=0,size=1.2,alpha=0.5)+
    geom_errorbar(aes(ymin=x.05,ymax=x.95,x=depth_cat_jitt,color=depth_range),width=0,size=0.4,alpha=0.5)+
    geom_line(aes(y=grand.mean,x=depth_cat_jitt,color=depth_range),alpha=0.5) +
    scale_x_reverse("Water Depth (m)") +
    scale_y_continuous(expression("Hake DNA (copies L"^-1*")"),limits = YLIM) +
    scale_color_npg(name="Bottom\nDepth (m)") +
    #scale_y_continuous(trans="log2") +
    coord_flip() +
    theme_bw()
# Key for water depth categories.
# water.depth.cat %in% c("w_2000_plus","w_1200_2000","w_750_1200") ~ "1000+",
# #water.depth.cat %in% c("w_750_1250") ~ "1000",
# water.depth.cat %in% c("w_400_750") ~ "500",
# water.depth.cat %in% c("w_210_400") ~ "300",
# water.depth.cat %in% c("w_125_210") ~ "150",
# water.depth.cat %in% c("w_75_125") ~ "100",
# water.depth.cat %in% c("w_0_75") ~ "50"))








# 
#   
#   coord_cartesian(ylim=c(500,0))+
#   scale_x_reverse() +
#   
 
# dat.SP.sum.by.depth.cat <- dat.SP %>% group_by(water.depth.cat.2,depth_cat) %>%
#   summarise(N=length(Mean),MEAN=mean(Mean),SD=sd(Mean),
#             MEDIAN=median(Mean),
#             q.25=quantile(Mean,probs=0.25),
#             q.75=quantile(Mean,probs=0.75))
# 
# MIN.bin <- 10
# dat.SP <- dat.SP %>% mutate(Mean.binned = ifelse(Mean<=MIN.bin,MIN.bin,Mean))
# LIM <- c(MIN.bin,max(dat.SP$Mean.binned))
# 
# depth.p2 <- ggplot(dat.SP) +
#   geom_point(aes(x=depth,y=Mean.binned,group=station,color=lat),alpha=0.5) +
#   geom_line(aes(x=depth,y=Mean.binned,group=station,color=lat),alpha=0.3) +
#   geom_point(data= dat.SP.sum.by.depth.cat,
#              aes(x=depth_cat,y=MEDIAN,group=water.depth.cat.2),color="red",shape=16,size=3) +
#   geom_point(data= dat.SP.sum.by.depth.cat,
#              aes(x=depth_cat,y=MEAN,group=water.depth.cat.2),color="red",shape=22,size=3) +
#   geom_line(data= dat.SP.sum.by.depth.cat,
#             aes(x=depth_cat,y=MEDIAN,group=water.depth.cat.2),color="red") +
#   geom_errorbar(data= dat.SP.sum.by.depth.cat,
#                 aes(x=depth_cat,ymin=q.25,ymax=q.75,group=water.depth.cat.2),color="red",width=0) +
#   geom_hline(yintercept=20, linetype="dashed") +
#   scale_color_viridis("Latitude",begin=0,end=0.8) +
#   scale_x_reverse() +
#   scale_y_continuous(limits=LIM,trans="log2") +
#   scale_shape(solid=FALSE) +
#   xlab("Depth") +
#   ylab(expression("Copies L"^-1)) +
#   coord_flip() +
#   facet_wrap(~water.depth.cat.2,ncol=2)+
#   theme_bw()
# print(depth.p2)
# 
# 
# # dat.SP.sum.by.depth.cat <- dat.SP.sum.by.depth.cat %>%
# #                               mutate(max.depth.cat.fact=as.factor(max.depth.cat))
# # LEV <- unique(dat.SP.sum.by.depth.cat$max.depth.cat.fact)
# # 
# # dat.SP.sum.by.depth.cat$max.depth.cat.fact <- factor(dat.SP.sum.by.depth.cat$max.depth.cat.fact,
# #                                                        levels=LEV)
# 
# LIM <- c(0,max(dat.SP.sum.by.depth.cat$q.75)*1.05)
# 
# depth.p3 <- ggplot(dat.SP.sum.by.depth.cat) +
#   geom_point(aes(x=-depth_cat,y=MEDIAN,color=water.depth.cat.2),
#              shape=16,size=3,alpha=0.8) +
#   # geom_point(aes(x=-depth_cat,y=MEAN,color=max.depth.cat.fact),
#   #           shape=16,size=3) +
#   geom_line(aes(x=-depth_cat,y=MEDIAN,color=water.depth.cat.2),alpha=0.6) + 
#   geom_ribbon(aes(x=-depth_cat,ymin=q.25,ymax=q.75,fill=water.depth.cat.2),
#               width=10,alpha=0.1) +
#   scale_color_viridis_d("Water Depth \nCategory(m)",begin=0,end=0.8) +
#   scale_fill_viridis_d("Water Depth \nCategory(m)",begin=0,end=0.8) +
#   #scale_y_continuous(label=LABEL,breaks=LABEL) +
#   scale_shape(solid=FALSE) +
#   xlab("Depth") +
#   ylab("Copies") +
#   coord_flip() +
#   scale_y_continuous(limits=LIM, expand=c(0,0)) +
#   theme_bw()
# print(depth.p3)
# 
# # depth.p4 <- depth.p3 +
# #   scale_y_continuous(trans="log2",label=LABEL,breaks=LABEL,expand=c(0,0.05)) 
# # 
# # print(depth.p4)
# 
# #######################
# 
# # Make a north-south transect by water depth
# 
# lat.by.depth.dat <- dat.summary %>% group_by(transect,depth_cat) %>% 
#   summarise(lat.mean = mean(lat), grand.mean = mean(Mean),grand.sd = sd(Mean))
# lat.by.depth.dat$depth_cat <-as.factor(lat.by.depth.dat$depth_cat )
# dat.summary2 <- dat.summary
# dat.summary2$depth_cat <-as.factor(dat.summary2$depth_cat )
# 
# p.lat.by.d <- ggplot() +
#   geom_point(data=lat.by.depth.dat,aes(y=grand.mean,x=lat.mean,color=depth_cat),alpha=0.5) +
#   geom_line(data=lat.by.depth.dat,aes(y=grand.mean,x=lat.mean,color=depth_cat),alpha=0.5) +
#   geom_point(data=dat.summary2,aes(y=Mean,x=lat,color=depth_cat),alpha=0.2) +
#   scale_y_log10() +
#   geom_hline(yintercept = 20,linetype="dashed",color="purple",alpha=0.5)+
#   #geom_hline(yintercept = 100,linetype="dashed",color="blue",alpha=0.5)+
#   theme_bw() +
#   scale_color_discrete("Depth") +
#   xlab("Latitude") +
#   ylab(expression("DNA copies L"^"-1"))
# 
# p.lat.by.d.facet <- ggplot() +
#   geom_point(data=lat.by.depth.dat,aes(y=grand.mean,x=lat.mean),alpha=0.5) +
#   geom_line(data=lat.by.depth.dat,aes(y=grand.mean,x=lat.mean),alpha=0.5) +
#   geom_point(data=dat.summary2,aes(y=Mean,x=lat),alpha=0.2) +
#   scale_y_log10() +
#   geom_hline(yintercept = 20,linetype="dashed",color="purple",alpha=0.5)+
#   #geom_hline(yintercept = 100,linetype="dashed",color="blue",alpha=0.5)+
#   theme_bw() +
#   #scale_color_discrete("Depth")
#   xlab("Latitude") +
#   ylab(expression("DNA copies L"^"-1")) + facet_wrap(~depth_cat)


#save output data to file.
# These can be used to merge with acoustic datasets.


# Combine the necessary data.frames into a list for use later.
Output.summary.qpcr <- list(
  # Data file that is helpful
  dat.obs.bin = dat.obs.bin,
  # Output (Standards)
  stand.plot = stand.plot,
  stand.plot.pres = stand.plot.pres,

  # Posterior summaries that are helpful.
  stanMod_summary_parts = Output.qpcr$stanMod_summary_parts,
  
  # Output (Field Summaries)
  tau_long_summary = tau_long_summary, 
  station_depth_out=station_depth_out,
  station_depth_out_liter=station_depth_out_liter,
  D_pred_combined = D_pred_combined,
  D_delta_out = D_delta_out,
  D_delta_out_liter = D_delta_out_liter,
  
  D_DNA_uncond_total = D_DNA_uncond_total,
  D_DNA_cum_sum = D_DNA_cum_sum,
  D_final_projected = D_final_projected,
  D_final_lat_1.0 = D_final_lat_1.0,
  D_final_lat_0.5 = D_final_lat_0.5,
  D_final_lat_equal = D_final_lat_equal,
  Resid = Resid,
  
  # summarized by MCMC sample for confidence intervals on 
  
  D_grid.cell_resample = D_grid.cell_resample, 
  D_1.0_resample = D_1.0_resample,
  D_0.5_resample = D_0.5_resample,
  
  D_smooth_summary = D_smooth_summary,
  
  #Marginal of D by depth 
  marginal_est_D_by_depth = marginal_est_D_by_depth,
  
  # Marginal smooth
  p_marginal_smooth = p_marginal_smooth,
  
  #Plots with maps
  p_DNA_lat_1.0 =p_DNA_lat_1.0,
  p_DNA_lat_0.5 = p_DNA_lat_0.5,
  p_DNA_lat_equal =p_DNA_lat_equal,
  p_DNA_lat_base =p_log_D_final_proj,
  
  p_D = p_D,
  p_D_legend_only = p_D_legend_only,
  
  p_D_CV = p_D_CV,
  p_D_CV_facet = p_D_CV_facet,
  p_D_CV_legend_only = p_D_CV_legend_only,
  
  # Variability plots 
  tau.sample.by.depth = tau.sample.by.depth,
  
  # Contamination and wash plots
  control.by.time = control.by.time ,
  control.hist = control.hist,
  control.hist1 = control.hist1,
  
  depth.hist.facet = depth.hist.facet,
  
  inhibit.plot = inhibit.plot ,
  inhibit.plot.by.depth  = inhibit.plot.by.depth ,
  inhibit.plot.by.depth.few = inhibit.plot.by.depth.few, 
  
  wash.plot = wash.plot ,
  wash.plot.by.depth = wash.plot.by.depth, 
  wash.param.hist = wash.param.hist,
  
  # Resid plots
  p_resid=p_resid,
  
  # Projections or projection helpers.
  N.POST = N.POST, # number of posterior samples used.
  dat_raster_fin = dat_raster_fin, 
  smooth.projections = smooth.projections,
  
  # Output(negative Controls, Contamination)
  field_neg_out=field_neg_out,
  field_neg_out_liter=field_neg_out_liter,
  sample_contam_total_out = sample_contam_total_out,
  sample_contam_total_out_liter = sample_contam_total_out_liter,
  STATION.DEPTH = STATION.DEPTH,
  delta_out=delta_out,
  mu_contam_out = mu_contam_out,
  sigma_contam_out =sigma_contam_out,
  lat.breaks = lat.breaks,
  
  sigma.pcr.param = sigma.pcr.param,
  sigma.stand.param =sigma.stand.param,
  base_map_trim_proj = base_map_trim_proj
)


setwd(results.dir)
#Output.summary.qpcr
save(Output.summary.qpcr,file=paste0("./_Summarized_Output/Qpcr_summaries_",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,"_NoSURF=",NO.SURFACE,".RData"))

#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
setwd(plot.dir)
###### Plots THAT MATTER
pdf(file=paste("Hake transect DNA stand.plot",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=6,width=7)
  print(stand.plot)
dev.off()
pdf(file=paste("Hake transect DNA stand.plot.pres",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=6,width=7)
  print(stand.plot.pres)
dev.off()

# Marginal Plots
pdf(file=paste("Hake transect DNA Marginals",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=4,width=6)
  print(wash.param.hist)
  print(p_marginal_smooth)
  print(depth_cat_plot)
  # Variability parameters.
  print(sd_param_plots)
dev.off()

# Dot plots for Niskins and Smooth Derived 
pdf(file=paste("Hake transect DNA D_delta inhibit",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=6,width=7)
  print(p2)
dev.off()
pdf(file=paste("Hake transect DNA D_delta wash",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=6,width=7)
  print(p3)
dev.off()
pdf(file=paste("Hake transect DNA D_delta red",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=10,width=7)
  print(p4)
dev.off()

# Spatial Smoothes
pdf(file=paste("Hake transect DNA p_D",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=6,width=5)
  print(p_D)
dev.off()
pdf(file=paste("Hake transect DNA p_D_facet",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=10,width=7)
  print(p_D_facet)
dev.off()
# pdf(file=paste("Hake transect DNA p_log_D_facet2",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=6,width=11)
#   print(p_log_D_facet2)
# dev.off()

pdf(file=paste(SPECIES," DNA_smoothes_",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf",sep=""),onefile = T,height=10,width=8)
  print(p_D_facet)
  print(p_D)
dev.off()

##### Summed to surface plots.
quartz(file=paste("Hake DNA final_proj 1_deg_lat_breaks",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),height=6,width=4,dpi=600,type="pdf")
  print(p_log_D_final_proj +
          geom_segment(data=lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                       linetype="dashed")+
          geom_text(data=lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2)
        )
dev.off()

quartz(file=paste("Hake DNA final_proj 0.5_deg_lat_breaks",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),height=6,width=4,dpi=600,type="pdf")
    print(p_log_D_final_proj +
        geom_segment(data=lat.breaks$lats.rounded.0.5,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                     linetype="dashed")+
        geom_text(data=lat.breaks$lats.rounded.0.5,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2)
    )
dev.off()

quartz(file=paste("Hake DNA final_proj equal_lat_breaks",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),height=6,width=4,dpi=600,type="pdf")
  print(p_log_D_final_proj +
        geom_segment(data=lat.breaks$lats.equal,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                     linetype="dashed")+
        geom_text(data=lat.breaks$lats.equal,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2)
)
dev.off()

##################################################
##################################################
##################################################


p_D
p_D_facet
p_log_D_SD
p_log_D_SD_facet

p_D_CV
p_D_CV_facet

p_log_D_smoothes_only_facet

# Latitudinal slices.
SP.transect.plots

# Factors for each depth
depth_cat_plot

# Variance parameters
print(sd_param_plots)


# Marginal Smoothes
p_marginal_smooth

# Residuals
p_log_resid_facet

# Contamination 
control.by.time
control.hist
control.hist1

depth.hist.facet

inhibit.plot 
inhibit.plot.by.depth 
inhibit.plot.by.depth.few 

wash.plot 
wash.plot.by.depth 
  
wash.param.hist