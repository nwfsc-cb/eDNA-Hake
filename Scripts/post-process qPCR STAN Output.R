###
library(tidyverse)
library(knitr)
library(reshape2)
library(viridis)
library(ggmap)
library(maps)
library(mapdata)
library(ggplot2)


script.dir <- "/Users/ole.shelton/Github/eDNA-Hake/Scripts"
# load and run in the acoustic data.
setwd(script.dir)
source("process acoustic data.R")
# get bathymetry data.
source("pull NOAA bathy for acoustic data.R")

# Read in the output from the Stan model
base.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Stan Model Fits"
setwd(base.dir)

# CHANGE THIS FOR switching between species.
#SPECIES <- "eulachon"

load(paste("qPCR 2019",SPECIES, "Fitted.RData"))

### Get some of the output on sample and station IDs from file.
dat.samp <- Output.qpcr$dat.obs.bin
dat.control.field.neg <- Output.qpcr$dat.control.field.neg %>% left_join(.,Output.qpcr$field_neg_out_liter)

dat.summary <- dat.samp %>%
  group_by(date,transect,station,station.depth,depth,station_depth_idx) %>%
  dplyr::summarise(N=length(unique(sample_idx)))
dat.summary <- left_join(dat.summary,Output.qpcr$station_depth_out_liter)
dat.summary <- dat.summary %>% mutate(depth_cat=case_when(depth < 25 ~ 0,
                                                          depth ==25 ~ 25,  
                                                          depth > 25  & depth <= 50  ~ 50,
                                                          depth > 50  & depth <= 100 ~ 100,
                                                          depth > 120 & depth <= 150 ~ 150,
                                                          depth > 151 & depth <= 200 ~ 200,
                                                          depth > 250 & depth <= 350 ~ 300,
                                                          depth > 400 & depth <= 500 ~ 500))

dat.inhibit <- Output.qpcr$dat.inhibit %>% mutate(depth_cat=case_when(depth < 25 ~ 0,
                                                                      depth == 25 ~ 25,            
                                                                      depth > 25  & depth <=50  ~ 50,
                                                                      depth > 50  & depth <=100 ~ 100,
                                                                      depth > 120 & depth <=150 ~ 150,
                                                                      depth > 151 & depth <=200 ~ 200,
                                                                      depth > 250 & depth <=350 ~ 300,
                                                                      depth > 400 & depth <=500 ~ 500))

##### Calculate some info about the field negative controls
field.neg.summ <- Output.qpcr$field_neg_out_liter 

N.station.run       <- length(unique(c(dat.inhibit$station,dat.samp$station)))
N.station.depth.run <- length(unique(c(dat.inhibit$station.depth,dat.samp$station.depth)))
N.sample.run        <- length(unique(c(dat.inhibit$sample,dat.samp$sample)))

###################################################################
###################################################################
###################################################################
###################################################################
###################################################################

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

dat.lat.lon.water <- dat.process %>% group_by(station) %>% summarise(lat=mean(lat),lon=mean(lon),water.depth=mean(water.depth)) %>%
                      mutate(water.depth.cat = case_when(water.depth >=2000 ~ "w_2000_plus",
                                                         water.depth < 2000 & water.depth >=1250 ~ "w_1250_2000",
                                                         water.depth < 1250 & water.depth >= 750 ~ "w_750_1250",
                                                         water.depth < 750 & water.depth >= 400 ~ "w_400_750",
                                                         water.depth < 400 & water.depth >= 250 ~ "w_250_400",
                                                         water.depth < 250  & water.depth >= 125 ~ "w_125_250",
                                                         water.depth < 125 & water.depth >= 75 ~ "w_75_125",
                                                         water.depth < 75  ~ "w_0_75"), 
                               water.depth.cat = case_when(station %in% c("26-3","27-8","49-7") ~ "w_2000_plus",
                                                         station %in% c("47-9") ~ "w_1250_2000",
                                                         TRUE ~ water.depth.cat )) %>%
                      mutate(water.depth.cat.2= water.depth.cat, 
                             water.depth.cat.2 = case_when(water.depth.cat %in% c("w_2000_plus","w_1250_2000","w_750_1250") ~ "1000+",
                                                          #water.depth.cat %in% c("w_750_1250") ~ "1000",
                                                          water.depth.cat %in% c("w_400_750") ~ "500",
                                                          water.depth.cat %in% c("w_250_400") ~ "300",
                                                          water.depth.cat %in% c("w_125_250") ~ "150",
                                                          #water.depth.cat %in% c("w_75_125") ~ "100",
                                                          water.depth.cat %in% c("w_0_75","w_75_125") ~ "50"))


#### Plotting information

# Plot for all samples
lat.lims <- c(min(dat.id$lat,na.rm=T),max(dat.id$lat,na.rm=T))
lon.lims <- c(min(dat.id$lon,na.rm=T),max(dat.id$lon,na.rm=T))

lat.lims.trim <- c(37.5,max(dat.id$lat,na.rm=T))
lon.lims.trim <- c(min(dat.id$lon,na.rm=T),-122.5)


# Call Base_map.R
setwd(script.dir)
source("Base_map.R")

sample.process <- base_map + 
  scale_color_manual(values=c(grey(0.6),"blue","red"))+
  geom_point(data=dat.process,aes(x=lon,y=lat,color=Level),size=0.8) 

sample.process.facet <- sample.process +
  facet_wrap(~Level)

#print(sample.process)


###################################################################
###################################################################
###################################################################
###################################################################
###################################################################


##### Calculate some inhibition characteristics.
inhibit.summ <- dat.inhibit %>% 
  group_by(date,transect,station,station.depth,depth,depth_cat) %>%
  summarise(N_pcr=length(unique(qPCR)),N_sample=length(unique(sample)),N_rep=length(sample)) %>%
  mutate(inhibited = "Y")
good.summ <- dat.summary %>%
  dplyr::select(date,transect,station,station.depth,depth,depth_cat,N_sample=N) %>%
  mutate(inhibited = "N")

g_i_dat <- bind_rows(inhibit.summ %>% dplyr::select(-N_pcr,-N_rep),good.summ) %>%
  group_by(station,depth_cat,station.depth,inhibited) %>%
  summarise(N_inhibited=length(inhibited)) %>%
  pivot_wider(.,names_from=c("inhibited"),values_from=c("N_inhibited")) %>%
  mutate(N=ifelse(is.na(N),0,N),
         Y=ifelse(is.na(Y),0,Y),
         Both=ifelse(N==1 & Y==1,1,0)) %>%
  group_by(depth_cat) %>% 
  summarise(Total=length(station.depth),
            Not_inhibited=Total-sum(Y),
            Partially_inhibited = sum(Both),
            Completely_inhibited = Total-sum(N)) %>%
  filter(is.na(depth_cat)==F)

g_i_dat <- bind_rows(g_i_dat %>% mutate(depth_cat=as.character(depth_cat)),
                     colSums(g_i_dat)%>%t()%>%as.data.frame()%>%mutate(depth_cat="TOTAL")) %>%
  rename("Depth Category (m)"=depth_cat,)

g_i_2 <- bind_rows(inhibit.summ %>% dplyr::select(-N_pcr,-N_rep),good.summ) %>%
  group_by(station,depth_cat,station.depth,inhibited) %>%
  dplyr::summarise(N_inhibited=length(inhibited)) %>%
  pivot_wider(.,names_from=c("inhibited"),values_from=c("N_inhibited")) %>%
  mutate(N=ifelse(is.na(N),0,N),
         Y=ifelse(is.na(Y),0,Y),
         Both=ifelse(N==1 & Y==1,1,0)) %>%
  left_join(.,dat.lat.lon.water)

#####
g_i_3 <- g_i_2 %>% mutate(inhib_cat=case_when(N==1 & Y==0 ~ "None",
                                              N==0 & Y==1 ~ "All",
                                              N==1 & Y==1 ~ "Partial"))
g_i_3$inhib_cat <- factor(g_i_3$inhib_cat,levels=c("All","Partial","None"))

### Make a spatial plot of inhibition

inhibit.plot <- base_map_trim +
  geom_point(data=g_i_3,aes(x=lon,y=lat,color=inhib_cat)) +
  scale_color_viridis_d(begin=0,end=0.8)+
  facet_wrap(~inhib_cat)

inhibit.plot.by.depth <- base_map_trim +
  geom_point(data=g_i_3 %>% filter(is.na(depth_cat)==F),aes(x=lon,y=lat,color=inhib_cat),alpha=0.5) +
  facet_wrap(~depth_cat,nrow = 2) +
  scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) 


inhibit.plot.by.depth.few <- base_map_trim +
  geom_point(data=g_i_3 %>% filter(depth_cat %in% c(0,50,100)) ,aes(x=lon,y=lat,color=inhib_cat),alpha=0.5) +
  facet_wrap(~depth_cat,nrow = 1) +
  scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) 


###################################################################
###################################################################
###----- Plots of standards and data
###################################################################
###################################################################

OFFSET <- Output.qpcr$OFFSET
pars <- Output.qpcr$samp
PCR <- Output.qpcr$PCR
dat.stand.pos <- Output.qpcr$dat.stand.pos
dat.stand.bin <- Output.qpcr$dat.stand.bin

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


lab.temp<- data.frame(variable=paste0("Y.",PCR$plate_idx), qPCR=PCR$qPCR)
STAND.REG <-left_join(STAND.REG,lab.temp)

stand.plot <- ggplot(dat.stand.pos) +
  geom_point(aes(x=log10_copies - OFFSET ,y=Ct,color=qPCR),alpha=0.5) +
  scale_shape_discrete(name="qPCR Plate") +
  theme_bw()
stand.plot <- stand.plot +
  geom_line(data=STAND.REG,aes(x=X,y=Y,color=qPCR),alpha=0.5) +
  scale_color_discrete(name="qPCR Plate") +
  ylab("PCR cycle")  +
  xlab("log10 copies DNA") 
#scale_x_continuous(labels = paste0("1e",LABS),limits = c(-2,4))
stand.plot



stand.plot.facet <- stand.plot +facet_wrap(~qPCR,ncol=2)

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

#######################################################
#######################################################
####### ---- Calculate Control quantities.
#######################################################
#######################################################

# Calculations for Field Negative Controls.
SAMPLES.CONTROL <- Output.qpcr$SAMPLES.CONTROL
dat.id.control <-  Output.qpcr$dat.id.control
field.neg.summ <- left_join(SAMPLES.CONTROL,field.neg.summ) %>% left_join(.,dat.id.control)

BREAKS <- seq(0,1000,by=10)
XLIM <- c(0,max(field.neg.summ$Mean)*1.1)
control.hist1 <- ggplot(field.neg.summ) +
  geom_histogram(aes(Mean,fill=field.negative.type),breaks=BREAKS) +
  xlab("Estimated copies per L") +
  scale_fill_discrete("Type")+
  ylab("Number of negative controls") +
  xlim(XLIM) +
  theme_bw()

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

#control.by.time.trim <- control.by.time + ylim(c(0,1000))
detect.lev = 20
N_control_positive <- field.neg.summ %>% filter(Mean > detect.lev ) %>% nrow(.)
N_control_run <- nrow(field.neg.summ)
N_control_total <- dat.id.control %>% filter(control=="field") %>% nrow(.)

#######################################################
#######################################################
####### ---- Compare distribution of field samples against negative controls.
#######################################################
#######################################################

dat.sum.v.control <- bind_rows(dat.summary %>% dplyr::select(date,
                                                             Mean.log,Log.Val.X5.,Log.Val.X95.,
                                                             Mean,Val.X5.,Val.X95.) %>% mutate(cat="field") %>% as.data.frame(),
                               field.neg.summ %>% dplyr::select(date,
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
  scale_fill_viridis_d("Category",begin=0,end=0.5,) +
  scale_x_log10() +
  theme_bw()
compare.hist2

#######################################################
#######################################################
####### ------- Make Maps
#######################################################
#######################################################


dat.SP <- left_join(dat.summary,dat.lat.lon.water) %>% filter(is.na(depth_cat)==F) %>% mutate(neg.depth = -depth)
LEV <- c("50","150","300","500","1000+") #,"1500+")
dat.SP$water.depth.cat.2 <- factor(dat.SP$water.depth.cat.2,levels=LEV) 
# drop sample sampled at 200m 
dat.SP <- dat.SP %>% filter(!depth==200)

BREAKS = c(-100,-2,0,1,2,3,4,5)
SP.p1 <- base_map_trim + 
  geom_point(data=dat.SP,aes(x=lon,y=lat,size=Mean.log),shape=21) +
  scale_size("Log Copies per L",breaks=BREAKS)+
  scale_shape("Log Copies",solid=FALSE)+
  facet_wrap(~depth_cat)
SP.p1


  summary(dat.SP$Mean)
  lower.lim <- 20
  LAB<-c(lower.lim,100,1000,10000,20000)
SP.p2 <- base_map_trim + 
    geom_point(data=dat.SP,aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
    scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA))+
    geom_point(data=dat.SP%>% filter(Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
    scale_shape("Copies / L",solid=FALSE)+
    facet_wrap(~depth_cat,nrow=2)
SP.p2
  
SP.p3 <- base_map_trim + 
  geom_point(data=dat.SP %>% filter(!depth_cat==25),aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
  scale_size("Copies / L",labels=LAB,breaks=LAB,range=c(0.01,10),limits=c(lower.lim,NA))+
  geom_point(data=dat.SP%>% filter(Mean<=lower.lim,!depth_cat==25),aes(x=lon,y=lat),shape="x",size=1,color="black") +
  scale_shape("Copies / L",solid=FALSE)+
  facet_wrap(~depth_cat,nrow=2)
SP.p3

######################## 
######################################################3
######################################################3
##### --------  Make Transect Cross-section plots.
######################################################3
######################################################3

## OK.  I'm hoping this is a way of plotting he overlap between the eDNA 
## and the acoustic data.

coast.ref.point <- dat.acoustic %>% filter(near.coast==1) %>% dplyr::select(transect, lat, lon,near.coast)
coast.ref.point$transect <- as.character(coast.ref.point$transect)


temp.all <- NULL
TRANS <-  sort(unique(dat.SP$transect))
for(i in 1:length(TRANS)){
  
  temp <- dat.SP %>% filter(transect == TRANS[i])
  ref  <- coast.ref.point %>%filter(transect == TRANS[i]) %>% select(lon,lat)
  D <- distGeo(p1=as.matrix(data.frame(lon=temp$lon,lat=temp$lat)),p2=ref)
  D <- D / 1000
  temp <- temp %>% as.data.frame() %>% mutate(dist.km = D) %>% arrange(dist.km) %>% mutate(id.numb = 1:nrow(temp))
  temp.all <- bind_rows(temp.all,temp)
}

dat.SP      <- left_join(dat.SP,temp.all)

dat.acoustic  <- dat.acoustic %>% mutate(mean.depth = ifelse(mean.depth==0,NA,mean.depth) )

dat.acoustic$transect <- as.character(dat.acoustic$transect)
dat.merge <- full_join(dat.SP %>% dplyr::select(transect,lat,lon,dist.km,neg.depth,Mean) ,
                       dat.acoustic %>% dplyr::select(transect,lat,lon,dist.km,mean.depth,biomass_mt,layer.thickness))

SP.transect.plots <- list()

#################
## OK. Make some plots
#################
if(SPECIES=="hake"){ 

ALP.red =0.5
ALP.blue = 0.3
ALP.blue.line = 0.8

lower.lim.copies = 20
BREAK.copies = c(20,100,500,1000,2500,5000,7500,10000,20000,30000)
LAB.copies= BREAK.copies

for(i in 1:length(TRANS)){
  nom <- as.name(paste0("Transect.",TRANS[i]))
  SP.transect.plots[[nom]] <- 
    ggplot(dat.merge%>% filter( transect== TRANS[i])) +
    geom_point(aes(x= dist.km, y=mean.depth,size=biomass_mt,shape="Biomass",color="Biomass"),alpha=ALP.blue.line) +  
    geom_point(aes(x= dist.km, y=depth,size=Mean,shape="Copies / L",color="Copies / L"),alpha=ALP.red) +
    geom_point(data=dat.merge %>% filter( transect== TRANS[i],Mean<lower.lim.copies),
                  aes(x= dist.km, y=depth),shape="x") +
    scale_size_continuous(name="Copies / L\nand\nBiomass (mt)",labels=LAB.copies,breaks=BREAK.copies,range=c(0.01,10),limits=c(lower.lim.copies,max(BREAK.copies))) +  
    # plot the bathymetry
    geom_polygon(data=bathy.transects %>% filter(transect==TRANS[i]),
                  aes(x=dist.km,y=-depth),fill="tan",alpha=ALP.red) +
    geom_ribbon(data=dat.acoustic %>% filter( transect== TRANS[i]),
               aes(x=dist.km,
                   ymin=mean.depth - 0.5*layer.thickness,
                   ymax=mean.depth + 0.5*layer.thickness),fill="blue",alpha=ALP.blue)+
    scale_shape_manual("Type",values=c(21,16),labels=c("Biomass","Copies / L")) +
    scale_color_manual("Type",values=c("blue","red"),labels=c("Biomass","Copies / L")) +
    xlab("Distance from start of transect (km)") +
    ylab("Depth(m)") +
    ggtitle(paste0("Transect ",TRANS[i],"; Latitude ", 
                  dat.bathy %>% filter(transect==TRANS[i]) %>% dplyr::select(lat_Mean) %>% round(.,1) %>% .$lat_Mean)) +
    coord_cartesian(ylim=c(0,500))+
    scale_x_reverse() +
    scale_y_reverse() +
    theme_bw()
}


setwd("../Plots and figures")
pdf(file="Hake transect profiles.pdf",onefile = T,height=6,width=7)
  print(SP.transect.plots)
dev.off()

}else{ #### THIS IS A SECTION FOR non-hake plots
  
  
  ALP.red =0.5
  ALP.blue = 0.3
  ALP.blue.line = 0.8
  
  lower.lim.copies = 20
  BREAK.copies = c(20,100,500,1000,2500,5000,7500,10000,20000,40000)
  LAB.copies= BREAK.copies
  
  for(i in 1:length(TRANS)){
    nom <- as.name(paste0("Transect.",TRANS[i]))
    SP.transect.plots[[nom]] <- 
      ggplot(dat.merge%>% filter( transect== TRANS[i])) +

      geom_point(aes(x= dist.km, y=depth,size=Mean),alpha=ALP.red,color="red",shape=16) +
      geom_point(data=dat.merge %>% filter( transect== TRANS[i],Mean<lower.lim.copies),
                 aes(x= dist.km, y=depth),shape="x") +
      scale_size_continuous(name="Copies / L",labels=LAB.copies,breaks=BREAK.copies,range=c(0.01,10),limits=c(lower.lim.copies,max(BREAK.copies))) +  
      # plot the bathymetry
      geom_polygon(data=bathy.transects %>% filter(transect==TRANS[i]),
                   aes(x=dist.km,y=-depth),fill="tan",alpha=ALP.red) +
      # scale_shape_manual("Copies / L",values=c(16),labels=c("Copies / L")) +
      # scale_color_manual("Copies / L",values=c("red"),labels=c("Copies / L")) +
      xlab("Distance from start of transect (km)") +
      ylab("Depth(m)") +
      ggtitle(paste0("Transect ",TRANS[i],"; Latitude ", 
                     dat.bathy %>% filter(transect==TRANS[i]) %>% dplyr::select(lat_Mean) %>% round(.,1) %>% .$lat_Mean)) +
      coord_cartesian(ylim=c(0,500))+
      scale_x_reverse() +
      scale_y_reverse() +
      theme_bw()
  }

  setwd("../Plots and figures")
  pdf(file=paste(SPECIES,"transect profiles.pdf"),onefile = T,height=6,width=7)
  print(SP.transect.plots)
  dev.off()
}



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

dat.SP.sum.by.depth.cat <- dat.SP %>% group_by(water.depth.cat.2,depth_cat) %>%
  summarise(N=length(Mean),MEAN=mean(Mean),SD=sd(Mean),
            MEDIAN=median(Mean),
            q.25=quantile(Mean,probs=0.25),
            q.75=quantile(Mean,probs=0.75))

MIN.bin <- 10
dat.SP <- dat.SP %>% mutate(Mean.binned = ifelse(Mean<=MIN.bin,MIN.bin,Mean))
LIM <- c(MIN.bin,max(dat.SP$Mean.binned))

depth.p2 <- ggplot(dat.SP) +
  geom_point(aes(x=depth,y=Mean.binned,group=station,color=lat),alpha=0.5) +
  geom_line(aes(x=depth,y=Mean.binned,group=station,color=lat),alpha=0.3) +
  geom_point(data= dat.SP.sum.by.depth.cat,
             aes(x=depth_cat,y=MEDIAN,group=water.depth.cat.2),color="red",shape=16,size=3) +
  geom_point(data= dat.SP.sum.by.depth.cat,
             aes(x=depth_cat,y=MEAN,group=water.depth.cat.2),color="red",shape=22,size=3) +
  geom_line(data= dat.SP.sum.by.depth.cat,
            aes(x=depth_cat,y=MEDIAN,group=water.depth.cat.2),color="red") +
  geom_errorbar(data= dat.SP.sum.by.depth.cat,
                aes(x=depth_cat,ymin=q.25,ymax=q.75,group=water.depth.cat.2),color="red",width=0) +
  geom_hline(yintercept=20, linetype="dashed") +
  scale_color_viridis("Latitude",begin=0,end=0.8) +
  scale_x_reverse() +
  scale_y_continuous(limits=LIM,trans="log2") +
  scale_shape(solid=FALSE) +
  xlab("Depth") +
  ylab(expression("Copies L"^-1)) +
  coord_flip() +
  facet_wrap(~water.depth.cat.2,ncol=2)+
  theme_bw()
print(depth.p2)









# dat.SP.sum.by.depth.cat <- dat.SP.sum.by.depth.cat %>%
#                               mutate(max.depth.cat.fact=as.factor(max.depth.cat))
# LEV <- unique(dat.SP.sum.by.depth.cat$max.depth.cat.fact)
# 
# dat.SP.sum.by.depth.cat$max.depth.cat.fact <- factor(dat.SP.sum.by.depth.cat$max.depth.cat.fact,
#                                                        levels=LEV)

LIM <- c(0,max(dat.SP.sum.by.depth.cat$q.75)*1.05)

depth.p3 <- ggplot(dat.SP.sum.by.depth.cat) +
  geom_point(aes(x=-depth_cat,y=MEDIAN,color=water.depth.cat.2),
             shape=16,size=3,alpha=0.8) +
  # geom_point(aes(x=-depth_cat,y=MEAN,color=max.depth.cat.fact),
  #           shape=16,size=3) +
  geom_line(aes(x=-depth_cat,y=MEDIAN,color=water.depth.cat.2),alpha=0.6) + 
  geom_ribbon(aes(x=-depth_cat,ymin=q.25,ymax=q.75,fill=water.depth.cat.2),
                width=10,alpha=0.1) +
  scale_color_viridis_d("Water Depth \nCategory(m)",begin=0,end=0.8) +
  scale_fill_viridis_d("Water Depth \nCategory(m)",begin=0,end=0.8) +
  #scale_y_continuous(label=LABEL,breaks=LABEL) +
  scale_shape(solid=FALSE) +
  xlab("Depth") +
  ylab("Copies") +
  coord_flip() +
  scale_y_continuous(limits=LIM, expand=c(0,0)) +
  theme_bw()
print(depth.p3)

depth.p4 <- depth.p3 +
  scale_y_continuous(trans="log2",label=LABEL,breaks=LABEL,expand=c(0,0.05)) 

print(depth.p4)




