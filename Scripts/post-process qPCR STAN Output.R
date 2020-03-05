###

library(tidyverse)
library(knitr)
library(reshape2)
library(viridis)
library(ggmap)
library(maps)
library(mapdata)
library(ggplot2)



# Read in the output from the Stan model
base.dir <- "/Users/ole.shelton/Github/eDNA-Hake/Stan Model Fits"
setwd(base.dir)
load("qPCR Hake 2019 hake Fitted.RData")

### Get some of the output on sample and station IDs from file.
dat.samp <- Output.qpcr$dat.obs.bin
dat.control.field.neg <- Output.qpcr$dat.control.field.neg %>% left_join(.,Output.qpcr$field_neg_out)

dat.summary <- dat.samp %>%
  group_by(date,transect,station,station.depth,depth,station_depth_idx,) %>%
  dplyr::summarise(N=length(unique(sample_idx)))
dat.summary <- left_join(dat.summary,Output.qpcr$station_depth_out)
dat.summary <- dat.summary %>% mutate(depth_cat=case_when(depth < 25~0,
                                                          depth == 25~ 25,  
                                                          depth > 25 & depth <=50 ~ 50,
                                                          depth > 50 & depth <=100 ~ 100,
                                                          depth > 120 & depth <=150 ~ 150,
                                                          depth > 151 & depth <=200 ~ 200,
                                                          depth > 250 & depth <=350 ~ 300,
                                                          depth > 400 & depth <=500 ~ 500))

dat.inhibit <- Output.qpcr$dat.inhibit %>% mutate(depth_cat=case_when(depth < 25~ 0,
                                                                      depth == 25~ 25,            
                                                                      depth > 25 & depth <=50 ~ 50,
                                                                      depth > 50 & depth <=100 ~ 100,
                                                                      depth > 120 & depth <=150 ~ 150,
                                                                      depth > 151 & depth <=200 ~ 200,
                                                                      depth > 250 & depth <=350 ~ 300,
                                                                      depth > 400 & depth <=500 ~ 500))

##### Calculate some info about the field negative controls
field.neg.summ <- Output.qpcr$field_neg_out 

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

dat.A <- dat.station.id.trim %>% dplyr::select(station,lat,lon) %>% mutate(Level="CTD Stations")
dat.B <- left_join(dat.sample.id,dat.id) %>% dplyr::select(station,lat,lon) %>%
  group_by(station,lat,lon) %>% summarise(N=length(station)) %>% filter(is.na(lat)==F) %>%
  mutate(Level="Sampled Stations")
dat.C <- dat.samp %>% dplyr::select(station,lat,lon) %>% group_by(station,lat,lon) %>% summarise(N=length(station)) %>%
  dplyr::select(-N) %>% mutate(Level="Processed Stations")


### THIS IS ALL OF THE CTD STATIONS.  IF YOU NEED A LAT LONG FOR EACH STATION, MERGE IN FROM HERE.
dat.process <- bind_rows(dat.A,dat.B,dat.C)
dat.process$Level <- factor(dat.process$Level, levels=c("CTD Stations","Sampled Stations","Processed Stations"))

dat.lat.lon <- dat.process %>% dplyr::select(station,lat,lon) %>% 
  group_by(station) %>% summarise(mean.lat=mean(lat),mean.lon=mean(lon)) %>% 
  filter(!station=="") %>% rename(lat=mean.lat,lon=mean.lon)


#### Plotting information
states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington"))

# Plot for all samples
lat.lims <- c(min(dat.id$lat,na.rm=T),max(dat.id$lat,na.rm=T))
lon.lims <- c(min(dat.id$lon,na.rm=T),max(dat.id$lon,na.rm=T))

base_map <-ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
  coord_fixed(xlim=lon.lims,ylim=lat.lims,ratio=1.3) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

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
  left_join(.,dat.lat.lon)

#####
g_i_3 <- g_i_2 %>% mutate(inhib_cat=case_when(N==1 & Y==0 ~ "None",
                                              N==0 & Y==1 ~ "All",
                                              N==1 & Y==1 ~ "Partial"))
g_i_3$inhib_cat <- factor(g_i_3$inhib_cat,levels=c("All","Partial","None"))

### Make a spatial plot of inhibition

inhibit.plot <- base_map +
  geom_point(data=g_i_3,aes(x=lon,y=lat,color=inhib_cat)) +
  scale_color_viridis_d(begin=0,end=0.8)+
  facet_wrap(~inhib_cat)

inhibit.plot.by.depth <- base_map +
  geom_point(data=g_i_3 %>% filter(is.na(depth_cat)==F),aes(x=lon,y=lat,color=inhib_cat),alpha=0.5) +
  facet_wrap(~depth_cat,nrow = 2) +
  scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) 


inhibit.plot.by.depth.few <- base_map +
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

BREAKS <- seq(0,30,by=0.1)
XLIM <- c(0,max(field.neg.summ$Mean)*1.1)
control.hist <- ggplot(field.neg.summ) +
  geom_histogram(aes(Mean),breaks=BREAKS) +
  xlab("Estimated copies per uL") +
  ylab("Number of negative controls") +
  xlim(XLIM) +
  theme_bw()

ALPHA = 0.5 
control.by.time <- ggplot(field.neg.summ) +
  geom_point(aes(x=date,y=Mean),alpha=ALPHA) +
  geom_errorbar(aes(x=date,ymin=Val.X5.,ymax=Val.X95.),alpha=ALPHA,width=0) +
  xlab("Date") +
  ylab("Copies per uL") +
  theme_bw()


N_control_positive <- field.neg.summ %>% filter(Mean>0.10) %>% nrow(.)
N_control_run <- nrow(field.neg.summ)
N_control_total <- nrow(dat.control.field.neg)

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



BREAKS <- 10^(c(-4,seq(-1,3,by=0.25)))
XLIM <- c(0,max(dat.sum.v.control$Mean)*1.01)
compare.hist <- ggplot(dat.sum.v.control) +
  geom_histogram(aes(Mean,fill=cat),breaks=BREAKS,alpha=0.5,position="identity",) +
  xlab("Estimated copies per uL") +
  ylab("Number of Samples") +
  xlim(XLIM) +
  scale_fill_viridis_d("Category",begin=0,end=0.5,) +
  scale_x_log10() +
  theme_bw()
compare.hist

#######################################################
#######################################################
####### ------- Make Maps
#######################################################
#######################################################


dat.hake <- left_join(dat.summary,dat.lat.lon) %>% filter(is.na(depth_cat)==F)

BREAKS = c(-100,-2,0,1,2,3,4,5)
hake.p1 <- base_map + 
  geom_point(data=dat.hake,aes(x=lon,y=lat,size=Mean.log),shape=21) +
  scale_size("Log Copies",breaks=BREAKS)+
  scale_shape("Log Copies",solid=FALSE)+
  facet_wrap(~depth_cat)
hake.p1

LAB<-c(0.1,1,10,50,100,200,400)
hake.p2 <- base_map + 
  geom_point(data=dat.hake,aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
  scale_size("Copies",labels=LAB,breaks=LAB,range=c(0.1,6))+
  scale_shape("Copies",solid=FALSE)+
  facet_wrap(~depth_cat,nrow=2)
#hake.p2

hake.p3 <- base_map + 
  geom_point(data=dat.hake%>%filter(depth_cat %in% c(0,50,150,300,500)),aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
  scale_size("Copies",labels=LAB,breaks=LAB,range=c(0.1,6))+
  scale_shape("Copies",solid=FALSE)+
  facet_wrap(~depth_cat,nrow=2)
#hake.p3

######################## 



#### NEED TO DIFFERENTIATE AMONG CTDs in very deep water and ~500m.
max.depth <- dat.hake %>% group_by(station) %>% summarise(max.depth.cat = max(depth_cat)) %>%
  mutate(max.depth.cat = ifelse(max.depth.cat<50,50,max.depth.cat))
dat.hake <- dat.hake %>% arrange(station, depth) %>% mutate(neg.depth = -depth) %>%
  left_join(.,max.depth)

dat.hake.sum.by.depth.cat <- dat.hake %>% group_by(max.depth.cat,depth_cat) %>%
  summarise(N=length(Mean),MEAN=mean(Mean),SD=sd(Mean),
            MEDIAN=median(Mean),
            q.25=quantile(Mean,probs=0.25),
            q.75=quantile(Mean,probs=0.75))

LAB<-c(0.1,1,10,50,100,200,400)
depth.p1 <- ggplot(dat.hake %>% filter(max.depth.cat %in% c(150,300))) +
  geom_point(aes(x=lat,y=neg.depth,size=Mean),alpha=0.5,color="red",shape=21) +
  #geom_line(aes(x=neg.depth,y=Mean,group=station,color=lat),alpha=0.3) +
  #scale_y_continuous(trans="log2") +
  scale_size("Copies",labels=LAB,breaks=LAB,range=c(0.1,6))+
  scale_shape("Copies",solid=FALSE) +
  #scale_scale_color_viridis(begin=0,end=0.8) +
  xlab("Latitude") +
  ylab("Depth") +
  theme_bw()
print(depth.p1)


DEPTH.MAX <- 500
LABEL = c(0.5,1,2,4,8,16,32,64,128,256) 
depth.p2 <- ggplot(dat.hake %>% filter(max.depth.cat %in% c(DEPTH.MAX))) +
  geom_point(aes(x=neg.depth,y=Mean,group=station,color=lat),alpha=0.5) +
  geom_line(aes(x=neg.depth,y=Mean,group=station,color=lat),alpha=0.3) +
  geom_point(data= dat.hake.sum.by.depth.cat %>% 
               filter(max.depth.cat %in% c(DEPTH.MAX)),
             aes(x=-depth_cat,y=MEDIAN),color="red",shape=16,size=3) +
  geom_point(data= dat.hake.sum.by.depth.cat %>% 
               filter(max.depth.cat %in% c(DEPTH.MAX)),
             aes(x=-depth_cat,y=MEAN),color="red",shape=22,size=3) +
  geom_line(data= dat.hake.sum.by.depth.cat %>% 
              filter(max.depth.cat %in% c(DEPTH.MAX)),
            aes(x=-depth_cat,y=MEDIAN),color="red") +
  geom_errorbar(data= dat.hake.sum.by.depth.cat %>% 
                  filter(max.depth.cat %in% c(DEPTH.MAX)),
                aes(x=-depth_cat,ymin=q.25,ymax=q.75),color="red",width=0) +
  scale_color_viridis("Latitude",begin=0,end=0.8) +
  scale_y_continuous(trans="log2",label=LABEL,breaks=LABEL) +
  scale_shape(solid=FALSE) +
  xlab("Depth") +
  ylab("Copies") +
  coord_flip() +
  theme_bw()
print(depth.p2)


# dat.hake.sum.by.depth.cat <- dat.hake.sum.by.depth.cat %>%
#                               mutate(max.depth.cat.fact=as.factor(max.depth.cat))
# LEV <- unique(dat.hake.sum.by.depth.cat$max.depth.cat.fact)
# 
# dat.hake.sum.by.depth.cat$max.depth.cat.fact <- factor(dat.hake.sum.by.depth.cat$max.depth.cat.fact,
#                                                        levels=LEV)

d.temp <- dat.hake.sum.by.depth.cat %>% 
  filter(!max.depth.cat %in% c(0,200)) %>% 
  arrange(max.depth.cat)
d.temp$max.depth.cat = factor(d.temp$max.depth.cat,levels=c("50","100","150","300","500"))

LIM <- c(0,max(d.temp$q.75)*1.05)
depth.p3 <- ggplot(d.temp) +
  geom_point(aes(x=-depth_cat,y=MEDIAN,color=max.depth.cat),
             shape=16,size=3,alpha=0.8) +
  # geom_point(aes(x=-depth_cat,y=MEAN,color=max.depth.cat.fact),
  #           shape=16,size=3) +
  geom_line(aes(x=-depth_cat,y=MEDIAN,color=max.depth.cat),alpha=0.5) + 
  geom_errorbar(aes(x=-depth_cat,ymin=q.25,ymax=q.75,color=max.depth.cat),
                width=10,alpha=0.5) +
  scale_color_viridis_d("Max \nDepth",begin=0,end=0.8) +
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




