# 
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(paletteer)
library(dplyr)

results.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Stan Model Fits/_Summarized_Output"
setwd(results.dir)
plot.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Plots and figures"


load("Acoustics 2019 lat.long.smooth 8_16_6_10_smooth_hurdle Base_Var Derived Q.RData")

D_acoustics_points    <- Acoustic.dat.figs$smooth.projections$D_acoustic_uncond_mt
D_acoustics_lat_equal <- Acoustic.dat.figs$smooth.projections$D_acoustic_uncond_mt_lat_equal
D_acoustics_lat_1.0   <- Acoustic.dat.figs$smooth.projections$D_acoustic_uncond_mt_lat_1.0
D_acoustics_lat_0.5   <- Acoustic.dat.figs$smooth.projections$D_acoustic_uncond_mt_lat_0.5

load("Qpcr_summaries_lat.long.smooth_Base_Var_5_10_fix_nu.RData")

D_DNA_points    <- Output.summary.qpcr$D_final_projected
D_DNA_lat_equal <- Output.summary.qpcr$D_final_lat_equal
D_DNA_lat_1.0   <- Output.summary.qpcr$D_final_lat_1.0
D_DNA_lat_0.5   <- Output.summary.qpcr$D_final_lat_0.5

### 
Both_points <- full_join(D_acoustics_points %>% rename(Mean_ac = Mean,
                                             Median_ac = Median,    
                                             SD_ac = SD,
                                             Q.0.01_ac = Q.0.01, 
                                             Q.0.025_ac =Q.0.025,
                                             Q.0.05_ac = Q.0.05,
                                             Q.0.25_ac = Q.0.25,
                                             Q.0.75_ac = Q.0.75,
                                             Q.0.95_ac = Q.0.95,
                                             Q.0.975_ac = Q.0.975,
                                             Q.0.99_ac = Q.0.99),
                         D_DNA_points %>% rename(Mean_dna = Mean,
                                                       Median_dna = Median,
                                                       SD_dna = SD,
                                                       Q.0.01_dna = Q.0.01, 
                                                       Q.0.025_dna =Q.0.025,
                                                       Q.0.05_dna = Q.0.05,
                                                       Q.0.25_dna = Q.0.25,
                                                       Q.0.75_dna = Q.0.75,
                                                       Q.0.95_dna = Q.0.95,
                                                       Q.0.975_dna = Q.0.975,
                                                       Q.0.99_dna = Q.0.99))



Both_equal <- full_join(D_acoustics_lat_equal %>% rename(Mean_ac = Mean,
                                                       Median_ac = Median,    
                                                       SD_ac = SD,
                                                       Q.0.01_ac = Q.0.01, 
                                                       Q.0.025_ac =Q.0.025,
                                                       Q.0.05_ac = Q.0.05,
                                                       Q.0.25_ac = Q.0.25,
                                                       Q.0.75_ac = Q.0.75,
                                                       Q.0.95_ac = Q.0.95,
                                                       Q.0.975_ac = Q.0.975,
                                                       Q.0.99_ac = Q.0.99),
                         D_DNA_lat_equal %>% rename(Mean_dna = Mean,
                                                 Median_dna = Median,
                                                 SD_dna = SD,
                                                 Q.0.01_dna = Q.0.01, 
                                                 Q.0.025_dna =Q.0.025,
                                                 Q.0.05_dna = Q.0.05,
                                                 Q.0.25_dna = Q.0.25,
                                                 Q.0.75_dna = Q.0.75,
                                                 Q.0.95_dna = Q.0.95,
                                                 Q.0.975_dna = Q.0.975,
                                                 Q.0.99_dna = Q.0.99))

Both_lat_0.5 <- full_join(D_acoustics_lat_0.5 %>% rename(Mean_ac = Mean,
                                                         Median_ac = Median,    
                                                         SD_ac = SD,
                                                         Q.0.01_ac = Q.0.01, 
                                                         Q.0.025_ac =Q.0.025,
                                                         Q.0.05_ac = Q.0.05,
                                                         Q.0.25_ac = Q.0.25,
                                                         Q.0.75_ac = Q.0.75,
                                                         Q.0.95_ac = Q.0.95,
                                                         Q.0.975_ac = Q.0.975,
                                                         Q.0.99_ac = Q.0.99),
                        D_DNA_lat_0.5 %>% rename(Mean_dna = Mean,
                                                   Median_dna = Median,
                                                   SD_dna = SD,
                                                   Q.0.01_dna = Q.0.01, 
                                                   Q.0.025_dna =Q.0.025,
                                                   Q.0.05_dna = Q.0.05,
                                                   Q.0.25_dna = Q.0.25,
                                                   Q.0.75_dna = Q.0.75,
                                                   Q.0.95_dna = Q.0.95,
                                                   Q.0.975_dna = Q.0.975,
                                                   Q.0.99_dna = Q.0.99))


Both_lat_1.0 <- full_join(D_acoustics_lat_1.0 %>% rename(Mean_ac = Mean,
                                                         Median_ac = Median,    
                                                         SD_ac = SD,
                                                         Q.0.01_ac = Q.0.01, 
                                                         Q.0.025_ac =Q.0.025,
                                                         Q.0.05_ac = Q.0.05,
                                                         Q.0.25_ac = Q.0.25,
                                                         Q.0.75_ac = Q.0.75,
                                                         Q.0.95_ac = Q.0.95,
                                                         Q.0.975_ac = Q.0.975,
                                                         Q.0.99_ac = Q.0.99),
                          D_DNA_lat_1.0 %>% rename(Mean_dna = Mean,
                                                   Median_dna = Median,
                                                   SD_dna = SD,
                                                   Q.0.01_dna = Q.0.01, 
                                                   Q.0.025_dna =Q.0.025,
                                                   Q.0.05_dna = Q.0.05,
                                                   Q.0.25_dna = Q.0.25,
                                                   Q.0.75_dna = Q.0.75,
                                                   Q.0.95_dna = Q.0.95,
                                                   Q.0.975_dna = Q.0.975,
                                                   Q.0.99_dna = Q.0.99))
Both_equal_1000 <- Both_equal
Both_lat_0.5_1000 <- Both_lat_0.5
Both_lat_1.0_1000 <- Both_lat_1.0
Both_equal_1000[,2:ncol(Both_equal_1000)] <-Both_equal_1000[,2:ncol(Both_equal_1000)] /1000
Both_lat_0.5_1000[,2:ncol(Both_lat_0.5_1000)] <- Both_lat_0.5_1000[,2:ncol(Both_lat_0.5_1000)] /1000
Both_lat_1.0_1000[,2:ncol(Both_lat_1.0_1000)] <- Both_lat_1.0_1000[,2:ncol(Both_lat_1.0_1000)] /1000

##############################################################################

ggplot(Both_equal) +
  geom_point(aes(x=Mean_ac,y=Mean_dna),alpha=0.2)+
  geom_errorbar(aes(x=Mean_ac,ymin=Q.0.05_dna,ymax=Q.0.95_dna),alpha=0.2)+
  geom_errorbarh(aes(y=Mean_dna,xmin=Q.0.05_ac,xmax=Q.0.95_ac),alpha=0.2)+
  theme_bw()

FILL.COL <- grey(0.8)

ggplot(Both_lat_0.5_1000) +
  geom_errorbar(aes(x=Mean_ac,ymin=Q.0.05_dna,ymax=Q.0.95_dna),alpha=1,width=0)+
  geom_errorbarh(aes(y=Mean_dna,xmin=Q.0.05_ac,xmax=Q.0.95_ac),alpha=1,height=0)+
  geom_point(aes(x=Mean_ac,y=Mean_dna),alpha=1,color=FILL.COL,size=2.5)+
  geom_text(aes(x=Mean_ac,y=Mean_dna,label=as.character(ID.lat.0.5)),alpha=1,size=2)+
  scale_x_continuous(expand=c(0,0),n.breaks=6) +
  scale_y_continuous(expand=c(0,0),n.breaks=6) +
  xlab("Acoustic Biomass (1000s mt)") +
  ylab("DNA Index (1000s)") +
  theme_classic()

p_cor_1.0 <- ggplot(Both_lat_1.0_1000) +
  geom_errorbar(aes(x=Mean_ac,ymin=Q.0.05_dna,ymax=Q.0.95_dna),alpha=1,width=0)+
  geom_errorbarh(aes(y=Mean_dna,xmin=Q.0.05_ac,xmax=Q.0.95_ac),alpha=1,height=0)+
  geom_point(aes(x=Mean_ac,y=Mean_dna),alpha=1,color=FILL.COL,size=2.5)+
  geom_text(aes(x=Mean_ac,y=Mean_dna,label=as.character(ID.lat.1.0)),alpha=1,size=2)+
  expand_limits(x=0,y=0) +
  scale_x_continuous(expand=c(0,0),n.breaks=6) +
  scale_y_continuous(expand=c(0,0),n.breaks=6) +
  xlab("Acoustic Biomass (1000s mt)") +
  ylab("DNA Index (1000s)") +
  theme_classic()

setwd(plot.dir)
quartz(file="Hake DNA-Acoustic bivariate correlation.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
print(p_cor_1.0)
dev.off()



cor.test(Both_equal$Mean_ac,Both_equal$Mean_dna)
cor.test(Both_lat_1.0$Mean_ac,Both_lat_1.0$Mean_dna)
cor.test(Both_lat_1.0$Mean_ac[Both_lat_1.0$ID.lat.1.0!=7],
         Both_lat_1.0$Mean_dna[Both_lat_1.0$ID.lat.1.0!=7])

cor.test(Both_lat_0.5$Mean_ac[Both_lat_0.5$ID.lat.0.5!=12],
         Both_lat_0.5$Mean_dna[Both_lat_0.5$ID.lat.0.5!=12])
cor.test(Both_lat_0.5$Mean_ac,
          Both_lat_0.5$Mean_dna)

cor.test(Both_points$Mean_ac,Both_points$Mean_dna)


# Make pairwise plot of acoustics and DNA with marginal densities.

max.dna = max(Both_points$Mean_dna) *1.01
max.ac = max(Both_points$Mean_ac) *1.01

g_base <- ggplot(Both_points) +
    geom_density2d_filled(aes(x=Mean_ac,y=Mean_dna),alpha=0.7,bins=8,shape=20) +  
    geom_point(aes(x=Mean_ac,y=Mean_dna),alpha=0.2)+
    xlab("Acoustic Biomass (mt)") +
    ylab("DNA Index") +
    scale_fill_brewer(palette = "Blues") +
    expand_limits(x=c(0,max.ac),y=c(-0.1,max.dna)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    theme(legend.position = "none")
    
g_top <- ggplot(Both_points) +
        geom_histogram(aes(Mean_ac),fill=grey(0.9),color="black",boundary=0) +
        theme_classic() +
        expand_limits(x=0,y=0.1) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        ylab("Frequency") +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())

g_right <- ggplot(Both_points) +
  geom_histogram(aes(y=Mean_dna),fill=grey(0.9),color="black",boundary=0) +
  theme_classic() +
  expand_limits(x=0,y=0.1) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Frequency") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

lay = rbind(c(1,1,NA),
            c(2,2,3),
            c(2,2,3))

setwd(plot.dir)
quartz(file="Hake bivariate point-level.jpeg",height=5.5,width=6,dpi=600,type="jpeg")
  grid.arrange(g_top,g_base,g_right,
      widths = c(1, 1, 0.5),
      heights=c(0.5,1,1),
      layout_matrix = lay)
dev.off()

################################################
### put DNA and Acoustics maps side-by-side.
################################################

transect_lines_summary <- Acoustic.dat.figs$dat.acoustic.bin %>%
                              group_by(transect) %>% 
                              summarise(lat= mean(lat),lon.min = min(lon),lon.max=max(lon))
CTD_locs <- Output.summary.qpcr$dat.obs.bin %>% 
                          group_by(station) %>% 
                          summarise(lat=mean(lat),lon=mean(lon))
  
survey.map <- Acoustic.dat.figs$base_map_trim_proj +
                  geom_segment(data=transect_lines_summary,aes(x=lon.min,xend=lon.max,y=lat,yend=lat)) +
                  geom_point(data=CTD_locs,aes(x=lon,y=lat),fill="transparent",color="red",shape=21)
  
survey.map              
  
  
size.regions = 3
y.axis.text <- 12

setwd(plot.dir)
quartz(file="Hake maps combined.jpeg",height=7.5,width=9,dpi=600,type="jpeg")
  grid.arrange(survey.map + xlab("") +
                geom_segment(data=lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                             linetype="dashed")+
                geom_text(data=lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
                 theme(plot.margin = unit(c(0.1,-6,0.1,-2), "lines"),
                       axis.text.x = element_text(size=y.axis.text),
                       axis.text.y = element_text(size=y.axis.text),
                       axis.title=element_text(size=y.axis.text+2)),
        
              Output.summary.qpcr$p_DNA_lat_base+
                geom_segment(data=lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                             linetype="dashed")+
                geom_text(data=lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
                ylab("") +
                theme( legend.position = c(1.07, .55),
                       legend.justification = c("right", "top"),
                       legend.box.just = "right",
                       plot.margin = unit(c(0.1,-4,0.1,-4), "lines"),
                       axis.text.x = element_text(size=y.axis.text),
                       axis.text.y = element_blank(),
                       axis.title=element_text(size=y.axis.text+2)),
        
              Acoustic.dat.figs$p_Acoustics_lat_base+
                geom_segment(data=lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                             linetype="dashed")+
                geom_text(data=lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
                ylab("") +
                xlab("") +
                theme( legend.position = c(1.10, .55),
                     legend.justification = c("right", "top"),
                     legend.box.just = "right",
                     plot.margin = unit(c(0.1,-2,0.1,-6), "lines"),
                     axis.text.x = element_text(size=y.axis.text),
                     axis.text.y = element_blank(),
                     axis.title=element_text(size=y.axis.text+2)),
               nrow=1)
  
dev.off()

################################################
# Spatial series.
################################################

acoustic.series <- ggplot(Both_lat_1.0_1000) +
  geom_point(aes(x=ID.lat.1.0,y=Mean_ac)) +
  geom_line(aes(x=ID.lat.1.0,y=Mean_ac),linetype="dashed") +
  geom_errorbar(aes(x=ID.lat.1.0,ymin=Q.0.25_ac,ymax=Q.0.75_ac),size=2,width=0) +
  geom_errorbar(aes(x=ID.lat.1.0,ymin=Q.0.05_ac,ymax=Q.0.95_ac),size=1,width=0) +
  geom_point(aes(x=ID.lat.1.0,y=Mean_ac),fill="white",size=3,shape=21) +
  expand_limits(x=c(0.5,11.5),y=c(0,max(Both_lat_1.0_1000$Q.0.99_ac))) +
  scale_x_continuous(expand=c(0,0),breaks=1:max(Both_lat_1.0_1000$ID.lat.1.0)) +
  scale_y_continuous(expand=c(0,NA))+
  ylab("Acoustic Biomass (1000s mt)") +
  xlab("Region") +
  theme_bw()
  

dna.series <- ggplot(Both_lat_1.0_1000) +
  geom_point(aes(x=ID.lat.1.0,y=Mean_dna)) +
  geom_line(aes(x=ID.lat.1.0,y=Mean_dna),linetype="dashed") +
  geom_errorbar(aes(x=ID.lat.1.0,ymin=Q.0.25_dna,ymax=Q.0.75_dna),size=2,width=0) +
  geom_errorbar(aes(x=ID.lat.1.0,ymin=Q.0.05_dna,ymax=Q.0.95_dna),size=1,width=0) +
  geom_point(aes(x=ID.lat.1.0,y=Mean_dna),fill="white",size=3,shape=21) +
  expand_limits(x=c(0.5,11.5),y=c(0,max(Both_lat_1.0_1000$Q.0.99_dna))) +
  scale_x_continuous(expand=c(0,0),breaks=1:max(Both_lat_1.0_1000$ID.lat.1.0)) +
  scale_y_continuous(expand=c(0,NA))+
  ylab("DNA Index (1000s)") +
  xlab("Region") +
  theme_bw()



gA <- ggplotGrob(acoustic.series)
gB <- ggplotGrob(dna.series)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1)

### Create dimensionless index for the two series by dividing by the mean value for the southernmost area.

mean_ac_1 <- Both_lat_1.0_1000 %>% filter(ID.lat.1.0==1) %>% pull(Mean_ac)
mean_dna_1 <- Both_lat_1.0_1000 %>% filter(ID.lat.1.0==1) %>% pull(Mean_dna)

Both_lat_1.0_1000_nondim <- Both_lat_1.0_1000 %>% dplyr::select(ID.lat.1.0, Mean_ac,
                                                                Q.0.05_ac,
                                                                Q.0.25_ac,
                                                                Q.0.75_ac,
                                                                Q.0.95_ac,
                                                                Mean_dna ,
                                                                Q.0.05_dna,
                                                                Q.0.25_dna,
                                                                Q.0.75_dna,
                                                                Q.0.95_dna) %>% 
                                                    mutate(Mean_ac=Mean_ac/mean_ac_1,
                                                           Q.0.05_ac= Q.0.05_ac/mean_ac_1,
                                                           Q.0.25_ac= Q.0.25_ac/ mean_ac_1,
                                                           Q.0.75_ac= Q.0.75_ac/ mean_ac_1,
                                                           Q.0.95_ac= Q.0.95_ac/ mean_ac_1,
                                                           Mean_dna=  Mean_dna / mean_dna_1,
                                                           Q.0.05_dna= Q.0.05_dna/mean_dna_1,
                                                           Q.0.25_dna= Q.0.25_dna/ mean_dna_1,
                                                           Q.0.75_dna= Q.0.75_dna/ mean_dna_1,
                                                           Q.0.95_dna= Q.0.95_dna/ mean_dna_1)
plus = 0.10
minus= -plus

both.series_nondim <- ggplot(Both_lat_1.0_1000_nondim) +
  #geom_point(aes(x=ID.lat.1.0+plus,y=Mean_dna)) +
  geom_line(aes(x=ID.lat.1.0+plus,y=Mean_dna),linetype="dashed") +
  geom_errorbar(aes(x=ID.lat.1.0+plus,ymin=Q.0.25_dna,ymax=Q.0.75_dna),size=2,width=0) +
  geom_errorbar(aes(x=ID.lat.1.0+plus,ymin=Q.0.05_dna,ymax=Q.0.95_dna),size=1,width=0) +
  geom_point(aes(x=ID.lat.1.0+plus,y=Mean_dna,shape="DNA",fill="white"),size=3) +
  
  geom_line(aes(x=ID.lat.1.0+minus,y=Mean_ac),linetype="dotted") +
  geom_errorbar(aes(x=ID.lat.1.0+minus,ymin=Q.0.25_ac,ymax=Q.0.75_ac),size=2,width=0) +
  geom_errorbar(aes(x=ID.lat.1.0+minus,ymin=Q.0.05_ac,ymax=Q.0.95_ac),size=1,width=0) +
  geom_point(aes(x=ID.lat.1.0+minus,y=Mean_ac,shape="Acoustic",fill="white"),size=3) +
  expand_limits(x=c(0.5,11.5),y=c(0,max(Both_lat_1.0_1000_nondim$Q.0.95_ac))*1.05) +
  scale_x_continuous(expand=c(0,0),breaks=1:max(Both_lat_1.0_1000$ID.lat.1.0)) +
  scale_y_continuous(expand=c(0,NA))+
  ylab("Relative Index") +
  xlab("Region") +
  scale_shape_manual("Survey",values = c("DNA"=21,"Acoustic"=22)) +
  scale_fill_manual("Survey",values = c("white"),labels=c("DNA","Acoustic")) +
  theme_bw() + 
  guides(fill=F) +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right")

quartz(file="Hake relative index.jpeg",height=4,width=6,dpi=600,type="jpeg")
  print(both.series_nondim)
dev.off()

######
### DO SOME CALCULATIONS ABOUT THE VARIABILITY IN THE SPATIAL FIELDS...
######
EXPAND <- 40
LOG.EXPAND <- log10(EXPAND)

D_pred_smooth <- Output.summary.qpcr$smooth.projections$D_pred_smooth_combined

ggplot(D_pred_smooth) +
  geom_point(aes(x=depth_cat_factor,y=Mean-LOG.EXPAND),alpha=0.5)










