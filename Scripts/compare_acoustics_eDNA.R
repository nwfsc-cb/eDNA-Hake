# 
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(paletteer)
library(dplyr)
library(cowplot)
library(viridis)

results.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Stan Model Fits/_Summarized_Output"
setwd(results.dir)
plot.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Plots and figures"

load("Acoustics 2019 lat.long.smooth 6_14_6_10_smooth_hurdle Base_Var Derived Q.RData")

D_acoustics_points    <- Acoustic.dat.figs$smooth.projections$D_acoustic_uncond_mt
D_acoustics_lat_equal <- Acoustic.dat.figs$smooth.projections$D_acoustic_uncond_mt_lat_equal
D_acoustics_lat_1.0   <- Acoustic.dat.figs$smooth.projections$D_acoustic_uncond_mt_lat_1.0
D_acoustics_lat_0.5   <- Acoustic.dat.figs$smooth.projections$D_acoustic_uncond_mt_lat_0.5

load("Qpcr_summaries_lat.long.smooth_Base_Var_4_10_fix_nu_T-FIN_NoSURF=FALSE.RData")

D_DNA_points    <- Output.summary.qpcr$D_final_projected
D_DNA_lat_equal <- Output.summary.qpcr$D_final_lat_equal
D_DNA_lat_1.0   <- Output.summary.qpcr$D_final_lat_1.0
D_DNA_lat_0.5   <- Output.summary.qpcr$D_final_lat_0.5

p_D    <- Output.summary.qpcr$p_D
p_D_CV <- Output.summary.qpcr$p_D_CV

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

##### Calculate pairwise correlation between eDNA and Acoustics using MCMC draws.
A.resamp.grid <- Acoustic.dat.figs$D_grid.cell_uncond_resample %>% rename(A.val = D)
DNA.resamp.grid <- Output.summary.qpcr$D_grid.cell_resample %>% rename(DNA.val = tot)

A.resamp.1.0 <-Acoustic.dat.figs$D_1.0_uncond_resample %>% rename(A.val = tot)
DNA.resamp.1.0 <- Output.summary.qpcr$D_1.0_resample %>% rename(DNA.val = tot)

A.resamp.0.5 <-Acoustic.dat.figs$D_0.5_uncond_resample %>% rename(A.val = tot)
DNA.resamp.0.5 <- Output.summary.qpcr$D_0.5_resample  %>% rename(DNA.val = tot)

# Do CI for grid.cell by grid.cell
resamp.grid.comb <- full_join(A.resamp.grid,DNA.resamp.grid)

resamp.grid.summary <- resamp.grid.comb %>% group_by(Gridcell_ID) %>% summarise(A.mean = mean(A.val),DNA.mean =mean(DNA.val))

cor.grid <- resamp.grid.comb %>% group_by(MCMC.rep) %>% summarise(COR =cor(A.val,DNA.val,method="pearson"))

ggplot(resamp.grid.comb) +
    geom_point()

# Get CI from 1.0 degree bins
resamp.1.0.comb <- full_join(A.resamp.1.0,DNA.resamp.1.0)
cor.1.0.vals <- resamp.1.0.comb %>% group_by(MCMC.rep) %>% summarise(COR =cor(A.val,DNA.val,method="pearson")) 
cor.1.0.vals.summary <- cor.1.0.vals %>% ungroup() %>% summarise(Mean = mean(COR), 
                            Median = median(COR),
                            Q.25 = quantile(COR,probs=0.25),
                            Q.75 = quantile(COR,probs=0.75),
                            Q.05 = quantile(COR,probs=0.05),
                            Q.95 = quantile(COR,probs=0.95),
                            Q.025 = quantile(COR,probs=0.025),
                            Q.975 = quantile(COR,probs=0.975))

# Get CI from 0.5 degree bins
resamp.0.5.comb <- full_join(A.resamp.0.5,DNA.resamp.0.5)

cor.0.5.vals <- resamp.0.5.comb %>% group_by(MCMC.rep) %>% summarise(COR =cor(A.val,DNA.val,method="pearson"))
cor.0.5.vals.summary <- cor.0.5.vals %>%
  ungroup() %>% summarise(Mean = mean(COR), 
                          Median = median(COR),
                          Q.25 = quantile(COR,probs=0.25),
                          Q.75 = quantile(COR,probs=0.75),
                          Q.05 = quantile(COR,probs=0.05),
                          Q.95 = quantile(COR,probs=0.95),
                          Q.025 = quantile(COR,probs=0.025),
                          Q.975 = quantile(COR,probs=0.975))
                            

#
##############################################################################

ggplot(Both_equal) +
  geom_point(aes(x=Mean_ac,y=Mean_dna),alpha=0.2)+
  geom_errorbar(aes(x=Mean_ac,ymin=Q.0.05_dna,ymax=Q.0.95_dna),alpha=0.2)+
  geom_errorbarh(aes(y=Mean_dna,xmin=Q.0.05_ac,xmax=Q.0.95_ac),alpha=0.2)+
  theme_bw()

FILL.COL <- grey(0.8)

p_cor_0.5 <- ggplot(Both_lat_0.5_1000) +
  geom_errorbar(aes(x=Mean_ac,ymin=Q.0.05_dna,ymax=Q.0.95_dna),alpha=1,width=0)+
  geom_errorbarh(aes(y=Mean_dna,xmin=Q.0.05_ac,xmax=Q.0.95_ac),alpha=1,height=0)+
  geom_point(aes(x=Mean_ac,y=Mean_dna),alpha=1,color=FILL.COL,size=2.5)+
  geom_text(aes(x=Mean_ac,y=Mean_dna,label=as.character(ID.lat.0.5)),alpha=1,size=2)+
  expand_limits(x=0,y=0) +
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
quartz(file="Hake_DNA-Acoustic_bivariate_correlation_1.0.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
print(p_cor_1.0)
dev.off()

setwd(plot.dir)
quartz(file="Hake_DNA-Acoustic_bivariate_correlation_0.5.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
print(p_cor_0.5)
dev.off()



cor.test(Both_equal$Mean_ac,Both_equal$Mean_dna,conf.level=0.90)
cor.test(Both_lat_1.0$Mean_ac,Both_lat_1.0$Mean_dna,conf.level=0.90)
# cor.test(Both_lat_1.0$Mean_ac[Both_lat_1.0$ID.lat.1.0!=7],
#          Both_lat_1.0$Mean_dna[Both_lat_1.0$ID.lat.1.0!=7])

# cor.test(Both_lat_0.5$Mean_ac[Both_lat_0.5$ID.lat.0.5!=12],
#          Both_lat_0.5$Mean_dna[Both_lat_0.5$ID.lat.0.5!=12])
cor.test(Both_lat_0.5$Mean_ac,
          Both_lat_0.5$Mean_dna)

cor.test(Both_points$Mean_ac,Both_points$Mean_dna,conf.level=0.90)

###########################################################################################
###########################################################################################
# Make pairwise plot of acoustics and DNA with marginal densities.

max.dna = max(Both_points$Mean_dna) *1.01
max.ac = max(Both_points$Mean_ac) *1.01

g_base <- ggplot(Both_points) +
    #geom_density2d_filled(aes(x=Mean_ac,y=Mean_dna),alpha=0.7,bins=8,shape=20) +  
    geom_point(aes(x=Mean_ac,y=Mean_dna),alpha=0.2)+
    xlab("Acoustic Biomass (mt)") +
    ylab("DNA Index") +
    #scale_fill_brewer(palette = "Blues") +
    expand_limits(x=c(0,max.ac),y=c(-0.1,max.dna)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    annotate(geom="text",x=120,y=6700,label="A") +
    theme(legend.position = "none")
    
g_top <- ggplot(Both_points) +
        geom_histogram(aes(Mean_ac),fill=grey(0.9),color="black",boundary=0) +
        theme_classic() +
        expand_limits(x=c(0,max.ac),y=c(0,1000)) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0),breaks=c(0,250,500,750,1000)) +
        ylab("Frequency") +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())

g_right <- ggplot(Both_points) +
  geom_histogram(aes(y=Mean_dna),fill=grey(0.9),color="black",boundary=0) +
  theme_classic() +
  expand_limits(x=0,y=0.1) +
  scale_x_continuous(expand=c(0,0),breaks=c(0,200,420)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Frequency") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

lay = rbind(c(1,1,NA),
            c(2,2,3),
            c(2,2,3))

setwd(plot.dir)
quartz(file="Hake_bivariate_point-level.jpeg",height=5.5,width=6,dpi=600,type="jpeg")
  grid.arrange(g_top,g_base,g_right,
      widths = c(1, 1, 0.5),
      heights=c(0.5,1,1),
      layout_matrix = lay)
dev.off()

# Combine point level and 1 degree plots into one figure.
setwd(plot.dir)
quartz(file="Hake_point-level_and_1_degree.pdf",height=4,width=8,dpi=900,type="pdf")
grid.arrange(
    g_top,
    g_base,
    g_right ,
    p_cor_1.0 +
      annotate(geom="text",x=15,y=1445,label="B"),
  widths=c(1,0.25,1),
  heights=c(0.25,1),
  layout_matrix=rbind(c(1,NA,NA),
                      c(2,3,4))
  )
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
y.axis.text <- 9

setwd(plot.dir)
quartz(file="Hake_maps_combined_to_surface.pdf",height=6,width=7.5,dpi=900,type="pdf")
  grid.arrange(survey.map + xlab("") +
                geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                             linetype="dashed")+
                geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
                 annotate(geom="text",x=-123.3,y=48.35,label="A") +
                 theme(plot.margin = unit(c(0.1,-6,0.1,-2), "lines"),
                       axis.text.x = element_text(size=y.axis.text),
                       axis.text.y = element_text(size=y.axis.text),
                       axis.title=element_text(size=y.axis.text+2)),
        
              Output.summary.qpcr$p_DNA_lat_base+
                geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                             linetype="dashed")+
                geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
                annotate(geom="text",x=-123.3,y=48.35,label="B") +
                ylab("") +
                theme( legend.position = c(1.15, .55),
                       legend.justification = c("right", "top"),
                       legend.box.just = "right",
                       plot.margin = unit(c(0.1,-4,0.1,-4), "lines"),
                       axis.text.x = element_text(size=y.axis.text),
                       axis.text.y = element_blank(),
                       axis.title=element_text(size=y.axis.text+2)),
        
              Acoustic.dat.figs$p_Acoustics_lat_base+
                geom_segment(data=Acoustic.dat.figs$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                             linetype="dashed")+
                geom_text(data=Acoustic.dat.figs$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
                ylab("") +
                xlab("") +
                annotate(geom="text",x=-123.3,y=48.35,label="C") +
                theme( legend.position = c(1.23, .55),
                     legend.justification = c("right", "top"),
                     legend.box.just = "right",
                     plot.margin = unit(c(0.1,-2,0.1,-6), "lines"),
                     axis.text.x = element_text(size=y.axis.text),
                     axis.text.y = element_blank(),
                     axis.title=element_text(size=y.axis.text+2)),
               nrow=1)
  
dev.off()


# Repeat with alternate gridding of coast
setwd(plot.dir)
quartz(file="Hake_maps_combined_to_surface_half_degree.jpeg",height=7.5,width=9,dpi=600,type="jpeg")
grid.arrange(survey.map + xlab("") +
               geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.0.5,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                            linetype="dashed")+
               geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.0.5,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
               theme(plot.margin = unit(c(0.1,-6,0.1,-2), "lines"),
                     axis.text.x = element_text(size=y.axis.text),
                     axis.text.y = element_text(size=y.axis.text),
                     axis.title=element_text(size=y.axis.text+2)),
             
             Output.summary.qpcr$p_DNA_lat_base+
               geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.0.5,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                            linetype="dashed")+
               geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.0.5,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
               ylab("") +
               theme( legend.position = c(1.07, .55),
                      legend.justification = c("right", "top"),
                      legend.box.just = "right",
                      plot.margin = unit(c(0.1,-4,0.1,-4), "lines"),
                      axis.text.x = element_text(size=y.axis.text),
                      axis.text.y = element_blank(),
                      axis.title=element_text(size=y.axis.text+2)),
             
             Acoustic.dat.figs$p_Acoustics_lat_base+
               geom_segment(data=Acoustic.dat.figs$lat.breaks$lats.rounded.0.5,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                            linetype="dashed")+
               geom_text(data=Acoustic.dat.figs$lat.breaks$lats.rounded.0.5,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
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

#####################################33
setwd(plot.dir)
#### Make Grid of MEAN Predicitons 

lay = rbind(c(1,2,3,4),
            c(5,6,7,8))


axis.text <- 9
title.text <-11
quartz(file="Hake_maps_Mean_depth_manual_facet.jpeg",height=9,width=8,dpi=600,type="jpeg")
  grid.arrange(survey.map + xlab("") + ggtitle("Survey")+
                 #geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                 #              linetype="dashed")+
                 #geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
                 theme(plot.margin = unit(c(0.1,-10,-0.8,-2), "lines"),
                       axis.text.x = element_text(size=axis.text),
                       axis.text.y = element_text(size=axis.text),
                       axis.title=element_text(size=title.text),
                       plot.title = element_text(size=title.text,vjust=-1)),
              p_D$x_0 + theme(legend.position = "none",
                               plot.margin = unit(c(0.1,-8,-0.8,-4), "lines"),
                               axis.text.x = element_text(size=axis.text,color="white"),
                               axis.text.y = element_text(size=axis.text,color="white"),
                               axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                               axis.title.x = element_text(size=title.text,color="white"),
                               plot.title = element_text(size=title.text,vjust=-1),
                                ),
               p_D$x_50 + theme(legend.position = "none",
                                plot.margin = unit(c(0.1,-6,-0.8,-6), "lines"),
                                axis.text.x = element_text(size=axis.text,color="white"),
                                axis.text.y = element_text(size=axis.text,color="white"),
                                axis.title.x = element_text(size=title.text,color="white"),
                                axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                                plot.title = element_text(size=title.text,vjust=-1),
                                   ),
  
               p_D$x_100 + theme(legend.position = "none",
                                 plot.margin = unit(c(0.1,-4,-0.8,-8), "lines"),
                                 axis.text.x = element_text(size=axis.text,color="white"),
                                 axis.text.y = element_text(size=axis.text,color="white"),
                                 axis.title.x = element_text(size=title.text,color="white"),
                                 axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                                 plot.title = element_text(size=title.text,vjust=-1),
                                 ),
              get_legend(p_D$x_0 + theme(legend.margin = margin(0,-70,0,70),
                                         legend.title = element_text(size=title.text))),
              
               p_D$x_150 + theme(legend.position = "none",
                                 plot.margin = unit(c(-0.8,-8,0.1,-4), "lines"),
                                 axis.text.x = element_text(size=axis.text),
                                 axis.text.y = element_text(size=axis.text),
                                 axis.title.y = element_text(size=title.text),
                                 axis.title.x = element_text(size=title.text,color="white"),
                                 plot.title = element_text(size=title.text,vjust=-1)
                                  ),
               p_D$x_300 + theme(legend.position = "none",
                                 plot.margin = unit(c(-0.8,-6,0.1,-6), "lines"),
                                 axis.text.x = element_text(size=axis.text),
                                 axis.text.y = element_text(size=axis.text,color="white"),
                                 axis.title.x = element_text(size=title.text),
                                 axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                                 plot.title = element_text(size=title.text,vjust=-1)
                                  ),
               p_D$x_500 + theme(legend.position = "none",
                                 plot.margin = unit(c(-0.8,-4,0.1,-8), "lines"),
                                 axis.text.x = element_text(size=axis.text),
                                 axis.text.y = element_text(size=axis.text,color="white"),
                                 axis.title.x = element_text(size=title.text,color="white"),
                                 axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                                 plot.title = element_text(size=title.text,vjust=-1),
                                     ),
               
               widths=c(1,1,1,1),
               layout_matrix = lay)
dev.off()                 






#### Make Grid of MEAN Predicitons VERSION 2

# Make marginal plot of CV vs. depth.  

CV.summary <- Output.summary.qpcr$D_pred_combined %>% group_by(depth_cat_factor) %>% 
        summarise(grand.mean.cv = mean(CV),
                  grand.median.cv =median(CV),
                  cv.q.025 = quantile(CV,probs=0.025),
                  cv.q.05 = quantile(CV,probs=0.05),
                  cv.q.10 = quantile(CV,probs=0.10),
                  cv.q.25 = quantile(CV,probs=0.25),
                  cv.q.75 = quantile(CV,probs=0.75),
                  cv.q.90 = quantile(CV,probs=0.90),
                  cv.q.95 = quantile(CV,probs=0.95),
                  cv.q.975 = quantile(CV,probs=0.975)
                ) %>%
        mutate(plot_depth=ifelse(depth_cat_factor==0,3,depth_cat_factor))
  
MAX.CV <- CV.summary$cv.q.90 %>% max()

# Total coastwide CV.
Output.summary.qpcr$D_DNA_uncond_total$SD / Output.summary.qpcr$D_DNA_uncond_total$Mean
Acoustic.dat.figs$D_acoustic_uncond_total_mt$SD / Acoustic.dat.figs$D_acoustic_uncond_total_mt$Mean

CV_marg_plot <- ggplot(CV.summary) +
    geom_errorbarh(aes(y=plot_depth,xmin=cv.q.10,xmax=cv.q.90),height=0) +
    geom_errorbarh(aes(y=plot_depth,xmin=cv.q.25,xmax=cv.q.75),height=0,size=2) +
    geom_point(aes(y=plot_depth,x=grand.mean.cv),shape=21,fill="white",size=1.5) +
    geom_point(aes(y=plot_depth,x=grand.median.cv),shape="|",size=5) +
    scale_y_reverse("Depth (m)")+
    scale_x_continuous("CV",expand=c(0.0,0.0),limits = c(0,2.75),
                       breaks=c(0,0.5,1,2,3), labels=c(0,0.5,1,2,3)) +
    theme_classic() +
    ggtitle("g)")


lay = rbind(c(1,2,3,4),
            c(5,6,7,8))

axis.text <- 9
title.text <-11
quartz(file="Hake_maps_Mean_depth_manual_facet_V2.pdf",height=9,width=8,dpi=900,type="pdf")
grid.arrange(
               # #geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
               # #              linetype="dashed")+
               # #geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
               # theme(plot.margin = unit(c(0.1,-10,-0.8,-2), "lines"),
               #       axis.text.x = element_text(size=axis.text),
               #       axis.text.y = element_text(size=axis.text),
               #       axis.title=element_text(size=title.text),
               #       plot.title = element_text(size=title.text,vjust=-1)),
             p_D$x_0 +  ggtitle("A 3m") +
                      theme(legend.position = "none",
                             plot.margin = unit(c(0.1,0,-0.8,0), "lines"),
                             axis.text.x = element_text(size=axis.text,color="white"),
                             axis.text.y = element_text(size=axis.text),
                             axis.title.y = element_text(size=title.text),
                             axis.title.x = element_text(size=title.text,color="white"),
                             plot.title = element_text(size=title.text,vjust=-1),
             ),
             p_D$x_50 + ggtitle("B 50m") + 
               theme(legend.position = "none",
                              plot.margin = unit(c(0.1,-2,-0.8,-6), "lines"),
                              axis.text.x = element_text(size=axis.text,color="white"),
                              axis.text.y = element_text(size=axis.text,color="white"),
                              axis.title.x = element_text(size=title.text,color="white"),
                              axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                              plot.title = element_text(size=title.text,vjust=-1),
             ),
             
             p_D$x_100 + ggtitle("C 100m") +
               theme(legend.position = "none",
                               plot.margin = unit(c(0.1,-2,-0.8,-10), "lines"),
                               axis.text.x = element_text(size=axis.text,color="white"),
                               axis.text.y = element_text(size=axis.text,color="white"),
                               axis.title.x = element_text(size=title.text,color="white"),
                               axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                               plot.title = element_text(size=title.text,vjust=-1),
             ),
             get_legend(p_D$x_0 + theme(legend.margin = margin(0,70,0,-70),
                                        legend.title = element_text(size=title.text))),
             
             p_D$x_150 + ggtitle("D 150m") +
               theme(legend.position = "none",
                               plot.margin = unit(c(-0.8,-0,0.1,0), "lines"),
                               axis.text.x = element_text(size=axis.text),
                               axis.text.y = element_text(size=axis.text),
                               axis.title.y = element_text(size=title.text),
                               axis.title.x = element_text(size=title.text,color="white"),
                               plot.title = element_text(size=title.text,vjust=-1)
             ),
             p_D$x_300 + ggtitle("E 300m") +
               theme(legend.position = "none",
                               plot.margin = unit(c(-0.8,-2,0.1,-6), "lines"),
                               axis.text.x = element_text(size=axis.text),
                               axis.text.y = element_text(size=axis.text,color="white"),
                               axis.title.x = element_text(size=title.text),
                               axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                               plot.title = element_text(size=title.text,vjust=-1)
             ),
             p_D$x_500  + ggtitle("F 500m") + 
               theme(legend.position = "none",
                               plot.margin = unit(c(-0.8,-2,0.1,-10), "lines"),
                               axis.text.x = element_text(size=axis.text),
                               axis.text.y = element_text(size=axis.text,color="white"),
                               axis.title.x = element_text(size=title.text,color="white"),
                               axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                               plot.title = element_text(size=title.text,vjust=-1),
             ),
             CV_marg_plot +
               ggtitle("G") + 
               theme(legend.position = "none",
                     plot.margin = unit(c(-0.8,2,0.1,-4), "lines"),
                     axis.text.x = element_text(size=axis.text),
                     axis.text.y = element_text(size=axis.text),
                     axis.title.x = element_text(size=title.text),
                     axis.title.y = element_text(size=title.text),
                     plot.title = element_text(size=title.text,vjust=-1,hjust=-0.2)),
             
             widths=c(1,1,1,1),
             layout_matrix = lay)
dev.off()                 

#### Make Grid of CV Predictions 
axis.text <- 9
title.text <-11
quartz(file="Hake_maps_CV_depth_manual_facet.jpeg",height=9,width=8,dpi=600,type="jpeg")
grid.arrange(
        survey.map + xlab("") + ggtitle("Survey")+
          #geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
          #               linetype="dashed")+
          #geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
            theme(plot.margin = unit(c(0.1,-10,-0.8,-2), "lines"),
              axis.text.x = element_text(size=axis.text),
              axis.text.y = element_text(size=axis.text),
              axis.title=element_text(size=title.text),
              plot.title = element_text(size=title.text,vjust=-1)),
        p_D_CV$x_0 + theme(legend.position = "none",
                              plot.margin = unit(c(0.1,-8,-0.8,-4), "lines"),
                              axis.text.x = element_text(size=axis.text,color="white"),
                              axis.text.y = element_text(size=axis.text,color="white"),
                              axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                              axis.title.x = element_text(size=title.text,color="white"),
                              plot.title = element_text(size=title.text,vjust=-1)),
            p_D_CV$x_50 + theme(legend.position = "none",
                 plot.margin = unit(c(0.1,-6,-0.8,-6), "lines"),
                 axis.text.x = element_text(size=axis.text,color="white"),
                 axis.text.y = element_text(size=axis.text,color="white"),
                 axis.title.x = element_text(size=title.text,color="white"),
                 axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                 plot.title = element_text(size=title.text,vjust=-1),
              ),
            p_D_CV$x_100 + theme(legend.position = "none",
                  plot.margin = unit(c(0.1,-4,-0.8,-8), "lines"),
                  axis.text.x = element_text(size=axis.text,color="white"),
                  axis.text.y = element_text(size=axis.text,color="white"),
                  axis.title.x = element_text(size=title.text,color="white"),
                  axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                  plot.title = element_text(size=title.text,vjust=-1),
              ),
            get_legend(p_D_CV$x_0 + theme(legend.margin = margin(0,-60,0,60),
                                   legend.title = element_text(size=title.text))),
        
            p_D_CV$x_150 + theme(legend.position = "none",
                  plot.margin = unit(c(-0.8,-8,0.1,-4), "lines"),
                  axis.text.x = element_text(size=axis.text),
                  axis.text.y = element_text(size=axis.text),
                  axis.title.y = element_text(size=title.text),
                  axis.title.x = element_text(size=title.text,color="white"),
                  plot.title = element_text(size=title.text,vjust=-1)
              ),
            p_D_CV$x_300 + theme(legend.position = "none",
                  plot.margin = unit(c(-0.8,-6,0.1,-6), "lines"),
                  axis.text.x = element_text(size=axis.text),
                  axis.text.y = element_text(size=axis.text,color="white"),
                  axis.title.x = element_text(size=title.text),
                  axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                  plot.title = element_text(size=title.text,vjust=-1)
              ),
            p_D_CV$x_500 + theme(legend.position = "none",
                  plot.margin = unit(c(-0.8,-4,0.1,-8), "lines"),
                  axis.text.x = element_text(size=axis.text),
                  axis.text.y = element_text(size=axis.text,color="white"),
                  axis.title.x = element_text(size=title.text,color="white"),
                  axis.title.y = element_text(size=title.text,color="white",hjust=-2),
                  plot.title = element_text(size=title.text,vjust=-1)),
widths=c(1,1,1,1),
layout_matrix = lay)
dev.off()                 

###############################################################
###############################################################
###############################################################

setwd(plot.dir)
quartz(file="Hake maps combined SD.jpeg",height=7.5,width=9,dpi=600,type="jpeg")
grid.arrange(survey.map + xlab("") +
               geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                            linetype="dashed")+
               geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
               theme(plot.margin = unit(c(0.1,-6,0.1,-2), "lines"),
                     axis.text.x = element_text(size=y.axis.text),
                     axis.text.y = element_text(size=y.axis.text),
                     axis.title=element_text(size=y.axis.text+2)),
             
             Output.summary.qpcr$p_DNA_lat_base+
               geom_segment(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                            linetype="dashed")+
               geom_text(data=Output.summary.qpcr$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
               ylab("") +
               theme( legend.position = c(1.07, .55),
                      legend.justification = c("right", "top"),
                      legend.box.just = "right",
                      plot.margin = unit(c(0.1,-4,0.1,-4), "lines"),
                      axis.text.x = element_text(size=y.axis.text),
                      axis.text.y = element_blank(),
                      axis.title=element_text(size=y.axis.text+2)),
             
             Acoustic.dat.figs$p_Acoustics_lat_base+
               geom_segment(data=Acoustic.dat.figs$lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                            linetype="dashed")+
               geom_text(data=Acoustic.dat.figs$lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=size.regions) +
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

quartz(file="Hake_relative_index.jpeg",height=4,width=6,dpi=600,type="jpeg")
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


####################################################################3
####################################################################
####################################################################
####################################################################
### Make some plots akin to Hovmuller plots
####################################################################
####################################################################
####################################################################
####################################################################

dat.project <- Output.summary.qpcr$smooth.projections$D_pred_combined %>% mutate(CV=SD/Mean)

# break into XXX degree latitudinal bins 
XXX = 0.25
lat.breaks <- seq(floor(min(dat.project$lat)),ceiling(max(dat.project$lat)),by=XXX)

lat.breaks.plot <- lat.breaks %>% as.data.frame() %>%
                      mutate(lat.group = lat.breaks + XXX/2) %>% filter(lat.group < max(lat.breaks)) %>%
                      mutate(id= 1:n())
# trim to provide just areas that are near 500m.
dat.500 <- dat.project %>% filter(bottom.depth.consensus < 750,bottom.depth.consensus > 400) %>% 
  mutate(lat.bin = cut(lat,breaks=lat.breaks))
lat.lab <- data.frame(lat.bin=levels(dat.500$lat.bin)) %>% as.data.frame() %>% mutate(id= 1:n()) %>% 
  left_join(.,lat.breaks.plot %>% dplyr::select(-"."))
dat.500 <- left_join(dat.500,lat.lab) %>%
  group_by(depth_cat_factor,lat.group) %>% 
  summarise(grand.mean=mean(Mean),grand.sd = sd(Mean),grand.cv=mean(CV))

# trim to provide just areas that are near 300m.
dat.300 <- dat.project %>% filter(bottom.depth.consensus < 400,bottom.depth.consensus > 225) %>% 
              mutate(lat.bin = cut(lat,breaks=lat.breaks))
lat.lab <- data.frame(lat.bin=levels(dat.300$lat.bin)) %>% as.data.frame() %>% mutate(id= 1:n()) %>% 
                left_join(.,lat.breaks.plot %>% dplyr::select(-"."),grand.cv=grand.sd/grand.mean)
dat.300 <- left_join(dat.300,lat.lab) %>%
              group_by(depth_cat_factor,lat.group) %>% 
              summarise(grand.mean=mean(Mean),grand.sd = sd(Mean),grand.cv=mean(CV))
# trim to provide just areas that are near 300m.
dat.150 <- dat.project %>% filter(bottom.depth.consensus < 225,bottom.depth.consensus > 125) %>% 
                  mutate(lat.bin = cut(lat,breaks=lat.breaks))
lat.lab <- data.frame(lat.bin=levels(dat.150$lat.bin)) %>% as.data.frame() %>% mutate(id= 1:n()) %>% 
                  left_join(.,lat.breaks.plot %>% dplyr::select(-"."))
dat.150 <- left_join(dat.150,lat.lab) %>%
              group_by(depth_cat_factor,lat.group) %>% 
              summarise(grand.mean=mean(Mean),grand.sd = sd(Mean),grand.cv=mean(CV))

dat.150$depth_cat_factor = factor(dat.150$depth_cat_factor,levels=c("500", "300","150","100","50","0"))
dat.300$depth_cat_factor = factor(dat.300$depth_cat_factor,levels=c("500","300","150","100","50","0"))
dat.500$depth_cat_factor = factor(dat.500$depth_cat_factor,levels=c("500","300","150","100","50","0"))

ggplot(dat.300)+
  geom_point(aes(x=lat.group,y=grand.mean))+
  facet_wrap(~depth_cat_factor)

H_150 <- ggplot(dat.150)+
  geom_tile(aes(x=lat.group,y=as.factor(depth_cat_factor),fill=grand.mean)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw()
H_150_CV <- ggplot(dat.150)+
  geom_tile(aes(x=lat.group,y=as.factor(depth_cat_factor),fill=grand.cv)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw()

H_300 <- ggplot(dat.300)+
  geom_tile(aes(x=lat.group,y=as.factor(depth_cat_factor),fill=grand.mean)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw()
H_300_CV <- ggplot(dat.300)+
  geom_tile(aes(x=lat.group,y=as.factor(depth_cat_factor),fill=grand.cv)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw()

H_500 <- ggplot(dat.500)+
  geom_tile(aes(x=lat.group,y=as.factor(depth_cat_factor),fill=grand.mean)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw()
H_500_CV <- ggplot(dat.500)+
  geom_tile(aes(x=lat.group,y=as.factor(depth_cat_factor),fill=grand.cv)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw()


setwd(plot.dir)
quartz(file="Hake DNA Hovmuller 150.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
  print(H_150)
dev.off()

setwd(plot.dir)
quartz(file="Hake DNA Hovmuller 300.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
print(H_300)
dev.off()

setwd(plot.dir)
quartz(file="Hake DNA Hovmuller 500.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
print(H_500)
dev.off()

setwd(plot.dir)
quartz(file="Hake DNA Hovmuller_CV 150.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
print(H_150_CV)
dev.off()

setwd(plot.dir)
quartz(file="Hake DNA Hovmuller_CV 300.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
print(H_300_CV)
dev.off()

setwd(plot.dir)
quartz(file="Hake DNA Hovmuller_CV 500.jpeg",height=4.5,width=5,dpi=600,type="jpeg")
print(H_500_CV)
dev.off()


####
## Make Distributional Plots
####
 
cdf_dist_plot <- ggplot() +
    geom_line(data=Acoustic.dat.figs$D_acoustic_uncond_cum_sum, 
              aes(x=lat,y=Mean,color="Acoustic",linetype="Acoustic"))+
    geom_ribbon(data=Acoustic.dat.figs$D_acoustic_uncond_cum_sum, 
                  aes(x=lat,ymin=Q.0.05,ymax=Q.0.95),alpha=0.3,fill=viridis(1,begin=0.0)) +
    geom_line(data=Output.summary.qpcr$D_DNA_cum_sum, aes(x=lat,y=Mean,color="eDNA",linetype="eDNA")) +
    geom_ribbon(data=Output.summary.qpcr$D_DNA_cum_sum, 
            aes(x=lat,ymin=Q.0.05,ymax=Q.0.95),alpha=0.3,fill=viridis(1,begin=0.4)) + 
    scale_y_continuous(expand=c(0.01,0.01)) +
    ylab("Cumulative Distribution") +
    xlab("Latitude") +
    scale_color_viridis_d("Survey",end=0.4)+
    scale_linetype_manual("Survey",values=c("solid","dashed"))+
    theme_bw()

both_dist_summary <- 
  data.frame(id="Acoustic",
             x.lon = -125.5,
             min.lat = unlist(Acoustic.dat.figs$D_acoustic_uncond_cum_sum %>% filter(Q.0.05 >= 0.5) %>% slice(which.min(lat)) %>% dplyr::select(lat)),
            max.lat  = unlist(Acoustic.dat.figs$D_acoustic_uncond_cum_sum %>% filter(Q.0.95 <= 0.5) %>% slice(which.max(lat)) %>% dplyr::select(lat)),
            mean.lat = unlist(Acoustic.dat.figs$D_acoustic_uncond_cum_sum %>% filter(Mean >= 0.5) %>% slice(which.min(lat)) %>% dplyr::select(lat))) %>%
  bind_rows(.,
    data.frame(id="eDNA",
            x.lon = -125,
             min.lat = unlist(Output.summary.qpcr$D_DNA_cum_sum %>% filter(Q.0.05 >= 0.5) %>% slice(which.min(lat)) %>% dplyr::select(lat)),
             max.lat  = unlist(Output.summary.qpcr$D_DNA_cum_sum %>% filter(Q.0.95 <= 0.5) %>% slice(which.max(lat)) %>% dplyr::select(lat)),
             mean.lat = unlist(Output.summary.qpcr$D_DNA_cum_sum %>% filter(Mean >= 0.5) %>% slice(which.min(lat)) %>% dplyr::select(lat)))
  )

cog_plot <-
  Acoustic.dat.figs$base_map_trim_proj +
    geom_point(data=both_dist_summary,aes(x=x.lon,y=mean.lat,color=id),size=2) +
    geom_errorbar(data=both_dist_summary,aes(x=x.lon,ymin=min.lat,ymax=max.lat,color=id),width=0,size=0.8)+
    scale_color_viridis_d("Center of Gravity\n   (Median)",end = 0.4)

quartz(file="Hake_compare_distribution.jpeg",height=4,width=8,dpi=600,type="jpeg")
  grid.arrange(
      cdf_dist_plot + theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"),
                            legend.position = c(0.75, .4),
                            legend.title=element_text(size=10)),
      cog_plot + 
        theme(plot.margin = unit(c(0.25,-4,0.25,-12), "lines"),
              legend.position = c(1.2, .5),
              legend.title=element_text(size=10)),
      nrow=1,
      widths=c(1.5,1))
dev.off()  

quartz(file="Hake_compare_distribution_plus_depth.pdf",height=6,width=6,dpi=900,type="pdf")
grid.arrange(
  cdf_dist_plot + labs(subtitle="A") +
                  theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"),
                        legend.position = c(0.75, .4),
                        legend.title=element_text(size=9),
                        plot.subtitle = element_text(vjust=-7,hjust=0.01,size=9)),
  cog_plot + labs(subtitle="B") + scale_x_continuous(breaks=c(-126,-124)) +
    theme(plot.margin = unit(c(0.25,-0,0.25,-6.5), "lines"),
          legend.position = c(1.4, .5),
          legend.title=element_text(size=9),
          plot.subtitle = element_text(vjust=-7,hjust=0.02,size=9)),
  Output.summary.qpcr$marginal_est_D_by_depth + labs(subtitle="C") +
    theme(plot.margin = unit(c(-1,0.25,1,0.4), "lines"),
          legend.position = c(1.2, .5),
          legend.title=element_text(size=9),
          plot.subtitle = element_text(vjust=-7,hjust=0.01,size=9)),
  nrow=2,
  widths=c(1.5,1))
dev.off()  
