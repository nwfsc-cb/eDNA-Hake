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
library(gtools)

results.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Stan Model Fits"
script.dir <- "/Users/ole.shelton/Github/eDNA-Hake/Scripts"
plot.dir   <- "/Users/ole.shelton/Github/eDNA-Hake/Plots and figures"
# load and run in the acoustic data.
setwd(script.dir)
#source("process acoustic data.R")

# get bathymetry data.
# source("pull NOAA bathy for acoustic data.R")

# Read in the output from the Stan model

setwd(results.dir)

# CHANGE THIS FOR switching between species.
SPECIES <- "hake" # eulachon, hake 

#load(paste("qPCR 2019",SPECIES, MOD, "7_12 Fitted.RData"))
load("Acoustics 2019 lat.long.smooth 6_14_6_10_smooth_hurdle Base_Var Fitted.RData")

#save(Output.qpcr,file=paste("qPCR 2019",SPECIES, MOD, "Fitted.RData"))

# Read in agreed upon dat_raster_fin that has been trimmed to 
dat_raster_fin <- readRDS(file="../Data/_projection_rds/dat_raster_fin.rds") 

# Read in a few breaks based on latitude for summarizing the projections on different spatial scales.
load(file="../Data/lat_breaks_for_projections.RData")

### Cacluate posterior summaries and predictive surfaces
setwd(script.dir)
source("summarize_stan_output_acoustics.R")

#### MAKE OBSERRVED v. PREDICTION Plots
pred_obs_bin_p1 <- ggplot(pred_obs_bin) + 
                   geom_point(aes(x=Mean,y=bin_weight_dens),
                              position=position_jitter(width=0,height=0.1),alpha=0.5) +
                  geom_smooth(aes(x=Mean,y=bin_weight_dens),
                              method = "glm", 
                              method.args = list(family = "binomial")) +
                              #formula = y ~splines::ns(x,6)) +
                  geom_abline(intercept=0,slope=1,color="red") +
                  theme_bw()

pred_obs_pos_p1 <- ggplot(pred_obs_pos) + 
                      geom_point(aes(x=Mean,y=weight_dens_mt_km2),alpha=0.5) +
                    geom_smooth(aes(x=Mean,y=weight_dens_mt_km2)) +
                    geom_abline(intercept=0,slope=1,color="red") +
                    theme_bw()

pred_obs_pos_p2 <- ggplot(pred_obs_pos) + 
                  geom_point(aes(x=Mean.log,y=log(weight_dens_mt_km2)),alpha=0.5) +
                  geom_smooth(aes(x=Mean.log,y=log(weight_dens_mt_km2))) +
                  geom_abline(intercept=0,slope=1,color="red") +
                  theme_bw()

pred_obs_pos_resid <- ggplot(pred_obs_pos) + 
                      geom_point(aes(x=Mean,y=resid),alpha=0.5) +
  geom_smooth(aes(x=log(Mean),y=log(weight_dens_mt_km2))) +
  #geom_abline(intercept=0,slope=1,color="red") +
  theme_bw()

#### MAP MAKING COMMANDS
# Plot Limits for all samples (see Base_Map.R)
# lat.lims <- c(min(dat.id$lat,na.rm=T),max(dat.id$lat,na.rm=T))
# lon.lims <- c(min(dat.id$lon,na.rm=T),max(dat.id$lon,na.rm=T))

#lat.lims.trim <- c(37.5,max(dat.id$lat,na.rm=T))
# Call Base_map.R
setwd(script.dir)
source("Base_map.R")


## Spatial Residuals

resid_pos_map <- base_map_trim + 
            geom_point(data=pred_obs_pos,aes(x=lon,y=lat,fill=resid,color=resid),alpha=0.5,size=4) +
            scale_color_gradient2(low="red",mid=grey(0.8),high="blue") + 
            scale_fill_gradient2(low="red",mid=grey(0.8),high="blue") +
            theme_bw()
                         
##########################################################
## Variance parameters.
##########################################################

kappa_plot <- ggplot() +
  geom_density(aes(pars$sigma),alpha=0.5,fill=viridis(1))+
  geom_vline(xintercept = mean(pars$sigma),linetype="dashed") +
  xlab(expression(kappa)) +
  theme_bw()

#######################################################
#######################################################
####### ------- Make Maps
#######################################################
#######################################################
# smooth.projections %>% names()

# 

# Raw Data projections
p_bin <- base_map_trim + 
  geom_point(data=dat.acoustic.bin,aes(x=lon,y=lat,color=as.factor(bin_weight_dens),
                                       fill=as.factor(bin_weight_dens),shape=as.factor(bin_weight_dens)),alpha=0.5,size=0.5) +
  scale_shape_manual("Occurrence",values=c(4,21)) +
  scale_color_viridis_d("Occurrence",begin=0,end=0.8) +
  scale_fill_viridis_d("Occurrence",begin=0,end=0.8) 
p_bin


dat.acoustic.pos$weight_dens_mt_km2_plot <- dat.acoustic.pos$weight_dens_mt_km2
dat.acoustic.pos$weight_dens_mt_km2_plot[dat.acoustic.pos$weight_dens_mt_km2_plot<1] <- 1

LAB.pos <- c(1,10,50,100,200,400,800)
p_pos <- base_map_trim + 
          geom_point(data=dat.acoustic.pos,aes(x=lon,y=lat,size=weight_dens_mt_km2_plot),alpha=0.9,shape=21,fill=NA,color="red") +
          scale_size(name=expression("Density (mt km"^-2*")"),limits = c(1,900), breaks=LAB.pos,range=c(0.01,5)) #trans="sqrt") 
p_pos

######################################################3
# Smooth projections.
######################################################3
SIZE = 1.3
STROKE = 0 

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

z.lim = c(log10(20),log10(2500)) 
z.lim.bin = c(0,1)
z.breaks <- log10(c(20,100,250,500,1000,2500))
z.bin.breaks <- seq(0,1,length.out = 6)

#DEPTH <- levels(STATION.DEPTH$depth_cat_factor)
p_log_D <- list()

#D_pred_background <- D_pred_log_combined %>% filter(depth_cat_factor == 0)

# Color scales for viridis
OPT = "plasma" # options are "viridis"(default), "magma", "plasma", "inferno"

  p_log_D[[as.name("binom")]] <-
    base_map_trim_proj +
    geom_point(data= D_pred_bin_combined ,
            aes(x=lon,y=lat,color=Mean,fill=Mean),alpha=1,size=SIZE,stroke=STROKE,shape=15) +
    scale_color_viridis_c(name="Probability of\n Occurrence",option=OPT,limits=z.lim.bin,breaks=z.bin.breaks) +
    scale_fill_viridis_c(name="Probability of\n Occurrence",option=OPT,trans="sqrt",limits=z.lim.bin,breaks=z.bin.breaks) +
    
    theme_bw() 
  p_log_D$binom
    
  sqrt_breaks <- c(0,10,50,100,seq(200,900,by=200))
  p_log_D[[as.name("pos")]] <- 
    base_map_trim_proj +
    geom_point(data= D_pred_pos_combined,
               aes(x=lon,y=lat,color=Mean,fill=Mean),alpha=1,size=SIZE,stroke=STROKE,shape=15) +
    scale_color_viridis_c(name=expression("Conditional\nDensity (mt km"^-2*")"),limits=c(0,900),option=OPT,trans="sqrt",breaks=sqrt_breaks,na.value=viridis(1,begin=1,end=1)) +
    scale_fill_viridis_c(name=expression("Conditional\nDensity (mt km"^-2*")"),limits=c(0,900),option=OPT,trans="sqrt",breaks=sqrt_breaks,na.value=viridis(1,begin=1,end=1)) +
    theme_bw() 
  p_log_D$pos
  
  sqrt_breaks <- c(0,100,500,seq(1000,5000,by=1000))
  z.lim = c(0,5200) 
  p_log_D[[as.name("uncond")]] <- 
    base_map_trim_proj +
    geom_point(data= D_pred_uncond_mt_combined ,
               aes(x=lon,y=lat,color=Mean,fill=Mean),alpha=1,size=SIZE,stroke=STROKE,shape=15) +
    scale_color_viridis_c(name=expression("Biomass (mt)"),option=OPT,trans="sqrt",breaks=sqrt_breaks,limits=z.lim) +
    scale_fill_viridis_c(name=expression("Biomass (mt)"),option=OPT,trans="sqrt",breaks=sqrt_breaks,limits=z.lim) +
    theme_bw() 
  
  #########
  
  p_Acoustics_lat_equal <- p_log_D$uncond +
          geom_segment(data=lat.breaks$lats.equal,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                       linetype="dashed")+
          geom_text(data=lat.breaks$lats.equal,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2)

  p_Acoustics_lat_1.0 <- p_log_D$uncond +
          geom_segment(data=lat.breaks$lats.rounded.1.0,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                       linetype="dashed")+
          geom_text(data=lat.breaks$lats.rounded.1.0,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2)

  p_Acoustics_lat_0.5 <- p_log_D$uncond +
          geom_segment(data=lat.breaks$lats.rounded.0.5,aes(x=lon.min,xend=lon.max,y=lat,yend=lat),
                       linetype="dashed")+
          geom_text(data=lat.breaks$lats.rounded.0.5,aes(x=lon.lab,y=mid.lat,label=ID),nudge_x = -0.1,size=2)
  
  
  
  
  ######################################################3
  ## SD in Space
  ######################################################3
  # SIZE = 1.3
  # STROKE = 0 
  # 
  # z.lim = c(0,1.25) 
  # z.breaks <- seq(z.lim[1],z.lim[2],length.out=5)
  # 
  # DEPTH <- levels(STATION.DEPTH$depth_cat_factor)
  # p_log_D_SD <- list()
  # 
  # D_pred_background <- D_pred_log_combined %>% filter(depth_cat_factor == 0)
  # 
  # for(i in 1:length(DEPTH)){
  #   p_log_D_SD[[as.name(paste0("x_",DEPTH[i]))]] <-
  #     base_map_trim +
  #     geom_point(data= D_pred_background,
  #                aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
  #     geom_point(data= D_pred_log_combined %>% filter(depth_cat_factor == DEPTH[i]),
  #                aes(x=lon,y=lat,color=SD),alpha=1,size=SIZE,stroke=STROKE) +
  #     geom_point(data= D_pred_log_combined %>% filter(SD > z.lim[2]),
  #                aes(x=lon,y=lat),alpha=1,size=SIZE,stroke=STROKE,color=viridis(1,begin=0.999,end=1)) +
  #     scale_color_viridis_c(option=OPT,limits=z.lim,breaks=z.breaks,labels=round(z.breaks,2),name=expression("SD DNA Copies L"^-1)) +
  #     ggtitle(paste(DEPTH[i],"m")) +
  #     theme_bw() 
  #   
  #   p_log_D_SD[[as.name(paste0("x_",DEPTH[i],"_dots"))]] <- 
  #     p_log_D_SD[[as.name(paste0("x_",DEPTH[i]))]] +
  #     geom_point(data=STATION.DEPTH %>% filter(depth_cat_factor == DEPTH[i]),aes(x=lon,y=lat),col="red")
  # }  
  # 
  # SIZE = 1.25
  # p_log_D_SD_facet <- base_map_trim +
  #   # geom_point(data= D_pred_background,
  #   #            aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
  #   geom_point(data= D_pred_log_combined ,
  #              aes(x=lon,y=lat,color=SD),alpha=1,size=SIZE,stroke=STROKE) +
  #   geom_point(data= D_pred_log_combined %>% filter(SD > z.lim[2]),
  #              aes(x=lon,y=lat),alpha=1,size=SIZE,stroke=STROKE,color=viridis(1,begin=0.999,end=1)) +
  #   scale_color_viridis_c(option=OPT,limits=z.lim,breaks=z.breaks,labels=round(z.breaks,2),name=expression("SD DNA Copies L"^-1)) +
  #   facet_wrap(~depth_cat_factor) +
  #   theme_bw() 


  
  p_marginal_smooth <- list()
  if(MODEL.TYPE=="lat.long.smooth"){
      p_marginal_smooth[[as.name("bin_bottom_depth")]] <- ggplot(D_pred_bottom_depth_bin_combined) +
        geom_line(aes(x=bathy.bottom.depth,y=Mean),alpha=0.8) +
        geom_ribbon(aes(x=bathy.bottom.depth,ymin=Val.95.,ymax=Val.5.),alpha=0.4) +
        theme_bw()
      p_marginal_smooth[[as.name("pos_bottom_depth")]] <- ggplot(D_pred_bottom_depth_pos_combined) +
        geom_line(aes(x=bathy.bottom.depth,y=Mean),alpha=0.8) +
        geom_ribbon(aes(x=bathy.bottom.depth,ymin=Val.95.,ymax=Val.5.),alpha=0.4) +
        theme_bw()
  }
  
 
  
  #################### 
  ## Smoothes Only
  ###################
   D_pred_smooth_pos_combined <-   smooth.projections$D_pred_smooth_pos_combined
  
   z.lim = c(min(D_pred_smooth_combined),max(D_pred_smooth_combined))
   z.breaks <- log10(c(20,100,250,500,1000,2500))
   #z.lim.labs <- 10^z.breaks
   
   
   SIZE = 1.3
   p_log_D_smoothes_only <- base_map_trim +
     # geom_point(data= D_pred_background,
     #            aes(x=lon,y=lat),size=SIZE,color=grey(0.6),alpha=1,stroke=STROKE) +
     geom_point(data= D_pred_smooth_pos_combined ,
                aes(x=lon,y=lat,color=Mean),alpha=0.75,size=SIZE,stroke=STROKE) +
     geom_point(data= D_pred_smooth_pos_combined %>% filter(Mean < z.lim[1]),
                aes(x=lon,y=lat),alpha=0.75,size=SIZE,stroke=STROKE,color=viridis(1,begin=0,end=0.001)) +
     geom_point(data= D_pred_smooth_pos_combined %>% filter(Mean > z.lim[2]),
                aes(x=lon,y=lat),alpha=0.75,size=SIZE,stroke=STROKE,color=viridis(1,begin=0.999,end=1)) +
     scale_color_viridis_c(option=OPT,
                           # limits=z.lim,
                           # breaks=z.breaks,
                           # labels=z.lim.labs,
                           name=expression("Biomass")) +
     #facet_wrap(~depth_cat_factor) +
     theme_bw() 


######################################################3
######################################################3

###### Plots for files

# Spatial Smoothes
setwd(plot.dir)
pdf(file=paste("Hake transect p_log_D Acoustics",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=6,width=5)
  print(p_bin)
  print(p_pos)
  print(p_log_D)
dev.off()

# Spatial smoothes with spatial strata
setwd(plot.dir)
quartz(file=paste("Hake Acoustics p_log_D lat_equal",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),height=6,width=4,dpi=600,type="pdf")
  print(p_Acoustics_lat_equal  )
dev.off()

quartz(file=paste("Hake Acoustics p_log_D lat_1.0_degree",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),height=6,width=4,dpi=600,type="pdf")
  print(p_Acoustics_lat_1.0  )
dev.off()

quartz(file=paste("Hake Acoustics p_log_D lat_0.5_degree",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),height=6,width=4,dpi=600,type="pdf")
  print(p_Acoustics_lat_0.5  )
dev.off()



# pdf(file=paste("Hake transect p_log_D_facet2",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf"),onefile = T,height=6,width=11)
#   print(p_log_D_facet2)
# dev.off()

#
pdf(file=paste("Acoustics_marginal_smoothes_",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf",sep=""),onefile = T,height=5,width=7)
  print(p_marginal_smooth)
dev.off()


pdf(file=paste("Acoustics_kappa_plot_",MODEL.TYPE,"_",MODEL.VAR,"_",MODEL.ID,".pdf",sep=""),onefile = T,height=3.5,width=5)
print(kappa_plot)
dev.off()


# p_log_D
# p_log_D_facet
# p_log_D_SD
# p_log_D_SD_facet
# 
# p_log_D_smoothes_only_facet

# Latitudinal slices.
# SP.transect.plots

# Factors for each depth
# depth_cat_plot

# Variance parameters
# sd_param_plots

# Marginal Smoothes
p_marginal_smooth

# Residuals
p_resid_facet

# Contamination 

##############################################3
# Save output (plots, data.frames) to file for comparison with the eDNA data.

Acoustic.dat.figs <- list(
  # Figures (Maps)
  p_bin = p_bin,
  p_pos = p_pos,
  p_log_D = p_log_D,
  p_marginal_smooth = p_marginal_smooth,

    # Figures (pred-obs)
  pred_obs_bin_p1 = pred_obs_bin_p1,
  pred_obs_pos_p1 = pred_obs_pos_p1,
  pred_obs_pos_p2 = pred_obs_pos_p2,
  resid_pos_map = resid_pos_map,
  
  # variance
  kappa_plot = kappa_plot,
  
  # Derived data files
  smooth.projections = smooth.projections,
  stanMod_summary_parts = Output.acoustic$stanMod_summary_parts,

  # Covariates, etc needed to make spatial predictions.
  smooth.dat.pred.bin = smooth.dat.pred.bin,
  smooth.dat.pred.pos = smooth.dat.pred.pos,
  new_data_trim = new_data_trim,
  
  # Raw Data
  dat.acoustic.bin = dat.acoustic.bin, 
  
  # Derived MCMC predictions.
  D_pred_bin = D_pred_bin,
  D_pred_pos = D_pred_pos,
  D_pred_uncond_mt_combined = D_pred_uncond_mt_combined,
  D_acoustic_uncond_cum_sum = D_acoustic_uncond_cum_sum,
  D_acoustic_uncond_total_mt = D_acoustic_uncond_total_mt,
  
  D_1.0_uncond_resample = D_1.0_uncond_resample, 
  D_0.5_uncond_resample = D_0.5_uncond_resample,
  D_grid.cell_uncond_resample = D_grid.cell_uncond_resample,
  
  # Map Plots
  p_Acoustics_lat_0.5 = p_Acoustics_lat_0.5,
  p_Acoustics_lat_equal = p_Acoustics_lat_equal,
  p_Acoustics_lat_base = p_log_D$uncond,
  lat.breaks = lat.breaks,
  base_map_trim_proj = base_map_trim_proj
)

setwd(results.dir)
save(Acoustic.dat.figs,file=paste("./_Summarized_Output/Acoustics 2019",MODEL.TYPE,MODEL.ID,MODEL.VAR,"Derived Q.RData"))
