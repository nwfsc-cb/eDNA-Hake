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
library(rgdal)
library(sp)
library(brms)
library(loo)
## 

# Working directories
base.dir <- "/Users/ole.shelton/Github/eDNA-Hake/"
data.dir <- paste0(base.dir,"Data/acoustics 2019")
script.dir <- paste0(base.dir,"/Scripts")
plot.dir <- paste0(base.dir,"Plots and figures")

###########################################################################
# Make a choice about the kind of model to run.
# Options are "Base", "lat.long.smooth", "lat.long.smooth.base"
MODEL.TYPE = "lat.long.smooth.base"
###########################################################################
DATA.INCLUDED = "raw" # options are "raw" or "combined"
###########################################################################
# identifier
MODEL.ID <- "10_20_10_20_smooth_hurdle"
###########################################################################
# variance scenario # options are "Base_Var", "Linear_Var"
MODEL.VAR <- "Base_Var" #
###########################################################################
#set.seed(111)
# Construct smoothes for each 
# define knots.
N.knots.lon.pos  <- 10
N.knots.lat.pos  <- 20
N.knots.lon.bin  <- 10
N.knots.lat.bin  <- 20
N.knots.bd <- 5

# load and run the acoustic data. This script is called by the qPCR script
setwd(script.dir)
source("process acoustic data for qPCR.R",local=T)

###################################################################
### HERE ADD SCRIPT FOR PULLING DEPTH FROM bathy
setwd(script.dir)
source("pull NOAA bathy for acoustic data.R",local=T)

temp <- marmap::get.depth(b,
                          x= dat.acoustic %>% dplyr::select(lon,lat),
                          locator = FALSE)  
dat.acoustic <- dat.acoustic %>% mutate(bathy.bottom.depth =  -temp$depth,
                                        bin_weight_dens = ifelse(weight_dens_kg_nm2>0,1,0),
                                        weight_dens_mt_nm2 = weight_dens_kg_nm2 / 1000)


### Subset the data to include only transects that were sampled for eDNA in 2019.
dat.acoustic <- dat.acoustic %>% filter(transect >=22,transect <=85)

### Convert half-nm increments into single mile increments.
dat.acoustic.N <- dat.acoustic %>% arrange(transect,VL_start) %>% group_by(transect) %>%
                summarise(N=length(transect),max.ID=max(ID))

all <- NULL
for(i in 1:nrow(dat.acoustic.N)){
  trans <- dat.acoustic.N$transect[i]
  x <- 1:dat.acoustic.N$max.ID[i]
  y <- 1:floor(dat.acoustic.N$max.ID[i]/2)
  y <- sort(c(y,y))
  if(length(y) < dat.acoustic.N$max.ID[i]){
    y = c(y,max(y)+1)
  }
  temp <- bind_cols(transect=rep(trans,dat.acoustic.N$max.ID[i]),
                    ID=x,new.ID=y)
  
  all <- bind_rows(all,temp) 
}

dat.acoustic <- left_join(dat.acoustic,all)

dat.acoustic.trim <- dat.acoustic %>% group_by(transect,new.ID) %>%
                      summarize(utm.lon = mean(utm.lon),
                                utm.lat = mean(utm.lat),
                                bathy.bottom.depth = mean(bathy.bottom.depth),
                                weight_dens_mt_nm2_new = mean(weight_dens_mt_nm2),
                                diff1 = max(weight_dens_mt_nm2)-weight_dens_mt_nm2_new,
                                diff2 = weight_dens_mt_nm2_new - min(weight_dens_mt_nm2),
                                sig_b = mean(sig_b),
                                mean.depth.m = mean(mean.depth.m),
                                N=length(new.ID)) %>%
                      mutate(bin_weight_dens = ifelse(weight_dens_mt_nm2_new>0,1,0)) %>%
                      rename(weight_dens_mt_nm2 = weight_dens_mt_nm2_new)
  
ggplot(dat.acoustic.trim %>% filter(N>1)) +
    geom_point(aes(y = diff1/weight_dens_mt_nm2,x= weight_dens_mt_nm2)) +
    geom_abline(intercept = 0, slope=1) +
    geom_abline(intercept = 0, slope=0.5,color="red") 

ggplot(dat.acoustic.trim%>% filter(N>1)) +
  geom_histogram(aes(diff1/weight_dens_mt_nm2)) 

dat.acoustic.trim <- dat.acoustic.trim %>% filter(N>1) %>% mutate(cv.coarse = diff1/weight_dens_mt_nm2)

summary(dat.acoustic.trim)

# Sort data so that all of the positive values are first
dat.acoustic <- dat.acoustic %>% arrange(desc(bin_weight_dens))
dat.acoustic.trim <- dat.acoustic.trim %>% arrange(desc(bin_weight_dens))

if(DATA.INCLUDED == "raw"){
  dat.acoustic.bin  <- dat.acoustic
}else if(DATA.INCLUDED == "combined"){
  dat.acoustic.bin  <- dat.acoustic.trim
}
N_obs_bin <- nrow(dat.acoustic.bin)
N_obs_pos <- sum(dat.acoustic.bin$bin_weight_dens)

dat.acoustic.pos  <- dat.acoustic.bin[1:N_obs_pos,]





#ggplot(dat.acoustic) + geom_point(aes(x=lon,y=lat,color=bathy.bottom.depth),alpha=0.25) +theme_bw()
###################################################################

### Go get the projected points for the coast from Blake's data
dat_raster_fin <- readRDS(file="../Data/_projection_rds/dat_raster_fin.rds")   

# Information from Chu on detection thresholds.

#  Email string from April 12,2021.
# Chu says:
#   In exporting the NASC from echoview, we set a threshold of Sv = -70 dB (RT - correct me if I am wrong), 
#   the resultant minimum NASC_min = 4*pi*1852^2*H*10^(-70/10) = 43.1 m^2/nmi^2 with H=10 m. 
#   here is an approximate estimate of minimum value of biomass density:
#   1) for adult hake, say 40-cm long: TS = 35.96 dB, sA = 4*pi*10^(TS/10) = 0.032  m^2
#   2) the biomass of the average hake of 40-cm is  wgt_sngle_hake40 = 0.45 kg
#   3) the minimum aerial biomass density =  NASC_min/sA*wgt_single_hake40 = 606 kg/nmi^2
#   
#   So the value of 606 kg/nmi^2 is the roughly estimated threshold of aerial biomass density.
#   
#   Rebecca.Thomas - NOAA Federal
#   Chu-
#     Interesting, I had never thought about it this way before. However, our threshold is -69 dB, not -70.  
#     Close :-)  I guess if you want to get technical, you could run to the calculation again with -69.
# 
#     I see what you are saying Chu, that the smallest number we could produce would be ~606 kg per nautical mile squared. 
#     I think in reality our cut off would be significantly higher than that, but since it would be difficult 
#     to put a number on what is essentially a subjective process, your number is probably the best one we can produce. 
#     At the scales you're talking about, we are just as likely to have false positives as well as real positives.  
#     Some of the zeros that we have almost certainly have hake in them, but we don't know when they would, of course.
# 
#   RT
# 
#   Hi RT,
# 
#   Thank RT for correcting me on the export threshold, i.e.  -69 dB instead of -17 dB.  
#   Using the same calculation in my previous email, this 1 dB difference will result in 10^0.1 = 1.2589 time difference in biomass estimate, 
#   i.e. the minimum aerial biomass density =  NASC_min/sA*wgt_single_hake40 = 1.2589*606 kg/nmi^2 = 763  kg/nmi^2
#   

#### Here are some data plots..
setwd(script.dir)
source("Base_map.R")
base_map_trim +
    geom_point(data= dat.acoustic %>% filter(weight_dens_mt_nm2>0),
               aes(x=lon,y=lat,
                 color=log(weight_dens_mt_nm2),
                 fill=log(weight_dens_mt_nm2)),alpha=0.5) +
  geom_point(data= dat.acoustic %>% filter(weight_dens_mt_nm2==0),
             aes(x=lon,y=lat),shape="x",alpha=0.2) +
  scale_color_viridis_c()

base_map_trim +
  geom_point(data= dat.acoustic,
             aes(x=lon,y=lat,
                 color=(bin_weight_dens),
                 alpha=0.5))+
  scale_color_viridis_c()



# Plot on log scale
ggplot(dat.acoustic %>% filter(transect<40)) +
  geom_point(aes(x=utm.lon,y=log(weight_dens_mt_nm2)))+
  facet_wrap(~transect,scales="free_x") +
  theme_bw()

ggplot(dat.acoustic.trim %>% filter(transect<40)) +
  geom_point(aes(x=utm.lon,y=log(weight_dens_mt_nm2)))+
  facet_wrap(~transect,scales="free_x") +
  theme_bw()





# Plot on identity scale
ggplot(dat.acoustic %>% filter(transect<40)) +
  geom_point(aes(x=utm.lon,y=weight_dens_mt_nm2))+
  facet_wrap(~transect,scales="free_x") +
  theme_bw()

ggplot(dat.acoustic.trim %>% filter(transect<40)) +
  geom_point(aes(x=utm.lon,y=log(weight_dens_mt_nm2)))+
  facet_wrap(~transect,scales="free_x") +
  theme_bw()

##################################
thresh_mt_nm2 <- 0.763 # This is the threshold determined from emails pasted above.
prob_val      <- 0.99  # This is the probability of observing presence at the threshold.

# cloglog calculations.
phi0 <- log(-log(1-0.99))
phi1 <- 2

x <- seq(-5,5,length.out=1000)
y <- 1- exp(-exp(phi0+phi1*x))
plot(y~x)

# Look at intecept with cloglog
phi_0_mean <- log(-log(1-prob_val))  
phi_0_fix <- phi_0_mean
phi_0_sd   <- 0.1
phi_0_rand <- rnorm(1e6,phi_0_mean,phi_0_sd)
test<-1- exp(-exp(phi_0_rand))
hist(test,breaks=1000)
mean(test)

# priors for phi_1 (normal prior)
phi_1_mean <- 20
phi_1_sd  <- 3

x <- seq(0,3,length.out=1000)
y1 <- 1- exp(-exp(phi0+0.9*x))
y2 <- 1- exp(-exp(phi0+1.5*x))
y3 <- 1- exp(-exp(phi0+2.1*x))

##### LOGIT CALCULATIONS
phi0 <- -log(1/prob_val - 1)
phi1 <- 20
offset = 0.76

x <- seq(0,3,length.out=1000)
y <- 1 / (1+exp(-(phi0+phi1*(x-offset))))

plot(y~x)
abline(h=prob_val,v=offset,col=2)

# Look at intecept with logit
phi_0_mean <- -log(1/prob_val - 1)
phi_0_fix  <- phi_0_mean


###################################################################
###################################################################
###################################################################
if(MODEL.TYPE == "lat.long.smooth"){
  # THis is a new version.  Derive basis function set up from the positive observations, make a second set of 
  # projections to the binomial components.
  
  ## FILL THIS IN TO HAVE A MODEL WITHOUT THE Bathymetry term.
  # LAT-LON TENSOR SMOOTHES
  utm.lon.lims <- c(min(dat.acoustic.bin$utm.lon), max(dat.acoustic.bin$utm.lon))
  utm.lat.lims <- c(min(dat.acoustic.bin$utm.lat), max(dat.acoustic.bin$utm.lat))
  knots.lon.pos    <- seq(utm.lon.lims[1],utm.lon.lims[2],length.out=N.knots.lon.pos)
  knots.lat.pos    <- seq(utm.lat.lims[1],utm.lat.lims[2],length.out=N.knots.lat.pos)
  
  # bottom.depth smooth
  #N.knots.bd <- 6
  bd.lims <- c(min(dat.acoustic.bin$bathy.bottom.depth), max(dat.acoustic.bin$bathy.bottom.depth))
  knots.bd <- seq(bd.lims[1],bd.lims[2],length.out=N.knots.bd)
  
  dat.acoustic.pos.temp <- dat.acoustic.pos %>% ungroup() %>% mutate(Y=rnorm(N_obs_pos,0,1))
  
  brms.object.pos <- brm(Y ~ t2(utm.lon,utm.lat,k=c(N.knots.lon.pos,N.knots.lat.pos),bs="cr") +
                         s(bathy.bottom.depth,k=N.knots.bd),
                         #knots=list(#utm.lon=knots.lon,
                         data=dat.acoustic.pos.temp,
                         chains=0)
  
  smooth.dat.pos <- standata(brms.object.pos)
  code.dat.pos   <- stancode(brms.object.pos)
  
  # Repeat for the bin
  dat.acoustic.bin.temp <- dat.acoustic.bin %>% ungroup() %>% mutate(Y=rnorm(N_obs_bin,0,1))
  brms.object.bin <- brm(Y ~ t2(utm.lon,utm.lat,k=c(N.knots.lon.bin,N.knots.lat.bin),bs="cr") +
                         s(bathy.bottom.depth,k=N.knots.bd),
                         #knots=list(#utm.lon=knots.lon,
                         data=dat.acoustic.bin.temp,
                         chains=0)
  smooth.dat.bin <- standata(brms.object.bin)
  code.dat.bin   <- stancode(brms.object.bin)
  
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
  # THis is a new version.  Derive basis function set up from the positive observations, make a second set of 
  # projections to the binomial components.
  
  ## FILL THIS IN TO HAVE A MODEL WITHOUT THE Bathymetry term.
  # LAT-LON TENSOR SMOOTHES
  utm.lon.lims <- c(min(dat.acoustic.bin$utm.lon), max(dat.acoustic.bin$utm.lon))
  utm.lat.lims <- c(min(dat.acoustic.bin$utm.lat), max(dat.acoustic.bin$utm.lat))
  knots.lon.pos    <- seq(utm.lon.lims[1],utm.lon.lims[2],length.out=N.knots.lon.pos)
  knots.lat.pos    <- seq(utm.lat.lims[1],utm.lat.lims[2],length.out=N.knots.lat.pos)
  
  # bottom.depth smooth
  #N.knots.bd <- 6
  bd.lims <- c(min(dat.acoustic.bin$bathy.bottom.depth), max(dat.acoustic.bin$bathy.bottom.depth))
  knots.bd <- seq(bd.lims[1],bd.lims[2],length.out=N.knots.bd)
  
  dat.acoustic.pos.temp <- dat.acoustic.pos %>% ungroup() %>% mutate(Y=rnorm(N_obs_pos,0,1))
  
  brms.object.pos <- brm(Y ~ t2(utm.lon,utm.lat,k=c(N.knots.lon.pos,N.knots.lat.pos),bs="cr"),
                     #s(bathy.bottom.depth,k=N.knots.bd),
                     #knots=list(#utm.lon=knots.lon,
                     #utm.lat=knots.lat,
                     #      bottom.depth = knots.bd),
                     data=dat.acoustic.pos.temp,
                     chains=0)
  
  smooth.dat.pos <- standata(brms.object.pos)
  code.dat.pos   <- stancode(brms.object.pos)
  
  # Repeat for the bin
  dat.acoustic.bin.temp <- dat.acoustic.bin %>% ungroup() %>% mutate(Y=rnorm(N_obs_bin,0,1))
  brms.object.bin <- brm(Y ~ t2(utm.lon,utm.lat,k=c(N.knots.lon.bin,N.knots.lat.bin),bs="cr"),
                         #s(bathy.bottom.depth,k=N.knots.bd),
                         #knots=list(#utm.lon=knots.lon,
                         #utm.lat=knots.lat,
                         #      bottom.depth = knots.bd),
                         data=dat.acoustic.bin.temp,
                         chains=0)
  smooth.dat.bin <- standata(brms.object.bin)
  code.dat.bin   <- stancode(brms.object.bin)
 
  # THIS IS THE MAGIC SAUCE FOR MAKING PREDICTIONS ON THE SAME SET OF BASIS FUNCTIONS. 
     # dat.acoustic.bin.temp <- dat.acoustic.bin %>% ungroup() %>% mutate(Y=rnorm(N_obs_bin,0,1))
     # smooth.dat.bin <- standata(brms.object.pos,newdata=dat.acoustic.bin.temp)

  # C <- brm(Y ~ s(bottom.depth,k=N.knots.bd),
  #          knots=list(bottom.depth = knots.bd),
  #          data=Q,
  #          chains=0)
  # 
  # brms.object.bin <-
  # smooth.datC <- standata(C)
  # code.datC   <- stancode(C)
}





#################################################################### 
stan_data = list(
  # Observations
  "bin_weight_dens" = dat.acoustic.bin$bin_weight_dens,
  "pos_weight_dens" = dat.acoustic.pos$weight_dens_mt_nm2,
  "N_obs_bin"  = N_obs_bin,
  "N_obs_pos"  = N_obs_pos,

  # Index for leave.one.out
  # "loo_pos_idx" = ,

  # SMOOTH COMPONENTS.
  # Data for linear effects
  # "K" = smooth.dat$K, # number of population-level effects
  # "X" = smooth.dat$X,  # population-level design matrix
  
  # data for splines (positive values)
  "Ks_pos" = smooth.dat.pos$Ks, # number of linear effects
  "Ks_bin" = smooth.dat.bin$Ks,
  "Xs_pos" = smooth.dat.pos$Xs, # design matrix for the linear effects
  "Xs_bin" = smooth.dat.bin$Xs, # design matrix for the linear effects
  
  # data for spline t2(utm.lon,utm.lat,k=c(N.knots.lon,N.knots.lat),bs="cr",by=depth_cat_factor)0
  "nb_1_pos" = smooth.dat.pos$nb_1,  # number of bases
  "knots_1_pos" = smooth.dat.pos$knots_1, # number of knots
  "nb_1_bin" = smooth.dat.bin$nb_1,  # number of bases
  "knots_1_bin" = smooth.dat.bin$knots_1, # number of knots
  # basis function matrices
  "Zs_1_1_pos" = smooth.dat.pos$Zs_1_1,
  "Zs_1_2_pos" = smooth.dat.pos$Zs_1_2,
  "Zs_1_3_pos" = smooth.dat.pos$Zs_1_3,

  "Zs_1_1_bin" = smooth.dat.bin$Zs_1_1,
  "Zs_1_2_bin" = smooth.dat.bin$Zs_1_2,
  "Zs_1_3_bin" = smooth.dat.bin$Zs_1_3,
  
  "thresh_mt_nm2" = thresh_mt_nm2, # This is the offset to make the phi_0 be the parameter 
                             # define the probability of observing
                             # positive kg nm^-2 at bin_offset kg nm^-2
  
  # Priors for cloglog parameters
  "phi_0_fix" = phi_0_fix,
  "phi_0_mean" = phi_0_mean,
  "phi_0_sd"   = phi_0_sd,
  "phi_1_mean" = phi_1_mean,
  "phi_1_sd"   = phi_1_sd
)

if(MODEL.TYPE == "lat.long.smooth"){
  stan_data = c(stan_data,
                list(# data for spline s(bottom.depth,k=N.knots.bd)
                  "nb_2"= smooth.dat$nb_2,  # number of bases
                  "knots_2"= smooth.dat$knots_2,  # number of knots
                  # basis function matrices
                  "Zs_2_1_pos"= smooth.dat.pos$Zs_2_1,
                  "Zs_2_1_bin"= smooth.dat.bin$Zs_2_1))
}


######################### STAN PARS #############################
stan_pars = c(
  #"phi_0",  # logit intercept for standards
  #"phi_1",  # logit slope for standard,
  #"sigma"     # variability among samples, given individual bottle, site, and month 
  #"log_lik"
)   

if(MODEL.TYPE=="lat.long.smooth"){
  stan_pars <- c(#stan_pars,
    # smooth and linear coefficients
    "Intercept_pos",
    "Intercept_bin",
    # "b",
    "bs_pos",
    "bs_bin",
    "s_1_1_pos","s_1_2_pos","s_1_3_pos",
    "s_1_1_bin","s_1_2_bin","s_1_3_bin",
    "s_2_1_pos","s_2_1_bin")
  # "s_3_1","s_3_2","s_3_3",
  # "s_4_1","s_4_2","s_4_3",
  # "s_5_1","s_5_2","s_5_3",
  # "s_6_1","s_6_2","s_6_3"
  #)   
  N_knots_all <- list("N.knots.lon.pos"=N.knots.lon.pos,
                      "N.knots.lat.pos"=N.knots.lat.pos,
                      "N.knots.lon.bin"=N.knots.lon.bin,
                      "N.knots.lat.bin"=N.knots.lat.bin,
                      "N.knots.bd"=N.knots.bd)
}
if(MODEL.TYPE=="lat.long.smooth.base"){
  stan_pars <- c(#stan_pars,
                 # smooth and linear coefficients
                  "Intercept_pos",
                  "Intercept_bin",
                 # "b",
                  "bs_pos",
                  "bs_bin",
                  "s_1_1_pos","s_1_2_pos","s_1_3_pos",
                  "s_1_1_bin","s_1_2_bin","s_1_3_bin")
                 # "s_2_1","s_2_2","s_2_3",
                 # "s_3_1","s_3_2","s_3_3",
                 # "s_4_1","s_4_2","s_4_3",
                 # "s_5_1","s_5_2","s_5_3",
                 # "s_6_1","s_6_2","s_6_3"
  #)   
  N_knots_all <- list("N.knots.lon.pos"=N.knots.lon.pos,
                      "N.knots.lat.pos"=N.knots.lat.pos,
                      "N.knots.lon.bin"=N.knots.lon.bin,
                      "N.knots.lat.bin"=N.knots.lat.bin)
}
if(MODEL.VAR=="Base_Var"){
  stan_pars <- c(stan_pars)
}


### INTIAL VALUES
stan_init_f2 <- function(n.chain){ 
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(

      # phi_0  = rnorm(1,phi_0_mean,phi_0_sd),
      Intercept = rnorm(1,3,0.5),
      # v_0 = runif(1,0.25,1),
      # v_1 = runif(1,-0.25,-0.05),
      # phi_1  = rnorm(1,20,1)
      sigma  = runif(1,0.01,0.1) 
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
Warm = 1000
Iter = 1000
Treedepth = 13
Adapt_delta = 0.80

LOC <- paste0(base.dir,"Scripts/Stan Files/")
setwd(LOC)

if(MODEL.TYPE=="lat.long.smooth"){
  if(MODEL.VAR=="Base_Var"){
  stanMod = stan(file = "acoustic_Hake_smoothes.stan" ,data = stan_data, 
                 verbose = FALSE, chains = N_CHAIN, thin = 1, 
                 warmup = Warm, iter = Warm + Iter, 
                 control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                 pars = stan_pars,
                 boost_lib = NULL,
                 sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,"_Acoustics.csv"),
                 init = stan_init_f2(n.chain=N_CHAIN)
                                     # phi_0_mean = phi_0_mean,
                                     # phi_0_sd   = phi_0_sd,
                                     # phi_1_mean = phi_1_mean,
                                     # phi_1_sd   = phi_1_sd)
  )
  }
}




if(MODEL.TYPE=="lat.long.smooth.base"){
  if(MODEL.VAR=="Base_Var"){
    stanMod = stan(file = "acoustic_Hake_smoothes_base.stan" ,data = stan_data, 
                   verbose = FALSE, chains = N_CHAIN, thin = 1, 
                   warmup = Warm, iter = Warm + Iter, 
                   control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                   pars = stan_pars,
                   boost_lib = NULL,
                   sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,"_Acoustics.csv"),
                   init = stan_init_f2(n.chain=N_CHAIN)
                   # phi_0_mean = phi_0_mean,
                   # phi_0_sd   = phi_0_sd,
                   # phi_1_mean = phi_1_mean,
                   # phi_1_sd   = phi_1_sd)
    )
  }
}






# Make diagnostic plots of fitting.
# Predicted Observed plots
# Sum to 1 nm bins.















# get_adaptation_info(stanMod)
# pars <- rstan::extract(stanMod, permuted = TRUE)
# log_lik_1 <- extract_log_lik(stanMod, merge_chains = FALSE)
# r_eff <- relative_eff(exp(log_lik_1))
# loo_val <- loo(log_lik_1, r_eff = r_eff)
# print(loo_val)
# plot(loo_val)
# 
# bad <- loo_val$pointwise %>% as.data.frame()
# bad$ID <- 1:nrow(bad)
# bad2 <- bad %>% filter(influence_pareto_k>1)
# 
# bad.dat <- dat.samp[bad2$ID,]
# dim(bad.dat)
# 
# ggplot()+
#   geom_histogram(data=bad.dat %>% filter(Ct>0),aes(Ct),fill="red",alpha=0.5) +
#   geom_histogram(data=dat.samp %>% filter(Ct>0),aes(Ct),fill="blue",alpha=0.5) + theme_bw()
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
  "phi_1",
  "sigma"))$summary     # variability among samples, given individual bottle, site, and month ))$summary
if(MODEL.TYPE == "lat.long.smooth"){
  stanMod_summary_parts[[as.name("smoothes")]] <- summary(stanMod, pars=c("bs",
                                                                          "s_1_1","s_1_2","s_1_3",  
                                                                          "s_2_1"
                                                                          ))$summary
  
}
if(MODEL.TYPE == "lat.long.smooth.base"){
  stanMod_summary_parts[[as.name("smoothes")]] <- summary(stanMod, pars=c("b", "bs",
                                                                          "s_1_1","s_1_2","s_1_3"
                                                                          ))$summary
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
  TRACE[[as.name("s_6")]] <- traceplot(stanMod,pars=c("s_6_1","s_6_2","s_6_3"),inc_warmup=FALSE)
  TRACE[[as.name("s_7")]] <- traceplot(stanMod,pars=c("s_7_1"),inc_warmup=FALSE)
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
  TRACE = TRACE,
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
  dat_raster_fin = dat_raster_fin
  #dat_raster_trim = dat_raster_trim
  
)

setwd(base.dir)
setwd("./Stan Model Fits/")
save(Output.qpcr,file=paste("qPCR 2019",SP,MODEL.TYPE,MODEL.ID,MODEL.VAR,"Fitted.RData"))
######




















