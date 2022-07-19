# This script is for reading in acoustic data, merging it with locations,
# estimating a stan model file, and writing results to file.

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
base.dir <- "/Users/ole.shelton/Github/eDNA-Hake-public/"
data.dir <- paste0(base.dir,"Data/acoustics 2019")
script.dir <- paste0(base.dir,"/Scripts")
plot.dir <- paste0(base.dir,"Plots and figures")

###########################################################################
# Make a choice about the kind of model to run.
# Options are 
  # "lat.long.smooth" - a model with  "lat.long.smooth.base"
MODEL.TYPE = "lat.long.smooth"
###########################################################################
DATA.INCLUDED = "raw" # options are "raw" or "combined"
###########################################################################
# identifier
MODEL.ID <- "6_14_6_10_smooth_hurdle"
###########################################################################
# variance scenario # options are "Base_Var", "Linear_Var"
MODEL.VAR <- "Base_Var" #
###########################################################################
#set.seed(111)
# Construct smoothes for each 
# define knots.
N.knots.lon.pos  <- 6
N.knots.lat.pos  <- 14
N.knots.lon.bin  <- 6
N.knots.lat.bin  <- 10
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
                                        weight_dens_mt_nm2 = weight_dens_kg_nm2 / 1000,
                                        weight_dens_mt_km2 = weight_dens_mt_nm2 / 3.43429)
# converted nm2 to km2

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
                                weight_dens_mt_km2_new = mean(weight_dens_mt_km2),
                                diff1 = max(weight_dens_mt_km2)-weight_dens_mt_km2_new,
                                diff2 = weight_dens_mt_km2_new - min(weight_dens_mt_km2),
                                sig_b = mean(sig_b),
                                mean.depth.m = mean(mean.depth.m),
                                N=length(new.ID)) %>%
                      mutate(bin_weight_dens = ifelse(weight_dens_mt_km2_new>0,1,0)) %>%
                      rename(weight_dens_mt_km2 = weight_dens_mt_km2_new)
  
ggplot(dat.acoustic.trim %>% filter(N>1)) +
    geom_point(aes(y = diff1/weight_dens_mt_km2,x= weight_dens_mt_km2)) +
    geom_abline(intercept = 0, slope=1) +
    geom_abline(intercept = 0, slope=0.5,color="red") 

ggplot(dat.acoustic.trim%>% filter(N>1)) +
  geom_histogram(aes(diff1/weight_dens_mt_km2)) 

dat.acoustic.trim <- dat.acoustic.trim %>% filter(N>1) %>% mutate(cv.coarse = diff1/weight_dens_mt_km2)

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

#### Here are some data plots..
setwd(script.dir)
source("Base_map.R")
base_map_trim +
    geom_point(data= dat.acoustic %>% filter(weight_dens_mt_km2>0),
               aes(x=lon,y=lat,
                 color=log(weight_dens_mt_km2),
                 fill=log(weight_dens_mt_km2)),alpha=0.5) +
  geom_point(data= dat.acoustic %>% filter(weight_dens_mt_km2==0),
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
  geom_point(aes(x=utm.lon,y=log(weight_dens_mt_km2)))+
  facet_wrap(~transect,scales="free_x") +
  theme_bw()

ggplot(dat.acoustic.trim %>% filter(transect<40)) +
  geom_point(aes(x=utm.lon,y=log(weight_dens_mt_km2)))+
  facet_wrap(~transect,scales="free_x") +
  theme_bw()

# Plot on identity scale
ggplot(dat.acoustic %>% filter(transect<40)) +
  geom_point(aes(x=utm.lon,y=weight_dens_mt_km2))+
  facet_wrap(~transect,scales="free_x") +
  theme_bw()

ggplot(dat.acoustic.trim %>% filter(transect<40)) +
  geom_point(aes(x=utm.lon,y=log(weight_dens_mt_km2)))+
  facet_wrap(~transect,scales="free_x") +
  theme_bw()

##################################
thresh_mt_km2 <- 0.763 / 3.43429 # This is the threshold determined from emails pasted above.
prob_val      <- 0.99  # This is the probability of observing presence at the threshold.

###################################################################
###################################################################
###################################################################
if(MODEL.TYPE == "lat.long.smooth"){
  # Derive basis function set up from the positive observations, make a second set of 
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
  # Derive basis function set up from the positive observations, make a second set of 
  # projections to the binomial components.
  
  ## FILL THIS IN TO HAVE A MODEL WITHOUT THE Bathymetry term.
  # LAT-LON TENSOR SMOOTHES
  utm.lon.lims <- c(min(dat.acoustic.bin$utm.lon), max(dat.acoustic.bin$utm.lon))
  utm.lat.lims <- c(min(dat.acoustic.bin$utm.lat), max(dat.acoustic.bin$utm.lat))
  knots.lon.pos    <- seq(utm.lon.lims[1],utm.lon.lims[2],length.out=N.knots.lon.pos)
  knots.lat.pos    <- seq(utm.lat.lims[1],utm.lat.lims[2],length.out=N.knots.lat.pos)
  
  # bottom.depth smooth
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
  "pos_weight_dens" = dat.acoustic.pos$weight_dens_mt_km2,
  "N_obs_bin"  = N_obs_bin,
  "N_obs_pos"  = N_obs_pos,

  # Index for leave.one.out
  # "loo_pos_idx" = ,

  # SMOOTH COMPONENTS.
  # Data for linear effects
  
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
  "nb_2" = smooth.dat.bin$nb_2,
  "knots_2" = smooth.dat.bin$knots_2,
  
  # basis function matrices
  "Zs_1_1_pos" = smooth.dat.pos$Zs_1_1,
  "Zs_1_2_pos" = smooth.dat.pos$Zs_1_2,
  "Zs_1_3_pos" = smooth.dat.pos$Zs_1_3,
  "Zs_2_1_pos" = smooth.dat.pos$Zs_2_1,
  
  "Zs_1_1_bin" = smooth.dat.bin$Zs_1_1,
  "Zs_1_2_bin" = smooth.dat.bin$Zs_1_2,
  "Zs_1_3_bin" = smooth.dat.bin$Zs_1_3,
  "Zs_2_1_bin" = smooth.dat.bin$Zs_2_1,
  
  "thresh_mt_km2" = thresh_mt_km2, # This is the offset to make the phi_0 be the parameter 
                             # define the probability of observing
                             # positive kg km^-2 at bin_offset kg km^-2
  
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
  "sigma",     # variability among samples, given individual bottle, site, and month 
  "log_lik_bin",
  "log_lik_pos"
)   

if(MODEL.TYPE=="lat.long.smooth"){
  stan_pars <- c(stan_pars,
    # smooth and linear coefficients
    "Intercept_pos",
    "Intercept_bin",
    
    "bs_pos",
    "bs_bin",
    "s_1_1_pos","s_1_2_pos","s_1_3_pos",
    "s_1_1_bin","s_1_2_bin","s_1_3_bin",
    "s_2_1_pos","s_2_1_bin",
    "theta_bin_pred",
    "D_pos_pred")
     
  N_knots_all <- list("N.knots.lon.pos"=N.knots.lon.pos,
                      "N.knots.lat.pos"=N.knots.lat.pos,
                      "N.knots.lon.bin"=N.knots.lon.bin,
                      "N.knots.lat.bin"=N.knots.lat.bin,
                      "N.knots.bd"=N.knots.bd)
}
if(MODEL.TYPE=="lat.long.smooth.base"){
  stan_pars <- c(stan_pars,
                 # smooth and linear coefficients
                  "Intercept_pos",
                  "Intercept_bin",
                  "bs_pos",
                  "bs_bin",
                  "s_1_1_pos","s_1_2_pos","s_1_3_pos",
                  "s_1_1_bin","s_1_2_bin","s_1_3_bin",
                 "theta_bin_pred",
                 "D_pos_pred")
                 
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
Warm = 1200
Iter = 3000
Treedepth = 12
Adapt_delta = 0.99

LOC <- paste0(base.dir,"Scripts/Stan Files/")
setwd(LOC)

if(MODEL.TYPE=="lat.long.smooth"){
  if(MODEL.VAR=="Base_Var"){
  stanMod = stan(file = "acoustic_Hake_smoothes.stan" ,data = stan_data, 
                 verbose = FALSE, chains = N_CHAIN, thin = 2, 
                 warmup = Warm, iter = Warm + Iter, 
                 control = list(max_treedepth=Treedepth,adapt_delta=Adapt_delta,metric="diag_e"),
                 pars = stan_pars,
                 boost_lib = NULL,
                 sample_file = paste0("./Output files/",MODEL.TYPE,"_",MODEL.ID,"_Acoustics.csv"),
                 init = stan_init_f2(n.chain=N_CHAIN)
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
              )
  }
}

samp_params <- get_sampler_params(stanMod)
stanMod_summary <- summary(stanMod)$summary
pars <- rstan::extract(stanMod, permuted = TRUE)
#########################################################################

# get_adaptation_info(stanMod)
pars <- rstan::extract(stanMod, permuted = TRUE)
log_lik_bin <- extract_log_lik(stanMod, parameter_name = "log_lik_bin", merge_chains = FALSE)
log_lik_pos <- extract_log_lik(stanMod, parameter_name = "log_lik_pos", merge_chains = FALSE)
r_eff_bin <- relative_eff(exp(log_lik_bin))
r_eff_pos <- relative_eff(exp(log_lik_pos))
loo_val_bin <- loo(log_lik_bin, r_eff = r_eff_bin)
loo_val_pos <- loo(log_lik_pos, r_eff = r_eff_pos)
print(loo_val_bin)
plot(loo_val_bin)
print(loo_val_pos)
plot(loo_val_pos)

#########################################################################

samp_params <- get_sampler_params(stanMod)
#samp_params 
stanMod_summary <- summary(stanMod)$summary
stanMod_summary_parts <- list()

stanMod_summary_parts[[as.name("param")]] <- summary(stanMod,pars=c(
  "sigma"))$summary     # variability among samples, given individual bottle, site, and month ))$summary
if(MODEL.TYPE == "lat.long.smooth"){
  stanMod_summary_parts[[as.name("smoothes")]] <- summary(stanMod, pars=c("Intercept_bin",
                                                                          "Intercept_pos",
                                                                          "bs_bin",
                                                                          "bs_pos",
                                                                          "s_1_1_bin","s_1_2_bin","s_1_3_bin",
                                                                          "s_2_1_bin",
                                                                          "s_1_1_pos","s_1_2_pos","s_1_3_pos",
                                                                          "s_2_1_pos"))$summary
}
if(MODEL.TYPE == "lat.long.smooth.base"){
  stanMod_summary_parts[[as.name("smoothes")]] <- summary(stanMod, pars=c("Intercept_bin",
                                                                          "Intercept_pos",
                                                                          "bs_bin",
                                                                          "bs_pos",
                                                                          "s_1_1_bin","s_1_2_bin","s_1_3_bin",
                                                                          "s_1_1_pos","s_1_2_pos","s_1_3_pos"
                                                                          ))$summary
  
}

# Make diagnostic plots of fitting.
# Predicted Observed plots
# Sum to 1 nm bins.


##### MAKE SOME DIAGNOSTIC PLOTS
TRACE <- list()
TRACE[[as.name("D_pos_pred")]] <- traceplot(stanMod,pars=c("lp__","D_pos_pred[12]","D_pos_pred[234]","D_pos_pred[456]","D_pos_pred[855]"),inc_warmup=FALSE)
TRACE[[as.name("theta_bin_pred")]] <- traceplot(stanMod,pars=c("lp__","theta_bin_pred[12]","theta_bin_pred[234]","theta_bin_pred[456]","theta_bin_pred[855]"),inc_warmup=FALSE)
#TRACE[[as.name("mu_smooth")]] <- traceplot(stanMod,pars=c("lp__","mu_smooth[12]","mu_smooth[234]","mu_smooth[456]","mu_smooth[878]"),inc_warmup=FALSE)

if(MODEL.TYPE=="lat.long.smooth"){
  TRACE[[as.name("Int")]] <- traceplot(stanMod,pars=c("Intercept_bin","Intercept_pos"),inc_warmup=FALSE)
  TRACE[[as.name("bs")]] <- traceplot(stanMod,pars=c("bs_bin","bs_pos"),inc_warmup=FALSE)
  TRACE[[as.name("s_1")]] <- traceplot(stanMod,pars=c("s_1_1_pos","s_1_2_pos","s_1_3_pos",
                                                      "s_1_1_bin","s_1_2_bin","s_1_3_bin"),inc_warmup=FALSE)
  TRACE[[as.name("s_2")]] <- traceplot(stanMod,pars=c("s_2_1_bin","s_2_1_pos"),inc_warmup=FALSE)
}else if(MODEL.TYPE=="lat.long.smooth.base"){
  TRACE[[as.name("Int")]] <- traceplot(stanMod,pars=c("Intercept_bin","Intercept_pos"),inc_warmup=FALSE)
  TRACE[[as.name("bs")]] <- traceplot(stanMod,pars=c("bs_bin","bs_pos"),inc_warmup=FALSE)
  TRACE[[as.name("s_1")]] <- traceplot(stanMod,pars=c("s_1_1_pos","s_1_2_pos","s_1_3_pos",
                                                      "s_1_1_bin","s_1_2_bin","s_1_3_bin"),inc_warmup=FALSE)
}

Output.acoustic <- list(
  # STAN MODEL ASSOCIATED THINGS
  stanMod = stanMod, 
  stanMod_summary = stanMod_summary,
  loo_val_pos = loo_val_pos, #loo summary
  log_lik_pos = log_lik_pos, #log_likelihood object for summary
  r_eff_pos = r_eff_pos,
  loo_val_bin = loo_val_bin, #loo summary
  log_lik_bin = log_lik_bin, #log_likelihood object for summary
  r_eff_bin = r_eff_bin,
  Iter=Iter,
    stanMod_summary_parts = stanMod_summary_parts,
  #stanMod_summary_D = stanMod_summary_D,
  samp = pars, samp_params=samp_params,
  TRACE = TRACE,
  SPECIES = "hake",
  MODEL.TYPE = MODEL.TYPE,
  MODEL.VAR = MODEL.VAR,
  MODEL.ID = MODEL.ID,
  # Input Data
  brms.object.pos =brms.object.pos, # for constructing the model.
  brms.object.bin =brms.object.bin, # for constructing the model.
  N_knots_all = N_knots_all,# for constructing the model.
  dat.acoustic.bin =dat.acoustic.bin,
  dat.acoustic.pos = dat.acoustic.pos,
  dat_raster_fin = dat_raster_fin
  #dat_raster_trim = dat_raster_trim
  
)

setwd(base.dir)
setwd("./Stan Model Fits/")
save(Output.acoustic,file=paste("Acoustics 2019",MODEL.TYPE,MODEL.ID,MODEL.VAR,"Fitted.RData"))
######
