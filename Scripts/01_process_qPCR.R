# This is the main file for reading in data, calling cleaning scripts for each year, 
#  and then combining the results

# This script is for reading in qPCR data, merging it with locations, and starting a set of plots.  
# It will likely be called in a rMarkdown file later.

rm(list=ls())
#setwd('./Github/eDNA-Hake')
# Libraries
library(tidyverse)
library(marmap)
library(ggplot2)
library(rstan)
library(lubridate)
library(reshape2)
library(gridExtra)
library(raster)
#library(rgdal)
library(sp)
library(brms)
library(sdmTMB)
library(loo)
library(here)

###########################################################################
## DECLARED SPECIES OF INTEREST
SP <- "hake" # options: hake, lamprey, eulachon

# Specify an inhibition limit for retaining samples.
INHIBIT.LIMIT <- 0.5

# load and run the acoustic data. this is needed to reference the offshore-ness of 
#setwd(script.dir)

source(here('Scripts',"process acoustic data for qPCR 2019 on.R"),local=TRUE)
# dat.acoustic and dat.acoustic.binned are the relevant data frames

# Read in 2019 data, clean, add some indexes.
source(here('Scripts',"process 2019 qPCR results.R"),local=TRUE)

# Read in 2021 data, clean, add some indexes.
source(here('Scripts',"process 2021 qPCR results.R"),local=TRUE)


# relevant checks on dimension and names for matrices
identical(names(dat.samp.2019),names(dat.samp.2021))
identical(names(dat.stand.2019),names(dat.stand.2021))
identical(names(PCR.2019),names(PCR.2021))
identical(names(dat.control.field.neg.2021),names(dat.control.field.neg.2019))
identical(names(dat.control.2021),names(dat.control.2019))
identical(names(dat.inhibit.2021),names(dat.inhibit.2019))

##### Combine outputs into single data.frame frame across years.
dat.samp <- rbind(dat.samp.2019,dat.samp.2021) %>%
                mutate(qPCR_year = paste0(qPCR,"_",year)) %>%
                dplyr::select(-plate_idx)
dat.stand <- rbind(dat.stand.2019%>% mutate(year=2019),dat.stand.2021%>% mutate(year=2021)) %>% 
                mutate(qPCR_year = paste0(qPCR,"_",year)) %>% arrange(year,qPCR_year) %>%
                mutate(ln_copies = log(copies_ul)) %>%
                dplyr::select(-plate_idx)
# Keep only the PCR plates that are in the field samples data.frame
PCR <- data.frame(qPCR_year = levels(as.factor(dat.samp$qPCR_year))) %>% mutate(plate_idx = 1:nrow(.))
n_plate <- nrow(PCR)

# Read in spatial data and data for bottom depth at each sample
# REPLACES THE SPATIAL DATA ANALYSIS DONE IN THE INDIVIDUAL YEAR SCRIPTS.
source(here('Scripts',"Process Spatial data.R"),local=TRUE)
# relevant files here are dat_raster_trim an the new depth data in dat.samp.

# Merge in plate index to samples and standards.
# check the the levels for the standards and the samples are identical:
dat.stand <- dat.stand %>% filter(qPCR_year %in% PCR$qPCR_year)
identical(sort(unique(dat.samp$qPCR_year)) ,sort(unique(dat.stand$qPCR_year)) )

dat.samp <- dat.samp %>% left_join(.,PCR) %>% arrange(year,sample)
# add transformation of bottom_depth
dat.samp <- dat.samp %>% mutate(ln.bottom.depth.NGDC = log(bottom.depth.NGDC))

dat.stand <- dat.stand %>% left_join(.,PCR) 
# ACCOUNT FOR 2ul of volume in standards reaction:
dat.stand <- dat.stand %>% mutate(log10_copies = log10(copies_ul) + log10(2),
                                  ln_copies = log(copies_ul) + log(2))

dat.control.field.neg <- rbind(dat.control.field.neg.2019,dat.control.field.neg.2021 )
dat.control <- rbind(dat.control.2019,dat.control.2021 )
dat.inhibit <- rbind(dat.inhibit.2019 %>% mutate(year=2019), dat.inhibit.2021 %>% mutate(year=2021) )

##################################################################
# Call raw data plotting scripts.
source(here('Scripts',"plot_raw_observations.R"),local=T)
##################################################################

##################################################################
# Make design matrices needed for fitting TMB model
# this includes:
# ------ making helper objects for the standards
# ------ spatial meshes and design matrices
# ------ fixed effect matrices
# ------ make offset vectors
# ------ smooth design matrices
# ------ random effect matrices 

##################################################################
# Define formulas for each model component
# Offsets:
# Always keep dilution, volume filtered as offsets.
# -- these are ln_vol_sample, ln_dilution.
# Add wash error estimate 
# -- wash_cov

# Fixed effects (basically, just year)
FORM.fixed <- copies_ul ~ year
#Smooth effects (allow a intercept for each depth category x year, add smooth effect of bottom depth.
FORM.smoothes <- "copies_ul ~ s(ln.bottom.depth.NGDC,k=3)" #s(depth_cat,k=3) + 
                                                        #+ s(bottom.depth.NGDC,by=year,k=4)"
# Smooth used in the weight matrix L.  basically only a factor of depth_cat
FORM.L <- "Y ~ s(depth_cat,k=3)"

# Define whether to use a Random bottle effect or not.
TRAD.RE.BOTTLE = FALSE

# Define whether to allow anisotropy (0=yes, 1=No)
anisotropy = 1

# Prediction set for depth categories.
depth_pred <- seq(0,500,by=50)
n_d_pred <- length(depth_pred)

# Define the number of spatial fields to use for each year
n_f <- 3
# Define the number of years
n_y <- dat.samp %>% distinct(year) %>% pull(year) %>% length()
# Define a vector controlling the number of spatial fields to use for each year
n_fy <- n_f*n_y

# Define bottle SD form
tau_has_smooths <- 0 # 0 for yes, 1 for no.
n_tau = ifelse(tau_has_smooths == 0,3,2)


# THIS IS THE CONVERSION FROM COPIES IN THE REACTION TO COPIES PER L
# Derived assuming a 2.5L water sample eluted in 100uL Longmire and 2 uL used in each PCR reaction.
# Each 2uL of sample corresponds to 2% (2 of 100) of the total sample and therefore is equivalent
# to the copies DNA contained in 50mL of water. (0.02 x 2500 ml = 50 ml)
EXPAND = 20

dat.samp$ln_expand_offset = log(1/EXPAND)

# Make design matrices needed for fitting TMB model.
source(here('Scripts',"Make_matrices.R"),local=T)



##################################################################
# Put the data in a form that can be fed to TMB
###################################################################

library(here)
library(TMB)
setwd("../src")

tmb_data <- list(#Y_i = Y_i, 
                 # Read in standards data
                 Y_stand = dat.stand$Ct ,
                 Z_stand = dat.stand$Ct_bin ,
                 Z_stand_int = dat.stand$Ct_bin ,
                 plate_stand_idx = dat.stand$plate_idx,
                 ln_copies = dat.stand$ln_copies,
                
                # Read in samples data.
                 Y_samp = dat.samp$Ct ,
                 Z_samp = dat.samp$Ct_bin ,
                 Z_samp_int = dat.samp$Ct_bin ,
                 plate_samp_idx = dat.samp$plate_idx,

                 # Read in counters
                 n_i = n_i,   # number of total observations
                 n_stand = n_stand, # number of observations in the standards data frame
                 n_f = n_f,  # number of factors per year for weight matrix L. (currently)
                 n_d = n_d,   # number of distinct depths in model
                 n_y = n_y,  # number of years in model
                 n_s = n_s,  # number of knot locations (constant across )
                 n_st = n_st, # number of unique station locations
                 n_fy = n_fy, # number of factor-year combinations
                 n_plate = n_plate, # number of qPCR plates.
                 n_bottle =n_bottle, # number of unique bottles.
                 n_d_pred = n_d_pred, # number of predicted depths
                 n_pred = n_pred, # number of prediction locations x depths total.
                
                # Read in offsets
                ln_dilution_offset = dat.samp$ln_dilution,
                ln_vol_offset = dat.samp$ln_vol_sample,
                ln_expand_offset = dat.samp$ln_expand_offset,
                
                # Wash offset for 2019 (0/1 vector that I want to keep separate from other covariates).
                wash_cov = dat.samp$wash_idx,

                 # Fixed and smoothes for the mean effect
                 Xf =  Xf,
                 Xs = SM$Xs,
                 Zs = SM$Zs,
                 n_smooth = n_smooth,
                 has_smooths = as.integer(has_smooths),
                 b_smooth_start = b_smooth_start,
                
                  # smaller matrices that make the predictions easier
                Xf_new =  Xf_new,
                # Xs_new = SM_new$Xs,
                # Zs_new = SM_new$Zs,
                
                # Random effect for bottle Design matrix.
                X_bottle = X_bottle,
                bottle_RE_idx = bottle_RE$RE_idx,
                
                # Smooth of depth_cat SD for bottle RE
                Xs_ln_tau = TAU$Xs %>% as.vector()  ,
                Zs_ln_tau  = TAU$Zs %>% as.vector() ,
                #n_tau_smooth = n_tau_smooth,
                tau_has_smooths = as.integer(tau_has_smooths),
                #b_tau_smooth_start = b_tau_smooth_start,
                Xs_ln_tau_sm = TAU_sm$Xs %>% as.vector(),
                Zs_ln_tau_sm  = TAU_sm$Zs %>% as.vector(),

                 #. Indexes for observations
                 year_idx = year_idx,
                 depth_idx = depth_idx,
                 year_pred_idx = year_pred_idx,
                 depth_pred_idx = depth_pred_idx, 

                 # Weight Matrices and helpers
                 b_L_smooth_start = b_L_smooth_start,
                 Z_L = Z_L$Z_L, # [L]ist [O]f (basis function matrices) [Matrices]
                 X_L = X_L$X_L, # smoother linear effect matrix
                 
                 #factor indexes for weight matrix 
                 F_L_idx = F_L_idx, #  // integer index for years (of length n_fy)
                 Y_L_idx = Y_L_idx,  #// integer index for years (of length n_fy)
                 FY_start_idx = FY_start_idx, # // integer index for years (of length n_fy)
                 
                 # Spatial smoothers
                 A_st =A_st, # projection from knot locations to observationss
                 A_spatial_index = spde$sdm_spatial_id - 1L,
                 spde = spde,
                 spde_aniso = spde_aniso,
                 anisotropy = anisotropy,
                 A_ID_idx = A_ID_idx, #station ID to match up observation and A_st output
                
                 # PREDICTION MATRICES # from dat.pred
                 Xf_pred  =  Xf_pred,   # 
                 Xs_pred  = SM_pred$Xs,
                 Zs_pred  = SM_pred$Zs,
                 Z_L_pred = Z_L_pred$Z_L, # [L]ist [O]f (basis function matrices) [Matrices]
                 X_L_pred = X_L_pred$X_L, # smoother linear effect matrix
                 A_pred   = A_pred,
                 n_pred_locs = n_pred_locs,
                 A_ID_pred_idx = A_ID_pred_idx #station ID to match up observation and A_st output
                
                # Priors
                # std_mu_prior = c(2, # mean for logit-scale intercept
                #            2, # mean for logit-scale slope
                #            40, # mean for positive intercept
                #            -1.35),# mean for positive slope
                # std_sigma_prior = c(2, # SD for logit-scale intercept
                #               2, # SD for logit-scale slope
                #               5, # SD for positive intercept
                #               0.1)# SD for positive slope
)

tmb_params <- list(# Regression terms
  
  # standards phi = binomial part, beta = positive part
  # hierarchical mean and sd in std_mu and std_sigma.
  std_mu = c(2,1.5,40,-1.35),
  ln_std_sigma = log(c(1,0.2,5,0.1)),
  
  std_phi_0 = rep(0.1,n_plate),
  std_phi_1 = rep(0,n_plate),
  std_beta_0 = rep(40,n_plate) ,
  std_beta_1 = rep(-1.3,n_plate),
  
  # observation variance term (Ct scale) intercept and slope
  ln_std_tau = c(1,-0.1),
  
  # fixed effect terms.
  betaf=c(rnorm(1,2,0.1),rnorm(ncol_beta-1,1,0.1)), 
  
  # wash_cov
  gamma_wash = 0,
  #random effect of bottle
  gamma_bottle = rep(0,n_bottle),
  ln_tau_bottle = c(-2,rep(0,n_tau-1)),
  #ln_tau_bottle_test =0,
  
  # smooth terms for covariates.
  bs=rep(0,n_bs),
  b_smooth = b_smooth,
  ln_smooth_sigma = runif(n_smooth,-0.1,0.1),
  
  # Weight terms
  bs_L = matrix(0,n_f,2), # linear components of weight matrix (intercept, slope)
  b_L_smooth = array(0,dim=c(n_f,ncol(Z_L$Z_L[[1]]))),
  ln_L_smooth_sigma = 0,#matrix(runif(n_f,-0.1,0.1),n_f,1),

  # Extra obs varaition,
  ln_sample_tau = -4,
  
  # SPDE terms
  ln_kappa = -4,
  
  # Latent field
  omega_s = matrix(rnorm(n_s*n_f,-0.1,0.1),n_s,n_fy)
  
)

if(anisotropy == 0){
  tmb_params <- c(tmb_params,
                  list(ln_H_input= matrix(0, nrow = 2L, ncol = 1)))
}

tmb_random <- c(#"bs",
                "b_smooth",
                "std_beta_0",
                "std_beta_1",
                "std_phi_0",
                "std_phi_1",
                "b_L_smooth",
                "omega_s",
                "gamma_bottle")

# Compile the model, load to make it accessible to the optimizer.
TMB::compile("Hake_full.cpp")
dyn.load(dynlib("Hake_full"))
# dyn.unload(dynlib("Hake_full")) 

# makethe objective function with the data
obj <- MakeADFun(data=tmb_data, 
                 parameters=tmb_params, 
                 random=tmb_random, 
                 DLL="Hake_full",
                 #profile="gamma_bottle",
                 hessian=TRUE)

opt <- nlminb(obj$par,obj$fn,obj$gr,
              control=list(eval.max=10000,iter.max=8000))
report <- obj$report() 
sd_rep <- sdreport(obj,getJointPrecision=TRUE)

report$jnll
report$range
report$betaf
report$ln_std_tau
report$tau_bottle_sm
report$ln_sample_tau

report$omega_s
report$gamma_bottle


opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian





#### WORK ON EXTRACTING ESTIMATES.
