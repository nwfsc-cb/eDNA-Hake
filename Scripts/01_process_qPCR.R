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
                mutate(qPCR_year = paste0(qPCR,"_",year))
dat.stand <- rbind(dat.stand.2019%>% mutate(year=2019),dat.stand.2021%>% mutate(year=2021)) %>% 
                mutate(qPCR_year = paste0(qPCR,"_",year)) %>% arrange(year,qPCR_year) %>%
                mutate(ln_copies = log(copies_ul))
# Keep only the samples that are in the field samples data.frame
PCR <- data.frame(qPCR_year = levels(as.factor(dat.samp$qPCR_year))) %>% mutate(plate_idx = 1:nrow(.))
# Merge in plate index to samples and standards.
# check the the levels for the standards and the samples are identical:
dat.stand <- dat.stand %>% filter(qPCR_year %in% PCR$qPCR_year)
identical(sort(unique(dat.samp$qPCR_year)) ,sort(unique(dat.stand$qPCR_year)) )

dat.samp <- dat.samp %>% left_join(.,PCR) %>% arrange(year,sample)
dat.stand <- dat.stand %>% left_join(.,PCR) 

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
# ------ spatial meshes and design matrices
# ------ fixed effect matrices
# ------ smooth design matrices
# ------ random effect matrices 

##################################################################
# Call raw data plotting scripts.
source(here('Scripts',"Make_matrices.R"),local=T)
##################################################################

###################################################################

library(here)
library(TMB)
setwd("../src")

tmb_data <- list(Y_i = Y_i, 
                 # Read in standards data
                 Y_stand = dat.stand$Ct ,
                 Z_stand = dat.stand$Ct_bin ,
                 plate_stand_idx = dat.stand$plate_idx,
                 ln_copies = dat.stand$ln_copies,
                
                # Read in samples data.
                 Y_samp = dat.samp$Ct ,
                 Z_samp = dat.samp$Ct_bin ,
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
                 
                 # Fixed and smoothes for the mean effect
                 Xf =  Xf,
                 Xs = SM$Xs,
                 Zs = SM$Zs,
                 n_smooth = n_smooth,
                 has_smooths = as.integer(has_smooths),
                 b_smooth_start = b_smooth_start,
                 
                 #. Indexes for observations
                 year_idx = year_idx,
                 depth_idx = depth_idx,
                 
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
                 A_ID_idx = A_ID_idx, #station ID to match up observation and A_st output
                
                # Priors
                std_mu = c(2, # mean for logit-scale intercept
                           2, # mean for logit-scale slope
                           40, # mean for positive intercept
                           -1.35),# mean for positive slope
                std_sigma = c(2, # SD for logit-scale intercept
                              2, # SD for logit-scale slope
                              5, # SD for positive intercept
                              0.1),# SD for positive slope
)

tmb_params <- list(# Regression terms
  
  # standards phi = binomial part, beta = positive part
  std_phi_0 = 0.1,
  std_phi_1 = -2,
  std_beta_0 = 30 ,
  std_beta_0 = -1.3,
  
  # fixed effect terms.
  betaf=rnorm(ncol_beta,0,0.1), 
  
  # smooth terms
  bs=rep(0,n_bs),
  b_smooth = b_smooth,
  ln_smooth_sigma = runif(n_smooth,-0.1,0.1),
  
  # Weight terms
  # bs_L = matrix(0,n_y,n_f),
  # b_L_smooth = array(0,dim=c(n_y,n_f,ncol(Z_L$Z_L[[1]]))),
  # ln_L_smooth_sigma = matrix(runif(n_y*n_f,-0.1,0.1),n_y,n_f),

  # SPDE terms
  ln_kappa = 0,
  
  # Latent field
  omega_s = matrix(rnorm(n_s*n_f,-0.1,0.1),n_s,n_fy),
  # observation variance terms.
  lnSigma = 0
)
  
tmb_random <- c("bs",
                "b_smooth",
                # "bs_L",
                # "b_L_smooth",
                "omega_s")

TMB::compile("RegSmooth_space.cpp")
dyn.load(dynlib("RegSmooth_space"))

obj <- MakeADFun(data=tmb_data, 
                 parameters=tmb_params, 
                 random=tmb_random, 
                 DLL="RegSmooth_space",
                 hessian=TRUE)

# dyn.unload(dynlib("RegSmooth_space"))

opt <- nlminb(obj$par,obj$fn,obj$gr)
summary(sdreport(obj))
rep <- obj$report()


opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)

#### WORK ON EXTRACTING ESTIMATES.

SM$Xs

A <- data.frame(depth_cat = dat.samp$depth_cat,year=dat.samp$year,bdc = dat.samp$bottom.depth.consensus)
tmp <-rep$eta_smooth_all %>% as.data.frame()
colnames(tmp) <- c("x1")
                   
                   ,"x2")"x3","x4")
B <- bind_cols(A,tmp)

ggplot(B) +
    geom_point(aes(y=x1,x=depth_cat))+
    facet_wrap(~year)

ggplot(B) +
  geom_point(aes(y=x2,x=depth_cat))+
  facet_wrap(~year)

ggplot(B) +
  geom_point(aes(y=x3,x=bdc))+
  facet_wrap(~year)

ggplot(B) +
  geom_point(aes(y=x4,x=bdc))+
  facet_wrap(~year)










## STANDARD CURVE IMPLEMENTATION.


if(stdcurve_flag == 1) {
  // Presence-Absence component of model.
  vector<Type> theta_stand(N_stand_bin);
  vector<Type> std_sigma = exp(ln_std_sigma);
  
  for(int i = 0; i < std_phi_0.size(); i++){
    jnll -= dnorm(std_phi_0(i), std_mu(0), std_sigma(0), true);//phi_0 ~ normal(2, 2)
    jnll -= dnorm(std_phi_1(i), std_mu(1), std_sigma(1), true);//phi_1 ~ normal(4, 2)
  }
  for(int i = 0; i < N_stand_bin; i++){ // likelihood
    theta_stand(i) = std_phi_0(pcr_stand_bin_idx(i)) + std_phi_1(pcr_stand_bin_idx(i)) * (D_bin_stand(i));// - stand_offset);
jnll -= dbinom_robust(bin_stand(i), Type(1), theta_stand(i), true); // likelihood
  }
  // Positive component of model.
  for(int i = 0; i < std_beta_0.size(); i++){
    jnll -= dnorm(std_beta_0(i), std_mu(2), std_sigma(2), true);//beta_0 ~ normal(40,5)
    jnll -= dnorm(std_beta_1(i), std_mu(3), std_sigma(3), true);//beta_1 ~ normal(-3.32,0.1)
  }
  vector<Type> kappa_stand(N_stand_pos);
  //sigma_all_stand = exp(log_sigma_all_stand);
  for(int i = 0; i < N_stand_pos; i++){
    kappa_stand(i) = std_beta_0(pcr_stand_pos_idx(i)) + std_beta_1(pcr_stand_pos_idx(i)) * (D_pos_stand(i));// - stand_offset);
jnll -= dnorm(pos_stand(i), kappa_stand(i), phi(0), true); // likelihood
  }
  
  REPORT(std_beta_0);
  REPORT(std_beta_1);
  REPORT(std_phi_0);
  REPORT(std_phi_1);
  REPORT(std_mu);
  REPORT(std_sigma);
  REPORT(ln_std_sigma);
  //REPORT(sigma_all_stand);
  ADREPORT(std_beta_0);
  ADREPORT(std_beta_1);
  ADREPORT(std_phi_0);
  ADREPORT(std_phi_1);
  ADREPORT(std_mu);
  ADREPORT(std_sigma);
  ADREPORT(ln_std_sigma);
  //ADREPORT(sigma_all_stand);
}

## Section 2

case stdcurve_family: {
  if(y_i(i,m) > 0) {
    tmp_ll = dnorm(y_i(i,m), std_beta_0(pcr_idx(i)) + std_beta_1(pcr_idx(i)) * mu_i(i,m), phi(m), true);
    tmp_ll += dbinom_robust(Type(1), size(i), std_phi_0(pcr_idx(i)) + std_phi_1(pcr_idx(i)) * mu_i(i,m), true);
  } else {
    tmp_ll = dbinom_robust(Type(0), size(i), std_phi_0(pcr_idx(i)) + std_phi_1(pcr_idx(i)) * mu_i(i,m), true);
  }
  SIMULATE{y_i(i,m) = rbinom(size(i), invlogit(std_phi_0(pcr_idx(i)) + std_phi_1(pcr_idx(i)) * mu_i(i,m))) * rnorm(std_beta_0(pcr_idx(i)) + std_beta_1(pcr_idx(i)) * mu_i(i,m), phi(m));}
  break;
}






