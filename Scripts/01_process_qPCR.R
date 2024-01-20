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
FORM.smoothes <- "copies_ul ~ s(depth_cat,by=year,k=4) + s(bottom.depth.NGDC,by=year,k=4)"
# Smooth used in the weight matrix L.  basically only a factor of depth_cat
FORM.L <- "Y ~ s(depth_cat,k=3)"

# Define whether to use a Random bottle effect or not.
RANDOM_BOTTLE = TRUE

# Define the number of spatial fields to use for each year
n_f <- 2
# Define the number of years
n_y <- dat.samp %>% distinct(year) %>% pull(year) %>% length()
# Define a vector controlling the number of spatial fields to use for each year
n_fy <- n_f*n_y

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
                
                # Read in offsets
                ln_dilution_offset = dat.samp$ln_dilution,
                ln_vol_offset = dat.samp$ln_vol_sample,
                
                # Wash offset for 2019 (0/1 vector that I want to keep separate from other covariates).
                wash_cov = dat.samp$wash_idx,

                 # Fixed and smoothes for the mean effect
                 Xf =  Xf,
                 Xs = SM$Xs,
                 Zs = SM$Zs,
                 n_smooth = n_smooth,
                 has_smooths = as.integer(has_smooths),
                 b_smooth_start = b_smooth_start,
                
                # Random effect for bottle Design matrix.
                X_bottle = X_bottle,
                 
                 #. Indexes for observations
                 year_idx = year_idx,
                 depth_idx = depth_idx,
                 
                 # Weight Matrices and helpers
                 b_L_smooth_start = b_L_smooth_start,
                 Z_L = Z_L$Z_L, # [L]ist [O]f (basis function matrices) [Matrices]
                 X_L = X_L$X_L, # smoother linear effect matrix
                 
                 # Design Matrix for the weight matrix.
                 L_design = L_design,
                
                 #factor indexes for weight matrix 
                 F_L_idx = F_L_idx, #  // integer index for years (of length n_fy)
                 Y_L_idx = Y_L_idx,  #// integer index for years (of length n_fy)
                 FY_start_idx = FY_start_idx, # // integer index for years (of length n_fy)
                 
                 # Spatial smoothers
                 A_st =A_st, # projection from knot locations to observationss
                 A_spatial_index = spde$sdm_spatial_id - 1L,
                 spde = spde,
                 A_ID_idx = A_ID_idx #station ID to match up observation and A_st output
                
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
  std_mu = c(2,2,40,-1.35),
  ln_std_sigma = log(c(2,2,5,0.1)),
  
  std_phi_0 = rep(0.1,n_plate),
  std_phi_1 = rep(-2,n_plate),
  std_beta_0 = rep(30,n_plate) ,
  std_beta_1 = rep(-1.3,n_plate),
  
  # observation variance term (Ct scale) intercept and slope
  ln_std_tau = c(0,-0.1),
  
  # fixed effect terms.
  betaf=rnorm(ncol_beta,0,0.1), 
  
  # wash_cov
  gamma_wash = 0,
  #random effect of bottle
  # gamma_bottle = rep(0,n_bottle),
  # ln_sigma_bottle = 0,
  # smooth terms
  bs=rep(0,n_bs),
  b_smooth = b_smooth,
  ln_smooth_sigma = runif(n_smooth,-0.1,0.1),
  
  # Weight terms
  bs_L = matrix(0,n_y,n_f),
  b_L_smooth = array(0,dim=c(n_y,n_f,ncol(Z_L$Z_L[[1]]))),
  ln_L_smooth_sigma = matrix(runif(n_y*n_f,-0.1,0.1),n_f,n_y),

  # SPDE terms
  ln_kappa = 0,
  
  # Latent field
  omega_s = matrix(rnorm(n_s*n_f,-0.1,0.1),n_s,n_fy)
  
)
  
tmb_random <- c("bs",
                "b_smooth",
                "std_beta_0",
                "std_beta_1",
                "std_phi_0",
                "std_phi_1",
                "bs_L",
                "b_L_smooth",
                "omega_s")#,
                #"gamma_bottle")

# Compile the model, load to make it accessible to the optimier.
TMB::compile("Hake_full.cpp")
dyn.load(dynlib("Hake_full"))
# dyn.unload(dynlib("Hake_full"))

# makethe objective function with the data
obj <- MakeADFun(data=tmb_data, 
                 parameters=tmb_params, 
                 random=tmb_random, 
                 DLL="Hake_full",
                 hessian=TRUE)

opt <- nlminb(obj$par,obj$fn,obj$gr, 
              control=list(eval.max=10000,iter.max=8000))
report <- obj$report()
summary(sdreport(obj))


report$omega_s


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


#3 Pred- Obs Plots
pred.samp <- data.frame(theta_samp =report$theta_samp,
                        kappa_samp = report$kappa_samp)

                   
pred.stand <- data.frame(theta_stand = report$theta_stand,
                          kappa_stand = report$kappa_stand)


pred.samp <- cbind(dat.samp,pred.samp)

pred.stand <- cbind(dat.stand,pred.stand)

pred.stand$sd_pred <- exp(report$ln_std_tau[1] + report$ln_std_tau[2] * pred.stand$ln_copies)
pred.stand <- pred.stand %>% mutate(kappa_plus = kappa_stand+sd_pred,
                                              kappa_minus = kappa_stand-sd_pred)


### MAKE SOME STANDARDS PLOTS
STAND.BREAKS = c(1,5,10,100,1000,10000,100000)

ggplot(pred.stand) +
  geom_jitter(aes(y=Ct_bin,x=exp(ln_copies)),width=0,alpha=0.5,height=0.05) +
  geom_point(aes(y=plogis(theta_stand),x=exp(ln_copies)),color="red") +
  geom_line(aes(y=plogis(theta_stand),x=exp(ln_copies)),color="red") +
  scale_x_continuous(trans="log",breaks = STAND.BREAKS)+
  facet_wrap(~qPCR_year) + 
  theme_bw()

ggplot(pred.stand %>% filter(Ct_bin ==1)) +
  geom_jitter(aes(y=Ct,x=exp(ln_copies)),width=0,alpha=0.5,height=0.05) +
  #geom_point(aes(y=kappa_stand,x=exp(ln_copies)),color="red") +
  geom_line(aes(y=kappa_stand,x=exp(ln_copies)),color="red") +
  geom_line(aes(y=kappa_plus,x=exp(ln_copies)),color="red",linetype="dashed") +
  geom_line(aes(y=kappa_minus,x=exp(ln_copies)),color="red",linetype="dashed") +
  scale_x_continuous(trans="log",breaks = STAND.BREAKS)+
  facet_wrap(~qPCR_year) + 
  theme_bw()

# Unknown samples.
ggplot(pred.samp) +
  geom_jitter(aes(y=Ct_bin,x=plogis(theta_samp)),width=0,alpha=0.5,height=0.05) +
  stat_smooth(aes(y=Ct_bin,x=plogis(theta_samp))) +
  geom_abline(intercept=0,slope = 1,linetype="dashed",color="red") +
  scale_x_continuous(limits=c(0,1)) +
  facet_grid(year~depth_cat)  +
  theme_bw()

 ggplot(pred.samp %>% filter(Ct_bin ==1)) +
  geom_jitter(aes(y=Ct,x=kappa_samp),width=0,alpha=0.5,height=0.05) +
  stat_smooth(aes(y=Ct,x=kappa_samp)) +
  geom_abline(intercept=0,slope = 1,linetype="dashed",color="red") +
  facet_grid(year~depth_cat)  +
   theme_bw()
#  scale_x_continuous(limits=c(0,1))

## Make some marginal plots of the fixed and smooth effects
 
 # First Make some plots of the fixed effects of depth_category by year.
    FORM.smoothes
  # 
    A <- dat.samp %>% dplyr::select(year,depth_cat) %>% bind_cols(val=rowSums(report$eta_smooth_i[,1:n_y]))
    B <- A %>% distinct(year,depth_cat,val) %>% mutate(year=as.factor(year))
 
    ggplot(B) +
      geom_point(aes(x=depth_cat,y=val,color=as.factor(year))) +
      geom_line(aes(x=depth_cat,y=val,color=as.factor(year))) +
      scale_y_continuous("Marginal effect of water depth") +
      scale_x_continuous("Depth(m)")+
      theme_bw()
    

    A <- dat.samp %>% dplyr::select(year,bottom.depth.NGDC) %>% bind_cols(val=rowSums(report$eta_smooth_i[,(n_y+1):(n_y*2)]))
    B <- A %>% mutate(year=as.factor(year))
    
    ggplot(B) +
      #geom_point(aes(x=bottom.depth.NGDC,y=val,color=as.factor(year))) +
      geom_line(aes(x=bottom.depth.NGDC,y=val,color=as.factor(year))) +
      scale_y_continuous("Marginal effect bottom depth") +
      scale_x_continuous("Depth(m)")+
      theme_bw()
    
#######
    ### Full spatial predictions at the observea locations.
#######
    D_pred  <- dat.samp %>% bind_cols(.,D_i=report$D_i) %>% 
                      distinct(year,station,sample,depth_cat,lon,lat,D_i) %>%
                      mutate(D_i_exp = exp(D_i),
                        D_mod = ifelse(D_i < -2, -2, D_i))
    
    ggplot(D_pred) +
        geom_point(aes())
    
    p_D_pred <- base_map_trim + 
      geom_point(data=D_pred,aes(x=lon,y=lat,fill=D_mod,color=D_mod),shape=21,alpha=0.5,size=3) +
      #scale_size("Copies / uL",labels=LAB,breaks=LAB,range=c(0.2,20),limits=c(0.5,NA))+
      # geom_point(data=SAMPLES%>% filter(Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
      # scale_shape("Copies / uL",solid=FALSE)+
      
      scale_color_viridis_c("Estimated\n log(DNA conc.)",option="plasma",begin=0,end=0.8) +
      scale_fill_viridis_c("Estimated\n log(DNA conc.)",option="plasma",begin=0,end=0.8) +
      facet_grid(year~depth_cat)
    
    ### Spatial Fields
    
    ST <- dat.samp %>% distinct(year,station,depth_cat,lon,lat,A_ID) 
    
    
    epsilon_st_2019 <- report$epsilon_st[,,1] %>% as.data.frame()
    colnames(epsilon_st_2019) <- uni.depths$depth_cat
    epsilon_st_2019 <-  epsilon_st_2019 %>% mutate(A_ID = 0:(nrow(.)-1))
    ep_A <- pivot_longer(epsilon_st_2019,values_to="val",names_to="depth_cat",cols = -A_ID)
    ep_A$pred_year <- 2019
    
    epsilon_st_2021 <- report$epsilon_st[,,2] %>% as.data.frame()
    colnames(epsilon_st_2021) <- uni.depths$depth_cat
    epsilon_st_2021 <- epsilon_st_2021 %>% mutate(A_ID = 0:(nrow(.)-1))
    ep_B <- pivot_longer(epsilon_st_2021,values_to="val",names_to="depth_cat",cols = -A_ID)
    ep_B$pred_year <- 2021
    
    ep_all <- bind_rows(ep_A,ep_B)
    
    ST <- ST %>% left_join(.,ep_all%>%mutate(depth_cat=as.numeric(depth_cat))) 
    
    p_spatial_eff <- base_map_trim + 
      geom_point(data=ST,aes(x=lon,y=lat,fill=val,color=val),shape=21,alpha=0.5,size=3) +
      #scale_size("Copies / uL",labels=LAB,breaks=LAB,range=c(0.2,20),limits=c(0.5,NA))+
      # geom_point(data=SAMPLES%>% filter(Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
      # scale_shape("Copies / uL",solid=FALSE)+
      
      scale_color_viridis_c("Spatial Fields",option="plasma",begin=0,end=0.8) +
      scale_fill_viridis_c("Spatial Fields",option="plasma",begin=0,end=0.8) +
      facet_grid(year~depth_cat)
 
    p_spatial_eff
 
 
 
