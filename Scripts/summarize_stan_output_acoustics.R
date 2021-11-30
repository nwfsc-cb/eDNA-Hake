#### Combine and Post-process stan results.
#####################################################3
#####################################################3
#####################################################3
#####################################################3
# Pull files from the list
pars <- Output.acoustic$samp
dat.acoustic.bin <- Output.acoustic$dat.acoustic.bin
dat.acoustic.pos <- Output.acoustic$dat.acoustic.pos

MODEL.TYPE <- Output.acoustic$MODEL.TYPE
MODEL.VAR <- Output.acoustic$MODEL.VAR
MODEL.ID <- Output.acoustic$MODEL.ID
#dat_raster_fin <- Output.acoustic$dat_raster_fin
TRACE <- Output.acoustic$TRACE
#dat_raster_fin <- Output.acoustic$dat_raster_fin

brms.object.pos <- Output.acoustic$brms.object.pos
brms.object.bin <- Output.acoustic$brms.object.bin

N.knots.lon.pos <- Output.acoustic$N_knots_all$N.knots.lon.pos
N.knots.lat.pos <- Output.acoustic$N_knots_all$N.knots.lat.pos
N.knots.lon.bin <- Output.acoustic$N_knots_all$N.knots.lon.bin
N.knots.lat.bin <- Output.acoustic$N_knots_all$N.knots.lat.bin
if(MODEL.TYPE =="lat.long.smooth"){
  N.knots.bd <- Output.acoustic$N_knots_all$N.knots.bd
}

# N posteriors to sample
N.POST <- 4000

### PLOTS OF MODEL ESTIMATES
#pairs(stanMod, pars = c(base_params), log = FALSE, las = 1)
Intercept_bin <- mean(pars$Intercept_bin)
Intercept_pos <- mean(pars$Intercept_pos)
sigma <-  mean(pars$sigma)


################################################################################################
################################################################################################
################################################################################################
################################################################################################
##### Extract data of interest, save to file for use elsewhere
################################################################################################
################################################################################################
################################################################################################
###############################################################################################

# Pres-Abs and Positive parts of the model
# Based on sampling volume of 2.5 L

PROBS <- c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975)
Log   <- paste0("log",PROBS)

theta_bin_out <- data.frame(trans_id= dat.acoustic.bin$trans_id, 
                        Mean.logit=apply(pars$theta_bin_pred,2,mean),
                        Sd.logit=apply(pars$theta_bin_pred,2,sd),
                        Logit.Val=data.frame(t(apply(pars$theta_bin_pred,2,quantile,probs=PROBS))),
                        Mean=apply(inv.logit(pars$theta_bin_pred),2,mean),
                        Sd=apply(inv.logit(pars$theta_bin_pred),2,sd),
                        Val=data.frame(t(apply(inv.logit(pars$theta_bin_pred),2,quantile,probs=PROBS))))
D_pos_out <- data.frame(trans_id= dat.acoustic.pos$trans_id, 
                                Mean.log=apply(pars$D_pos_pred,2,mean),
                                Sd.log=apply(pars$D_pos_pred,2,sd),
                                Log.Val=data.frame(t(apply(pars$D_pos_pred,2,quantile,probs=PROBS))),
                                Mean=apply(exp(pars$D_pos_pred),2,mean),
                                Sd=apply(exp(pars$D_pos_pred),2,sd),
                                Val=data.frame(t(apply(exp(pars$D_pos_pred),2,quantile,probs=PROBS))))

##############################
## Project Surfaces for lat.long.smooth and lat.long.smooth.base MODEL.TYPE
##############################
smooth.projections <- list()
if(MODEL.TYPE =="lat.long.smooth" | MODEL.TYPE=="lat.long.smooth.base"){
  
  # make a projection matrix for the new data
  new_data  <- dat_raster_fin  %>% rename(utm.lon=x,utm.lat=y,bathy.bottom.depth=depth_m) %>%
                    filter(bathy.bottom.depth > 40)
  
  if(MODEL.TYPE=="lat.long.smooth.base"){
    new_data_trim <- new_data %>% mutate(Y=rnorm(nrow(new_data))) %>% 
                        dplyr::select(-lon,-lat,-Gridcell_ID,-WM_depth_m,-bathy.bottom.depth)
  }else if(MODEL.TYPE=="lat.long.smooth"){
    new_data_trim <- new_data %>% mutate(Y=rnorm(nrow(new_data))) %>% 
                        dplyr::select(-lon,-lat,-Gridcell_ID,-WM_depth_m)
  }
  # create the basis function math needed to project from the posterior
  smooth.dat.pred.bin <-  standata(brms.object.bin,newdata = new_data_trim)
  smooth.dat.pred.pos <-  standata(brms.object.pos,newdata = new_data_trim)
  
  ### Extract the relevant items from smooth.dat.pred, make predictions for posterior samples.
  N.ID = seq(1,length(pars$Intercept_bin),length.out=N.POST) %>% round()
  D_pred_bin <- matrix(0,nrow(smooth.dat.pred.bin$Xs),N.POST)
  D_pred_pos <- D_pred_bin
  D_pred_bottom_depth_bin <- D_pred_bin
  D_pred_bottom_depth_pos <- D_pred_bin
  D_smooth_bin = D_pred_bin
  D_smooth_pos = D_pred_pos
  ID <- data.frame(id=1:ncol(D_pred_bin),N.ID)
  sigma.vec <- pars$sigma[N.ID]
  
  for(i in N.ID){
    D_pred_bin[,ID$id[ID$N.ID==i]] <- 
      pars$Intercept_bin[i] +
      smooth.dat.pred.bin$Xs %*% pars$bs_bin[i,] +
      smooth.dat.pred.bin$Zs_1_1 %*% pars$s_1_1_bin[i,] + 
      smooth.dat.pred.bin$Zs_1_2 %*% pars$s_1_2_bin[i,] + 
      smooth.dat.pred.bin$Zs_1_3 %*% pars$s_1_3_bin[i,] 
   
    D_pred_pos[,ID$id[ID$N.ID==i]] <- 
      pars$Intercept_pos[i] +
      smooth.dat.pred.pos$Xs %*% pars$bs_pos[i,] +
      smooth.dat.pred.pos$Zs_1_1 %*% pars$s_1_1_pos[i,] + 
      smooth.dat.pred.pos$Zs_1_2 %*% pars$s_1_2_pos[i,] + 
      smooth.dat.pred.pos$Zs_1_3 %*% pars$s_1_3_pos[i,]
    
    if(MODEL.TYPE=="lat.long.smooth"){
      D_pred_bin[,ID$id[ID$N.ID==i]] <- D_pred_bin[,ID$id[ID$N.ID==i]] + smooth.dat.pred.bin$Zs_2_1 %*% pars$s_2_1_bin[i,]
      D_pred_pos[,ID$id[ID$N.ID==i]] <- D_pred_pos[,ID$id[ID$N.ID==i]] + smooth.dat.pred.pos$Zs_2_1 %*% pars$s_2_1_pos[i,]
      
      # Calculate the conditional effect for bottom.depth
      D_pred_bottom_depth_bin[,ID$id[ID$N.ID==i]] <-
        smooth.dat.pred.bin$Xs[,ncol(smooth.dat.pred.bin$Xs)] * pars$bs_bin[i,ncol(smooth.dat.pred.bin$Xs)] +
        smooth.dat.pred.bin$Zs_2_1 %*% pars$s_2_1_bin[i,]
      D_pred_bottom_depth_pos[,ID$id[ID$N.ID==i]] <-
        smooth.dat.pred.pos$Xs[,ncol(smooth.dat.pred.pos$Xs)] * pars$bs_pos[i,ncol(smooth.dat.pred.pos$Xs)] +
        smooth.dat.pred.pos$Zs_2_1 %*% pars$s_2_1_pos[i,]
    }

    D_smooth_bin[,ID$id[ID$N.ID==i]] <- smooth.dat.pred.bin$Xs[,(1:ncol(smooth.dat.pred.bin$Xs)-1)] %*% pars$bs_bin[i,(1:ncol(smooth.dat.pred.bin$Xs)-1)] +
                                        smooth.dat.pred.bin$Zs_1_1 %*% pars$s_1_1_bin[i,] + 
                                        smooth.dat.pred.bin$Zs_1_2 %*% pars$s_1_2_bin[i,] + 
                                        smooth.dat.pred.bin$Zs_1_3 %*% pars$s_1_3_bin[i,]
    
    D_smooth_pos[,ID$id[ID$N.ID==i]] <- smooth.dat.pred.pos$Xs[,(1:ncol(smooth.dat.pred.pos$Xs)-1)] %*% pars$bs_pos[i,(1:ncol(smooth.dat.pred.pos$Xs)-1)] +
                                        smooth.dat.pred.pos$Zs_1_1 %*% pars$s_1_1_pos[i,] + 
                                        smooth.dat.pred.pos$Zs_1_2 %*% pars$s_1_2_pos[i,] + 
                                        smooth.dat.pred.pos$Zs_1_3 %*% pars$s_1_3_pos[i,]
  }    
  
  D_pred_bin_logit_summary <- data.frame(Mean=rowMeans(D_pred_bin),
                                    SD=apply(D_pred_bin,1,sd),
                                    Val=apply((D_pred_bin),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  D_pred_bin_summary <- data.frame(Mean=rowMeans(inv.logit(D_pred_bin)),
                                    SD=apply(inv.logit(D_pred_bin),1,sd),
                                    Val=apply(inv.logit(D_pred_bin),1,quantile,probs=PROBS) %>% t() %>% as.matrix())

  D_pred_pos_log_summary <- data.frame(Mean=rowMeans(D_pred_pos),
                                    SD=apply(D_pred_pos,1,sd),
                                    Val=apply((D_pred_pos),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  D_pred_pos_summary <- data.frame(Mean=rowMeans(exp(D_pred_pos)),
                                    SD=apply(exp(D_pred_pos),1,sd),
                                    Val=apply(exp(D_pred_pos),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  
  if(MODEL.TYPE=="lat.long.smooth"){
    D_pred_bottom_depth_bin_summary <- data.frame(Mean=rowMeans(D_pred_bottom_depth_bin),
                                            SD=apply((D_pred_bottom_depth_bin),1,sd),
                                            Val=apply((D_pred_bottom_depth_bin),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
    D_pred_bottom_depth_pos_summary <- data.frame(Mean=rowMeans(D_pred_bottom_depth_pos),
                                            SD=apply((D_pred_bottom_depth_pos),1,sd),
                                            Val=apply((D_pred_bottom_depth_pos),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  }
  
  D_pred_smooth_bin <- data.frame(Mean=rowMeans(D_smooth_bin),
                              SD=apply(D_smooth_bin,1,sd),
                              Val=apply(D_smooth_bin,1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  D_pred_smooth_pos <- data.frame(Mean=rowMeans(D_smooth_pos),
                                  SD=apply(D_smooth_pos,1,sd),
                                  Val=apply(D_smooth_pos,1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  
  ## Make unconditional predictions
  KM_2 <- 25 # each grid cell is 25km^2
  
   D_pred_pos_real <- D_pred_pos * 0
   for(i in 1:N.POST){
      D_pred_pos_real[,i] <- #exp(D_pred_pos[,i]) 
            rlnorm(nrow(D_pred_pos), D_pred_pos[,i] - 0.5*pars$sigma[N.ID[i]]^2, pars$sigma[N.ID[i]])
    }
   D_pred_uncond <- inv.logit(D_pred_bin) * D_pred_pos_real
  
   D_pred_uncond_summary <- data.frame(Mean=rowMeans(D_pred_uncond),
                                       SD=apply(D_pred_uncond,1,sd),
                                       Val=apply(D_pred_uncond,1,quantile,probs=PROBS) %>% t() %>% as.matrix())
   D_pred_uncond_mt_summary <- data.frame(Mean=rowMeans(D_pred_uncond*KM_2),
                                       SD=apply(D_pred_uncond*KM_2,1,sd),
                                       Val=apply(D_pred_uncond*KM_2,1,quantile,probs=PROBS) %>% t() %>% as.matrix())
   
   D_pred_uncond_combined <- bind_cols(D_pred_uncond_summary,new_data_trim %>% dplyr::select(-Y)) %>%
     left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID)) 
   D_pred_uncond_mt_combined <- bind_cols(D_pred_uncond_mt_summary,new_data_trim %>% dplyr::select(-Y)) %>%
     left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID)) 
   
   # Group unconditional predictions into latitudinal bins.
   d.all <-  bind_cols(as.data.frame(D_pred_uncond*KM_2),new_data_trim %>% dplyr::select(-Y)) %>%
     left_join(.,new_data %>% dplyr::select(utm.lon,utm.lat,Gridcell_ID)) %>%
     dplyr::select(-utm.lon,-utm.lat,-bathy.bottom.depth)
   
   #make long-form
   d.all <- d.all %>% pivot_longer(starts_with("V"),names_to="MCMC.rep",values_to = "w.dens") %>%
     dplyr::select(Gridcell_ID, MCMC.rep, D = w.dens)
   
   # Merge in the different ID categories for dividing the coast up into spatial strata.
   # I'm too tired to do this is plyr so here are some loops.
   new_dat_temp <- new_data %>% dplyr::select(Gridcell_ID,lon,lat) %>% distinct()
   
   new_dat_strata <- NULL
   for(i in 1:nrow(lat.breaks$lats.equal)){
     temp <- new_dat_temp %>% 
       filter(lat < lat.breaks$lats.equal$lat.max[i], lat > lat.breaks$lats.equal$lat[i] ) %>%
       mutate( ID.equal = lat.breaks$lats.equal$ID[i])
     new_dat_strata <- bind_rows(new_dat_strata,temp)
   }
   
   new_dat_temp <- new_dat_strata
   new_dat_strata <- NULL
   for(i in 1:nrow(lat.breaks$lats.rounded.0.5)){
     temp <- new_dat_temp %>% 
       filter(lat < lat.breaks$lats.rounded.0.5$lat.max[i], lat > lat.breaks$lats.rounded.0.5$lat[i] ) %>%
       mutate( ID.lat.0.5 = lat.breaks$lats.rounded.0.5$ID[i])
     new_dat_strata <- bind_rows(new_dat_strata,temp)
   }
   
   new_dat_temp <- new_dat_strata
   new_dat_strata <- NULL
   for(i in 1:nrow(lat.breaks$lats.rounded.1.0)){
     temp <- new_dat_temp %>% 
       filter(lat < lat.breaks$lats.rounded.1.0$lat.max[i], lat > lat.breaks$lats.rounded.1.0$lat[i] ) %>%
       mutate( ID.lat.1.0 = lat.breaks$lats.rounded.1.0$ID[i])
     new_dat_strata <- bind_rows(new_dat_strata,temp)
   }
   #new_dat_strata is the data.frame of interest
   
   # Merge new_dat_strata into d.all 
   d.all <- left_join(d.all,new_dat_strata)
   
   D_acoustic_uncond_mt <- d.all %>% group_by(Gridcell_ID,MCMC.rep) %>% 
     summarise(tot = sum(D)) %>% group_by(Gridcell_ID) %>%
     summarise(Mean=mean(tot), # These are summaries across MCMC for each cell.
               Median=median(tot),
               SD=sd(tot),
               Q.0.01 = quantile(tot,probs=c(0.01)),
               Q.0.025 = quantile(tot,probs=c(0.025)),
               Q.0.05 = quantile(tot,probs=c(0.05)),
               Q.0.25 = quantile(tot,probs=c(0.25)),
               Q.0.75 = quantile(tot,probs=c(0.75)),
               Q.0.95 = quantile(tot,probs=c(0.95)),
               Q.0.975 = quantile(tot,probs=c(0.975)),
               Q.0.99 = quantile(tot,probs=c(0.99)))
   
   D_acoustic_uncond_lat_equal <- d.all %>% group_by(MCMC.rep,ID.equal) %>% 
                                    summarise(tot = sum(D)) %>% group_by(ID.equal) %>%
           summarise(Mean=mean(tot), # These are summaries across MCMC for each cell.
                 Median=median(tot),
                  SD=sd(tot),
                  Q.0.01 = quantile(tot,probs=c(0.01)),
                  Q.0.025 = quantile(tot,probs=c(0.025)),
                  Q.0.05 = quantile(tot,probs=c(0.05)),
                  Q.0.25 = quantile(tot,probs=c(0.25)),
                  Q.0.75 = quantile(tot,probs=c(0.75)),
                  Q.0.95 = quantile(tot,probs=c(0.95)),
                  Q.0.975 = quantile(tot,probs=c(0.975)),
                  Q.0.99 = quantile(tot,probs=c(0.99)))
   
   D_acoustic_uncond_lat_1.0 <- d.all %>% group_by(MCMC.rep,ID.lat.1.0) %>% 
                                summarise(tot = sum(D)) %>% group_by(ID.lat.1.0) %>%
          summarise(Mean=mean(tot), # These are summaries across MCMC for each cell.
               Median=median(tot),
               SD=sd(tot),
               Q.0.01 = quantile(tot,probs=c(0.01)),
               Q.0.025 = quantile(tot,probs=c(0.025)),
               Q.0.05 = quantile(tot,probs=c(0.05)),
               Q.0.25 = quantile(tot,probs=c(0.25)),
               Q.0.75 = quantile(tot,probs=c(0.75)),
               Q.0.95 = quantile(tot,probs=c(0.95)),
               Q.0.975 = quantile(tot,probs=c(0.975)),
               Q.0.99 = quantile(tot,probs=c(0.99)))
   
   D_acoustic_uncond_lat_0.5 <- d.all %>% group_by(MCMC.rep,ID.lat.0.5) %>% 
                                summarise(tot = sum(D)) %>% group_by(ID.lat.0.5) %>%
          summarise(Mean=mean(tot), # These are summaries across MCMC for each cell.
               Median=median(tot),
               SD=sd(tot),
               Q.0.01 = quantile(tot,probs=c(0.01)),
               Q.0.025 = quantile(tot,probs=c(0.025)),
               Q.0.05 = quantile(tot,probs=c(0.05)),
               Q.0.25 = quantile(tot,probs=c(0.25)),
               Q.0.75 = quantile(tot,probs=c(0.75)),
               Q.0.95 = quantile(tot,probs=c(0.95)),
               Q.0.975 = quantile(tot,probs=c(0.975)),
               Q.0.99 = quantile(tot,probs=c(0.99)))
   
   D_acoustic_uncond_total_mt <- d.all %>% group_by(Gridcell_ID,MCMC.rep) %>% 
     summarise(tot = sum(D)) %>% 
     left_join(.,new_dat_strata %>% dplyr::select(Gridcell_ID,lat)) %>%
     arrange(lat,MCMC.rep) %>%
     group_by(MCMC.rep) %>%
     mutate(cum_sum = cumsum(tot),max_cum_sum = max(cum_sum)) %>%
     mutate(cum_sum_prob = cum_sum / max_cum_sum)
   
  D_acoustic_uncond_cum_sum <-   D_acoustic_uncond_total_mt %>%    
     group_by(lat) %>%
     summarise(Mean=mean(cum_sum_prob), # These are summaries across MCMC for each cell.
               Median=median(cum_sum_prob),
               SD=sd(cum_sum_prob),
               Q.0.01 = quantile(cum_sum_prob,probs=c(0.01)),
               Q.0.025 = quantile(cum_sum_prob,probs=c(0.025)),
               Q.0.05 = quantile(cum_sum_prob,probs=c(0.05)),
               Q.0.25 = quantile(cum_sum_prob,probs=c(0.25)),
               Q.0.75 = quantile(cum_sum_prob,probs=c(0.75)),
               Q.0.95 = quantile(cum_sum_prob,probs=c(0.95)),
               Q.0.975 = quantile(cum_sum_prob,probs=c(0.975)),
               Q.0.99 = quantile(cum_sum_prob,probs=c(0.99)))
   
   D_acoustic_uncond_total_mt <- D_acoustic_uncond_total_mt %>% ungroup() %>%
      summarise(Mean=mean(max_cum_sum), # These are summaries across MCMC for each cell.
               Median=median(max_cum_sum),
               SD=sd(max_cum_sum),
               Q.0.01 = quantile(max_cum_sum,probs=c(0.01)),
               Q.0.025 = quantile(max_cum_sum,probs=c(0.025)),
               Q.0.05 = quantile(max_cum_sum,probs=c(0.05)),
               Q.0.25 = quantile(max_cum_sum,probs=c(0.25)),
               Q.0.75 = quantile(max_cum_sum,probs=c(0.75)),
               Q.0.95 = quantile(max_cum_sum,probs=c(0.95)),
               Q.0.975 = quantile(max_cum_sum,probs=c(0.975)),
               Q.0.99 = quantile(max_cum_sum,probs=c(0.99)))
   
   # Summarize but maintain MCMC to do uncertainty in the correlation between eDNA and acoustics.
   
   D_1.0_uncond_resample <- d.all %>% group_by(MCMC.rep,ID.lat.1.0) %>% 
                              summarise(tot = sum(D))
   
   D_0.5_uncond_resample <- d.all %>% group_by(MCMC.rep,ID.lat.0.5) %>% 
                              summarise(tot = sum(D))
   
   D_grid.cell_uncond_resample <- d.all %>% dplyr::select(Gridcell_ID,MCMC.rep,D)
   
   
   
   #####
   
  D_pred_bin_combined <- bind_cols(D_pred_bin_summary,new_data_trim %>% dplyr::select(-Y)) %>%
              left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID))
  D_pred_pos_log_combined <- bind_cols(D_pred_pos_log_summary,new_data_trim %>% dplyr::select(-Y)) %>%
              left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID))
  D_pred_pos_combined <- bind_cols(D_pred_pos_summary,new_data_trim %>% dplyr::select(-Y)) %>%
              left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID))

  D_pred_smooth_bin_combined <- bind_cols(D_pred_smooth_bin,new_data_trim %>% dplyr::select(-Y)) %>%
    left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID))
  D_pred_smooth_pos_combined <- bind_cols(D_pred_smooth_pos,new_data_trim %>% dplyr::select(-Y)) %>%
    left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID))
  
  # 
  D_pred_bottom_depth_bin_combined <- bind_cols(D_pred_bottom_depth_bin_summary,new_data_trim %>% dplyr::select(-Y)) %>%
      left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID))
  D_pred_bottom_depth_pos_combined <- bind_cols(D_pred_bottom_depth_pos_summary,new_data_trim %>% dplyr::select(-Y)) %>%
    left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID))
  # 
  # 
  if(MODEL.TYPE=="lat.long.smooth"){
    smooth.projections <- list(
                              D_acoustic_uncond_mt = D_acoustic_uncond_mt,
                              D_acoustic_uncond_mt_lat_equal = D_acoustic_uncond_lat_equal,
                              D_acoustic_uncond_mt_lat_0.5 = D_acoustic_uncond_lat_0.5,                              
                              D_acoustic_uncond_mt_lat_1.0 = D_acoustic_uncond_lat_1.0,
                               D_pred_uncond_combined = D_pred_uncond_combined,
                               D_pred_uncond_mt_combined =  D_pred_uncond_mt_combined,
                               D_pred_bin_combined = D_pred_bin_combined,
                               D_pred_pos_combined = D_pred_pos_combined,
                               D_pred_bottom_depth_bin_combined = D_pred_bottom_depth_bin_combined,
                               D_pred_bottom_depth_pos_combined = D_pred_bottom_depth_pos_combined,
                               D_pred_smooth_bin_combined = D_pred_smooth_bin_combined,
                               D_pred_smooth_pos_combined = D_pred_smooth_pos_combined
                               )
  }else if(MODEL.TYPE=="lat.long.smooth.base"){
    smooth.projections <- list(
                              D_acoustic_uncond_mt = D_acoustic_uncond_mt,
                              D_acoustic_uncond_mt_lat_equal = D_acoustic_uncond_lat_equal,
                              D_acoustic_uncond_mt_lat_0.5 = D_acoustic_uncond_lat_0.5,                              
                              D_acoustic_uncond_mt_lat_1.0 = D_acoustic_uncond_lat_1.0,
                              D_pred_uncond_combined = D_pred_uncond_combined,
                              D_pred_uncond_mt_combined =  D_pred_uncond_mt_combined,
                              D_pred_bin_combined = D_pred_bin_combined,
                              D_pred_pos_combined = D_pred_pos_combined,
                              D_pred_smooth_bin_combined = D_pred_smooth_bin_combined,
                              D_pred_smooth_pos_combined = D_pred_smooth_pos_combined
    )
  }
}

# Calculation for Observed-Predicted and Residual
pred_obs_bin <- bind_cols(dat.acoustic.bin %>% 
                           dplyr::select(transect,lat,lon,utm.lat,utm.lon,bin_weight_dens,weight_dens_mt_km2,bathy.bottom.depth),
                          theta_bin_out)

pred_obs_pos <- bind_cols(dat.acoustic.pos %>% 
                            dplyr::select(transect,lat,lon,utm.lat,utm.lon,bin_weight_dens,weight_dens_mt_km2,bathy.bottom.depth),
                          D_pos_out) %>%
                          mutate(resid = weight_dens_mt_km2 - Mean) #log.resid=log(resid))

###  
# Combine the necessary data.frames into a list for use later.
Output.summary <- list(

  D_pos_out = D_pos_out,
  theta_bin_out = theta_bin_out,
  pred_obs_bin =pred_obs_bin,
  pred_obs_pos =pred_obs_pos,

  # Projections or projection helpers.
  N.POST = N.POST, # number of posterior samples used.
  dat_raster_fin = dat_raster_fin, 
  smooth.projections = smooth.projections,

  D_pred_uncond = D_pred_uncond 
)
 
