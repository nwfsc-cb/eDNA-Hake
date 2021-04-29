#### Combine and Post-process stan results.
#####################################################3
#####################################################3
#####################################################3
#####################################################3
# Pull files from the list
pars <- Output.qpcr$samp
dat.stand.bin <- Output.qpcr$dat.stand.bin
dat.stand.pos <- Output.qpcr$dat.stand.pos
dat.id <- Output.qpcr$dat.id
PCR <- Output.qpcr$PCR
STATION.DEPTH <- Output.qpcr$STATION.DEPTH
SAMPLES <- Output.qpcr$SAMPLES
MODEL.TYPE <- Output.qpcr$MODEL.TYPE
MODEL.VAR <- Output.qpcr$MODEL.VAR
dat_raster_fin <- Output.qpcr$dat_raster_fin
depth.fact <- sort(unique(STATION.DEPTH$depth_cat_factor))
TRACE <- Output.qpcr$TRACE

brms.object <- Output.qpcr$brms.object

N.knots.lon <- Output.qpcr$N_knots_all[[1]]
N.knots.lat <- Output.qpcr$N_knots_all[[2]]
if(MODEL.TYPE =="lat.long.smooth"){
  N.knots.bd <- Output.qpcr$N_knots_all[[3]]
}

# N posteriors to sample
N.POST <- 2000

### PLOTS OF STANDARDS and associated REGRESSIONS
#pairs(stanMod, pars = c(base_params), log = FALSE, las = 1)
B1 <- apply(pars$beta_0,2,mean)
B2 <- apply(pars$beta_1,2,mean)

P0 <- apply(pars$phi_0,2,mean)
P1 <- apply(pars$phi_1,2,mean)

V0 <-  mean(pars$sigma_stand_int)
#V1 <- mean(pars$sigma_stand_slope)
#V2 <- mean(pars$sigma_stand_slope2)

# Plot regression against Standard
X <- seq(0,5,length.out=1000) 
Y <- t(B1 + B2 %*% t(X ))

STAND.REG <- data.frame(X=X,Y=Y)
STAND.REG <- melt(STAND.REG,id.vars="X",value.name="Y")

x.lim=c(min(X),max(X))
y.lim=c(20,40)
for(i in 1:ncol(Y)){
  plot(Y[,i]~X,xlim=x.lim,ylim=y.lim,type="l",col=2)
  par(new=T)
}
plot(dat.stand.pos$Ct~ dat.stand.pos$log10_copies,xlim=x.lim,ylim=y.lim)

lab.temp<- data.frame(variable=paste0("Y.",PCR$plate_idx), qPCR=PCR$qPCR)
STAND.REG <-left_join(STAND.REG,lab.temp)

stand.plot <- ggplot(dat.stand.pos) +
  geom_point(aes(x=log10_copies ,y=Ct,color=qPCR),alpha=0.75) +
  scale_shape_discrete(name="qPCR Plate") +
  theme_bw()
stand.plot <- stand.plot +
  geom_line(data=STAND.REG,aes(x=X,y=Y,color=qPCR)) +
  scale_color_discrete(name="qPCR Plate") +
  ylab("PCR cycle")  +
  xlab("log10 copies DNA") +
  facet_wrap(~qPCR)
#scale_x_continuous(labels = paste0("1e",LABS),limits = c(-2,4))
stand.plot

# Plot occurrence of standard
Y <- t(P0 + P1 %*% t(X ))
LOGIT <- data.frame(X=X,Y=plogis(Y))
LOGIT <- melt(LOGIT,id.vars="X",value.name="Y")

LOGIT <-left_join(LOGIT,lab.temp)

stand.plot.pres <- ggplot(dat.stand.bin) +
  geom_jitter(aes(x=log10_copies ,y=Ct_bin,color=qPCR),alpha=0.75,width=0,height=0.05) +
  geom_line(data=LOGIT,aes(y=Y,x=X,color=qPCR)) +
  theme_bw() +
  scale_shape_discrete(name="qPCR Plate") +
  scale_color_discrete(name="qPCR plate") +
  xlab("log10 copies DNA") +
  ylab("Amplification success") +
  scale_y_continuous(breaks=c(0,1),labels = c("No","Yes"))
#scale_x_continuous(breaks=BREAKS,labels = paste0("1e",LABS),limits = c(-2.5,4))

stand.plot.pres

################################################################################################
################################################################################################
################################################################################################
################################################################################################
##### Extract data of interest, save to file for use elsewhere
################################################################################################
################################################################################################
################################################################################################
###############################################################################################

# THIS IS THE FACTOR FOR EXPANDING FROM COPIES PER ul to COPIES PER L 
# Based on sampling volume of 2.5 L

EXPAND <- 40
LOG.EXPAND <- log10(EXPAND)

PROBS <- c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975)
Log   <- paste0("log",PROBS)
station_depth_out <- data.frame(station_depth_idx= 1:ncol(pars$D), 
                                Mean.log=apply(pars$D,2,mean),
                                Sd.log=apply(pars$D,2,sd),
                                Log.Val=data.frame(t(apply(pars$D,2,quantile,probs=PROBS))),
                                Mean=apply(10^pars$D,2,mean),
                                Sd=apply(10^pars$D,2,sd),
                                Val=data.frame(t(apply(10^pars$D,2,quantile,probs=PROBS))))

station_depth_out_liter <- data.frame(station_depth_idx= 1:ncol(pars$D), 
                                      Mean.log=apply(pars$D+LOG.EXPAND,2,mean),
                                      Sd.log=apply(pars$D+LOG.EXPAND,2,sd),
                                      Log.Val=data.frame(t(apply(pars$D+LOG.EXPAND,2,quantile,probs=PROBS))),
                                      Mean=apply(10^(pars$D+LOG.EXPAND),2,mean),
                                      Sd=apply(10^(pars$D+LOG.EXPAND),2,sd),
                                      Val=data.frame(t(apply(10^(pars$D+LOG.EXPAND),2,quantile,probs=PROBS))))

D_delta_out <- data.frame(sample_idx= 1:ncol(pars$D_delta), 
                          Mean.log=apply(pars$D_delta,2,mean),
                          Sd.log=apply(pars$D_delta,2,sd),
                          Log.Val=data.frame(t(apply(pars$D_delta,2,quantile,probs=PROBS))),
                          Mean=apply(10^pars$D_delta,2,mean),
                          Sd=apply(10^pars$D_delta,2,sd),
                          Val=data.frame(t(apply(10^pars$D_delta,2,quantile,probs=PROBS))))

D_delta_out_liter <- data.frame(sample_idx= 1:ncol(pars$D_delta), 
                                Mean.log=apply(pars$D_delta+LOG.EXPAND,2,mean),
                                Sd.log=apply(pars$D_delta+LOG.EXPAND,2,sd),
                                Log.Val=data.frame(t(apply(pars$D_delta+LOG.EXPAND,2,quantile,probs=PROBS))),
                                Mean=apply(10^(pars$D_delta+LOG.EXPAND),2,mean),
                                Sd=apply(10^(pars$D_delta+LOG.EXPAND),2,sd),
                                Val=data.frame(t(apply(10^(pars$D_delta+LOG.EXPAND),2,quantile,probs=PROBS))))

sample_contam_total_out <- data.frame(sample_idx= 1:ncol(pars$D_contam), 
                                      Mean.log=apply(pars$D_contam,2,mean),
                                      Sd.log=apply(pars$D_contam,2,sd),
                                      Log.Val=data.frame(t(apply(pars$D_contam,2,quantile,probs=PROBS))),
                                      Mean=apply(10^pars$D_contam,2,mean),
                                      Sd=apply(10^pars$D_contam,2,sd),
                                      Val=data.frame(t(apply(10^pars$D_contam,2,quantile,probs=PROBS))))

sample_contam_total_out_liter <- data.frame(sample_idx= 1:ncol(pars$D_contam), 
                                            Mean.log=apply(pars$D_contam+LOG.EXPAND,2,mean),
                                            Sd.log=apply(pars$D_contam+LOG.EXPAND,2,sd),
                                            Log.Val=data.frame(t(apply(pars$D_contam+LOG.EXPAND,2,quantile,probs=PROBS))),
                                            Mean=apply(10^(pars$D_contam+LOG.EXPAND),2,mean),
                                            Sd=apply(10^(pars$D_contam+LOG.EXPAND),2,sd),
                                            Val=data.frame(t(apply(10^(pars$D_contam+LOG.EXPAND),2,quantile,probs=PROBS))))

delta_out <- data.frame(sample_idx= 1:ncol(pars$delta), 
                        Mean.log=apply(pars$delta,2,mean) ,
                        Sd.log=apply(pars$delta,2,sd),
                        Log.Val=data.frame(t(apply(pars$delta,2,quantile,probs=PROBS))),
                        Mean=apply(10^pars$delta,2,mean),
                        Sd=apply(10^pars$delta,2,sd),
                        Val=data.frame(t(apply(10^pars$delta,2,quantile,probs=PROBS))))

field_neg_out   <-      data.frame(sample_control_idx= 1:ncol(pars$D_control), 
                                   Mean.log=apply(pars$D_control,2,mean),
                                   Sd.log=apply(pars$D_control,2,sd),
                                   Log.Val=data.frame(t(apply(pars$D_control,2,quantile,probs=PROBS))),
                                   Mean=apply(10^pars$D_control,2,mean),
                                   Sd=apply(10^pars$D_control,2,sd),
                                   Val=data.frame(t(apply(10^pars$D_control,2,quantile,probs=PROBS))))

field_neg_out_liter <-  data.frame(sample_control_idx= 1:ncol(pars$D_control), 
                                   Mean.log=apply(pars$D_control+LOG.EXPAND,2,mean),
                                   Sd.log=apply(pars$D_control+LOG.EXPAND,2,sd),
                                   Log.Val=data.frame(t(apply(pars$D_control+LOG.EXPAND,2,quantile,probs=PROBS))),
                                   Mean=apply(10^(pars$D_control+LOG.EXPAND),2,mean),
                                   Sd=apply(10^(pars$D_control+LOG.EXPAND),2,sd),
                                   Val=data.frame(t(apply(10^(pars$D_control+LOG.EXPAND),2,quantile,probs=PROBS))))


mu_contam_out <- data.frame(type = c("raw","liter"),
                            Mean = c(mean(pars$mu_contam),mean(pars$mu_contam+log(EXPAND))),
                            Sd = sd(pars$mu_contam))

sigma_contam_out <- data.frame(Mean = mean(pars$sigma_contam), Sd = sd(pars$sigma_contam))

##############################
### Depth factor fixed effects
##############################
depth_id<- levels(STATION.DEPTH$depth_cat_factor) %>% as.numeric(as.character())
      depth_fact_raw <- cbind(pars$b[,1], pars$b[,1]+pars$b[,2:length(depth_id)]) %>% as.data.frame()
      depth_fact_summary <- data.frame(depth_id = depth_id, 
                                          Mean.log=apply(depth_fact_raw+LOG.EXPAND,2,mean),
                                          Sd.log=apply(depth_fact_raw+LOG.EXPAND,2,sd),
                                          Log.Val=data.frame(t(apply(depth_fact_raw+LOG.EXPAND,2,quantile,probs=PROBS))),
                                          Mean=apply(10^(depth_fact_raw+LOG.EXPAND),2,mean),
                                          Sd=apply(10^(depth_fact_raw+LOG.EXPAND),2,sd),
                                          Val=data.frame(t(apply(10^(depth_fact_raw+LOG.EXPAND),2,quantile,probs=PROBS))))

##############################
## Project Surfaces for lat.long.smooth and lat.long.smooth.base MODEL.TYPE
##############################

smooth.projections <- list()
if(MODEL.TYPE =="lat.long.smooth" | MODEL.TYPE=="lat.long.smooth.base"){
  
  # make a projection matrix
  #depth.fact <- sort(rep(STATION.DEPTH$depth_cat_factor %>% unique() %>% sort(),nrow(dat_raster_fin)))
  d_cat <- sort(unique(STATION.DEPTH$depth_cat))
  n.rep <- length(unique(depth.fact))
  
  new_data <- NULL
  orig_data <- NULL
  for(i in 1:n.rep){
    temp <- dat_raster_fin %>% filter(depth_m > (d_cat[i]-10))
    temp$depth_cat_factor <- d_cat[i]
    new_data <- rbind(new_data,temp)
  }
  #new_data$depth_cat_factor <- depth.fact
  new_data  <- new_data %>% rename(utm.lon=x,utm.lat=y,bottom.depth.consensus=depth_m)
  orig_data <- STATION.DEPTH %>% dplyr::select(utm.lon,utm.lat,bottom.depth.consensus,depth_cat_factor,lon,lat,station_depth_idx) 
                      
  if(MODEL.TYPE=="lat.long.smooth.base"){
    new_data_trim <- new_data %>% mutate(Y=rnorm(nrow(new_data))) %>% 
                        dplyr::select(-lon,-lat,-Gridcell_ID,-WM_depth_m,-bottom.depth.consensus)
    orig_data_trim <- orig_data %>% mutate(Y=rnorm(nrow(orig_data))) %>%
                          dplyr::select(-lon,-lat,-station_depth_idx,-bottom.depth.consensus)
  }else if(MODEL.TYPE=="lat.long.smooth"){
    new_data_trim <- new_data %>% mutate(Y=rnorm(nrow(new_data))) %>% 
                        dplyr::select(-lon,-lat,-Gridcell_ID,-WM_depth_m)
    orig_data_trim <- orig_data %>% mutate(Y=rnorm(nrow(orig_data))) %>% 
                        dplyr::select(-lon,-lat,-station_depth_idx)
  }
  # create the basis function math needed to project from the posterior
  smooth.dat.pred <-  standata(brms.object,newdata = new_data_trim)
  orig.dat.pred   <-  standata(brms.object,newdata = orig_data_trim)
  
  ### Extract the relevant items from smooth.dat.pred, make predictions for posterior samples.
  
  N.ID = seq(1,nrow(pars$b),length.out=N.POST) %>% round()
  D_pred <- matrix(0,nrow(smooth.dat.pred$X),N.POST)
  D_pred_bottom_depth <- D_pred
  D_smooth = D_pred
  ID <- data.frame(id=1:ncol(D_pred),N.ID)
  
  D_orig_pred <- matrix(0,nrow(orig.dat.pred$X),N.POST)
  D_orig_pred_bottom_depth <- D_pred
  D_orig_smooth = D_orig_pred
  Depth_cat_fact <- matrix(0,nrow(orig.dat.pred$X),N.POST)
  
  for(i in N.ID){
    D_pred[,ID$id[ID$N.ID==i]] <- smooth.dat.pred$X %*% pars$b[i,] +
      smooth.dat.pred$Xs %*% pars$bs[i,] +
      smooth.dat.pred$Zs_1_1 %*% pars$s_1_1[i,] + smooth.dat.pred$Zs_1_2 %*% pars$s_1_2[i,] + smooth.dat.pred$Zs_1_3 %*% pars$s_1_3[i,] + 
      smooth.dat.pred$Zs_2_1 %*% pars$s_2_1[i,] + smooth.dat.pred$Zs_2_2 %*% pars$s_2_2[i,] + smooth.dat.pred$Zs_2_3 %*% pars$s_2_3[i,] + 
      smooth.dat.pred$Zs_3_1 %*% pars$s_3_1[i,] + smooth.dat.pred$Zs_3_2 %*% pars$s_3_2[i,] + smooth.dat.pred$Zs_3_3 %*% pars$s_3_3[i,] + 
      smooth.dat.pred$Zs_4_1 %*% pars$s_4_1[i,] + smooth.dat.pred$Zs_4_2 %*% pars$s_4_2[i,] + smooth.dat.pred$Zs_4_3 %*% pars$s_4_3[i,] + 
      smooth.dat.pred$Zs_5_1 %*% pars$s_5_1[i,] + smooth.dat.pred$Zs_5_2 %*% pars$s_5_2[i,] + smooth.dat.pred$Zs_5_3 %*% pars$s_5_3[i,] + 
      smooth.dat.pred$Zs_6_1 %*% pars$s_6_1[i,] + smooth.dat.pred$Zs_6_2 %*% pars$s_6_2[i,] + smooth.dat.pred$Zs_6_3 %*% pars$s_6_3[i,] 
    #smooth.dat.pred$Zs_7_1 %*% pars$s_7_1[i,]
    D_orig_pred[,ID$id[ID$N.ID==i]] <- orig.dat.pred$X %*% pars$b[i,] +
      orig.dat.pred$Xs %*% pars$bs[i,] +
      orig.dat.pred$Zs_1_1 %*% pars$s_1_1[i,] + orig.dat.pred$Zs_1_2 %*% pars$s_1_2[i,] + orig.dat.pred$Zs_1_3 %*% pars$s_1_3[i,] + 
      orig.dat.pred$Zs_2_1 %*% pars$s_2_1[i,] + orig.dat.pred$Zs_2_2 %*% pars$s_2_2[i,] + orig.dat.pred$Zs_2_3 %*% pars$s_2_3[i,] + 
      orig.dat.pred$Zs_3_1 %*% pars$s_3_1[i,] + orig.dat.pred$Zs_3_2 %*% pars$s_3_2[i,] + orig.dat.pred$Zs_3_3 %*% pars$s_3_3[i,] + 
      orig.dat.pred$Zs_4_1 %*% pars$s_4_1[i,] + orig.dat.pred$Zs_4_2 %*% pars$s_4_2[i,] + orig.dat.pred$Zs_4_3 %*% pars$s_4_3[i,] + 
      orig.dat.pred$Zs_5_1 %*% pars$s_5_1[i,] + orig.dat.pred$Zs_5_2 %*% pars$s_5_2[i,] + orig.dat.pred$Zs_5_3 %*% pars$s_5_3[i,] + 
      orig.dat.pred$Zs_6_1 %*% pars$s_6_1[i,] + orig.dat.pred$Zs_6_2 %*% pars$s_6_2[i,] + orig.dat.pred$Zs_6_3 %*% pars$s_6_3[i,] 
    
    if(MODEL.TYPE=="lat.long.smooth"){
      D_pred[,ID$id[ID$N.ID==i]] <- D_pred[,ID$id[ID$N.ID==i]] + smooth.dat.pred$Zs_7_1 %*% pars$s_7_1[i,]
      D_orig_pred[,ID$id[ID$N.ID==i]] <- D_orig_pred[,ID$id[ID$N.ID==i]] + orig.dat.pred$Zs_7_1 %*% pars$s_7_1[i,]
      
      # Calculate the conditional effect for bottom.depth
      D_pred_bottom_depth[,ID$id[ID$N.ID==i]] <-
        smooth.dat.pred$Xs[,ncol(smooth.dat.pred$Xs)] * pars$bs[i,ncol(smooth.dat.pred$Xs)] +
        smooth.dat.pred$Zs_7_1 %*% pars$s_7_1[i,]
    }
    
    
    D_smooth[,ID$id[ID$N.ID==i]] <-         smooth.dat.pred$Xs[,(1:ncol(smooth.dat.pred$Xs)-1)] %*% pars$bs[i,(1:ncol(smooth.dat.pred$Xs)-1)] +
      smooth.dat.pred$Zs_1_1 %*% pars$s_1_1[i,] + smooth.dat.pred$Zs_1_2 %*% pars$s_1_2[i,] + smooth.dat.pred$Zs_1_3 %*% pars$s_1_3[i,] + 
      smooth.dat.pred$Zs_2_1 %*% pars$s_2_1[i,] + smooth.dat.pred$Zs_2_2 %*% pars$s_2_2[i,] + smooth.dat.pred$Zs_2_3 %*% pars$s_2_3[i,] + 
      smooth.dat.pred$Zs_3_1 %*% pars$s_3_1[i,] + smooth.dat.pred$Zs_3_2 %*% pars$s_3_2[i,] + smooth.dat.pred$Zs_3_3 %*% pars$s_3_3[i,] + 
      smooth.dat.pred$Zs_4_1 %*% pars$s_4_1[i,] + smooth.dat.pred$Zs_4_2 %*% pars$s_4_2[i,] + smooth.dat.pred$Zs_4_3 %*% pars$s_4_3[i,] + 
      smooth.dat.pred$Zs_5_1 %*% pars$s_5_1[i,] + smooth.dat.pred$Zs_5_2 %*% pars$s_5_2[i,] + smooth.dat.pred$Zs_5_3 %*% pars$s_5_3[i,] + 
      smooth.dat.pred$Zs_6_1 %*% pars$s_6_1[i,] + smooth.dat.pred$Zs_6_2 %*% pars$s_6_2[i,] + smooth.dat.pred$Zs_6_3 %*% pars$s_6_3[i,] 
    
  }    
  
  D_pred_summary <- data.frame(Mean=rowMeans(10^(D_pred+LOG.EXPAND)),
                               SD=apply(10^(D_pred+LOG.EXPAND),1,sd),
                               Val=apply(10^(D_pred+LOG.EXPAND),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  D_pred_log_summary <- data.frame(Mean=rowMeans(D_pred+LOG.EXPAND),
                                   SD=apply((D_pred+LOG.EXPAND),1,sd),
                                   Val=apply((D_pred+LOG.EXPAND),1,quantile,probs=PROBS) %>% t() %>% as.matrix())

  D_orig_pred_summary <- data.frame(Mean=rowMeans(10^(D_orig_pred+LOG.EXPAND)),
                               SD=apply(10^(D_orig_pred+LOG.EXPAND),1,sd),
                               Val=apply(10^(D_orig_pred+LOG.EXPAND),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  D_orig_pred_log_summary <- data.frame(Mean=rowMeans(D_orig_pred+LOG.EXPAND),
                                   SD=apply((D_orig_pred+LOG.EXPAND),1,sd),
                                   Val=apply((D_orig_pred+LOG.EXPAND),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  
  if(MODEL.TYPE=="lat.long.smooth"){
    D_pred_bottom_depth_summary <- data.frame(Mean=rowMeans(D_pred_bottom_depth),
                                            SD=apply((D_pred_bottom_depth),1,sd),
                                            Val=apply((D_pred_bottom_depth),1,quantile,probs=PROBS) %>% t() %>% as.matrix())
    D_pred_bottom_depth_combined <- bind_cols(D_pred_bottom_depth_summary,new_data_trim %>% dplyr::select(-Y)) %>%
      left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID,depth_cat_factor))
  }
  
  D_pred_smooth <- data.frame(Mean=rowMeans(D_smooth+LOG.EXPAND),
                              SD=apply(D_smooth+LOG.EXPAND,1,sd),
                              Val=apply(D_smooth+LOG.EXPAND,1,quantile,probs=PROBS) %>% t() %>% as.matrix())
  
  D_pred_combined <- bind_cols(D_pred_summary,new_data_trim %>% dplyr::select(-Y)) %>%
        left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID,depth_cat_factor))
  D_pred_log_combined <- bind_cols(D_pred_log_summary,new_data_trim %>% dplyr::select(-Y)) %>%
        left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID,depth_cat_factor))
  
  D_orig_pred_combined <- bind_cols(D_orig_pred_summary,
                                    orig_data_trim %>% dplyr::select(-Y),
                                    orig_data %>% dplyr::select(station_depth_idx)) %>%
                        left_join(.,orig_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,depth_cat_factor,station_depth_idx))
  D_orig_pred_log_combined <- bind_cols(D_orig_pred_log_summary,
                                        orig_data_trim %>% dplyr::select(-Y),
                                        orig_data %>% dplyr::select(station_depth_idx))  %>%
                        left_join(.,orig_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,depth_cat_factor,station_depth_idx))
  
  D_pred_smooth_combined <- bind_cols(D_pred_smooth,new_data_trim %>% dplyr::select(-Y)) %>%
    left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID,depth_cat_factor))
  
  if(MODEL.TYPE=="lat.long.smooth"){
    smooth.projections <- list(D_pred_combined = D_pred_combined,
                               D_pred_log_combined = D_pred_log_combined,
                               D_pred_bottom_depth_combined = D_pred_bottom_depth_combined,
                               D_pred_smooth_combined = D_pred_smooth_combined,
                               D_orig_pred_combined = D_orig_pred_combined,
                               D_orig_pred_log_combined = D_orig_pred_log_combined
                               )
  }else if(MODEL.TYPE=="lat.long.smooth.base"){
    smooth.projections <- list(D_pred_combined = D_pred_combined,
                               D_pred_log_combined= D_pred_log_combined,
                               D_pred_smooth_combined=D_pred_smooth_combined,
                               D_orig_pred_combined = D_orig_pred_combined,
                               D_orig_pred_log_combined = D_orig_pred_log_combined
                                )
  }
}


# Calculate some residual 
Resid <-D_delta_out_liter %>% dplyr::select(sample_idx,Mean.log,Mean) %>% 
          left_join(.,
             SAMPLES %>% 
               dplyr::select(sample_idx,samp_station_depth_idx) %>% 
               rename(station_depth_idx=samp_station_depth_idx)) %>%
          left_join(.,
              station_depth_out_liter %>% 
                dplyr::select(station_depth_idx,Mean.log,Mean) %>% 
                rename(sd.mean.log = Mean.log,sd.mean = Mean))

Resid <- Resid %>% mutate(Resid.log = Mean.log - sd.mean.log,Resid = Mean - sd.mean) %>%
              left_join(., 
                        STATION.DEPTH %>% 
                          dplyr::select(utm.lon,utm.lat,
                                        lon,lat,
                                        station_depth_idx,depth_cat_factor))

# Combine the necessary data.frames into a list for use later.
Output.summary <- list(
  # Output (Standards)
  stand.plot = stand.plot,
  stand.plot.pres = stand.plot.pres,
  # Output (Field Summaries)
  station_depth_out=station_depth_out,
  station_depth_out_liter=station_depth_out_liter,
  D_delta_out = D_delta_out,
  D_delta_out_liter = D_delta_out_liter,
  Resid = Resid,
  
  # Projections or projection helpers.
  N.POST = N.POST, # number of posterior samples used.
  dat_raster_fin = dat_raster_fin, 
  smooth.projections = smooth.projections,
  
  # Output(negative Controls, Contamination)
  field_neg_out=field_neg_out,
  field_neg_out_liter=field_neg_out_liter,
  sample_contam_total_out = sample_contam_total_out,
  sample_contam_total_out_liter = sample_contam_total_out_liter,
  delta_out=delta_out,
  mu_contam_out = mu_contam_out,
  sigma_contam_out =sigma_contam_out
)
 
