#### Combine and Post-process stan results.
#####################################################3
#####################################################3
#####################################################3
#####################################################3
# Pull files from the list
pars <- Output.qpcr$samp
dat.stand.bin <- Output.qpcr$dat.stand.bin
dat.stand.pos <- Output.qpcr$dat.stand.pos
dat.obs.bin <- Output.qpcr$dat.obs.bin
dat.obs.pos <- Output.qpcr$dat.obs.pos
dat.id <- Output.qpcr$dat.id
PCR <- Output.qpcr$PCR
STATION.DEPTH <- Output.qpcr$STATION.DEPTH
SAMPLES <- Output.qpcr$SAMPLES
MODEL.TYPE <- Output.qpcr$MODEL.TYPE
MODEL.VAR <- Output.qpcr$MODEL.VAR
MODEL.ID <- Output.qpcr$MODEL.ID
NO.SURFACE = Output.qpcr$NO.SURFACE
#dat_raster_fin <- Output.qpcr$dat_raster_fin
depth.fact <- sort(unique(STATION.DEPTH$depth_cat_factor))
TRACE <- Output.qpcr$TRACE

brms.object <- Output.qpcr$brms.object

N.knots.lon <- Output.qpcr$N_knots_all[[1]]
N.knots.lat <- Output.qpcr$N_knots_all[[2]]
if(MODEL.TYPE =="lat.long.smooth"){
  N.knots.bd <- Output.qpcr$N_knots_all[[3]]
}

# N posteriors to sample
N.POST <- 4000

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


# THIS IS THE FACTOR FOR EXPANDING FROM COPIES PER ul to COPIES PER L 
# Based on sampling volume of 2.5 L

EXPAND.CONTROL <- 40
LOG.EXPAND.CONTROL <- log10(EXPAND.CONTROL)

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
                                            Mean.log=apply(pars$D_contam+LOG.EXPAND.CONTROL,2,mean),
                                            Sd.log=apply(pars$D_contam+LOG.EXPAND.CONTROL,2,sd),
                                            Log.Val=data.frame(t(apply(pars$D_contam+LOG.EXPAND.CONTROL,2,quantile,probs=PROBS))),
                                            Mean=apply(10^(pars$D_contam+LOG.EXPAND.CONTROL),2,mean),
                                            Sd=apply(10^(pars$D_contam+LOG.EXPAND.CONTROL),2,sd),
                                            Val=data.frame(t(apply(10^(pars$D_contam+LOG.EXPAND.CONTROL),2,quantile,probs=PROBS))))

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
  
  # Drop the 300m observation for Gridcell 133469 (deep spot right near Cape Flattery)
  new_data_1 <- new_data %>% 
                filter(Gridcell_ID!=133469)
  new_data <- bind_rows(new_data_1,
                        new_data %>% 
                          filter(Gridcell_ID==133469,
                                 depth_cat_factor!=300))
  
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
      smooth.dat.pred$Zs_5_1 %*% pars$s_5_1[i,] + smooth.dat.pred$Zs_5_2 %*% pars$s_5_2[i,] + smooth.dat.pred$Zs_5_3 %*% pars$s_5_3[i,] 
      
    #smooth.dat.pred$Zs_7_1 %*% pars$s_7_1[i,]
    D_orig_pred[,ID$id[ID$N.ID==i]] <- orig.dat.pred$X %*% pars$b[i,] +
      orig.dat.pred$Xs %*% pars$bs[i,] +
      orig.dat.pred$Zs_1_1 %*% pars$s_1_1[i,] + orig.dat.pred$Zs_1_2 %*% pars$s_1_2[i,] + orig.dat.pred$Zs_1_3 %*% pars$s_1_3[i,] + 
      orig.dat.pred$Zs_2_1 %*% pars$s_2_1[i,] + orig.dat.pred$Zs_2_2 %*% pars$s_2_2[i,] + orig.dat.pred$Zs_2_3 %*% pars$s_2_3[i,] + 
      orig.dat.pred$Zs_3_1 %*% pars$s_3_1[i,] + orig.dat.pred$Zs_3_2 %*% pars$s_3_2[i,] + orig.dat.pred$Zs_3_3 %*% pars$s_3_3[i,] + 
      orig.dat.pred$Zs_4_1 %*% pars$s_4_1[i,] + orig.dat.pred$Zs_4_2 %*% pars$s_4_2[i,] + orig.dat.pred$Zs_4_3 %*% pars$s_4_3[i,] + 
      orig.dat.pred$Zs_5_1 %*% pars$s_5_1[i,] + orig.dat.pred$Zs_5_2 %*% pars$s_5_2[i,] + orig.dat.pred$Zs_5_3 %*% pars$s_5_3[i,] 
      
    
    if(NO.SURFACE!="TRUE"){
      D_pred[,ID$id[ID$N.ID==i]] <- D_pred[,ID$id[ID$N.ID==i]] +
          smooth.dat.pred$Zs_6_1 %*% pars$s_6_1[i,] + smooth.dat.pred$Zs_6_2 %*% pars$s_6_2[i,] + smooth.dat.pred$Zs_6_3 %*% pars$s_6_3[i,] 
      D_orig_pred[,ID$id[ID$N.ID==i]] <- D_orig_pred[,ID$id[ID$N.ID==i]] +
          orig.dat.pred$Zs_6_1 %*% pars$s_6_1[i,] + orig.dat.pred$Zs_6_2 %*% pars$s_6_2[i,] + orig.dat.pred$Zs_6_3 %*% pars$s_6_3[i,] 
    }  
    if(MODEL.TYPE=="lat.long.smooth" & NO.SURFACE!="TRUE"){
      D_pred[,ID$id[ID$N.ID==i]] <- D_pred[,ID$id[ID$N.ID==i]] + smooth.dat.pred$Zs_7_1 %*% pars$s_7_1[i,]
      D_orig_pred[,ID$id[ID$N.ID==i]] <- D_orig_pred[,ID$id[ID$N.ID==i]] + orig.dat.pred$Zs_7_1 %*% pars$s_7_1[i,]
      
      # Calculate the conditional effect for bottom.depth
      D_pred_bottom_depth[,ID$id[ID$N.ID==i]] <-
        smooth.dat.pred$Xs[,ncol(smooth.dat.pred$Xs)] * pars$bs[i,ncol(smooth.dat.pred$Xs)] +
        smooth.dat.pred$Zs_7_1 %*% pars$s_7_1[i,]
    }
    
    if(MODEL.TYPE=="lat.long.smooth" & NO.SURFACE=="TRUE"){
      D_pred[,ID$id[ID$N.ID==i]] <- D_pred[,ID$id[ID$N.ID==i]] + smooth.dat.pred$Zs_6_1 %*% pars$s_6_1[i,]
      D_orig_pred[,ID$id[ID$N.ID==i]] <- D_orig_pred[,ID$id[ID$N.ID==i]] + orig.dat.pred$Zs_6_1 %*% pars$s_6_1[i,]
      
      # Calculate the conditional effect for bottom.depth
      D_pred_bottom_depth[,ID$id[ID$N.ID==i]] <-
        smooth.dat.pred$Xs[,ncol(smooth.dat.pred$Xs)] * pars$bs[i,ncol(smooth.dat.pred$Xs)] +
        smooth.dat.pred$Zs_6_1 %*% pars$s_6_1[i,]
    }
    
    D_smooth[,ID$id[ID$N.ID==i]] <-         smooth.dat.pred$Xs[,(1:(ncol(smooth.dat.pred$Xs)-1))] %*% pars$bs[i,(1:(ncol(smooth.dat.pred$Xs)-1))] +
      smooth.dat.pred$Zs_1_1 %*% pars$s_1_1[i,] + smooth.dat.pred$Zs_1_2 %*% pars$s_1_2[i,] + smooth.dat.pred$Zs_1_3 %*% pars$s_1_3[i,] +
      smooth.dat.pred$Zs_2_1 %*% pars$s_2_1[i,] + smooth.dat.pred$Zs_2_2 %*% pars$s_2_2[i,] + smooth.dat.pred$Zs_2_3 %*% pars$s_2_3[i,] +
      smooth.dat.pred$Zs_3_1 %*% pars$s_3_1[i,] + smooth.dat.pred$Zs_3_2 %*% pars$s_3_2[i,] + smooth.dat.pred$Zs_3_3 %*% pars$s_3_3[i,] +
      smooth.dat.pred$Zs_4_1 %*% pars$s_4_1[i,] + smooth.dat.pred$Zs_4_2 %*% pars$s_4_2[i,] + smooth.dat.pred$Zs_4_3 %*% pars$s_4_3[i,] +
      smooth.dat.pred$Zs_5_1 %*% pars$s_5_1[i,] + smooth.dat.pred$Zs_5_2 %*% pars$s_5_2[i,] + smooth.dat.pred$Zs_5_3 %*% pars$s_5_3[i,] 
      
    if(NO.SURFACE!="TRUE"){ 
      D_smooth[,ID$id[ID$N.ID==i]] <-  D_smooth[,ID$id[ID$N.ID==i]] +
        smooth.dat.pred$Zs_6_1 %*% pars$s_6_1[i,] + smooth.dat.pred$Zs_6_2 %*% pars$s_6_2[i,] + smooth.dat.pred$Zs_6_3 %*% pars$s_6_3[i,]
    }
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

  D_smooth_summary <- bind_cols(as.data.frame(D_smooth),new_data_trim %>% dplyr::select(-Y)) %>%
                              left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID,depth_cat_factor)) %>%
                              pivot_longer(starts_with("V"),names_to="MCMC.rep",values_to = "dens") %>%
                              group_by(MCMC.rep,depth_cat_factor) %>%
                              summarise(SD = sd(dens)) %>%
                              group_by(depth_cat_factor) %>%
                              summarise(Mean_SD = mean(SD), SD_SD = sd(SD),
                                        Q.0.025= quantile(SD,probs=c(0.025)),
                                          Q.0.05= quantile(SD,probs=c(0.05)),
                                          Q.0.25= quantile(SD,probs=c(0.25)),
                                          Median = median(SD),
                                          Q.0.75=quantile(SD,probs=c(0.75)),
                                          Q.0.95= quantile(SD,probs=c(0.95)),
                                          Q.0.975=quantile(SD,probs=c(0.975)))

    
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


# Calculate some residuals and pred-observed for the NISKIN-MEAN FIELD COMPARISONS.  
Resid <- D_delta_out_liter %>% dplyr::select(sample_idx,Mean.log,Mean) %>% 
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

Resid <- Resid %>% rename(pred.mean.log=sd.mean.log,pred.mean=sd.mean,
                          obs.mean=Mean,obs.mean.log=Mean.log)

# Make Pred-Obs plots
p_resid <- list()

p_resid$p1 <- ggplot(Resid) +
    geom_point(aes(y=obs.mean.log,x=pred.mean.log),alpha=0.3) +
    geom_abline(intercept=0,slope=1,color="red") +
    theme_bw()

Resid <- Resid %>% mutate(depth_cat_factor2 = depth_cat_factor,
                          depth_cat_factor2 = as.character(depth_cat_factor2),
                          depth_cat_factor2=ifelse(depth_cat_factor2=="0","3",depth_cat_factor2),
                          depth_cat_factor2 = as.factor(depth_cat_factor2))
Resid$depth_cat_factor2 <- factor(Resid$depth_cat_factor2,levels=sort(as.numeric(as.character(unique(Resid$depth_cat_factor2)))))

p_resid$p2 <- ggplot(Resid) +
  geom_point(aes(y=obs.mean.log,x=pred.mean.log),alpha=0.3) +
  geom_abline(intercept=0,slope=1,color="red") +
  xlab("prediction (log_10)") +
  ylab("observation (log_10)") + 
  theme_bw() +
  facet_wrap(~depth_cat_factor2)

p_resid$p3 <- ggplot(Resid) +
  geom_point(aes(y=obs.mean,x=pred.mean),alpha=0.3) +
  geom_abline(intercept=0,slope=1,color="red") +
  theme_bw() +
  facet_wrap(~depth_cat_factor2)


p_resid$p4 <- ggplot(Resid) +
  geom_point(aes(y=Resid.log,x=pred.mean.log),alpha=0.3) +
  geom_abline(intercept=0,slope=0,color="red") +
  theme_bw() +
  facet_wrap(~depth_cat_factor2)

p_resid$p5 <- ggplot(Resid) +
  geom_point(aes(y=Resid.log,x=lat),alpha=0.3) +
  geom_abline(intercept=0,slope=0,color="red") +
  theme_bw() +
  facet_wrap(~depth_cat_factor)

## Make a second set of resid for the actual qPCR count numbers predicted vs. observed.


# Calculate some predicted - observed for the raw, observed Ct 
##    Both some presence-absence and postive Ct.

p_pred_bin <- 1 / (1+exp(-Output.qpcr$pred_bin))

dat.temp.bin <- dat.obs.bin %>% dplyr::select(sample_idx,qPCR,sample,Ct,Ct_bin,Niskin, 
                       depth,date,transect,lat,lon,bottom.depth,station.depth,
                       depth_cat_factor) %>% bind_cols(.,p_pred_bin)

p_resid$p_raw_bin_resid1 <- ggplot(dat.temp.bin) +
    geom_jitter(aes(y=Ct_bin,x=Mean),alpha=0.2,width=0,height=0.1)+
    geom_abline(intercept=0,slope=1,color="red") +
    theme_bw()



dat.temp.pos <- dat.obs.pos %>% dplyr::select(sample_idx,qPCR,sample,Ct,Ct_bin,Niskin, 
                                        depth,date,transect,lat,lon,bottom.depth,station.depth,
                                        depth_cat_factor) %>% bind_cols(.,Output.qpcr$pred_pos)

# This is for the raw observations
p_resid$p_raw_pos_resid1 <-  ggplot(dat.temp.pos) +
  geom_jitter(aes(y=Ct,x=Mean),alpha=0.2,width=0,height=0.1)+
  geom_abline(intercept=0,slope=1,color="red") +
  theme_bw()

dat.temp.pos <- dat.obs.pos %>% dplyr::select(sample_idx,qPCR,sample,Ct,Ct_bin,Niskin, 
                                              depth,date,transect,lat,lon,bottom.depth,station.depth,
                                              depth_cat_factor) %>% bind_cols(.,Output.qpcr$pred_pos)

dat.temp.pos.sum <- dat.temp.pos %>% group_by(sample_idx,qPCR,sample,Niskin, 
                                              depth,date,transect,lat,lon,bottom.depth,station.depth,
                                              depth_cat_factor) %>% 
                                    summarise(Ct.mean = mean(Ct),pred_Ct = mean(Mean))

p_resid$p_raw_pos_resid2 <-ggplot(dat.temp.pos.sum) +
  geom_jitter(aes(y=Ct.mean,x=pred_Ct),alpha=0.2,width=0,height=0.1)+
  geom_abline(intercept=0,slope=1,color="red") +
  theme_bw()

#############################################
#############################################
#############################################
#############################################
#### ----- From the projected smoothes, interpolate to new depth strata
#############################################
#############################################
#############################################
#############################################

# We use D_pred + LOG.EXPAND from above.

D_pred_interp <- bind_cols(as.data.frame(10^(D_pred+LOG.EXPAND)),new_data_trim %>% dplyr::select(-Y)) %>%
  left_join(.,new_data %>% dplyr::select(lon,lat,utm.lon,utm.lat,Gridcell_ID,depth_cat_factor,bottom.depth.consensus))

# Linear interpolate to these depths (DEPTHS)
d.step = 50 # depth interval in meters to project to
depth.obs <- sort(unique(D_pred_interp$depth_cat_factor))
DEPTHS <- seq(0,max(depth.obs),by= d.step) %>% as.data.frame() %>% rename(d='.')
DEPTHS <- DEPTHS %>% filter(!d %in% depth.obs) 

# Derive depths to define interpolation:
for(i in 1:nrow(DEPTHS)){
  DEPTHS$above[i] <- max(depth.obs[which(depth.obs < DEPTHS$d[i])])
  DEPTHS$below[i] <- min(depth.obs[which(depth.obs > DEPTHS$d[i])])
}
# Define weights based on depth distance:
DEPTHS <- DEPTHS %>% mutate(w.above = (d-above) / (d-above + below-d),
                            w.below = 1-w.above,
                            d.cutoff = d-10)
# loop over interpolation depths:
d.all <- NULL
for(i in 1:nrow(DEPTHS)){
  temp.a <- D_pred_interp %>% filter(depth_cat_factor == DEPTHS$above[i]) %>% 
                              mutate(weight = DEPTHS$w.above[i])
  temp.b <- D_pred_interp %>% filter(depth_cat_factor == DEPTHS$below[i]) %>%
                              mutate(weight = DEPTHS$w.below[i])
  temp <- bind_rows(temp.a,temp.b)
  # Pick out the instances where there is only one depth prediction (e.g. 150 but no 300)
  temp2 <- temp %>% group_by(Gridcell_ID) %>% summarise(N=length(Gridcell_ID)) %>% filter(N==1)
  # change weights for observations with only one adjacent prediction
  temp <- temp %>% mutate(weight = ifelse(Gridcell_ID %in% temp2$Gridcell_ID,1,weight))
  # get rid of of predictions if the depth is not within 10m of the depth (less that d.cutoff) 
  temp <- temp %>% filter(bottom.depth.consensus > DEPTHS$d.cutoff[i])
  
  # Check to make sure that all the weights sum to 1
  check <- temp %>% group_by(Gridcell_ID) %>% summarise(SUM=sum(weight))
  if(min(check$SUM)<1 | max(check$SUM)>1){ print("ERROR, ERROR"); break}
  
  temp[,1:N.POST] <- temp[,1:N.POST] * temp$weight
  temp <- temp %>% dplyr::select(-utm.lat,-utm.lon,-bottom.depth.consensus, 
                                 -depth_cat_factor,-lon,-lat,-weight) 
  temp.long <- temp %>% pivot_longer(!Gridcell_ID,names_to="MCMC.rep",values_to = "w.dens")

  d.weighted <- temp.long %>% group_by(Gridcell_ID,MCMC.rep) %>% summarise(D = sum(w.dens)) %>% 
                  mutate(depth_cat_factor= DEPTHS$d[i] )
  
  d.all <- bind_rows(d.all,d.weighted)
}


# add on the depths that we already have predictions for here.

  d.measured <-  bind_cols(as.data.frame(10^(D_pred+LOG.EXPAND)),new_data_trim %>% dplyr::select(-Y)) %>%
        left_join(.,new_data %>% dplyr::select(utm.lon,utm.lat,Gridcell_ID,depth_cat_factor)) %>%
        dplyr::select(-utm.lon,-utm.lat,-bottom.depth.consensus)

  d.measured <- d.measured %>% pivot_longer(starts_with("V"),names_to="MCMC.rep",values_to = "w.dens") %>%
                  dplyr::select(Gridcell_ID, MCMC.rep, D = w.dens, depth_cat_factor)

  d.all <- bind_rows(d.all,d.measured)

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

  # Summarize concentration across depths.
  d.all.no.depth <- d.all %>% filter(depth_cat_factor>0) %>% 
    group_by(Gridcell_ID, MCMC.rep,ID.equal,ID.lat.0.5,ID.lat.1.0) %>%
    summarise(tot = sum(D))  # This sums the projections across depths
  

  # Produce different summaries at different spatial aggregations.
  # This is the summary for each cell.
  D_final_projected <- d.all.no.depth %>%
                group_by(Gridcell_ID) %>% 
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
 
  # This is the summary for one way of breaking up the coast (equal latitudinal divisions)
  D_final_lat_equal <- d.all.no.depth %>% rename(Tot = tot) %>%
              group_by(MCMC.rep,ID.equal) %>% 
              summarise(tot = sum(Tot)) %>%  # This sums the projections across strata for each MCMC
              group_by(ID.equal) %>% 
              summarise(Mean=mean(tot), # This summarizes the averages within each strata
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
  
  D_final_lat_0.5 <- d.all.no.depth %>% rename(Tot = tot) %>%
    group_by(MCMC.rep,ID.lat.0.5) %>% 
    summarise(tot = sum(Tot)) %>%  # This sums the projections across strata for each MCMC
    group_by(ID.lat.0.5) %>% 
    summarise(Mean=mean(tot), # This summarizes the averages within each strata
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
  
  D_final_lat_1.0 <- d.all.no.depth %>% rename(Tot = tot) %>%
    group_by(MCMC.rep,ID.lat.1.0) %>% 
    summarise(tot = sum(Tot)) %>%  # This sums the projections across strata for each MCMC
    group_by(ID.lat.1.0) %>% 
    summarise(Mean=mean(tot), # This summarizes the averages within each strata
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

  D_DNA_uncond_total <- d.all.no.depth %>% 
    left_join(.,new_dat_strata %>% dplyr::select(Gridcell_ID,lat)) %>%
    arrange(lat,MCMC.rep) %>%
    group_by(MCMC.rep) %>%
    mutate(cum_sum = cumsum(tot),max_cum_sum = max(cum_sum)) %>%
    mutate(cum_sum_prob = cum_sum / max_cum_sum)

  D_DNA_cum_sum <-   D_DNA_uncond_total %>%    
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
  
  D_DNA_uncond_total <- D_DNA_uncond_total %>% 
    ungroup() %>%
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
  
  # repeat the summaries but calculate an non-dimension index referenced to the southern most area
  # D_final_lat_1.0_NON_DIM <- d.all.no.depth %>% rename(Tot = tot) %>%
  #   group_by(MCMC.rep,ID.lat.1.0) %>% 
  #   summarise(tot = sum(Tot))
  # 
  # D_first_area <- D_final_lat_1.0_NON_DIM
  # 
  # 
  
  D_final_projected <- left_join(D_final_projected, 
                                 new_data %>% dplyr::select(Gridcell_ID,utm.lat,utm.lon,lon,lat) %>% distinct())

  ##### 
  
  D_grid.cell_resample <- d.all.no.depth %>% dplyr::select(-ID.equal,-ID.lat.0.5,-ID.lat.1.0)
                      
  D_1.0_resample <- d.all.no.depth %>%
                      rename(Tot = tot) %>%
                      group_by(MCMC.rep,ID.lat.1.0) %>% 
                      summarise(tot = sum(Tot))
  D_0.5_resample <- d.all.no.depth %>%
                      rename(Tot = tot) %>%
                      group_by(MCMC.rep,ID.lat.0.5) %>% 
                      summarise(tot = sum(Tot))
  
  
  