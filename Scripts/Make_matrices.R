# Prep spatial matrices for hake 
###################################################################
# Call function needed by INLA section below
source(here("Scripts","mesh_functions.R"),local=TRUE)
####

## Dealing with standards.

# Make fixed effect indices for standards
# make indexes for the standards and the samples.
dat.stand$plate_idx = dat.stand$plate_idx - 1
dat.samp$plate_idx = dat.samp$plate_idx - 1

n_plate <- nrow(PCR)
n_stand= nrow(dat.stand)
###################################################################
## Make offsets 
###################################################################

# Dilution offset
dat.samp <- dat.samp %>% mutate(ln_dilution = log(dilution))

# Volume sampled offset
dat.samp <- dat.samp %>% mutate(ln_vol_sample = log(vol.standard))

###################################################################
##3 MAKE BOTTLE INDEXES
###################################################################

bottle_all <- dat.samp %>% group_by(year,station,depth_cat,sample) %>% 
  summarise(N=length(sample))
bottle_all$bottle_idx <- 1:nrow(bottle_all)

bottle_id <- bottle_all  %>%
  ungroup() %>% 
  group_by(year,station,depth_cat) %>%
  summarise(N_bottle=length(sample))%>%
  mutate(drop_RE=ifelse(N_bottle==1,1,0)) %>%
  ungroup() %>%
  mutate(station_idx = 1:length(year)) # add 1 or 2 for each bottle_pair

  bottle_all <- left_join(bottle_all,bottle_id) %>%
      arrange(station_idx,bottle_idx) %>%
      group_by(year,station_idx) %>%
      mutate(bot_indicator = 1:length(year)) %>% # add 1 or 2 for each bottle_pair
      mutate(bot_indicator=ifelse(bot_indicator ==2,-1,bot_indicator))
  
  bottle_trim <- bottle_all %>% filter(drop_RE==0)
  bottle_trim$bottle_trim_idx <- 1:nrow(bottle_trim)
  
  st_trim <- bottle_trim %>% distinct(year,station,depth_cat,station_idx)
  st_trim$station_trim_idx <- 1:nrow(st_trim)
  
  bottle_all <- bottle_all %>% left_join(.,bottle_trim %>% dplyr::select(bottle_idx,bottle_trim_idx)) %>%
                        left_join(.,st_trim)
  bottle_all <- bottle_all %>% 
                      mutate(bottle_trim_idx = ifelse(is.na(bottle_trim_idx),-99,bottle_trim_idx)) %>%
                      mutate(station_trim_idx = ifelse(is.na(station_trim_idx),-99,station_trim_idx))

  bottle_all$bottle_trim_idx <- factor(bottle_all$bottle_trim_idx,levels=c(-99,bottle_trim$bottle_trim_idx))
  bottle_all$station_trim_idx <- factor(bottle_all$station_trim_idx,levels=c(-99,st_trim$station_trim_idx))

  ########
  # Check for no more than 2 samples per station.
    A <- left_join(dat.samp,bottle_all)
    A <- A %>% group_by(year,station,depth_cat) %>%
            mutate(max.bot = max(bot_indicator)) 
    B <- A %>% filter(max.bot>1)
    B %>% distinct(year,station,depth_cat)
  #########
  
  dat.samp <- left_join(dat.samp,bottle_all)
  
###################################################################
###################################################################
# Fixed effects  
dat.samp <- dat.samp %>% mutate(year = as.factor(year)) %>%
                mutate(copies_ul = ifelse(is.na(copies_ul),0,copies_ul))

# dat.temp <- dat.samp 
# dat.temp <- dat.temp %>% mutate(copies_ul = ifelse(is.na(copies_ul),0,copies_ul))

# start with year terms and include wash_idx as  0/1 but continuous
dat.samp$year <- as.factor(dat.samp$year)
model_frame   <- model.frame(FORM.fixed, dat.samp)
# So This is the fixed effects matrix.
Xf <- model.matrix(FORM.fixed, model_frame)
#attr(X_f,"names") <- list(names=colnames(X_f))
ncol_beta = ncol(Xf)

# Make some indexes for year and depth categories.
#Y_i <- dat.samp$copies_ul 
n_i <- nrow(dat.samp)

uni.year <- dat.samp %>% distinct(year) %>% arrange(year)
year_idx <- as.numeric(factor(dat.samp$year,levels = uni.year$year)) - 1L

### Make a small prediction matrix for unique fixed effect combinations.
fixed.new <- model_frame %>% mutate(copies_ul = 1) %>% distinct() 
Xf_new <- model.matrix(FORM.fixed, fixed.new)

# this is the unique depth categories that are contained in the observations.
uni.depths <- dat.samp %>% distinct(depth_cat) %>% arrange(depth_cat)
depth_idx <- as.numeric(factor(dat.samp$depth_cat,levels = uni.depths$depth_cat)) - 1L

n_d = nrow(uni.depths)
###################################################################
###################################################################
###################################################################
# Don't forget about wash mistake that affected only 2019 samples.
# make single vector that can be multiplied by a single parameter

###################################################################
###################################################################
###################################################################
# Smooth effects for observations.
# This section defines smooths that act on all of the observations 
# Xs are the linear effects of the smooth
# Zs are the spline components of the smooth (a list)
source(here("src","smoothers.R"),local=TRUE)

# Model form for smoothes
#FORM.smoothes <- "copies_ul ~ s(bottom.depth.consensus,by=year,k=4)"
SM <- parse_smoothers(eval(FORM.smoothes) ,data=dat.samp)

# Files needed in TMB
n_bs     <- ncol(SM$Xs)
b_smooth_start <- SM$b_smooth_start
n_smooth <- length(b_smooth_start)
b_smooth <- if (SM$has_smooths) rep(0,sum(SM$sm_dims)) else array(0) 
has_smooths <- SM$has_smooths

# This is for making predictions to new data.
new_smooth_pred <- parse_smoothers(eval(FORM.smoothes),data=dat.samp,
                                    #newdata= NEWDATA,
                                    basis_prev = SM$basis_out)

###################################################################
## Random effect matrix for bottles.  
###################################################################
tau_has_smooths <- 1 # 0 for yes, 1 for no.
n_tau = ifelse(tau_has_smooths == 0,3,2)

# First calculate depth-varying SD.
FORM.bot.tau <- Ct ~ 0+ s(depth_cat,k=3)
# This is for making the random effect of bottle a function of depth_category.
TAU <- parse_smoothers(eval(FORM.bot.tau) ,data=bottle_trim)

# Files needed in TMB
n_tau_bs     <- ncol(TAU$Xs)
b_tau_smooth_start <- TAU$b_tau_smooth_start
n_tau_smooth <- length(b_tau_smooth_start)
b_tau_smooth <- if (TAU$has_smooths) rep(0,sum(TAU$sm_dims)) else array(0) 

TAU$Zs <- as.matrix(TAU$Zs[[1]])


# This is for making precitions to new data.
new_tau_smooth_pred <- parse_smoothers(eval(FORM.bot.tau),data=bottle_trim,
                                       #newdata= NEWDATA,
                                       basis_prev = TAU$basis_out)

# make a small version of the tau_bottle matrices for giggles.
TAU_sm <- parse_smoothers(eval(FORM.bot.tau) ,data=bottle_trim,
                          newdata=uni.depths,
                          basis_prev = TAU$basis_out)
TAU_sm$Zs <- as.matrix(TAU_sm$Zs[[1]])

### Make design matrices for random bottle effects
### Only include stations that have replicate bottles.
# go through and figure out which of the samples are singletons 
# (only one niskin included instead of two at each location)

if(TRAD.RE.BOTTLE == TRUE){

  FORM.bot <- Ct ~ 0+ as.factor(bottle_trim_idx)

  model_frame   <- model.frame(FORM.bot, dat.samp)
  # So This is the fixed effects matrix.
  X_bottle <- model.matrix(FORM.bot, model_frame)
  # get rid of first column (-99 factor)
  X_bottle <- X_bottle[,-c(1)]

  n_bottle = ncol(X_bottle)

}else if(TRAD.RE.BOTTLE == FALSE){ # This is the stronger, sum to zero constraint design matrix
    # THIS ONLY WORKS If YOU HAVE EXACTLY 2 SAMPLES PER 
  
  FORM.bot <- Ct ~ 0 + as.factor(station_trim_idx)
  
  model_frame   <- model.frame(FORM.bot, dat.samp)
  # So This is the fixed effects matrix.
  X_bottle <- model.matrix(FORM.bot, model_frame)
  # get rid of first column (-99 factor)
  X_bottle <- X_bottle[,-c(1)]
  
  # multiply each column of X_bottle by the column in the dat.samp that is 1 or -1
  X_bottle <- X_bottle * dat.samp$bot_indicator
  
  n_bottle = ncol(X_bottle)
}

# Fix indexing so it starts at 0
bottle_trim$bottle_trim_idx <- bottle_trim$bottle_trim_idx -1L
bottle_all$bottle_trim_idx <- as.numeric(as.character(bottle_all$bottle_trim_idx)) -1

st_trim$station_trim_idx <-st_trim$station_trim_idx -1L
bottle_all$station_trim_idx <- as.numeric(as.character(bottle_all$station_trim_idx)) -1

if(TRAD.RE.BOTTLE == TRUE){ bottle_RE <- bottle_trim %>% mutate(RE_idx = bottle_trim_idx)}
if(TRAD.RE.BOTTLE == FALSE){ bottle_RE <- st_trim %>% mutate(RE_idx = station_trim_idx)}

###################################################################
###################################################################
####################################################################
#### SPATIAL FIELD
#### OK. Make a nice mesh for the observed sample locations.
locs <- dat.samp[,c("year","station","utm.lon","utm.lat")] %>% 
              distinct(year,station,utm.lon,utm.lat) %>% 
              mutate(A_ID = 1:nrow(.)-1L)
inp <- sf::st_as_sf(locs,coords=3:4)
domain <-fmesher::fm_nonconvex_hull(inp,
                              concave = -0.025,
                              convex = -0.025)
  
inla_mesh <- fmesher::fm_mesh_2d_inla(
  #loc=locs[,c("utm.lon","utm.lat")],
  loc.domain = domain, # coordinates
  boundary=domain,
  max.edge = c(40, 1000), # max triangle edge length; inner and outer meshes
  offset = c(30, 80),  # inner and outer border widths
  #max.n.strict=100,#,
  cutoff = 62 , # minimum triangle edge length
  min.angle=20
)
mesh <- make_mesh(dat.samp, c("utm.lon", "utm.lat"), mesh = inla_mesh)
mesh$mesh$n
plot(mesh)

mesh <- mesh
mesh$loc_xy <- as.matrix(locs[, mesh$xy_cols, drop = FALSE])
mesh$A_st <- fmesher::fm_basis(mesh$mesh, loc = mesh$loc_xy)
mesh$sdm_spatial_id <- locs$A_ID

A_st = mesh$A_st
A_ID = mesh$sdm_spatial_id 

spde = get_spde_matrices(mesh)
spde_aniso = make_anisotropy_spde(mesh)

# merge location ID back into data file.
dat.samp <- left_join(dat.samp, locs %>% dplyr::select(year,station,A_ID))
A_ID_idx <- dat.samp$A_ID

n_s  <- nrow(mesh$mesh$loc) # number of knot locations
n_st <- nrow(A_st) # number of unique station locations among years

################################################################################
################################################################################
# Smooth effects for the weight matrices.
# This section defines smooths that act on the weights for the factor analysis.
# Xs are the linear effects of the smooth
# Zs are the spline components of the smooth (a list)

# For now, use a single basis function set for all years. Make sure 
depth_dat <- list()
for(f in 1:n_f){
  f.end <- length(uni.depths$depth_cat)
  depth_dat[[f]] <- data.frame( Y = rep(1,f.end-f+1),depth_cat=c(uni.depths$depth_cat)[f:f.end])
}
# This will become a list with one value entry for each factor
X_L <- list()
Z_L <- list()
b_L_smooth_start <- rep(0,n_y)

L_basis_start <-parse_smoothers(eval(FORM.L), data=depth_dat[[1]])

for(f in 1:n_f){
    L_basis  <-parse_smoothers(eval(FORM.L), data=depth_dat[[1]],
                               newdata=depth_dat[[f]],### CHANGE THIS TO f if you want a different matrix (first few rows ==0)
                               basis_prev = L_basis_start$basis_out)
    X_L$X_L[[f]] <- cbind(rep(1,nrow(L_basis$Xs)),L_basis$Xs)
    Z_L$Z_L[[f]] <- L_basis$Zs[[1]]
   if(f>1){X_L$X_L[[f]] =rbind(matrix(0,f-1,ncol(X_L$X_L[[f]])),X_L$X_L[[f]])
           Z_L$Z_L[[f]] =rbind(matrix(0,f-1,ncol(Z_L$Z_L[[f]])),Z_L$Z_L[[f]])}
}


if(i < n_y){b_L_smooth_start[i+1]  <- b_L_smooth_start[i] + ncol(as.data.frame(L_basis$Zs))} 

  for(f in 1:n_f){
    X_L$col_dims[f] <- ncol(X_L$X_L[[f]])
    Z_L$col_dims[f] <- ncol(Z_L$Z_L[[f]])
  }

n_bs_L    <- sum(X_L$col_dims)
b_L_smooth <- rep(0,sum(Z_L$col_dims)) 

#// factor indexes for weight matrix 
F_L_idx <- rep(1:n_f,n_y) -1L
Y_L_idx <- sort(rep(1:n_y,n_f)) -1L

FY_start_idx <- rep(0,n_y)
  for(i in 1:n_y){
    if(i==1){ FY_start_idx[i]= 0}
    if(i > 1){FY_start_idx[i] = FY_start_idx[i-1] + length(Y_L_idx[Y_L_idx ==(i-1)])}
  }

# This is how you make a prediction set.
newdepths_L <- parse_smoothers(eval(FORM.L), data=depth_dat[[1]],
                           newdata=data.frame(Y=1:3,depth_cat=c(50,100,150)),
                           basis_prev = L_basis_start$basis_out)

##########################################################################
##########################################################################
##########################################################################
##########################################################################
### MAKE PREDICITON SET FOR EACH COMPONENT.
##########################################################################
##########################################################################
##########################################################################
##########################################################################

# Start with the locations. to predict to.
# Make a projection matrix to all of the points included in the domain.
plot(x=dat_raster_trim$utm.lon,y=dat_raster_trim$utm.lat)
mesh_pred <-  make_mesh(dat_raster_trim, c("utm.lon", "utm.lat"), mesh = inla_mesh)
A_pred <- mesh_pred$A_st
n_pred_locs <- nrow(A_pred)


### Make and index trimming the years at different points.
lat.lim <- dat.samp %>% group_by(year) %>% summarise(min.lat=min(utm.lat),max.lat=max(utm.lat))

lat.lim.s <- max(lat.lim$min.lat) - 10


#########################################################################
# Make a set of prediction depths for each of those predicted locations.

dat.pred <- dat_raster_trim %>% mutate(A_ID_pred_idx = (1:nrow(.))-1) 
dat.pred <- expand_grid(dat.pred, depth_cat = depth_pred) %>%
              mutate(bottom.depth.NGDC=depth_m,
                     ln.bottom.depth.NGDC=log(bottom.depth.NGDC),
                     copies_ul = 1)

# Get rid of depth locations that are deeper than the bottom
dat.pred <- dat.pred %>% filter(bottom.depth.NGDC > depth_pred)

# Add flag for 2019- and 2021-specific locations
dat.pred <- dat.pred %>% left_join(.,ac.lims[["2019"]] %>% mutate(x2019 = 1) %>% dplyr::select(Gridcell_ID,x2019))
dat.pred <- dat.pred %>% left_join(.,ac.lims[["2021"]] %>% mutate(x2021 = 1) %>% dplyr::select(Gridcell_ID,x2021))

                    
# cross with year factor
FORM.fixed
yr <- dat.samp %>% distinct(year) %>% mutate(year=as.numeric(as.character(year)))

dat.pred1 <- NULL
for(y in 1:length(yr$year)){
  tmp <- dat.pred 
  tmp$year = yr$year[y]
  dat.pred1 <- rbind(dat.pred1,tmp)
}

dat.pred <- dat.pred1 %>% mutate(year = as.factor(year))
rm(dat.pred1)
################################################################################
# Fixed effects for prediction
model_frame   <- model.frame(FORM.fixed, dat.pred)
# So This is the fixed effects matrix.
Xf_pred <- model.matrix(FORM.fixed, model_frame)

# Smooth main effects for prediction
# This section defines smooths that act on all of the observations 
# Xs are the linear effects of the smooth
# Zs are the spline components of the smooth (a list)

# Smooth effects for predictions
# This is for making precitions to new data.
SM_pred <- parse_smoothers(eval(FORM.smoothes),data=dat.samp,
                                   newdata= dat.pred,
                                   basis_prev = SM$basis_out)
# Files needed in TMB :::: same as the smoothes above

##### Weight Matrix on predictions.
depth_dat_pred <- list()
for(f in 1:n_f){
  f.end <- length(uni.depths$depth_cat)
  depth_dat_pred[[f]] <- data.frame( Y = rep(1,length(depth_pred)-f+1),depth_cat=c(depth_pred)[f:length(depth_pred)])
}

# This will become a list with one value entry for each factor
X_L_pred <- list()
Z_L_pred <- list()

b_L_smooth_start <- rep(0,n_y)

# This is how you make a prediction set.
L_basis_pred <- parse_smoothers(eval(FORM.L), data=depth_dat[[1]],
                               newdata=data.frame(Y=rep(1,length(depth_pred)),depth_cat=depth_pred),
                               basis_prev = L_basis_start$basis_out)

for(f in 1:n_f){
  L_basis_pred <- parse_smoothers(eval(FORM.L), data=depth_dat[[1]],
                                  newdata=depth_dat_pred[[f]],
                                  basis_prev = L_basis_start$basis_out)  
  X_L_pred$X_L[[f]] <- cbind(rep(1,nrow(L_basis_pred$Xs)),L_basis_pred$Xs)
  Z_L_pred$Z_L[[f]] <- L_basis_pred$Zs[[1]]
  if(f>1){
    X_L_pred$X_L[[f]] =rbind(matrix(0,f-1,ncol(X_L_pred$X_L[[f]])),X_L_pred$X_L[[f]])
    Z_L_pred$Z_L[[f]] =rbind(matrix(0,f-1,ncol(Z_L_pred$Z_L[[f]])),Z_L_pred$Z_L[[f]])
  }
}

n_pred <- nrow(dat.pred)

# Add flag for minimum latitude
dat.pred <- dat.pred %>% mutate(shared_locs = ifelse(utm.lat>lat.lim.s,1,0))



# Add flag for 2021-specific locations

###########
###########
## ADD INDEX FOR DEPTH AND YEAR and LOCATION ID for prediction set
###########
###########

A_ID_pred_idx = dat.pred$A_ID_pred_idx

# this is the unique depth categories that are contained in the the prediction set (same levels as those in dat.samp)
year_pred_idx <- as.numeric(factor(dat.pred$year,levels = uni.year$year)) - 1L

# this is the unique depth categories that are contained in the the prediction set
depth_pred_idx <- as.numeric(factor(dat.pred$depth_cat,levels = depth_pred)) - 1L
unique(depth_pred_idx)

