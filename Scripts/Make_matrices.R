# Prep spatial matrices for hake 
###################################################################
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
  mutate(drop_RE=ifelse(N_bottle==1,1,0))



bottle_all <- left_join(bottle_all,bottle_id) %>%
                arrange(bottle_idx) 

bottle_trim <- bottle_all %>% filter(drop_RE==0)
bottle_trim$bottle_trim_idx <- 1:nrow(bottle_trim)

bottle_all <- bottle_all %>% left_join(.,bottle_trim %>% dplyr::select(bottle_idx,bottle_trim_idx))
bottle_all <- bottle_all %>% mutate(bottle_trim_idx = ifelse(is.na(bottle_trim_idx),-99,bottle_trim_idx))

bottle_all$bottle_trim_idx <- factor(bottle_all$bottle_trim_idx,levels=c(-99,bottle_trim$bottle_trim_idx))

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

# This is for making precitions to new data.
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

FORM.bot <- Ct ~ 0+ as.factor(bottle_trim_idx)

model_frame   <- model.frame(FORM.bot, dat.samp)
# So This is the fixed effects matrix.
X_bottle <- model.matrix(FORM.bot, model_frame)
# get rid of first column (-99 factor)
X_bottle <- X_bottle[,-c(1)]

# go through and figure out which of the samples are singletons 
# (only one niskin included instead of two at each location)
n_bottle = ncol(X_bottle)

# Fix indexing so it starts at 0
bottle_trim$bottle_trim_idx <- bottle_trim$bottle_trim_idx -1L
bottle_all$bottle_trim_idx <- as.numeric(as.character(bottle_all$bottle_trim_idx)) -1
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
  max.edge = c(80, 1000), # max triangle edge length; inner and outer meshes
  offset = c(25, 20),  # inner and outer border widths
  #max.n.strict=100,#,
  cutoff =45, # minimum triangle edge length
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

get_spde_matrices <- function(x) {
  x <- x$spde[c("c0", "g1", "g2")]
  names(x) <- c("M0", "M1", "M2") # legacy INLA names needed!
  x
}
spde = get_spde_matrices(mesh)

# merge location ID back into data file.
dat.samp <- left_join(dat.samp, locs %>% dplyr::select(year,station,A_ID))

A_ID_idx <- dat.samp$A_ID

n_s  <- nrow(mesh$mesh$loc) # number of knot locations
n_st <- nrow(A_st) # number of unique station locations among years

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
                               newdata=depth_dat[[f]],
                               basis_prev = L_basis_start$basis_out)
    X_L$X_L[[f]] <- L_basis$Xs
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

# Design matrix for L

L_design <- matrix(1,n_d,n_f)
if(n_f>1){
  for(i in 2:ncol(L_design)){
    L_design[1:(i-1),i] <- 0
  }
}

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
