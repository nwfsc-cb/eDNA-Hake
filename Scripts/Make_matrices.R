# Prep spatial matrices for hake 
###################################################################
## Dealing with standards.

# Make fixed effect matrix for standards

n_qpcr <- nrow(PCR)
plate_idx = PCR$plade_idx -1L

dat.stand$plate_idx = dat.stand$plate_idx -1
dat.samp$plate_idx = dat.samp$plate_idx -1






























###################################################################
###################################################################
# Fixed effects  

dat.samp <- dat.samp %>% mutate(year = as.factor(year)) %>%
                mutate(copies_ul = ifelse(is.na(copies_ul),0,copies_ul))

dat.temp <- dat.samp 
dat.temp <- dat.temp %>% mutate(copies_ul = ifelse(is.na(copies_ul),0,copies_ul))

# start with year terms.
dat.temp$year <- as.factor(dat.temp$year)
FORM.fixed <- copies_ul ~ year 
model_frame   <- model.frame(FORM.fixed, dat.temp)
# So This is the fixed effects matrix.
Xf <- model.matrix(FORM.fixed, model_frame)
#attr(X_f,"names") <- list(names=colnames(X_f))
ncol_beta = ncol(Xf)

# Make some indexes for year and depth categories.
Y_i <- dat.temp$copies_ul 

n_i <- length(Y_i)

uni.year <- dat.samp %>% distinct(year) %>% arrange(year)
year_idx <- as.numeric(factor(dat.samp$year,levels = uni.year$year)) - 1L

# this is the unique depth categories that are contained in the observations.
uni.depths <- dat.samp %>% distinct(depth_cat) %>% arrange(depth_cat)
depth_idx <- as.numeric(factor(dat.samp$depth_cat,levels = uni.depths$depth_cat)) - 1L

n_d = nrow(uni.depths)
###################################################################
###################################################################
###################################################################
# Collect Offsets.

# Dilution offset
dat.samp$log_dil <- log(dat.samp$dilution)
# Wash Offset (only affects 2019 samples)
# see dat.samp$wash.indicator

###################################################################
###################################################################
###################################################################
# Smooth effects for observations.
# This section defines smooths that act on all of the observations 
# Xs are the linear effects of the smooth
# Zs are the spline components of the smooth (a list)
source(here("src","smoothers.R"),local=TRUE)

# Model form for smoothes
FORM.smoothes <- "copies_ul ~ s(depth_cat,by=year,k=4)" # + s(bottom.depth.consensus,by=year,k=4)"
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
  max.edge = c(80, 2000), # max triangle edge length; inner and outer meshes
  offset = c(20, 20),  # inner and outer border widths
  #max.n.strict=100,#,
  cutoff =25, # minimum triangle edge length
  min.angle=21
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
# Define the number of spatial fields to use for each year
n_f <- 1
# Define the number of years
n_y <- dat.samp %>% distinct(year) %>% pull(year) %>% length()
# Define a vector controlling the number of spatial fields to use for each year
n_fy <- n_f*n_y

################################################################################
# Smooth effects for the weight matrices.
# This section defines smooths that act on the weights for the factor analysis.
# Xs are the linear effects of the smooth
# Zs are the spline components of the smooth (a list)

FORM.L <- "Y ~ s(depth_cat,k=4)"

# For now, use a single basis function set for all years. Make sure 
depth_dat <- list()
depth_dat[["all"]] <- data.frame( Y = rep(1,length(uni.depths$depth_cat)),depth_cat=c(uni.depths$depth_cat))

X_L <- list()
Z_L <- list()
b_L_smooth_start <- rep(0,n_y)

temp <-parse_smoothers(eval(FORM.L), data=depth_dat[["all"]])

  for(i in 1:n_y){
    X_L$X_L[[i]] <- temp$Xs
    Z_L$Z_L[[i]] <- temp$Zs[[1]]
    if(i < n_y){b_L_smooth_start[i+1]  <- b_L_smooth_start[i] + ncol(as.data.frame(temp$Zs))} 
  }
  for(i in 1:n_y){
    X_L$col_dims[i] <- ncol(X_L$X_L[[i]])
    Z_L$col_dims[i] <- ncol(Z_L$Z_L[[i]])
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
newdepths_L <- parse_smoothers(eval(FORM.L), data=depth_dat[["all"]],
                           newdata=data.frame(Y=1:3,depth_cat=c(0,25,50)),
                           basis_prev = temp$basis_out)

# depth_dat[["2019"]] <- data.frame( Y = rep(1,length(uni.depths$depth_cat)),depth_cat=uni.depths$depth_cat)
# depth_dat[["2021"]] <- data.frame( Y = rep(1,length(uni.depths$depth_cat)),depth_cat=uni.depths$depth_cat)
# 
# X_L <- list()
# Z_L <- list()
# b_L_smooth_start <- rep(0,n_y)
# 
# for(i in 1:n_y){
#   temp <-parse_smoothers(Y ~ s(depth_cat,k=4), data=depth_dat[[i]])
#   X_L[[i]] <- temp$Xs
#   Z_L[[i]] <- temp$Zs
#   if(i < n_y){b_L_smooth_start[i+1]  <- b_L_smooth_start[i] + ncol(as.data.frame(temp$Zs))} 
# }





















