library(tidyverse)
library(here)
library(sdmTMB)
library(sf)
library(viridis)
library(rnaturalearth)
library(marmap)

# Build a mesh to implement the SPDE approach:
# mesh <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 20)
# fit <- sdmTMB(
#   density ~ s(depth),
#   data = pcod_2011, mesh = mesh,
#   family = tweedie(link = "log")
# )

# Use the same model as above, but pass in a 2nd dataframe for std curve estimation
# fit2 <- sdmTMB( # should throw error, because wrong family
#   density ~ s(depth),
#   data = pcod_2011, mesh = mesh,
#   family = tweedie(link = "log"),
#   control = sdmTMBcontrol(stdcurve = TRUE, stdcurve_df = d)
# )
#
# fit3 <- sdmTMB(
#   density ~ s(depth),
#   data = pcod_2011, mesh = mesh,
#   family = delta_gamma(),
#   control = sdmTMBcontrol(stdcurve = TRUE, stdcurve_df = d)
# )

# load in observed data
# standards
d <- read_rds(here('Data',"eulachon qPCR 2019 and 2021 standards clean.rds")) %>% 
  mutate(Ct=replace_na(Ct,0))# delta-models in sdmTMB need this

d_obs <- read_rds(here('Data','eulachon qPCR 2019 and 2021 samples clean.rds')) %>% 
  mutate(Ct=replace_na(Ct,0)) %>% # delta-models in sdmTMB need this
  mutate(utm.lon.km=utm.lon.m/1000,
         utm.lat.km=utm.lat.m/1000)

mesh <- make_mesh(d_obs, c("utm.lon.km", "utm.lat.km"), cutoff = 20)

# using only data <=150m depth category
d_obs_filt <- d_obs %>% filter(depth_cat %in% c(0,50,150))
mesh_filt <- make_mesh(d_obs_filt,c("utm.lon.km", "utm.lat.km"), cutoff = 20)

# Fitting the model with no standard curve (intercept only)
# fit4 <- sdmTMB(
#   Ct ~ 1,
#   data = d_obs, mesh = mesh,
#   family = delta_gaussian()
# )
# 
# # Fit everything together # 15213.08
# fit5 <- sdmTMB(
#   Ct ~ 1,
#   data = d_obs, mesh = mesh,
#   family = delta_gaussian(),
#   control = sdmTMBcontrol(stdcurve_df = d)
# )

# this is interesting, runs for > 1 hr. I think there's some confounding
# between spatiotemporal effects and betas?
# fit6 <- sdmTMB(
#   Ct ~ year,
#   data = d_obs, mesh = mesh,
#   time = "year",
#   spatiotemporal = "iid",
#   family = delta_gaussian(),
#   control = sdmTMBcontrol(stdcurve = TRUE, stdcurve_df = d, multiphase = TRUE)
# )


# Add depth as smoother
# fit8 <- sdmTMB( # AIC 15105.79
#   Ct ~ s(depth_cat,k=3),
#   data = d_obs, mesh = mesh,
#   family = delta_gaussian(),
#   control = sdmTMBcontrol(stdcurve_df = d)
# )
# sanity(fit8)

# try depth as spatially varying coefficient
fit9 <- sdmTMB( # AIC 14313.59
  Ct ~ 1,
  spatial_varying = ~ depth_cat,
  data = d_obs_filt, 
  mesh = mesh_filt,
  # time="year",
  # spatiotemporal = "IID",
  family = stdcurve(),
  control = sdmTMBcontrol(stdcurve_df = d)
)

# use time as depth -- doesn't converge
# d_obs$depth_bin <- as.numeric(as.factor(d_obs$depth_cat))
# fit10 <- sdmTMB( # AIC 14313.59
#   Ct ~ 1,
#   time = "depth_bin",
#   spatiotemporal = "iid",
#   data = d_obs, mesh = mesh,
#   family = delta_gaussian(),
#   control = sdmTMBcontrol(stdcurve_df = d)
# )
# 
# # depth as SVC fields
# m <- model.matrix(Ct ~ as.factor(depth_cat), data = d_obs)
# d_obs$d1 <- m[,1]
# d_obs$d2 <- m[,2]
# d_obs$d3 <- m[,3]
# d_obs$d4 <- m[,4]
# d_obs$d5 <- m[,5]
# d_obs$d6 <- m[,6]
# d_obs$d7 <- m[,7]
# d_obs$d8 <- m[,8]

# estimate each depth bin with a separate spatial field / separate variance
# fit11 <- sdmTMB( # this version has a hard time converging
#   Ct ~ 1,
#   spatial_varying = ~ d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8,
#   data = d_obs, mesh = mesh,
#   family = delta_gaussian(),
#   control = sdmTMBcontrol(stdcurve_df = d)
# )

## Random effects
# fit5$sd_report$par.random with SEs in fit5$sd_report$diag.cov.random

# parameter estimates
tidy(fit9,"ran_pars") # range 48.3km +- 7.71km

## Prediction
# MAKE GRID
pred.crs <- terra::rast(here('Data','raster_grid_blake','fivekm_grid.tif')) %>% st_crs()
grid.pred <- read_rds(here('Data','prediction_grid_5km_sdmTMB.rds'))
# add depth categories
grid.pred0 <- grid.pred %>% mutate(depth_cat=0)
grid.pred50 <- grid.pred %>% filter(bathy.bottom.depth<=-50) %>% mutate(depth_cat=50)
grid.pred150 <- grid.pred %>% filter(bathy.bottom.depth<=-100) %>% mutate(depth_cat=150)

grid.predz <- bind_rows(grid.pred0,grid.pred50,grid.pred150)

#background map
coast <- ne_states(country='United States of America',returnclass = 'sf') %>% 
  filter(name %in% c('California','Oregon','Washington','Nevada')) %>%
  st_transform(crs = pred.crs)

# bathy
limits.for.map <- d_obs_filt %>% 
  st_as_sf(coords=c('lon','lat'),crs=4326) %>% 
  st_bbox()+c(-1,-1,1,1) #add an extra degree on each side for buffer


b = getNOAA.bathy(lon1 = limits.for.map["xmin"],
                  lon2 = limits.for.map[ "xmax" ],
                  lat1 = limits.for.map["ymin"],
                  lat2 = limits.for.map["ymax"],
                  resolution = 1,keep = TRUE)
plot(b, image=TRUE, deep=-1000, shallow=0, step=250)
# make into a dataframe of contours for ggplotting
library(contoureR)
bdf <- fortify(b) %>% 
  contoureR::getContourLines(levels=c(0, -100, -200, -500, -1000,-1500)) %>% 
  st_as_sf(coords=c("x","y"),crs=4326) %>% 
  st_transform(pred.crs) %>%
  mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2]) %>% 
  st_set_geometry(NULL)

bathy.only.map <- ggplot(bdf,aes(x,y,group=Group,colour=z)) + geom_path()
bathy.only.map+coord_equal()


# Predict/Project
fit9.pred <- predict(fit9,newdata=grid.predz) %>% 
  mutate(utm.lon.m=utm.lon.km*1000,
       utm.lat.m=utm.lat.km*1000)

# distribution of fitted value
ggplot(fit9.pred,aes(est))+geom_density()
fivenum(exp(fit9.pred$est))       

# Map of predictions
pred.map <- ggplot()+
  geom_path(data=bdf,aes(x,y,group=Group),color='gray20')+
  geom_raster(data=fit9.pred,aes(utm.lon.m,utm.lat.m,fill=est),alpha=0.8)+
  geom_sf(data=coast,fill="gray50")+
  scale_fill_viridis(option='turbo')+
  labs(x="",y="",fill="log(eDNA)")+
  facet_wrap(~depth_cat)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(fill=NA,color='black'))
  # facet_grid(year~depth_cat)
# pred.map
pred.map.zoom <- pred.map+ylim(1e6,2e6)+xlim(7e5,10e5)
pred.map.zoom

#SVC of depth
zeta_z_pred.map <- fit9.pred %>% 
  distinct(zeta_s_depth_cat,.keep_all = T) %>% 
  ggplot()+
  geom_raster(data=fit9.pred,aes(utm.lon.m,utm.lat.m,fill=zeta_s_depth_cat))+
  geom_sf(data=coast,fill="gray50")+
  scale_fill_gradient2()+
  ylim(1e6,2e6)+xlim(7e5,10e5)+
  labs(x="",y="",fill="SVC\nDepth",title="")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(fill=NA,color='black'),
        axis.text = element_blank())
zeta_z_pred.map

# Omega_S (spatial random effects)
#SVC of depth
omega_s_pred.map <- fit9.pred %>% 
  distinct(omega_s,.keep_all = T) %>% 
  ggplot()+
  geom_raster(data=fit9.pred,aes(utm.lon.m,utm.lat.m,fill=omega_s))+
  geom_sf(data=coast,fill="gray50")+
  scale_fill_viridis(option='turbo')+
  ylim(1e6,2e6)+xlim(7e5,10e5)+
  labs(x="",y="",fill=expression(omega),title="Spatial REs")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(fill=NA,color='black'),
        axis.text = element_blank())
omega_s_pred.map

# depth distributions?
lat.breaks <- seq(1e6,2e6,by=1e5)
depth.dist <- fit9.pred %>% 
  mutate(lat.bin=cut(utm.lat.m,lat.breaks)) %>% 
  filter(!is.na(lat.bin)) %>% 
  filter(bathy.bottom.depth>-300) %>% 
  ggplot(aes(bathy.bottom.depth,exp(est)))+
  geom_point()+
  geom_smooth()

depth.dist

# probability of positive?
# phis and betas
plateREs <- tibble(plate=d$plate,platenum=as.numeric(as.factor(d$plate))) %>% 
  distinct() %>% 
  mutate(phi0=fit9$sd_report$par.random[427:491],phi1=fit9$sd_report$par.random[492:556],
         beta0=fit9$sd_report$par.random[557:621],
         beta1=fit9$sd_report$par.random[622:686])
                   
ctphi0 <- fit9$sd_report$par.random[427:491]

# by month?
d_obs_filt %>% 
  mutate(monthname=month.name[month]) %>% 
  mutate(monthname=factor(monthname,levels=c("July","August","September"))) %>% 
  ggplot(aes(utm.lon.km,utm.lat.km))+
  geom_point()+
  coord_equal()+
  labs(x="Eastings",y="Northings")+
  facet_grid(year~monthname)
