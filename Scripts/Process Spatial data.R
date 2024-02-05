# ---- PULL IN SPATIAL DATA TO MAKE GRID ESTIMATE ON AND 
# ---- TO PROJECT TO 

#str_name <- paste0(base.dir,"Data/raster_grid_blake/fivekm_grid.tif")
dat_raster = raster(here("Data","raster_grid_blake","fivekm_grid.tif"))
dat_raster_extracted <- rasterToPoints(dat_raster)

# Get depth information.  This is for the broad 5km x 5 km grid.
raster_5km_depth <- read.csv(here("Data","raster_grid_blake","weighted_mean_NGDC_depths_for_5km_gridcells.csv"))
raster_5km_depth$depth_m <-  - raster_5km_depth$WM_depth_m

# Get fine scale depth information for specific points (100m-ish scale)
raster_CTD_depth <- read.csv(here("Data","raster_grid_blake","CTD_locations_2019to21_with_depth.csv"))
raster_CTD_depth$depth_m <-  - raster_CTD_depth$NGDC_m

# MERGE IN THESE DEPTHS TO THE OBSERVED SAMPLES
dat.samp <- dat.samp %>% left_join(.,
                                   raster_CTD_depth %>% dplyr::select(year,month,day,station,transect,
                                                                     bottom.depth.NGDC=depth_m))



# Merge dat_raster_extracted and associated depths (5km grid)
dat_raster_all <- dat_raster_extracted %>% as.data.frame() %>%
                      rename(Gridcell_ID = fivekm_grid) %>%
                      left_join(.,raster_depth) %>% 
                      mutate(utm.lon = x/1000,utm.lat=y/1000)
# Trim the location included to include only reasonable-ish locations based on observed sample locations.
# First make a box around the locations in the observed sample locations.

BUF <- 10
lon_trim <- c(dat.samp$utm.lon %>% min() - BUF,dat.samp$utm.lon %>% max() + BUF)
lat_trim <- c(dat.samp$utm.lat %>% min() - BUF,dat.samp$utm.lat %>% max() + BUF)
dat_raster_trim <- dat_raster_all %>% 
                      filter(utm.lon>lon_trim[1],utm.lon<lon_trim[2]) %>%
                      filter(utm.lat>lat_trim[1],utm.lat<lat_trim[2])

# add latitude and logitude back into this trimmed data frame

dat_raster_A <- NULL
dat_raster_trim <- dat_raster_trim %>% as.data.frame()
dat_raster_trim$y.km <- dat_raster_trim$y/1000
dat_raster_trim$x.km <- dat_raster_trim$x/1000

# Keep the point in bounds for each year and then make a master list of points included for all years.
ac.lims <- list()
yrs <- unique(acoustic.lims$year)

for(j in 1:length(yrs)){
  NOM <- as.name(yrs[j])
  ac.lims[[NOM]] <- NULL
  ac.tmp <- acoustic.lims %>% filter(year == yrs[j])
  
  
  for(i in 2:nrow(ac.tmp)){
    filt.temp <-  ac.tmp[(i-1):i,] %>% arrange(mean.lat.utm)
    temp <- dat_raster_trim %>% filter(y.km >= filt.temp$mean.lat.utm[1], 
                                          y.km < filt.temp$mean.lat.utm[2])
    uni.lat <- unique(temp$y.km) 
  
    SLOPE.min <-    (filt.temp$min.lon.utm[2] - filt.temp$min.lon.utm[1]) /
      (filt.temp$mean.lat.utm[2] - filt.temp$mean.lat.utm[1])
    SLOPE.max <-    (filt.temp$max.lon.utm[2] - filt.temp$max.lon.utm[1]) /
      (filt.temp$mean.lat.utm[2] - filt.temp$mean.lat.utm[1])
    
    uni.lat <- uni.lat %>% as.data.frame() %>% rename(y.km = ".") %>% 
      mutate(lon.min = SLOPE.min * (y.km -  filt.temp$mean.lat.utm[1]) + filt.temp$min.lon.utm[1]) %>%
      mutate(lon.max = SLOPE.max * (y.km -  filt.temp$mean.lat.utm[1]) + filt.temp$max.lon.utm[1])
  
    temp <- left_join(temp,uni.lat,by="y.km") %>% filter(x.km > lon.min,x.km < lon.max)
    ac.lims[[NOM]] <- bind_rows(ac.lims[[NOM]],temp)
  }
}
  
dat_raster_A <- ac.lims[[1]]
for(i in 2:length(yrs)){
  dat_raster_A <- rbind(dat_raster_A,ac.lims[[i]])
}

dat_raster_trim <- dat_raster_A %>% distinct(Gridcell_ID) %>% 
            left_join(.,dat_raster_trim)

ggplot(acoustic.lims) + 
  geom_point(aes(x=min.lon.utm,y=mean.lat.utm)) +
  geom_point(aes(x=max.lon.utm,y=mean.lat.utm),col=2) +
  geom_point(data=dat_raster_trim, aes(x=x.km,y=y.km),col="blue") +
  facet_wrap(~year)


# Add lat and lon to dat_raster for help with plotting later
PROJ.txt <- dat_raster@srs
proj      <- SpatialPointsDataFrame(coords = dat_raster_trim %>% ungroup() %>% dplyr::select(x,y),
                                    data = dat_raster_trim,
                                    proj4string = CRS(as.character(PROJ.txt)))
proj.latlon <- sp::spTransform(proj, CRSobj = CRS("+proj=longlat"))

A <-  data.frame(Gridcell_ID= proj.latlon$Gridcell_ID,
                 data.frame(proj.latlon@coords))
colnames(A)[2:3] <- c("lon","lat")
dat_raster_trim <- left_join(dat_raster_trim, A) %>% dplyr::select(-y.km,-x.km)

ggplot(acoustic.lims) + 
  geom_point(aes(x=min.lon.utm,y=mean.lat.utm)) +
  geom_point(aes(x=max.lon.utm,y=mean.lat.utm),col=2) +
  geom_point(data=dat_raster_trim, aes(x=x.km,y=y.km),col="blue") +
  facet_wrap(~year)

