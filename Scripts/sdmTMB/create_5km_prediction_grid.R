library(tidyverse)
library(sf)
library(terra)
library(here)

# Blake's raster with just grid cell ID numbers
dat_raster=rast(here('Data','raster_grid_blake','fivekm_grid.tif'))
dat_samps <- read_rds(here('Data','eulachon qPCR 2019 and 2021 samples clean.rds'))

# Depth associated with grid cell IDs
raster_depth <- read_csv(here('Data','raster_grid_blake','weighted_mean_NGDC_depths_for_5km_gridcells.csv')) %>% 
  rename(depth_m=WM_depth_m)

# join depths to cell numbers
cellnums <- tibble(rastID=as.numeric(values(dat_raster))) %>% 
  mutate(rastcell=row_number()) %>% 
  left_join(raster_depth,by=c("rastID"="Gridcell_ID"))

# make bathy raster
rast_depth <- setValues(dat_raster,cellnums$depth_m)
plot(rast_depth)

# as a tibble for prediction
grid.pred <- as.data.frame(rast_depth,xy=T) %>% 
  mutate(utm.lat.km=y/1000,utm.lon.km=x/1000,bathy.bottom.depth=fivekm_grid)
# reduce to a bounding box that matches the mesh
predbbox <- dat_samps %>% 
  dplyr::select(utm.lon.m,utm.lat.m) %>% 
  st_as_sf(coords=c("utm.lon.m","utm.lat.m"),crs=crs(dat_raster)) %>% 
  st_bbox()
grid.pred <- grid.pred %>% 
  st_as_sf(coords=c("x","y"),crs=crs(dat_raster),remove = F) %>% 
  st_crop(predbbox) %>% 
  st_set_geometry(NULL)

ggplot(grid.pred)+
  geom_point(aes(utm.lon.km,utm.lat.km),size=0.5)+
  coord_equal()
