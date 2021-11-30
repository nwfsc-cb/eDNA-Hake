### Go get bathymetric data from NOAA to overlay on the transects.

limits.for.map <- dat.acoustic %>% 
                    summarise_at(c("lat", "lon"), .funs = list( Max = ~ round( max(.x) + 1,0) ,
                                                                Min = ~ round( min(.x) - 1,0))) 

b = getNOAA.bathy(lon1 = limits.for.map["lon_Max"],
          lon2 = limits.for.map[ "lon_Min" ],
          lat1 = limits.for.map["lat_Min"],
          lat2 = limits.for.map["lat_Max"],
          resolution = 1,keep = TRUE)
#b <- fortify(b)

# find the first and last lat-lons from dat.acoustic data frame for each transect.
dat.bathy <- dat.acoustic %>% group_by(transect) %>% filter(lon==min(lon) | lon==max(lon)) %>% 
                    dplyr::select(transect,lon,lat) %>% arrange(transect,desc(lon)) %>%
                    group_by(transect) %>%
                    summarise_at(c("lon","lat"), .funs = list( Max = ~ max(.x)  ,
                                             Min = ~ min(.x) ,
                                             Mean = ~ mean(.x))) %>% group_by(transect)

bathy.transects <- NULL
for(i in 1:nrow(dat.bathy)){
  dat.trans <- get.transect(b,
                          x1=dat.bathy$lon_Max[i],x2=dat.bathy$lon_Min[i],
                          y1=dat.bathy$lat_Mean[i],y2=dat.bathy$lat_Mean[i],distance=TRUE)
  bathy.transects <- dat.trans %>% mutate(transect = dat.bathy$transect[i]) %>% bind_rows(bathy.transects,.)
  # Add a final point to each transect so that the geom_polygong command works without weirdness in the other scripts
  bathy.transects <- bind_rows(bathy.transects, 
                              data.frame(lon=-99,lat=-99,dist.km=0,
                                         depth=min(dat.trans$depth),transect=dat.bathy$transect[i]))
  }

###############

