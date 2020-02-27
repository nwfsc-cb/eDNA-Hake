library(sf)
library(raster)
library(dplyr)
library(spData)
#library(spDataLarge)
library(tmap)
library(ggplot2)

library(rnaturalearth)
library(rnaturalearthhires)
library(rnaturalearthdata)

dat.ctd.raw   <- read.csv("./Data/CTD_hake_data_10-2019.csv")
dat.water <- read.csv("./Data/eDNA hake water samples.csv")

# Clean up lat-longs for plotting.
#unique(nchar(as.character(dat.ctd$SciGPS.Lon)))

dat.ctd.raw$lat <- as.numeric(substr(dat.ctd.raw$SciGPS.Lat,1,2)) +
                      as.numeric(substr(dat.ctd.raw$SciGPS.Lat,3,9))/60 
dat.ctd.raw$lon <- - (as.numeric(substr(dat.ctd.raw$SciGPS.Lon,1,3)) +
                      as.numeric(substr(dat.ctd.raw$SciGPS.Lon,4,10))/60)

dat.ctd <- dat.ctd.raw %>% dplyr::select(Date,Time,Button,Station=Station..,lat,lon) %>% 
            mutate(Station=as.character(Station)) %>%
            filter(Button=="CTD at Depth") 
# This adds in location for 
dat.ctd <- dat.ctd.raw %>% dplyr::select(Date,Time,Button,Station=Station..,lat,lon) %>% 
              mutate(Station=as.character(Station)) %>%
              filter(Station=="49-9") %>% bind_rows(.,dat.ctd) 

# Manually add back in two stations not associated with CTD casts.

dat.ctd.extra <- data.frame(Station= c("78-MT506", "77-MT505"),
                            lat= c(47.13849,47.0382), 
                            lon=c(-124.54172, -124.57494))
dat.ctd <- full_join(dat.ctd,dat.ctd.extra)


dat.water <- dat.water %>% rename(Station=CTD.cast) %>% 
                mutate(depth.mod = depth,depth.mod = ifelse(depth.mod=="sfc",0,depth.mod))

dat.loc <- left_join(dat.water,dat.ctd)

dat.loc %>% filter(is.na(lon))

### NEED to figure out locations for "49-9", "78-MT506", "77-MT505"

### Make a simple plot of the locations that were sampled for eDNA.

world <- ne_countries(scale="large",returnclass="sf")
sf::st_crs(world)
plot.loc <- ggplot(data=world) +
  
  # geom_contour(data = b %>%  filter (z < 0), 
  #              aes(x=x, y=y, z=z),
  #              breaks=c( -100, -500),
  #              size=c(0.3),
  #              colour="grey")  +
  
  # geom_raster(data = b %>%  filter (z < 0), 
  #             aes(x=x, y=y, fill=z)) +
  geom_raster(data = slope,
              aes(x = x, y = y , fill = slope)) +
  geom_sf() +
  # scale_fill_gradient2(low = "darkblue",
  #                      mid = "white",
  #                      high = "yellow",
  #                      midpoint = -500) +
  
  geom_text_contour (data = b %>%  filter (z < 0, y < 47), 
                     aes(x=x, y=y,   z = z),
                     breaks=c(-100, -500),
                     colour = "black") +
  
  coord_sf(xlim=c(-127,-122),ylim=c(36,50),expand=F) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  geom_point(data=dat.loc,aes(y=lat,x=lon),col="red",alpha=0.1) +
  theme(axis.text.x = element_text(angle=90,vjust= 0.5))

plot.loc
quartz(file="./Plots and Figures/CTD Locations.pdf",height=8,width=3,type="pdf",dpi=300)
  print(plot.loc)
dev.off()


tm_shape(us_states) + tm_polygons() +
tm_shape(rivers) + tm_lines(col="blue") +
tm_shape(land) + tm_raster("elevation",palette= terrain.colors(10))


plot(dat.loc$lon,dat.loc$lat)