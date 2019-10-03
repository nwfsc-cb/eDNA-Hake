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

dat.ctd <- dat.ctd.raw %>% dplyr::select(Date,Time,Button,Station=Station..,lat,lon) %>% mutate(Station=as.character(Station)) %>%
            filter(Button=="CTD at Depth")
dat.water <- dat.water %>% rename(Station=CTD.cast) %>% 
                mutate(depth.mod = depth,depth.mod = ifelse(depth.mod=="sfc",0,depth.mod))

dat.loc <- left_join(dat.water,dat.ctd)



dat.loc %>% filter(is.na(lon)==T)



### Make a simple plot of the locations that were sampled for eDNA.

world <- ne_countries(scale="large",returnclass="sf")

plot.loc <- ggplot(data=world) +
  geom_sf() +
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