# Base Map script.

#### Plotting information
states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington"))

# lat.lims <- c(min(dat.id$lat,na.rm=T),max(dat.id$lat,na.rm=T))
# lon.lims <- c(min(dat.id$lon,na.rm=T),max(dat.id$lon,na.rm=T))
lat.lims.trim <- c(38.25,48)
lon.lims.trim <- c(-126.5,-122.5)


# base_map <-ggplot(data = west_coast) + 
#   geom_polygon(aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
#   coord_fixed(xlim=lon.lims,ylim=lat.lims,ratio=1.3) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   theme_bw()

base_map_trim <-ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
  coord_fixed(xlim=lon.lims.trim,ylim=lat.lims.trim,ratio=1.2) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()



lat.lims.trim.proj <- c(38.6,48.1)
lon.lims.trim.proj <- c(-126.55,-122.75)

base_map_trim_proj <-ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
  coord_fixed(xlim=lon.lims.trim.proj,ylim=lat.lims.trim.proj,ratio=1.2) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

base_map_trim_proj