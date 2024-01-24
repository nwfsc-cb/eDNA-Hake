# Base Map script.

#### Plotting information
states <- map_data("state")

canada  <- map_data("world","Canada")
mexico  <- map_data("world","Mexico")
#canada <- subset(world, region %in% c("canada"))

west_coast <- subset(states, region %in% c("california", "oregon", "washington","nevada"))

# lat.lims <- c(min(dat.id$lat,na.rm=T),max(dat.id$lat,na.rm=T))
# lon.lims <- c(min(dat.id$lon,na.rm=T),max(dat.id$lon,na.rm=T))
lat.lims.full.trim <- c(34,54)
lon.lims.full.trim <- c(-134,-120)

lat.lims.trim <- c(35,48)
lon.lims.trim <- c(-126.5,-120.5)

lat.lims.trim.2022 <- c(38.25,48)
lon.lims.trim.2022 <- c(-126.5,-122.5)

# base_map <-ggplot(data = west_coast) + 
#   geom_polygon(aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
#   coord_fixed(xlim=lon.lims,ylim=lat.lims,ratio=1.3) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   theme_bw()

base_map_whole_coast <-ggplot(data = west_coast) + 
  geom_polygon(data = canada,aes(x = long, y = lat, group=group), fill = grey(0.7), color = "black")+
  geom_polygon(data = mexico,aes(x = long, y = lat, group=group), fill = grey(0.7), color = "black")+
  geom_polygon(data = west_coast,aes(x = long, y = lat, group=group), fill = grey(0.7), color = "black")+
  coord_fixed(xlim=lon.lims.full.trim,ylim=lat.lims.full.trim,ratio=1.2) +  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

base_map_trim <-ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
  geom_polygon(data = canada,aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
  coord_fixed(xlim=lon.lims.trim,ylim=lat.lims.trim,ratio=1.2) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()


lat.lims.trim.proj <- c(38,48.1)
lon.lims.trim.proj <- c(-126.55,-122)
lon.lims.trim.proj.hake <- c(-126.55,-120.5)

base_map_trim_proj <-ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
  geom_polygon(data = canada,aes(x = long, y = lat, group=group), fill = grey(0.5), color = "black")+
  coord_fixed(xlim=lon.lims.trim.proj,ylim=lat.lims.trim.proj,ratio=1.2) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

base_map_trim_proj

lat.lims.trim.2022 <- c(33.2,48.6)
lon.lims.trim.2022 <- c(-126.55,-117.0)
#-118.4

base_map_trim_2022 <-ggplot() + 
  geom_polygon(data = west_coast,aes(x = long, y = lat, group=group), fill = grey(0.7), color = "black")+
  geom_polygon(data = canada,aes(x = long, y = lat, group=group), fill = grey(0.7), color = "black")+
  geom_polygon(data = mexico,aes(x = long, y = lat, group=group), fill = grey(0.7), color = "black")+
  coord_fixed(xlim=lon.lims.trim.2022,ylim=lat.lims.trim.2022,ratio=1.5) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

base_map_trim_2022





