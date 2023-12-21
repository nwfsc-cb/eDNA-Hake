# Plots of raw Hake data


## MAKE PLOTS OF DATA.

dat.samp.mod <- dat.samp %>% mutate(copies_plus = copies_ul / (dilution*vol.standard)) %>% 
  group_by(year,qPCR,transect,sample,depth_cat,vol.standard,
           bathy.bottom.depth,bottom.depth.consensus,
           transect.dist.km,
           dilution,station,lon, lat) %>%
  summarise(Mean= mean(copies_plus))


# Simple Mean plots
p1_hake_by_year <- base_map_trim + 
  geom_point(data=dat.samp.mod,aes(x=lon,y=lat),shape=21,color="black",size=0.1)+
  geom_point(data=dat.samp.mod,aes(x=lon,y=lat,size=Mean),shape=21,color="red") +
  scale_size("Copies / uL",labels=LAB,breaks=LAB,range=c(0.2,20),limits=c(0.5,NA))+
  # geom_point(data=SAMPLES%>% filter(Mean<=lower.lim),aes(x=lon,y=lat),shape="x",size=1,color="black") +
  # scale_shape("Copies / uL",solid=FALSE)+
  # scale_color_viridis_d("Inhibition\nCategory",begin=0,end=0.8) +
  facet_grid(year~depth_cat)


quartz(file = here("Plots and figures","2019, 2021, Hake raw copies.pdf"),  
       type = "pdf",dpi=300,width = 11,height=8.5)
print(p1_hake_by_year)
dev.off()

quartz(file = here("Plots and figures","2019, 2021, Hake acoustic ouput.pdf"),  
       type = "pdf",dpi=300,width = 6,height=6)
print(hake.acoustic.p1)
dev.off()

### Standards plots
dat.stand <- dat.stand %>% mutate(loge_copies = log(copies_ul))

ggplot(dat.stand %>% filter(Ct>0)) +
  geom_point(aes(x=log10_copies,y=Ct))

ggplot(dat.stand %>% filter(Ct>0)) +
  geom_point(aes(x=loge_copies,y=Ct))
