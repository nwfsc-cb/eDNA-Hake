# JOIN MULTIPLE YEARS OF PROCESSED QPCR DATA, INCLUDING UNKNOWN (FIELD) SAMPLES, STANDARDS
# right now, for eulachon

library(tidyverse)
library(here)

# import samples data
dat19 <- read_rds(here('Data','eulachon qPCR 2019 joined cleaned 11_15_2023.rds'))
dat21 <- read_rds(here('Data','eulachon qPCR 2021 joined cleaned.rds'))

# import standards data
stand19 <- read_rds(here('Data','eulachon qPCR 2019 standards cleaned 11_20_2023.rds'))
stand21 <- read_rds(here('Data','eulachon qPCR 2021 standards cleaned 11_20_2023.rds'))

## Thin, clean, filter, bind samples data
## For now, we include only the necessary identifiers and minimal covariates

# 2019 data
glimpse(dat19)
# samples by depth category
dat19 %>% count(depth_cat)
stations_depths_count_2019 <- dat19 %>% count(station,depth_cat) # there should be a lot of 6s (2 bio reps x 3 tech. reps)
# Some depth/station combinations have more, because they went through dilution tests etc. For now, we will keep all uninhibited samples (se filtering below)

dat19_thin <- dat19 %>%
  filter(is.na(Zymo)) %>% 
  dplyr::select(date,year,month,day,time,station,lat,lon,utm.lon.m,utm.lat.m,qPCR,inhibition_rate,dilution,Ct,copies_ul,depth_cat,bathy.bottom.depth,transect_dist_m)
dat19_thin %>% count(dilution)

# 2021 data
glimpse(dat21)
# samples by depth category
dat21 %>% count(depth_cat)
stations_depths_count_2021 <- dat21 %>% count(station,depth_cat)
dat21_thin <- dat21 %>% 
  dplyr::select(date,year,month,day,time,station,lat,lon,utm.lon.m,utm.lat.m,qPCR,inhibition_rate,dilution,Ct,copies_ul,depth_cat,bathy.bottom.depth,transect_dist_m) %>% 
  mutate(inhibition_rate=as.numeric(inhibition_rate))

datjoin <- dat19_thin %>% 
  bind_rows(dat21_thin)
# 16886 samples

# count of samples by year, station, and depth category. Note that this still includes inhibited samples and various controls and tests
stations_depths_count_all_raw <- datjoin %>% count(year,station,depth_cat)

# Filter (but keep track of what we filter out for later discussion)
# first filter= only keep non-inhibited samples
datjoin <- datjoin %>%
  filter(inhibition_rate<0.5)
# 14498 samples

# count NAs across vars
NA_tracker <- datjoin %>% summarise(across(everything(), ~sum(is.na(.))))
NA_tracker
# most vars have 2084 rows with NA- look at these
NADate <- datjoin %>% filter(is.na(date))# 2084 obs
# a lot of these are non-template controls or other types of methodological controls, but not 100% sure
datjoin <- datjoin %>% 
  filter(!is.na(date))
# 12414 samples

NA_tracker2 <- datjoin %>% summarise(across(everything(), ~sum(is.na(.))))
NA_tracker2

# depth category is NA for some other controls
NAdepthcat <- datjoin %>% filter(is.na(depth_cat))

datjoin <- datjoin %>% 
  filter(!is.na(inhibition_rate),!is.na(depth_cat))
# 12293 samples

NA_tracker3 <- datjoin %>% summarise(across(everything(), ~sum(is.na(.))))
NA_tracker3

# dilution is NA for many samples, but NAs should be 1. But keep a record of them in case we need for later QA/QC
NAdilut <- datjoin %>% filter(is.na(dilution))

datjoin <- datjoin %>% 
  mutate(dilution=ifelse(is.na(dilution),1,dilution))
NA_tracker4 <- datjoin %>% summarise(across(everything(), ~sum(is.na(.))))
NA_tracker4

stations_depths_count_all_filt<- datjoin %>% count(year,station,depth_cat)

# Set some final variable names
datjoin_final <- datjoin %>% 
  rename(plate=qPCR) %>% 
  mutate(Ct=na_if(Ct,-99)) %>% 
  dplyr::select(-inhibition_rate)

# Thin, clean, bind standards data
stand19_thin <- stand19 %>% 
  dplyr::select(qPCR,Ct,known_conc_ul=copies_ul) %>% 
  mutate(year=2019)
stand21_thin <- stand21 %>% 
  dplyr::select(qPCR,Ct,known_conc_ul=copies_ul) %>% 
  mutate(year=2021)
standjoin <- stand19_thin %>% 
  bind_rows(stand21_thin) %>% 
  rename(plate=qPCR) %>% 
  mutate(Ct=na_if(Ct,-99))
glimpse(standjoin)

# some quick visualizations
plot_theme <- theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA))
theme_set(plot_theme)

# Locations
all_samp_locs <- datjoin %>% 
  ggplot(aes(utm.lon.m/1000,utm.lat.m/1000))+
  geom_point(size=0.25)+
  coord_equal()+
  labs(x="Eastings",y="Northings")+
  facet_wrap(~year)+
  theme(axis.text.x=element_text(angle=45))
all_samp_locs

# Depths
samps_depths <- datjoin_final %>%
  group_by(year,depth_cat) %>% 
  summarise(num_samp=n(),num_pos=sum(Ct>0,na.rm=T)) %>% 
  ungroup()
samps_depths

stand.curves.all <- standjoin %>% 
  ggplot(aes(log10(known_conc_ul),Ct,color=plate))+
  geom_point()+
  guides(color='none')+
  geom_smooth(method='lm',se=F)+
  facet_wrap(~year)+
  labs(x="Log 10 (Conc.)",y="Ct",title="Standard Curves")
stand.curves.all

# Save!
write_rds(standjoin,here('Data','eulachon qPCR 2019 and 2021 standards clean.rds'))
write_rds(datjoin_final,here('Data','eulachon qPCR 2019 and 2021 samples clean.rds'))
