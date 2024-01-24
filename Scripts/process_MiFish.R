
rm(list=ls())
#setwd('./Github/eDNA-Hake')
# Libraries
library(tidyverse)
library(marmap)
library(ggplot2)
library(rstan)
library(lubridate)
library(reshape2)
library(gridExtra)
library(raster)
#library(rgdal)
library(sp)
library(here)
library(tidyverse)
library(scatterpie)
library(PNWColors)

#Read in Taxon Tables.

tabs <- dir(here("Data","MiFish","csv"))

dat.long <- NULL

for(i in 1:length(tabs)){
  A <- read.csv(here("Data","MiFish","csv",tabs[i]))
  if(colnames(A)[1] =="X"){
    A <- A %>% dplyr::select(-X)
  }
  seq.run <- sub("\\_.*", "", tabs[i])
  
  B <- pivot_longer(A,cols= - Species ,values_to="count",names_to="ID")
  B <- B %>% mutate(seq.run= seq.run)
  
  dat.long <- rbind(dat.long,B) 
}

dat.long <- dat.long %>% mutate(project.ID = substr(ID,1,9),
                                sample = substr(ID,11,nchar(ID)),
                                sample = gsub("\\..*", "", sample) ) %>%
                        filter(!sample == c(""),!sample =="ve") %>%# Drop controls.
                        mutate(seq.samp=paste0(seq.run,"_",sample))
#Make sure that all species x run x sample are padded into the data.frame

# Remove a couple duplicates that were irritating to find

dat.long %>% filter(sample==1542) %>% as.data.frame()

dat.long <- dat.long %>% filter(!ID == "MFU.52193.1485.d1") %>%
                        filter(!ID == "MFU.52193.1542.d5") %>%
                        filter(!ID == "MFU.52193.1836.d1") %>%
                        filter(!ID == "MFU.52193.1843.d5") %>%
                        filter(!ID == "MFU.52193.1867.d1") %>%
                        filter(!ID == "MFU.52193.1879.d1")

all_comb <- expand_grid(dat.long %>% distinct(sample,seq.run),
                        Species=unique(dat.long$Species)) 
dat.long <- dat.long %>% left_join(all_comb,.) %>% 
                    mutate(count=ifelse(is.na(count),0,count))


# Summarize in a few ways

dat.long <- dat.long %>% group_by(seq.run,sample) %>% 
                  mutate(tot.reads = sum(count)) %>%
                  ungroup() %>%
                  mutate(prop = count / tot.reads)

# Make a tile plot by seq.run

p_tile_plot <- ggplot(dat.long) +
      geom_tile(aes(y=Species,x=seq.run,color=count))


quartz(file = here("Plots and figures","Tile_plot.jpeg"),  
       type = "jpeg",dpi=300,width = 6,height=8.5)
      print(p_tile_plot)
dev.off()



### Pull in the various auxiliary location, depth, other datasets.

# This is the main file for reading in data, calling cleaning scripts for each year, 
#  and then combining the results

# This script is for reading in qPCR data, merging it with locations, and starting a set of plots.  
# It will likely be called in a rMarkdown file later.

###########################################################################
## DECLARED SPECIES OF INTEREST
SP <- "hake" # options: hake, lamprey, eulachon

# Specify an inhibition limit for retaining samples.
INHIBIT.LIMIT <- 0.5

# load and run the acoustic data. this is needed to reference the offshore-ness of 
#setwd(script.dir)

source(here('Scripts',"process acoustic data for qPCR 2019 on.R"),local=TRUE)
# dat.acoustic and dat.acoustic.binned are the relevant data frames

# Read in 2019 data, clean, add some indexes.
source(here('Scripts',"process 2019 qPCR results.R"),local=TRUE)

# Read in 2021 data, clean, add some indexes.
source(here('Scripts',"process 2021 qPCR results.R"),local=TRUE)

##### Combine outputs into single data.frame frame across years.
dat.samp <- rbind(dat.samp.2019,dat.samp.2021) %>%
  mutate(qPCR_year = paste0(qPCR,"_",year)) %>%
  dplyr::select(-plate_idx)
dat.stand <- rbind(dat.stand.2019%>% mutate(year=2019),dat.stand.2021%>% mutate(year=2021)) %>% 
  mutate(qPCR_year = paste0(qPCR,"_",year)) %>% arrange(year,qPCR_year) %>%
  mutate(ln_copies = log(copies_ul)) %>%
  dplyr::select(-plate_idx)

# Read in spatial data and data for bottom depth at each sample
# REPLACES THE SPATIAL DATA ANALYSIS DONE IN THE INDIVIDUAL YEAR SCRIPTS.
source(here('Scripts',"Process Spatial data.R"),local=TRUE)
# relevant files here are dat_raster_trim an the new depth data in dat.samp.



##############################################################################
# OK Merge in the sample data with the Mifish
##############################################################################
# CPS sub set

CPS <- c("Engraulis mordax", "Sardinops sagax", "Trachurus symmetricus",
         "Scomber japonicus","Clupea pallasi")

dat.samp.mifish.all <- dat.samp %>% distinct(year,sample,dilution, month, day,
                                              station,station.depth,Niskin,depth_cat,
                                              utm.lon,utm.lat,lon,lat,wash.indicator,
                                              bottom.depth.NGDC) 
                                

dat.samp.mifish <- left_join(dat.long,dat.samp.mifish.all %>%filter(year==2019))


# Where the full samples are in 2019 and 2021

p_mifish_loc <- base_map_trim +
    geom_point(data=dat.samp.mifish.all,aes(x=lon,y=lat)) + 
    geom_point(data=dat.samp.mifish,aes(x=lon,y=lat),col="red") + 
    theme_bw() +
    facet_grid(year~depth_cat)

quartz(file = here("Plots and figures","2019, 2021, mifish location.jpeg"),  
       type = "jpeg",dpi=300,width = 11,height=8.5)
print(p_mifish_loc)
dev.off()

### Read Depths.

A <- dat.samp.mifish %>% group_by(sample) %>% summarise(Total_reads = mean(tot.reads)+1)

ggplot(A) +
    geom_histogram(aes(Total_reads)) +
    scale_x_continuous("Total reads + 1 (log scale)",trans="log10",breaks=c(0,1,10,100,1000,10000,100000,1000000)) +
    theme_bw()

# Filter out those with more than 1000 reads

these <- A %>% filter(Total_reads>1000)

dat.samp.mifish %>% distinct(sample,seq.run) %>% as.data.frame()

dat.samp.mifish <- dat.samp.mifish %>% filter(sample %in% these$sample)
dat.samp.mifish %>% distinct(sample,seq.run) %>% as.data.frame()

##################
# Make Presence-Absence maps for all Species.
p_PA <- list()

SP <- dat.samp.mifish %>% distinct(Species) %>% arrange(Species)
dat.samp.mifish <- dat.samp.mifish %>% mutate(bin = ifelse(count>1,1,0), bin = as.factor(bin))

for(i in 1:length(SP$Species)){
  p_PA[[SP$Species[i]]] <- base_map_trim +
                  geom_point(data=dat.samp.mifish %>% filter(Species == SP$Species[i]),
                             aes(x=lon,y=lat,fill=bin,color=bin), alpha=0.5,shape=21) +
                        scale_fill_manual("Occurrence",values=c("black","red")) +
                        scale_color_manual("Occurrence",values=c("black","red")) +
                        ggtitle(SP$Species[i])+
                        facet_grid(year~depth_cat)
}

pdf(file = here("Plots and figures","Mifish all Pres-Abs.pdf"),onefile=TRUE,
      width = 11,height=7)
  print(p_PA)
dev.off()

##############################
# Make Plots of proportions for each species.
##############################

CPS <- c("Engraulis mordax", "Sardinops sagax", "Trachurus symmetricus",
         "Scomber japonicus","Clupea pallasii")

# subset data for the 5 CPS species

dat.CPS <- dat.samp.mifish %>% filter(Species %in% CPS) %>% 
                  group_by(sample,seq.run) %>% 
                  mutate(cps.total = sum(count)) %>%
                  ungroup() %>%                            
                  mutate(prop.sp.cps = count / cps.total)  %>% 
                  mutate(prop.all.cps = cps.total / tot.reads,
                         prop.sp.cps = ifelse(is.nan(prop.sp.cps),0,prop.sp.cps))
          

dat.CPS.sum <- dat.CPS %>%
                      group_by(year,depth_cat,sample,seq.run,tot.reads,lon,lat) %>%
                      summarise(prop.all.cps = mean(prop.all.cps))




# Proportion of DNA from the 5 main CPS species out of all reads

p_CPS_prop_5sp <- base_map_trim_proj +
  geom_point(data=dat.CPS.sum,
             aes(x=lon,y=lat,size=prop.all.cps), alpha=0.5,color="red") +
  scale_size_continuous("Proportion", limits = c(0.001,1),breaks=c(0.01,0.1,0.25,0.5,0.75,1.0)) +
  geom_point(data=dat.CPS.sum %>% filter(prop.all.cps==0),
             aes(x=lon,y=lat),shape="x", alpha=1) +
  ggtitle("Prop. DNA from 5 CPS")+
  facet_grid(year~depth_cat)


quartz(file = here("Plots and figures","Prop 5 CPS of total DNA.jpeg"),type="jpeg", dpi=300,
       width = 11,height=7)
  print(p_CPS_prop_5sp)
dev.off()
  
# Proportion of DNA from the hake out of all reads  
dat.hake <- dat.samp.mifish %>% filter(Species =="Merluccius productus")

p_hake_prop <- base_map_trim_proj +
  geom_point(data=dat.hake,
             aes(x=lon,y=lat,size=prop), alpha=0.5,color="blue") +
  scale_size_continuous("Proportion", limits = c(0.001,1),breaks=c(0.01,0.1,0.25,0.5,0.75,1.0)) +
  geom_point(data=dat.hake %>% filter(prop==0),
             aes(x=lon,y=lat),shape="x", alpha=1) +
  ggtitle("Prop. DNA from hake")+
  facet_grid(year~depth_cat)

quartz(file = here("Plots and figures","Prop hake of total DNA.jpeg"),type="jpeg", dpi=300,
    width = 11,height=7)
print(p_hake_prop)
dev.off()


# Make scatterpie plot for CPS - only 

PAL <-pnw_palette("Bay",5)


dat.CPS.wide <- dat.CPS %>% pivot_wider(.,id_cols=c("year","sample","depth_cat","seq.run","lon","lat","prop.all.cps","cps.total","tot.reads"),
                             names_from = "Species",values_from = "prop.sp.cps")


p_pies_all_depth <-   base_map_trim_proj +
  geom_scatterpie(aes(x=lon, y=lat,r = 0.2 * prop.all.cps), alpha=0.75,
                data = dat.CPS.wide, cols = CPS) +
  scale_fill_manual("Species",values=PAL) +
  geom_point(data=dat.CPS.wide %>% filter(prop.all.cps ==0),
             aes(x=lon,y=lat),shape="x", alpha=1) +
  ggtitle("Prop. DNA among CPS")+
  facet_grid(year~depth_cat)

p_pies_0_50 <-   base_map_trim_proj +
  geom_scatterpie(aes(x=lon, y=lat,r = 0.2 * prop.all.cps), alpha=0.75,
                  data = dat.CPS.wide %>% filter(depth_cat %in% c(0,50)), cols = CPS) +
  scale_fill_manual("Species",values=PAL) +
  geom_point(data=dat.CPS.wide %>% filter(prop.all.cps ==0,depth_cat %in% c(0,50)),
             aes(x=lon,y=lat),shape="x", alpha=1) +
  ggtitle("Prop. DNA among CPS")+
  facet_grid(year~depth_cat)


quartz(file = here("Plots and figures","CPS pie all depth.jpeg"),type="jpeg", dpi=300,
       width = 11,height=7)
  print(p_pies_all_depth)
dev.off()


quartz(file = here("Plots and figures","CPS pie 0, 50m depth.jpeg"),type="jpeg", dpi=300,
       width = 8,height=7)
  print(p_pies_0_50)
dev.off()



