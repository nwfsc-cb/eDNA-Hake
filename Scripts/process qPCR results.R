# This script is for reading in qPCR data, merging it with locations, and starting a set of plots.  
# It will likely be called in a rMarkdown file later.

# Libraries
library(tidyverse)
library(marmap)
library(ggplot2)

# Working directories
base.dir <- "/Users/ole.shelton/Github/eDNA-Hake/"
data.dir <- paste0(base.dir,"Data")
plot.dir <- paste0(base.dir,"Plots and figures")

# Pull in qPCR data, qPCR standards, sample id information
setwd(data.dir)
dat.all <- read.csv("./qPCR/Hake_eDNA_2019_qPCR_results.csv")
dat.stand <- read.csv("./qPCR/Hake_eDNA_2019_qPCR_standards.csv")
dat.station.id <- read.csv("./cleaned_CTD_stations.csv")
dat.sample.id  <- read.csv("./2019 Hake- Shimada cruise - eDNA water sampling.csv")

# Modify a few of the identifiers.
dat.station.id <- dat.station.id %>% mutate(transect = as.numeric(substr(Station,1,2)))

dat.sample.id  <- dat.sample.id %>% dplyr::select(sample=Tube..,
                                                  Station=CTD.cast,
                                                  Niskin,
                                                  depth,
                                                  volume = water.filtered..L.,
                                                  Fluor,
                                                  Zymo=Zymo.columns)

dat.id <- full_join(dat.sample.id,dat.station.id) %>% rename(station=Station) %>%
              mutate(control = case_when(station=="" ~ "extraction",
                                        station=="-" ~ "field",
                                        TRUE ~ "no")) %>%
              mutate(depth=as.character(depth)) %>%
              mutate(depth = case_when(depth=="-" ~ "-99",
                                       depth=="sfc" ~ "0",
                                       depth=="300/150" ~ "300",
                                       depth=="" ~ "0",
                                       TRUE ~depth))

dat.id.control <- dat.id %>% filter(!control == "no")
dat.id.samp    <- dat.id %>% filter(control == "no") %>% mutate(depth=as.numeric(as.character(depth)))

###########################################################################
## DECLARED SPECIES OF INTEREST
SP <- "hake" # options: hake, lamprey, eulachon
###########################################################################

# STANDARDS 
THESE <- c("qPCR","well","sample","type","IPC_Ct","inhibition_rate")
cut.plates <- c("H5","H6")

dat.stand <- dat.stand %>% dplyr::select(THESE,grep(SP,colnames(dat.stand))) %>% 
                      rename(task=paste0(SP,"_task"),Ct=paste0(SP,"_Ct"),copies_ul = paste0(SP,"_copies_ul")) %>%
                      filter(task=="STANDARD",!qPCR%in%cut.plates)

    # Make sure there are sufficient samples for each qPCR
    stand.summary <- dat.stand %>% group_by(qPCR,copies_ul) %>% summarise(N= length(copies_ul)) %>% as.data.frame()
    PCR <- data.frame(qPCR=unique(dat.stand$qPCR),plate_idx= 1:length(unique(dat.stand$qPCR)))
    dat.stand <- left_join(dat.stand,PCR)


# SAMPLES, ntc and control
dat.samp <- dat.all %>% dplyr::select(THESE,useful,Zymo,grep(SP,colnames(dat.all))) %>%
                  rename(Ct=paste0(SP,"_Ct"),copies_ul = paste0(SP,"_copies_ul")) %>%
                  dplyr::select(-grep("copies_l",colnames(.))) %>%
                  filter(!qPCR %in% cut.plates) %>%
                  mutate(Ct = as.character(Ct),
                         Ct=ifelse(Ct=="",-99,Ct),
                         Ct=ifelse(Ct=="Undetermined",-99,Ct),
                         Ct=ifelse(is.na(Ct)==T,-99,Ct), 
                         Ct=as.numeric(Ct)) %>%
                  mutate(IPC_Ct = as.character(IPC_Ct),
                        IPC_Ct=ifelse(IPC_Ct=="",-99,IPC_Ct),
                        IPC_Ct=ifelse(IPC_Ct=="Undetermined",-99,IPC_Ct),
                        IPC_Ct=ifelse(is.na(IPC_Ct)==T,-99,IPC_Ct), 
                        IPC_Ct=as.numeric(IPC_Ct)) %>%
                  filter(type=="unknowns")
dat.samp <- left_join(dat.samp,dat.id.samp,by="sample")


dat.ntc <- dat.all %>% filter(type=="ntc") %>% mutate(IPC_Ct = as.numeric(as.character(IPC_Ct))) %>%
                      group_by(qPCR) %>% 
                      dplyr::summarise(mean.ntc = mean(IPC_Ct),sd.ntc=sd(IPC_Ct))

dat.control <- dat.all %>% filter(grepl("neg",type),type=="extc")

# Merge in ntc to data to flag inhibited samples 
dat.samp <- left_join(dat.samp,dat.ntc) %>% mutate(mean.ntc = as.numeric(as.character(mean.ntc))) %>%
                  mutate(inhibit.val = IPC_Ct-mean.ntc,
                         inhibit.bin=ifelse(inhibit.val < 1 & inhibit.val > -1,0,1))

