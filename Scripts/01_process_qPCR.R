# This is the main file for reading in data, calling cleaning scripts for each year, 
#  and then combining the results

# This script is for reading in qPCR data, merging it with locations, and starting a set of plots.  
# It will likely be called in a rMarkdown file later.

rm(list=ls())

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
library(brms)
library(sdmTMB)
library(loo)
library(here)

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

# relevant checks on dimension and names for matrices
identical(names(dat.samp.2019),names(dat.samp.2021))
identical(names(dat.stand.2019),names(dat.stand.2021))
identical(names(PCR.2019),names(PCR.2021))
identical(names(dat.control.field.neg.2021),names(dat.control.field.neg.2019))
identical(names(dat.control.2021),names(dat.control.2019))
identical(names(dat.inhibit.2021),names(dat.inhibit.2019))

##### Combine outputs into single data.frame frame across years.
dat.samp <- rbind(dat.samp.2019,dat.samp.2021)
dat.stand <- rbind(dat.stand.2019%>% mutate(year=2019),dat.stand.2021%>% mutate(year=2019))
PCR <- rbind(PCR.2019 %>% mutate(year=2019),PCR.2021 %>% mutate(year=2021))
dat.control.field.neg <- rbind(dat.control.field.neg.2019,dat.control.field.neg.2021 )
dat.control <- rbind(dat.control.2019,dat.control.2021 )
dat.inhibit <- rbind(dat.inhibit.2019 %>% mutate(year=2019), dat.inhibit.2021 %>% mutate(year=2021) )

##################################################################
# Call raw data plotting scripts.
source(here('Scripts',"plot_raw_observations.R"),local=T)
##################################################################

##################################################################
# Make design matrices needed for fitting TMB model
# this includes:
# ------ spatial meshes and design matrices
# ------ fixed effect matrices
# ------ smooth design matrices
# ------ random effect matrices 












