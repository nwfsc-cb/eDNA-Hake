library(tidyverse)
library(rstan)
library(bayesplot)
library(here)
options(mc.cores = parallel::detectCores())
rstan_options("auto_write" = TRUE)

#library prep data, w PCR pipetting info
lib_mifish <- read.csv(here("Data/ole_data_calcofi/mifish_library_prep.csv"))
lib_mifish$New_name <-  paste0(   #need to match format of other dataset
  substr(lib_mifish$New_name, 1, 4), "_",
  substr(lib_mifish$New_name, 8, nchar(as.character(lib_mifish$New_name)))
)

# library prep data, w PCR pipetting info
lib_sebastes <- read.csv(here("Data/ole_data_calcofi/sebastes_library_prep.csv"))
lib_sebastes$New_name <-  paste0(   #need to match format of other dataset
  substr(lib_sebastes$New_name, 1, 4), "_",
  substr(lib_sebastes$New_name, 8, nchar(as.character(lib_sebastes$New_name)))
)

dat.mifish <- readRDS(here("Data/ole_data_calcofi/mifish_tech_nReads.RDS")) %>% 
  filter(!ID_mifish %in% c("MISSING", ""),
         !is.na(mifish_reads)) %>% 
  dplyr::select(Sample, ID_mifish, mifish_reads, "station_id", "ext_rep", "tech_rep", ) %>% 
  distinct() %>% 
  unite(c("station_id", "ext_rep", "tech_rep"), col = "New_name", sep = "_", remove = F) %>% #to join to library-prep info
  filter(ext_rep == 1) %>% 
  left_join(lib_mifish) %>% 
  mutate(name = ID_mifish, nReads = mifish_reads) %>%   #for compatibility w earlier versions
  group_by(Sample) %>% 
  mutate(tot = sum(mifish_reads)) %>% 
  filter(!is.na(tot)) %>%   #omit failed replicates
  filter(!is.na(name)) %>% 
  na.omit() %>% 
  group_by(ID_mifish, Sample, station_id, tech_rep, Percent.of.PCR.added) %>% 
  dplyr::summarise(nReads = sum(nReads)) %>% 
  distinct() %>% 
  mutate(log_r_mifish = log(Percent.of.PCR.added)) %>% 
  ungroup() %>% 
  dplyr::select(-Percent.of.PCR.added)

#dat.mifish$ID_mifish[dat.mifish$ID_mifish == "Lampanyctus tenuiformis"] <- "Nannobrachium ritteri"
dat.mifish <- dat.mifish %>% group_by(ID_mifish, Sample, tech_rep, station_id, log_r_mifish) %>% summarise(nReads = sum(nReads))

common_sp <- dat.mifish %>%
  group_by(ID_mifish) %>%
  dplyr::summarise(tot = sum(nReads)) %>%
  arrange(desc(tot)) %>%
  top_n(50) %>%
  pull(ID_mifish)

#this is a very high percentage of all reads
dat.mifish %>%  filter (ID_mifish %in% common_sp) %>% ungroup() %>% summarize(sum(nReads)) / 
  dat.mifish %>% ungroup() %>% summarize(sum(nReads))

#make a new file dat.mifish.init that has been filtered.  Do not overwrite dat.mifish
dat.mifish.init <- dat.mifish %>% 
  filter(ID_mifish %in% common_sp) #take just common species, for now

dat.sebastes <- readRDS(here("Data/ole_data_calcofi/sebastes_tech_nReads.RDS")) %>% 
  filter(!ID_sebastes %in% c("MISSING", "", "Unidentified"),
         !is.na(sebastes_reads)) %>% 
  dplyr::select(Sample, ID_sebastes, sebastes_reads, "station_id", "ext_rep", "tech_rep", ) %>% 
  distinct() %>% 
  unite(c("station_id", "ext_rep", "tech_rep"), col = "New_name", sep = "_", remove = F) %>% #to join to library-prep info
  filter(ext_rep == 1) %>% 
  left_join(lib_sebastes) %>% 
  mutate(name = ID_sebastes, nReads = sebastes_reads) %>%   #for compatibility w earlier versions
  group_by(Sample) %>% 
  mutate(tot = sum(sebastes_reads)) %>% 
  filter(!is.na(tot)) %>%   #omit failed replicates
  filter(!is.na(name)) %>% 
  na.omit() %>% 
  group_by(ID_sebastes, Sample, station_id, tech_rep, Percent.of.PCR.added) %>% 
  dplyr::summarise(nReads = sum(nReads)) %>% 
  distinct() %>% 
  mutate(log_r_sebastes = log(Percent.of.PCR.added)) %>% 
  ungroup() %>% 
  dplyr::select(-Percent.of.PCR.added)
dat.sebastes$ID_sebastes[dat.sebastes$ID_sebastes == "Sebastes"] <- "Sebastes AANOTHER"  #disambiguate from other, known rockfish species, in a way that preserves alpha order

# Find any sebastes that has zero reads and drop from data.  
drop.these.sebastes <- dat.sebastes %>% group_by(ID_sebastes) %>% summarise(tot = sum(nReads)) %>% filter(tot==0) %>% pull(ID_sebastes)
dat.sebastes <- dat.sebastes %>% filter(!ID_sebastes %in% drop.these.sebastes)

dat.count <- readRDS(here("Data/ole_data_calcofi/microscopy_tech_nReads.RDS")) %>% 
  filter(!ID_microscopy %in% c("MISSING", "", "Disintegrated fish larvae")) %>% 
  filter(station_id %in% unique(dat.mifish$station_id)) %>% 
  #dat.count$larval_counts[dat.count$ID_microscopy == "Vinciguerria sp."] <- 0  #DATA problem; check w Zack
  # dat.count <- dat.count %>% 
  #   mutate(ID_microscopy = ifelse(grepl("Sebastes", ID_microscopy), "Sebastes", ID_microscopy),
  #          ID_mifish = ifelse(grepl("Sebastes", ID_mifish), "Sebastes", ID_mifish),
  #          Unique_ID = ifelse(grepl("Sebastes", Unique_ID), "Sebastes", Unique_ID)
  #          ) %>%   ##TEMPORARY -- synonymize all Sebastes
  group_by(ID_microscopy, station_id, standard.haul.factor, volume.filtered, proportion.sorted) %>% 
  dplyr::summarise(larval_counts = sum(larval_counts)) 

#synonymize all Sebastes in counts
dat.count$ID_microscopy[grepl("Sebastes", dat.count$ID_microscopy)] <- "Sebastes"
dat.count <- dat.count %>% group_by(ID_microscopy, station_id, standard.haul.factor, volume.filtered, proportion.sorted) %>% dplyr::summarise(larval_counts = sum(larval_counts))

#import info to set up prior on Sebastes counts
#rockfish <- read.csv(here("Data/ole_data_calcofi/Thompson_2017_supp.csv", header = F)

#require data from all of the datasets at a given station
keepstations <- intersect(dat.count$station_id, dat.sebastes$station_id) %>% intersect(dat.mifish$station_id)
dat.count <- dat.count %>% filter(station_id %in% keepstations)
dat.sebastes <- dat.sebastes %>% filter(station_id %in% keepstations)
dat.mifish <- dat.mifish %>% filter(station_id %in% keepstations)

common_count <- dat.count %>%
  group_by(ID_microscopy) %>%
  dplyr::summarise(tot = sum(larval_counts)) %>%
  arrange(desc(tot)) %>% 
  top_n(25) %>%
  pull(ID_microscopy)

#this is a very high percentage of all counts
dat.count %>%  filter (ID_microscopy %in% common_count) %>% ungroup() %>% dplyr::summarize(sum(larval_counts)) / 
  dat.count %>% ungroup() %>% dplyr::summarize(sum(larval_counts))
dat.count %>%  filter(., is.na(larval_counts)) %>% View()

#make a new object dat.count.init that has been filtered.  Do not overwrite dat.count
dat.count.init <- dat.count %>%
  filter(ID_microscopy %in% common_count)  #take just common species, for now

# pull appropriate species from microscopy counts
species_mapping <- read.csv(here("Data/ole_data_calcofi/20210622_species_mapping_file.csv")) %>% 
  filter(ID_mifish %in% dat.mifish.init$ID_mifish  | 
           ID_microscopy %in% dat.count.init$ID_microscopy | 
           ID_sebastes %in% dat.sebastes$ID_sebastes) %>% 
  filter(!ID_sebastes %in% drop.these.sebastes) %>%
  arrange(Unique_ID)
names(species_mapping) <- c("ID_master", "ID_mifish", "ID_count", "ID_sebastes")
species_mapping <- species_mapping %>% arrange(ID_master,ID_count,ID_mifish)
species_mapping[species_mapping == ""] <- "MISSING"
species_mapping <- species_mapping %>% filter(species_mapping$ID_master != "Clupeiformes")
species_mapping$ID_count[grepl("Sebastes", species_mapping$ID_count)] <- "Sebastes"   #synonymizing sebastes in counts
#species_mapping <- species_mapping[-which(species_mapping$ID_master == "Sebastes diploproa"),] #species doesn't occur in dataset

# Re-filter dat.count to include only the species found in species_mapping, as a safety check
dat.count <- dat.count %>% 
  ungroup() %>% 
  dplyr::rename(ID_count = ID_microscopy) %>%  
  dplyr::rename(Abund = larval_counts) %>%
  filter(ID_count %in% species_mapping$ID_count)

MAX.VOL <- max(dat.count$volume.filtered)
dat.count <- dat.count %>% ungroup %>%  mutate(vol.ratio = volume.filtered / MAX.VOL,
                                               log_V=log(vol.ratio),
                                               log_P=log(proportion.sorted))
# Re-filter amplicon data too
dat.mifish <- dat.mifish %>% 
  filter(ID_mifish %in% species_mapping$ID_mifish) #This pulls species that were in ID_count or ID_sebastes but not in the top group of mifish.

dat.sebastes <- dat.sebastes %>% 
  filter(ID_sebastes %in% species_mapping$ID_sebastes) #This doesn't do much because we didn't filter the sebastes.

### ONLY KEEP INSTANCES WHERE EITHER READS or COUNTS for one of the observations > 0  
pos.count <- dat.count %>% filter(Abund > 0) %>% dplyr::select(ID_count,station_id) %>%
  left_join(.,species_mapping)  #note number of rows will increase here, because many distinct values of ID_sebastes map to a single value of ID_count
pos.mifish <- dat.mifish %>% ungroup %>% filter(nReads > 0) %>% dplyr::select(ID_mifish,station_id) %>%
  distinct() %>% left_join(.,species_mapping)
pos.sebastes <- dat.sebastes %>% filter(nReads > 0) %>% dplyr::select(ID_sebastes,station_id) %>%
  distinct() %>% left_join(.,species_mapping)

pos.all <- pos.count %>% dplyr::select(ID_master,station_id) %>%
  bind_rows(.,pos.mifish %>% dplyr::select(ID_master,station_id)) %>%
  bind_rows(.,pos.sebastes %>% dplyr::select(ID_master,station_id)) %>%
  distinct() %>% arrange(ID_master,station_id) %>%
  full_join(.,species_mapping %>% dplyr::select(ID_master,ID_count,ID_mifish, ID_sebastes)) 

  # This should weed out any samples with all zeros...  
  dat.count <- semi_join(dat.count,pos.count)
  dat.mifish <- semi_join(dat.mifish,pos.mifish) #CULLING FURTHER ZEROS
  dat.sebastes <- semi_join(dat.sebastes,pos.sebastes) #NOTE treated differently, because of the problem of having no counts for species-specific rockfish
  
  #create species indices
  species_mapping$sp_index_master <- match(species_mapping$ID_master, unique(species_mapping$ID_master))
  species_mapping$sp_index_mifish <-   data.frame(ID_mifish = unique(dat.mifish$ID_mifish)) %>% 
                                            mutate(sp_index_mifish = 1:n()) %>% 
                                            right_join(species_mapping) %>% 
                                            arrange(sp_index_master) %>% 
                                            pull(sp_index_mifish)
  species_mapping$sp_index_count <- data.frame(ID_count = unique(dat.count$ID_count)) %>% 
                                        mutate(sp_index_count = 1:n()) %>% 
                                        right_join(species_mapping) %>% 
                                        arrange(sp_index_master) %>% 
                                        pull(sp_index_count)
  species_mapping$sp_index_sebastes <- data.frame(ID_sebastes = unique(dat.sebastes$ID_sebastes)) %>% 
                                          mutate(sp_index_sebastes = 1:n()) %>% 
                                          right_join(species_mapping) %>% 
                                          arrange(sp_index_master, ID_sebastes) %>% 
                                          pull(sp_index_sebastes)

  
  ## Count number of species/entities 
  I_master <- length(species_mapping$sp_index_master)
  
    SP_count <- species_mapping$ID_count[species_mapping$ID_count!="MISSING"] %>% unique()
    I_count <- length(SP_count)
  
    SP_mifish <- species_mapping$ID_mifish[species_mapping$ID_mifish!="MISSING"] %>% unique()
    I_mifish <- length(SP_mifish)

    SP_sebastes <- species_mapping$ID_sebastes[species_mapping$ID_sebastes!="MISSING"] %>% unique()
    I_sebastes <- length(SP_sebastes)
    
    
      
    ## make a matrix from the master list to each primer / data 
      M_to_count <- matrix(0,I_count,I_master)
      M_to_mifish<- matrix(0,I_mifish,I_master)
      M_to_sebastes<- matrix(0,I_sebastes,I_master)
      
      for(i in 1:I_count){
        #M_to_count[i,which(SP_count[i]==species_mapping$ID_count)] <- 1
        
        M_to_count[species_mapping$sp_index_count[SP_count[i] == species_mapping$ID_count],
                   species_mapping$sp_index_master[SP_count[i] == species_mapping$ID_count]] <- 1
      }
      for(i in 1:I_mifish){
        #M_to_mifish[i,which(SP_mifish[i]==species_mapping$ID_mifish)] <- 1
        M_to_mifish[species_mapping$sp_index_mifish[SP_mifish[i] == species_mapping$ID_mifish],
                    species_mapping$sp_index_master[SP_mifish[i] == species_mapping$ID_mifish]] <- 1
      }
      for(i in 1:I_sebastes){
        #M_to_sebastes[i,which(SP_sebastes[i]==species_mapping$ID_sebastes)] <- 1
        M_to_sebastes[species_mapping$sp_index_sebastes[SP_sebastes[i] == species_mapping$ID_sebastes],
                   species_mapping$sp_index_master[SP_sebastes[i] == species_mapping$ID_sebastes]] <- 1
              }
      
      species_mapping %>%  View()
      colnames(M_to_count)  <- species_mapping$ID_master
      colnames(M_to_mifish) <- species_mapping$ID_master
      colnames(M_to_sebastes) <- species_mapping$ID_master
      rownames(M_to_count)   <- species_mapping %>% filter(!ID_count=="MISSING") %>% dplyr::select(ID_count,sp_index_count) %>%
                                  distinct() %>% pull(ID_count)
      rownames(M_to_mifish)  <- species_mapping %>% filter(!ID_mifish=="MISSING") %>% dplyr::select(ID_mifish,sp_index_mifish) %>%
                                  distinct() %>% pull(ID_mifish)
      rownames(M_to_sebastes)  <- species_mapping %>% filter(!ID_sebastes=="MISSING") %>% dplyr::select(ID_sebastes,sp_index_sebastes) %>%
                                    distinct() %>% pull(ID_sebastes)
      M_to_count %>% View()
      M_to_mifish %>%  View()
      M_to_sebastes %>%  View()
      #write these to file, for checking accuracy
      write.csv(M_to_count, "M_to_count.csv")
      write.csv(M_to_mifish, "M_to_mifish.csv")
      write.csv(M_to_sebastes, "M_to_sebastes.csv")
      
      
    ## CHECK THESE MATRICES. CHARACTERISTICS SHOULD BE:
        # COLSUMS all are 0 or 1.  Most should be 1.
        # ROWSUMS 0, 1, >1 all are possible.
      # does each ID_master map to (at most) a single entity in each dataset?
      MAP.ERROR <- max(colSums(M_to_count),colSums(M_to_mifish),colSums(M_to_sebastes))
      names(which(colSums(M_to_count)>1))  -> names
      species_mapping %>%  filter(., ID_master %in% names)
      
      # does each ID_master map to at least one entity across all datasets?
      MIN.ERROR <- M_to_count %>% 
        rbind(M_to_mifish) %>% 
        rbind(M_to_sebastes) %>% 
        colSums() %>% 
        `==`(0) %>% 
        sum()
        
  ##Technical PCR replicate count and index
    J_all_count <- dat.count %>% 
      group_by(station_id, ID_count) %>% 
      summarise(J=length(station_id)) %>% 
      group_by(station_id) %>% 
      summarise(MIN=min(J),MAX=max(J)) %>% 
      filter(MIN == MAX) %>% 
      dplyr::select(station_id,N_rep = MIN) %>% 
      rename(N_count=N_rep)
    J_count  = max(J_all_count$N_count)
    
    J_all_mifish <-  dat.mifish %>% 
      group_by(station_id,ID_mifish) %>% 
      summarise(J=length(station_id)) %>% # HOw many replicates per fish and station - this should always be
      group_by(station_id) %>% 
      summarise(MIN=min(J),MAX=max(J)) %>% 
      filter(MIN == MAX) %>% # all of them?
      dplyr::select(station_id,N_rep = MIN) %>% 
      rename(N_mifish=N_rep)
    J_mifish  = max(J_all_mifish$N_mifish)
  
    J_all_sebastes <-  dat.sebastes %>% 
      group_by(station_id,ID_sebastes) %>% 
      summarise(J=length(station_id)) %>% 
      group_by(station_id) %>% 
      summarise(MIN=min(J),MAX=max(J)) %>% 
      filter(MIN == MAX) %>% 
      dplyr::select(station_id,N_rep = MIN) %>% 
      rename(N_sebastes=N_rep)
    J_sebastes  = max(J_all_sebastes$N_sebastes)
    
  
  ##Create station/location index
  STATION_ID <- sort(unique(dat.count$station_id)) %>% as.data.frame() %>%
                    rename(station_id=".")
  STATION_ID$station_idx <- 1:nrow(STATION_ID)

  ##Station/location count  
  #note because we required data from all 3 datasets from all sites included, the values of these should all be the same
  K_master = max(nrow(J_all_mifish), nrow(J_all_sebastes))  
  K_mifish = nrow(J_all_mifish)
    #K_mifish_station_vec = J_all_mifish %>% left_join(STATION_ID) %>% pull(station_idx)
  K_count = nrow(J_all_count)
    #K_count_station_vec = J_all_count %>% left_join(STATION_ID) %>% pull(station_idx)
  K_sebastes = nrow(J_all_sebastes)
    #K_sebastes_station_vec = J_all_sebastes %>% left_join(STATION_ID) %>% pull(station_idx)
  
  
  #create index of unique PCR reactions for mifish (station-replicates)
  D_mifish <- dat.mifish %>% 
    ungroup() %>% 
    left_join(STATION_ID) %>% 
    left_join(species_mapping %>% dplyr::select(ID_mifish, sp_index_mifish) %>% distinct()) %>% 
    unite(station_idx, tech_rep, sep = "_", remove = F, col = "temp") %>% 
    mutate(mifish_station_rep_idx = match(temp, unique(temp))) %>% 
    dplyr::select(-temp)
  
  #create index of unique PCR reactions for sebastes (station-replicates)
  D_sebastes <- dat.sebastes %>% 
    ungroup() %>% 
    left_join(STATION_ID) %>% 
    left_join(species_mapping %>% dplyr::select(ID_sebastes, sp_index_sebastes) %>% distinct()) %>% 
    unite(station_idx, tech_rep, sep = "_", remove = F, col = "temp") %>% 
    mutate(sebastes_station_rep_idx = match(temp, unique(temp))) %>% 
    dplyr::select(-temp)

  #add station_index and species index to dat.count
  D_count <- dat.count %>% 
    left_join(STATION_ID) %>% 
    left_join(species_mapping %>% dplyr::select(ID_count, sp_index_count) %>% distinct())
 
  N_mifish_obs <- nrow(D_mifish)
  N_sebastes_obs <- nrow(D_sebastes)
  N_count_obs  <- nrow(D_count)
 
    
  #### Checks that all of the gymnastics didn't get messed up.  Each of these should = 1
  D_count$sp_index_count %>% unique() %>% sort() %>% length()/D_count$sp_index_count %>% max()
  D_mifish$sp_index_mifish %>% unique() %>% sort() %>% length()/D_mifish$sp_index_mifish %>% max()
  D_sebastes$sp_index_sebastes %>% unique() %>% sort() %>% length()/D_sebastes$sp_index_sebastes %>% max()
  D_count$station_idx %>% max() / D_count$station_idx %>% unique() %>% sort() %>% length() 
  D_mifish$station_idx %>% max() / D_mifish$station_idx %>% unique() %>% sort() %>% length() 
  D_sebastes$station_idx %>% max() / D_sebastes$station_idx %>% unique() %>% sort() %>% length()    
  
  # make prior for b_sebastes
  # names(rockfish) <- c("ID_sebastes", "abundance", "larvae")
  # sebastes_prior <- rockfish %>% right_join(species_mapping %>% filter(ID_sebastes != "MISSING")) %>% 
  #   mutate(abundance = ifelse(is.na(abundance), 1, abundance),
  #          prop = abundance / sum(abundance),
  #          prior = log(prop * 300)) %>% 
  #   dplyr::select(ID_sebastes, prior) %>% 
  #   left_join(species_mapping %>% dplyr::select(ID_sebastes, ID_master)) %>% 
  #   dplyr::select(-ID_sebastes)
  # 
  #we want to keep all stations/species for which there is at least one observation in one dataset (count, mifish, sebastes)
    Obs <- full_join(pos.all %>% dplyr::select(ID_master,station_id),
                   species_mapping %>% dplyr::select(ID_master,sp_index_master)) %>%
          left_join(.,STATION_ID) 
          #left_join(sebastes_prior) %>% 
          #mutate(prior = ifelse(is.na(prior), 0, prior))
  
  N_station_species_master <- nrow(Obs)
  
  
  
  N_pcr_mifish = 39
  N_pcr_sebastes = 49


###############################    
  
  
  stan_data <- list(
    "N_mifish_obs" = N_mifish_obs,
    "N_count_obs"  = N_count_obs,
    "N_sebastes_obs"= N_sebastes_obs,
    
    "D_count_obs"  = as.vector(D_count$Abund),
    "D_mifish_obs" = D_mifish$nReads,
    "D_sebastes_obs" = D_sebastes$nReads,
    
    # log of fraction of amplicons observed for each replicate-community pair
    "log_r_mifish" = D_mifish %>% dplyr::select(mifish_station_rep_idx, log_r_mifish) %>% 
                      distinct() %>% pull(log_r_mifish), 
    
    "log_r_sebastes" = D_sebastes %>% dplyr::select(sebastes_station_rep_idx, log_r_sebastes) %>% 
      distinct() %>% pull(log_r_sebastes), 
    
    ## Species count
    "I_master" = I_master ,
    "I_count"  = I_count,
    "I_mifish" = I_mifish,
    "I_sebastes" = I_sebastes,
    
    ##Technical PCR replicate count and index
    "J_mifish" = J_mifish,
    "J_mifish_vec" = J_all_mifish$N_mifish,
    "J_sebastes" = J_sebastes,
    "J_sebastes_vec" = J_all_sebastes$N_sebastes,
    
    ##Station / Site count
    "K_master" = K_master,
    "K_mifish" = K_mifish,
    "K_sebastes" = K_sebastes ,
    "K_count" = K_count,
    
    "N_pcr_mifish" = N_pcr_mifish, # Number of PCR cycles (assumed constant across all PCRs)
    "N_pcr_sebastes" = N_pcr_sebastes,
    
    # Master species list indexes and counters
    "N_station_species_master" = N_station_species_master,
    "master_station_idx" = Obs$station_idx,
    "master_sp_idx" = Obs$sp_index_master,
    
    # Indexes and covariate for mifish data
    "mifish_station_idx" = D_mifish$station_idx, # station index for mifish data.
    "mifish_sp_idx" = D_mifish$sp_index_mifish, # Species index for mifish data
    "mifish_station_rep_idx" = D_mifish$mifish_station_rep_idx,
    "N_mifish_station_rep_idx" = length(unique(D_mifish$mifish_station_rep_idx)),
    
    # Indexes and covariate for sebastes data
    "sebastes_station_idx" = D_sebastes$station_idx, # station index for sebastes data.
    "sebastes_sp_idx" = D_sebastes$sp_index_sebastes, # Species index for sebastes data
    "sebastes_station_rep_idx" = D_sebastes$sebastes_station_rep_idx,
    "N_sebastes_station_rep_idx" = length(unique(D_sebastes$sebastes_station_rep_idx)),
    
    # Indexes and covariate for count data
    "N_count_obs" = N_count_obs,
    "log_V" = D_count$log_V, # Volume of water sampled relative to the reference volume.
    "log_P" = D_count$log_P, # Fraction of larvae sorted.
    "count_station_idx" = D_count$station_idx, # station index for count data.
    "count_sp_idx" = D_count$sp_index_count, # Species index for count data
    
    # Identity matrices relating identifications for species among datasets
    "M_to_count" = M_to_count,
    "M_to_mifish" = M_to_mifish,
    "M_to_sebastes" = M_to_sebastes,
    
    #mean and SD of distribution from which eta is drawn. Eta is the expected fraction of amplicons that are sampled for sequencing. This is on a log scale, such that a mean of -4 represents exp(-4) == 0.018 
    "log_eta_prior" = c(-4,0.5),  #this parameter is shared for the mifish and sebastes datasets, reflecting the assumption that the same fraction of amplicons are sequenced in the two datasets
    "beta_prior_mifish" = c(1,1),  # priors for the beta distribution from which values of mifish amp efficiencies are drawn
    "a_sebastes_sd" = 0.05   # standard deviation for values of sebastes amp efficiencies, which will share a common beta distribution. Here, setting the SD and letting the model estimate the mean from the data. 
    )

  
  stan_pars = c(
    "b_master",
    "b_master_grid",
    "b_mifish",
    "b_sebastes",
    "b_count", 
    "a_mifish",
    "a_sebastes_all",
    "tau_mifish",
    "tau_sebastes",
    "log_eta_mifish",
    "log_eta_sebastes",
    "log_eta_mifish_mean",
    "log_eta_sebastes_mean"
    #"log_b_mifish_sardine",
    #"sim_obs"
    #"nu"
  )   
  
  stan_init_f2 <- function(n.chain,I_mifish,N_station_rep_idx){#J_seb,K_seb){ 
    A <- list()
    for(i in 1:n.chain){
      A[[i]] <- list(
        #log_eta_mifish = runif(1, -4, 0),
        #log_eta_sebastes = runif(1,-4, 0),
        a_mifish = runif(I_mifish,0.4,0.6),
        a_sebastes_all = .5,
        tau_mifish = runif(1,0.5,1.5),
        tau_sebastes = runif(1,0.5,1.5) 
      )
    }
    return(A)
  }

  
  # Define the MCMC, and run
  N_CHAIN = 3
  Warm = 100
  Iter = 200
  Treedepth = 12
  Adapt_delta = 0.95
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  write_rds(stan_data, here("Data/ole_data_calcofi/Calcofi_STAN_data.rds"))
  
  
  if(MAP.ERROR > 1 | MIN.ERROR > 0){print("YOU HAVE A MAPPING MATRIX ERROR, SPECIES ARE BEING DOUBLE COUNTED")
                  print("DON'T RUN THE STAN CODE UNTIL FIXED!")
                  print("YOUR ARMPITS SMELL OF ELDERBERRIES.")
                  print("YOUR PAL, OLE")
                  }else{
                    
  stanMod = stan(file = here("Data/ole_data_calcofi/CalCOFI_20210617_randomEta.stan") ,
                 data = stan_data, 
                 verbose = FALSE, chains = N_CHAIN, thin = 1, 
                 warmup = Warm, iter = Warm + Iter, 
                 control = list(max_treedepth=Treedepth,
                                stepsize=0.01,
                                adapt_delta=Adapt_delta,
                                metric="diag_e"),
                 pars = stan_pars,
                 refresh = 10,
                 boost_lib = NULL,
                 #sample_file = here("Data/ole_data_calcofi/twoTau_tmp.csv"),
                 init = stan_init_f2(n.chain=N_CHAIN,
                                     I_mifish = I_mifish,
                                     N_station_rep_idx=stan_data$N_station_rep_idx)
  )
  
  
  
  # get_adaptation_info(stanMod)
  pars <- rstan::extract(stanMod, permuted = TRUE)
  samp_params <- get_sampler_params(stanMod)
  stanMod_summary <- summary(stanMod)$summary
  
  TRACE <- list()
  TRACE[[as.name("Var")]] <- traceplot(stanMod,pars=c("lp__","tau_sebastes", "tau_mifish", "log_eta_mifish_mean", "log_eta_sebastes_mean"),inc_warmup=FALSE)
  TRACE[[as.name("a_mifish")]] <- traceplot(stanMod,pars=c("lp__","a_mifish"),inc_warmup=FALSE)
  TRACE[[as.name("b_master")]] <- traceplot(stanMod,pars=c("lp__","b_master[1]","b_master[20]","b_master[32]",
                                                           "b_master[64]","b_master[85]","b_master[123]","b_master[278]",
                                                           "b_master[448]","b_master[669]","b_master[810]","b_master[911]",
                                                           "b_master[1233]","b_master[1234]"))
  TRACE[[as.name("b_master_grid")]] <- traceplot(stanMod,pars=c("lp__","b_master_grid[1,1]","b_master_grid[2,1]","b_master_grid[28,10]",
                                                           "b_master_grid[32,14]","b_master_grid[41,10]","b_master_grid[55,40]"))
  
  stanMod_summary %>% as.data.frame %>% dplyr::select(Rhat) %>% rownames_to_column("name") %>% 
    filter(grepl("b_master_grid", name))
  
  stanMod_summary %>% as.data.frame %>% dplyr::select(Rhat) %>% rownames_to_column("name") %>% arrange(desc(Rhat)) %>% head(20)
  traceplot(stanMod, pars = c("b_mifish[38,29]","b_mifish[38,25]", "a_mifish[25]", "b_master[3277]"))
  

  
  # Parse model output.
  
  
  Output <- list(
    Nspecies = I_master,
    Ncommunities = K_master,
    Nreplicates = J_mifish,
    N_pcr_mifish = N_pcr_mifish,
    N_pcr_sebastes = N_pcr_sebastes,
    
    species_mapping = species_mapping,
    
    stan_data = stan_data,
    stan_pars = stan_pars,
    stanMod = stanMod,
    
    D_count = D_count,
   
    D_mifish =D_mifish,
    D_sebastes = D_sebastes,
    
    Stations = STATION_ID
    
  )
  if(file.exists("Stan_Fit_combined_last.RData")){
    file.rename(from = "Stan_Fit_combined_last.RData",
                to = paste0("Stan_Fit_combined_",format(Sys.time(), "%Y%m%d_%H%M"),".RData"))}
  save(Output,file="Stan_Fit_combined_last.RData")
  
  } # End Catch for Mapping errors
  

  
  
  
  