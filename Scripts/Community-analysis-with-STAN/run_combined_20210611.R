library(tidyverse)
library(rstan)
library(bayesplot)
options(mc.cores = parallel::detectCores())
rstan_options("auto_write" = TRUE)

  #library prep data, w PCR pipetting info
  lib_mifish <- read.csv("./data/mifish_library_prep.csv")
      lib_mifish$New_name <-  paste0(   #need to match format of other dataset
        substr(lib_mifish$New_name, 1, 4), "_",
        substr(lib_mifish$New_name, 8, nchar(as.character(lib_mifish$New_name)))
      )
  
  # library prep data, w PCR pipetting info
  lib_sebastes <- read.csv("./data/sebastes_library_prep.csv")
  lib_sebastes$New_name <-  paste0(   #need to match format of other dataset
      substr(lib_sebastes$New_name, 1, 4), "_",
      substr(lib_sebastes$New_name, 8, nchar(as.character(lib_sebastes$New_name)))
    )

  dat.mifish <- readRDS("./data/mifish_tech_nReads.RDS") %>% 
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
  
  dat.sebastes <- readRDS("./data/sebastes_tech_nReads.RDS") %>% 
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
    
  
  dat.count <- readRDS("./data/microscopy_tech_nReads.RDS") %>% 
    filter(!ID_microscopy %in% c("MISSING", "", "Disintegrated fish larvae")) %>% 
    filter(station_id %in% unique(dat.mifish$station_id)) %>% 
  #dat.count$larval_counts[dat.count$ID_microscopy == "Vinciguerria sp."] <- 0  #DATA problem; check w Zack
  # dat.count <- dat.count %>% 
  #   mutate(ID_microscopy = ifelse(grepl("Sebastes", ID_microscopy), "Sebastes", ID_microscopy),
  #          ID_mifish = ifelse(grepl("Sebastes", ID_mifish), "Sebastes", ID_mifish),
  #          Unique_ID = ifelse(grepl("Sebastes", Unique_ID), "Sebastes", Unique_ID)
  #          ) %>%   ##TEMPORARY -- synonymize all Sebastes
    group_by(ID_microscopy, station_id, standard.haul.factor, volume.filtered, proportion.sorted) %>% 
    dplyr::summarise(larval_counts = sum(larval_counts)) %>% 
    mutate(., ID_microscopy = recode(ID_microscopy, `Oxyjulis californica` = "Halichoeres californica")) #recode this
  
  #synonymize all Sebastes in counts
  dat.count$ID_microscopy[grepl("Sebastes", dat.count$ID_microscopy)] <- "Sebastes"
  dat.count <- dat.count %>% group_by(ID_microscopy, station_id, standard.haul.factor, volume.filtered, proportion.sorted) %>% dplyr::summarise(larval_counts = sum(larval_counts))
  
  #import info to set up prior on Sebastes counts
  #rockfish <- read.csv("./data/Thompson_2017_supp.csv", header = F)
  
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
    dat.count %>% ungroup() %>% summarize(sum(larval_counts))
  
  #make a new object dat.count.init that has been filtered.  Do not overwrite dat.count
  dat.count.init <- dat.count %>%
    filter(ID_microscopy %in% common_count)  #take just common species, for now

    # pull appropriate species from microscopy counts
  species_mapping <- read.csv("./data/20210610_species_mapping_file.csv") %>% 
    filter(ID_mifish %in% dat.mifish.init$ID_mifish  | 
             ID_microscopy %in% dat.count.init$ID_microscopy | 
             ID_sebastes %in% dat.sebastes$ID_sebastes) %>% 
    filter(!ID_sebastes %in% drop.these.sebastes) %>%
    arrange(Unique_ID)
  names(species_mapping) <- c("ID_master", "ID_mifish", "ID_count", "ID_sebastes")
  species_mapping <- species_mapping %>% arrange(ID_master,ID_count,ID_mifish)
  species_mapping[species_mapping == ""] <- "MISSING"

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
  pos.mifish <- dat.mifish %>% filter(nReads > 0) %>% dplyr::select(ID_mifish,station_id) %>%
    distinct() %>% left_join(.,species_mapping)
  pos.sebastes <- dat.sebastes %>% filter(nReads > 0) %>% dplyr::select(ID_sebastes,station_id) %>%
    distinct() %>% left_join(.,species_mapping)
  
  pos.all <- pos.count %>% dplyr::select(ID_master,station_id) %>%
    bind_rows(.,pos.mifish %>% dplyr::select(ID_master,station_id)) %>%
    bind_rows(.,pos.sebastes %>% dplyr::select(ID_master,station_id)) %>%
    distinct() %>% arrange(ID_master,station_id) %>%
    full_join(.,species_mapping %>% dplyr::select(ID_master,ID_count,ID_mifish, ID_sebastes)) 
  
  # This should weed out any samples with all zeros...  
  dat.count <- semi_join(dat.count,pos.all)
  dat.mifish <- semi_join(dat.mifish,pos.all)
  dat.sebastes <- semi_join(dat.sebastes,pos.all)
  
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
      
      
      colnames(M_to_count)  <- species_mapping$ID_master
      colnames(M_to_mifish) <- species_mapping$ID_master
      colnames(M_to_sebastes) <- species_mapping$ID_master
      rownames(M_to_count)   <- species_mapping %>% filter(!ID_count=="MISSING") %>% dplyr::select(ID_count,sp_index_count) %>%
                                  distinct() %>% pull(ID_count)
      rownames(M_to_mifish)  <- species_mapping %>% filter(!ID_mifish=="MISSING") %>% dplyr::select(ID_mifish,sp_index_mifish) %>%
                                  distinct() %>% pull(ID_mifish)
      rownames(M_to_sebastes)  <- species_mapping %>% filter(!ID_sebastes=="MISSING") %>% dplyr::select(ID_sebastes,sp_index_sebastes) %>%
                                    distinct() %>% pull(ID_sebastes)
      
      #write these to file, for checking accuracy
      write.csv(M_to_count, "M_to_count.csv")
      write.csv(M_to_mifish, "M_to_mifish.csv")
      write.csv(M_to_sebastes, "M_to_sebastes.csv")
      
      
    ## CHECK THESE MATRICES. CHARACTERISTICS SHOULD BE:
        # COLSUMS all are 0 or 1.  Most should be 1.
        # ROWSUMS 0, 1, >1 all are possible.
        MAP.ERROR <- max(colSums(M_to_count),colSums(M_to_mifish),colSums(M_to_sebastes))
        names(which(colSums(M_to_count)>1))  -> names
        species_mapping %>%  filter(., ID_master %in% names)
        #M_to_count[,names] %>%  View()
        
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
      summarise(J=length(station_id)) %>% 
      group_by(station_id) %>% 
      summarise(MIN=min(J),MAX=max(J)) %>% 
      filter(MIN == MAX) %>% 
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
  
  
  
  N_pcr = 39


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
    
    "N_pcr" = N_pcr, # Number of PCR cycles (assumed constant across all PCRs)
    
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
    #"mc_sp_idx" =D_count_obs$sp_idx,  
    "count_sp_idx" = D_count$sp_index_count, # Species index for count data
    
    # Identity matrices relating identifications for species among datasets
    "M_to_count" = M_to_count,
    "M_to_mifish" = M_to_mifish,
    "M_to_sebastes" = M_to_sebastes,
    
    #mean and SD of distribution from which eta is drawn. Eta is the expected fraction of amplicons that are sampled for sequencing. This is on a log scale, such that a mean of -4 represents exp(-4) == 0.018 
    "log_eta_prior" = c(-4,0.5),  #this parameter is shared for the mifish and sebastes datasets, reflecting the assumption that the same fraction of amplicons are sequenced in the two datasets
    
    
    "beta_prior_mifish" = c(1,1),  # priors for the beta distribution from which values of mifish amp efficiencies are drawn
    "a_sebastes_sd" = 0.05   # standard deviation for values of sebastes amp efficiencies, which will share a common beta distribution. Here, setting the SD and letting the model estimate the mean from the data. 
    #"b_master_prior_log_mean" = Obs$prior  # prior for inferred true value of counts, influenced by previous work on sebastes larvae in Thompson et al. 2017
    )
  # if(FIX.ANCHOVY==TRUE){  
  #     stan_data <- c(stan_data,
  #              list("a_mifish_fix" = a_mifish_fix$a_mifish_fix,
  #              "a_mifish_fix_idx" = a_mifish_fix$sp_index_mifish,
  #               "N_a_mifish_fix" = N_a_mifish_fix,
  #              "a_mifish_est_idx" = a_mifish_est$sp_index_mifish))
  # }
  
  stan_pars = c(
    "b_master",
    "b_master_grid",
    "b_mifish",
    "b_sebastes",
    "b_count", 
    "a_mifish",
    "a_sebastes_all",
    #"a_sebastes",
    #"phi",
    "tau_mifish",
    "tau_sebastes",
    #"tau_sebastes",
    #"tau_sebastes",
    #"a_mifish_sd",
    "log_eta"
    #"eta_mifish",
    #"pi_master_grid",
    #"pi_mifish_grid",
    #"pi_sebastes_grid"
    #"log_r_mifish_real"
  )   
  
  stan_init_f2 <- function(n.chain,I_mifish,N_station_rep_idx){#J_seb,K_seb){ 
    A <- list()
    for(i in 1:n.chain){
      A[[i]] <- list(
        log_eta_mifish = -4,
        log_eta_sebastes = -4,
        a_mifish = runif(I_mifish,0.4,0.6),
        a_sebastes_all = .5,
        #a_sebastes = runif(I_sebastes,0.6,0.65),
        #epsilon = rnorm(N_station_rep_idx, -1.609438, 0.1),
        tau = runif(1,0.01,0.3)
        #tau_sebastes = runif(1,0.01,0.3) 
      )
    }
    return(A)
  }

  
  # Define the MCMC, and run
  N_CHAIN = 3
  Warm = 200
  Iter = 300
  Treedepth = 12
  #Adapt_delta = 0.90
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  # if(FIX.ANCHOVY==TRUE){FILE = "CalCOFI_DNA_species_mifish_v0513.stan"}
  # if(FIX.ANCHOVY==FALSE){FILE = "CalCOFI_DNA_species_mifish_v0513_a_unknown.stan"}
  
  
  
  if(MAP.ERROR>1){print("YOU HAVE A MAPPING MATRIX ERROR, SPECIES ARE BEING DOUBLE COUNTED")
                  print("DON'T RUN THE STAN CODE UNTIL FIXED!")
                  print("YOUR ARMPITS SMELL OF ELDERBERRIES.")
                  print("YOUR PAL, OLE")
                  }else{
                    
  stanMod = stan(file = "CalCOFI_20210611_twoTau.stan" ,data = stan_data, 
                 verbose = FALSE, chains = N_CHAIN, thin = 1, 
                 warmup = Warm, iter = Warm + Iter, 
                 control = list(max_treedepth=Treedepth,
                                stepsize=0.01,
                                #adapt_delta=Adapt_delta,
                                metric="diag_e"),
                 pars = stan_pars,
                 refresh = 10,
                 boost_lib = NULL,
                 sample_file = paste0("./OutputFiles/test_Ole.csv"),
                 init = stan_init_f2(n.chain=N_CHAIN,
                                     I_mifish = I_mifish,
                                     N_station_rep_idx=stan_data$N_station_rep_idx)
  )
  
  
  
  # get_adaptation_info(stanMod)
  pars <- rstan::extract(stanMod, permuted = TRUE)
  samp_params <- get_sampler_params(stanMod)
  stanMod_summary <- summary(stanMod)$summary
  
  TRACE <- list()
  TRACE[[as.name("Var")]] <- traceplot(stanMod,pars=c("lp__","tau_mifish","tau_sebastes","log_eta"),inc_warmup=FALSE)
  TRACE[[as.name("a_mifish")]] <- traceplot(stanMod,pars=c("lp__","a_mifish"),inc_warmup=FALSE)
  #TRACE[[as.name("a_sebastes")]] <- traceplot(stanMod,pars=c("lp__","a_sebastes"),inc_warmup=FALSE)
  TRACE[[as.name("b_master")]] <- traceplot(stanMod,pars=c("lp__","b_master[1]","b_master[20]","b_master[32]",
                                                           "b_master[64]","b_master[85]","b_master[123]","b_master[278]",
                                                           "b_master[448]","b_master[669]","b_master[810]","b_master[911]",
                                                           "b_master[1233]","b_master[1234]"))
  TRACE[[as.name("b_master_grid")]] <- traceplot(stanMod,pars=c("lp__","b_master_grid[1,1]","b_master_grid[2,1]","b_master_grid[28,10]",
                                                           "b_master_grid[32,14]","b_master_grid[41,10]","b_master_grid[55,40]"))
  
  stanMod_summary %>% as.data.frame %>% dplyr::select(Rhat) %>% rownames_to_column("name") %>% 
    filter(grepl("b_master_grid", name))
  
  #stanMod_summary %>% as.data.frame %>% dplyr::select(Rhat) %>% rownames_to_column("name") %>% 
   # filter(grepl("tau_mifsh", name))
  
  stanMod_summary %>% as.data.frame %>% dplyr::select(Rhat) %>% rownames_to_column("name") %>% arrange(desc(Rhat)) %>% head(20)
  traceplot(stanMod, pars = c("b_mifish[8,46]","b_mifish[8,42]", "a_mifish[18]", "b_master[2386]"))
  
  ## NEXT:
  ### COLLECT ALL OF THE MODEL OUTPUT, write to file
  ### MAKE PREDICTED-OBSERVED PLOTS TO MAKE SURE WE ARE NOT IN CRAZY TOWN.
  ### LOOK AT SUMMARY OF PROPORTIONAL and COUNT estimates for master.
  ### LOOK AT SUMMARY OF PROPORTIONAL estimates for mifish
  ### LOOK AT VARIABILITY IN ESTIMATES OF a_mifish
  

  
  # Parse model output.
  
  
  Output <- list(
    Nspecies = I_master,
    Ncommunities = K_master,
    Nreplicates = J_mifish,
    N_pcr = N_pcr,
    
    species_mapping = species_mapping,
    
    stan_data = stan_data,
    stan_pars = stan_pars,
    stanMod = stanMod,
    
    D_count = D_count,
   
    D_mifish =D_mifish,
    D_sebastes = D_sebastes,
    
    Stations = STATION_ID
    
    #true_a = true_a,
    #true_b = true_b,
    #logit_eta_true = logit_eta_true,
    #eta_true = eta_true
  )
  
  save(Output,file=paste0("Stan_Fit_combined_",format(Sys.time(), "%Y%m%d_%H%M"),".RData"))
  
  } # End Catch for Mapping errors
  
  # 
  # 
  # 
  # 
  # #how to observed amplicons look, relative to posterior?
  # #mifish
  # 
  # species_mapping$sp_index_master[species_mapping$ID_mifish=="Engraulis mordax"] #get index
  # y <- D_mifish %>% filter(ID_mifish == "Engraulis mordax" & station_idx == 1) %>% pull(nReads) /
  #   D_mifish %>% filter(station_idx == 1) %>% summarize(sum(nReads)) %>% as.numeric()
  # yrep <- unlist(extract(stanMod, pars = "pi_mifish_grid[1,16]"))
  # sapply(1:length(y), function(x){sum(y[x] > yrep)/length(yrep)})  #observed right in the middle of posterior
  # 
  # species_mapping$sp_index_master[species_mapping$ID_mifish=="Sardinops sagax"] #get index
  # y <- D_mifish %>% filter(ID_mifish == "Sardinops sagax" & station_idx == 1) %>% pull(nReads) /
  #   D_mifish %>% filter(station_idx == 1) %>% summarize(sum(nReads)) %>% as.numeric()
  # yrep <- unlist(extract(stanMod, pars = "pi_mifish_grid[1,37]"))
  # sapply(1:length(y), function(x){sum(y[x] > yrep)/length(yrep)})  #observed right in the middle of posterior
  # 
  # species_mapping$sp_index_master[species_mapping$ID_mifish=="Merluccius productus"] #get index
  # y <- D_mifish %>% filter(ID_mifish == "Merluccius productus" & station_idx == 1) %>% pull(nReads) /
  #   D_mifish %>% filter(station_idx == 1) %>% summarize(sum(nReads)) %>% as.numeric()
  # yrep <- unlist(extract(stanMod, pars = "pi_mifish_grid[1,25]"))
  # sapply(1:length(y), function(x){sum(y[x] > yrep)/length(yrep)})  #observed right in the middle of posterior
  # 
  # 
  # #how to observed amplicons look, relative to posterior?
  # #sebastes
  # 
  # species_mapping$sp_index_sebastes[species_mapping$ID_sebastes=="Sebastes jordani"] #get index
  # y <- D_sebastes %>% filter(ID_sebastes == "Sebastes jordani" & station_idx == 2) %>% pull(nReads) /
  #   D_sebastes %>% filter(station_idx == 2) %>% summarize(sum(nReads)) %>% as.numeric()
  # yrep <- unlist(extract(stanMod, pars = "pi_sebastes_grid[2,11]"))
  # sapply(1:length(y), function(x){sum(y[x] > yrep)/length(yrep)})  #observed relative to posterior
  # 
  # species_mapping$sp_index_sebastes[species_mapping$ID_sebastes=="Sebastes levis"] #get index
  # y <- D_sebastes %>% filter(ID_sebastes == "Sebastes levis" & station_idx == 2) %>% pull(nReads) /
  #   D_sebastes %>% filter(station_idx == 2) %>% summarize(sum(nReads)) %>% as.numeric()
  # yrep <- unlist(extract(stanMod, pars = "pi_sebastes_grid[2,12]"))
  # sapply(1:length(y), function(x){sum(y[x] > yrep)/length(yrep)})  #observed relative to posterior
  # 
  # 
  
  
  
  # look at implied biomass/counts
  
  # b <- extract(stanMod, pars = "b_master_grid")$b_master_grid
  # str(b)
  # 
  # d <- data.frame(STATION_ID,
  #                 Anchovy_bMeans = colMeans(b[,,which(Output$species_mapping$ID_master == "Engraulis mordax")]),
  #                 Sardine_bMeans = colMeans(b[,,which(Output$species_mapping$ID_master == "Sardinops sagax")]),
  #                 Bocaccio_bMeans = colMeans(b[,,which(Output$species_mapping$ID_master == "Sebastes paucispinis")])
  #                 ) %>%
  #   mutate(year = substr(station_id, 1,4),
  #          year = as.numeric(year),
  #          station = str_replace(station_id, "[0-9]+_", ""))
  # 
  # d %>%
  #   pivot_longer(cols = c(Sardine_bMeans, Anchovy_bMeans, Bocaccio_bMeans)) %>%
  #   ggplot(aes(x = year, y = value, color = name)) +
  #   geom_point() +
  #   geom_line() +
  #   facet_grid(station~.) +
  #   ylab("Log Posterior Count Equivalent")
  # 
  # d %>% dplyr::select(station_id, Sardine_bMeans) %>% 
  #   left_join(D_count %>% filter(ID_count == "Sardinops sagax") ) %>% 
  #   ggplot(aes(x = log(Abund), y = Sardine_bMeans)) + 
  #   geom_point() 
    
  
   
   
# 
#   rockfish <- which(grepl("Sebastes", species_mapping$ID_master))
#   
#   rockfish_means <- sapply(1:length(rockfish), FUN = function(x){colMeans(b[,,x])})  
#     colnames(rockfish_means) <- species_mapping$ID_master[rockfish]
#     rownames(rockfish_means) <- STATION_ID$station_id
#   
#     rockfish_means %>% 
#       as.data.frame() %>% 
#       rownames_to_column("station_id") %>% 
#       mutate(year = substr(station_id, 1,4),
#              year = as.numeric(year),
#              station = str_replace(station_id, "[0-9]+_", "")) %>% 
#       dplyr::select(-station_id) %>% 
#       pivot_longer(cols = c(-year, -station), names_to = "species") %>% 
#       ggplot(aes(x = year, y = value, color = species)) +
#         geom_point() +
#         geom_line() +
#         facet_grid(station~.)
#     
    
    
    
    
    
  # ## plotting
  # 
  # bayesplot::ppc_intervals(
  #   y = true_a[order(true_a)],
  #   yrep = extract(stanMod, par = "a_mifish")$a_mifish[,order(true_a)],
  #   x = true_a[order(true_a)],
  #   prob = 0.75
  # ) +
  #   ggtitle("Amplification Efficiency")
  # 
  # 
  # extract(stanMod, par = "pi_master_grid")$pi_master_grid
  # 
  # pi_master_grid_est <- apply(pars$pi_master_grid,c(2,3),mean) 
  # pi_master_grid_sd <- apply(pars$pi_master_grid,c(2,3),sd) 
  # pi_master_grid_q05 <- apply(pars$pi_master_grid,c(2,3),quantile,probs=0.05) 
  # pi_master_grid_q95 <- apply(pars$pi_master_grid,c(2,3),quantile,probs=0.95)
  # 
  # pi_master_out <- NULL
  # for(i in 1:ncol(pi_master_grid_est)){
  #   temp <- data.frame(Mean = pi_master_grid_est[,i],
  #                      SD = pi_master_grid_sd[,i],
  #                      q05 = pi_master_grid_q05[,i],
  #                      q95 = pi_master_grid_q95[,i],
  #                      species_idx = i ,
  #                      community_idx = 1:nrow(pi_master_grid_est))
  #   pi_master_out <- bind_rows(pi_master_out,temp)
  # }
  # 
  # b_true_est <- left_join(true_b,pi_master_out)
  # 
  # predicted.v.true.B <- 
  #   ggplot(b_true_est) +
  #   geom_point(aes(x=b_proportion,y=Mean),alpha=0.5,size=2)+
  #   geom_errorbar(aes(x=b_proportion,ymin=q05,ymax=q95),alpha=0.5) +
  #   geom_abline(intercept=0,slope=1,alpha=0.5,color="red") +
  #   theme_bw() +
  #   ggtitle("Species Composition")
  # 
  # predicted.v.true.B
  # predicted.v.true.B + facet_wrap(~species_idx)
  # predicted.v.true.B + facet_wrap(~community_idx)
  # 
  # 
  # #
  # bayesplot::ppc_intervals(
  #   y = eta_true,
  #   yrep = extract(stanMod, par = "eta_mifish") %>% unlist() %>% as.matrix(),
  #   x = eta_true,
  #   prob = 0.75
  # ) +
  #   ggtitle("Amplification Efficiency")
  
  
  
  
  
  