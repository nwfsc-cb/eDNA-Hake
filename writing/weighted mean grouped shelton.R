# calculate "weighted mean" (WM) depth value for each 5km grid cell (Fivekm_g_full),
# using value (raster) attribute table (VAT) from composite NGDC depth grid (Comp_bath_0)
# depth values are in decimeters (dm) so need to divide WM depth by 10 to convert to meters
# general equation for calculating WM for each Fivekm_g_full grid cell is sum-product(depths vector and grid cell counts vector)/total grid cell count/10
# depths = Comp_bath_0 (in dm)
# grid cell counts = Count
# 5km grid cell ID# = Fivekm_g_full

# required packages
library(tidyverse)
library(magrittr)

# load .csv text file of VAT. Column types are integer, integer, integer, integer
grid_5km_raw <- read_csv("~/Documents/GitHub/Weighted-Means/src/5km_grid_combined_with_ngdc_dm_FULL_vat.csv",col_types='iiii')

# calculate weighted mean for each grid cell and divide by 10 to convert depth from decimeters to meters
weighted_mean <- grid_5km_raw %>%
  group_by(Fivekm_g_full) %>% 
  mutate(bath_output = weighted.mean(Comp_bath_0, Count)/10)

# remove duplicate rows in weighted_mean
dupe_remove <- weighted_mean[!duplicated(weighted_mean[ , c("Fivekm_g_full","bath_output")]),]

# create df with renamed 5km grid cell ID and weighted mean columns
dupe.sub <- dupe_remove %>% select(Fivekm_g_full, bath_output) %>%
  set_colnames(c("Gridcell_ID", "wm_ngdc_depth_m"))

# save to .csv file
write.csv(dupe.sub, "~/Documents/Projects/Ecosystem Science/Whale Entanglement/Regions/weighted_mean_NGDC_depths_for_5km_gridcells.csv",row.names = FALSE)