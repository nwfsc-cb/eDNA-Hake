# CTD plots Station 67-1

all.CTDS.to.keep <- read_csv(here(
  "Data","cleaned_CTD_stations.csv"))

Good_THermocline <- all.CTDS.to.keep %>% filter (Station %in% c("67-1", "66-1"))

files.inCTD <- tibble (file.name = list.files(path = here("Data/1_converted/"), pattern = "*.csv"),
                       path = list.files(path = here( "Data/1_converted"), pattern = "*.csv", full.names = T))


Good_THermocline %>% 
  mutate(file.name = paste0(CTD_Cast, ".csv")) %>% 
  left_join(files.inCTD) %>% 
  group_by(CTD_Cast) %>% 
  nest() %>%
  mutate(CTD_Data = map (data,  ~ read_csv(file =.x$path ))) -> all.CTDS.to.keep

all.CTDS.to.keep %>%  # Takes forever
  mutate (CTD_Data = map (CTD_Data, ~.x %>% mutate_at (c("Temperature (degC)", "Salinity (psu)", "Oxygen (ml_per_l)", "Fluorescence (ug_per_l)"), .funs = despike))) -> all.CTDS.to.keep
  
# Keep only ascent

all.CTDS.to.keep %>%  # Takes forever
  mutate(CTD_Data = map (CTD_Data, function(.x){
    maxdepth <- max(.x$`Depth (m)`)
    .x %>% rownames_to_column("step") %>% 
      filter (`Depth (m)`== maxdepth) %>%
    slice(1) %>%
    pull(step) %>% map(as.numeric) -> todelete
  .x %>% slice(todelete[[1]]:nrow(.x)) 
  })) -> all.CTDS.to.keep

# plot
all.CTDS.to.keep %>% 
  mutate(plot = map2(data, CTD_Data, function(.x, .y){
    .y %>% 
      
      select(Depth = `Depth (m)`, starts_with("Temp")) %>%
      rename(Sensor1 = 2, Sensor2 = 3) %>%
      pivot_longer(-Depth, names_to = "Sensor", values_to = "Temperature") %>%
      mutate(Station = .x$Station) %>% 
      filter(Depth>5) %>% 
      ggplot(aes(x = Temperature, y = Depth)) +
      geom_line(aes(group = Sensor, color = Sensor)) +
      scale_y_reverse() +
      scale_color_brewer(type = "div",palette = "Dark2") +
      facet_wrap(~Station) +
      guides(color = "none")+
      theme_linedraw()
  })) %>% pull(plot) -> plots

plots[[1]] + plots[[2]]
ggsave(here("Presentations", "images", "CTD.plot.png"), height = 3.5, width = 7)
