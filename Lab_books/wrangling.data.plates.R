all.plates <- read_csv("/Users/ramon.gallegosimon/Downloads/Hake_eDNA_2019 plates_layout.csv")

all.plates %>% 
  filter (str_detect(X5, "Hake")) %>% 
  select(X5) %>% 
  separate (X5, into = c(NA, "plate"), sep = "plate ") -> plates.left

all.plates %>% 
  filter (!is.na(X4)) %>%
  select(4:16) %>% 
  rename (row = 1) %>% 
  rename_with(~1:12, .cols= 2:13) %>% 
  mutate (group = rep(1:10, each = 8)) %>% 
  nest(-group) %>% 
  bind_cols(plates.left) %>% 
  select(data, plate) -> lefts
  
 
all.plates %>% 
  filter (str_detect(`inhibted - to purify`, "Hake")) %>% 
  select(`inhibted - to purify`) %>% 
  separate (`inhibted - to purify`, into = c(NA, "plate"), sep = "plate ") -> plates.right

all.plates %>% 
  filter (!is.na(X19)) %>%
  select(19:31) %>% 
  rename (row = 1) %>% 
  rename_with(~1:12, .cols= 2:13) %>% 
  mutate (group = rep(1:11, each = 8)) %>% 
  nest(-group) %>% 
  bind_cols(plates.right) %>% 
  select(data, plate) -> righta

bind_rows(lefts, righta) %>% 
  mutate (data = map(data , ~.x %>% 
                       mutate_all(as.character) %>% pivot_longer( -row, names_to="col", values_to="sample"))) %>% 
  unnest(data) %>% 
  filter (!is.na(sample)) %>% 
  filter (!str_detect(sample, "extc|empty")) %>% 
  write_csv(here("Data", "Samples_To_Sequence","Hake_eDNA_2019_plates_layout_long.csv"))
