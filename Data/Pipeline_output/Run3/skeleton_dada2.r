## If you don't need the report, and just want to execute the dad2.Rmd script 
## with a params file 
## Run this script with 
## Rscript skeleton_dada2.r params.txt
## For the params file, make a copy of the /params.example.txt

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library (tidyverse)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else{params <- read_csv(args[1])}

params 
params <- params %>% pivot_wider(names_from = Argument, values_from=value)
params %>% mutate(across(.cols= starts_with("trimming"), .fn = as.numeric)) %>% 
  muatte*