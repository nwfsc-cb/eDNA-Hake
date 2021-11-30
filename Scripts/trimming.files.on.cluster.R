library (tidyverse)
library (dada2)
library (digest)
library (seqinr)
library (knitr)
library (kableExtra)
library (patchwork)

sample.metadata <- read_csv("/home/rgallegosimon/data/noprimers/output.metadata.csv")
filt_path <- "/home/rgallegosimon/data/noprimers/filtered"
getwd()
filt_path
output.folder <- "/home/rgallegosimon/data/noprimers/output"
# Check if output directory exists - if not create as a subfolder of input dir

# trimming.length.Read1 <- 255
# trimming.length.Read2 <- 65

all.possibilities<-  list ("265_60" = c(265, 60),
      "250_75" = c(250, 75),
      "225_100"= c(225, 100),
      "200_100" = c(200, 120))

trimmer.test <- function(list.of.trims){

sample.metadata %>%
  separate(file1, into = "basename", sep= "_L001_R1", remove = F) %>% 
  mutate(filtF1 = file.path("filtered", paste0(basename, "_F1_filt.fastq.gz")),
         filtR1 = file.path("filtered", paste0(basename, "_R1_filt.fastq.gz"))) %>%
  select(-basename) %>% 
  mutate (outFs = pmap(.l= list (file1, filtF1, file2, filtR1),
                       .f = function(a, b, c, d) {
                         filterAndTrim(a,b,c,d,
                                       # truncLen=c(params$trimming.length.Read1,params$trimming.length.Read2),
                                       truncLen=c(list.of.trims[1],list.of.trims[2]),
                                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                       compress=TRUE, multithread=TRUE )
                       } )) -> output.dada2


  sample.metadata %>%
    slice(1:3) %>% 
    separate(file1, into = "basename", sep= "_L001_R1", remove = F) %>% 
    mutate(filtF1 = file.path("filtered", paste0(basename, "_F1_filt.fastq.gz")),
           filtR1 = file.path("filtered", paste0(basename, "_R1_filt.fastq.gz"))) %>%
    select(-basename) %>% 
    mutate (outFs = pmap_dfr(.l= list (file1, filtF1, file2, filtR1),
                         .f = function(a, b, c, d) {
                           filterAndTrim(a,b,c,d,
                                         # truncLen=c(params$trimming.length.Read1,params$trimming.length.Read2),
                                         truncLen=c(263,60),
                                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                         compress=TRUE, multithread=TRUE ) %>% 
                             as.data.frame()
                         } ))-> temp 
    
  temp %>% 
 filter (outFs$reads.out > 0)
 
  temp %>% 
    mutate(outFs = map_dfr(outFs, as.data.frame))
  
output.dada2 %>% 
   pull(outFs) %>% map(~ as.data.frame(.x) %>% 
                                    rownames_to_column) %>%  bind_rows() %>%  
  mutate(pass = reads.out/reads.in,
         sample = fct_reorder(sample, pass),
         fail = 1-pass) %>% 
  pivot_longer(cols = c(fail,pass), names_to = "process", values_to = "nreads") -> all.reads



ggplot(all.reads,aes(x = sample, y = nreads)) +
  geom_col(aes(fill = process), position="fill") +
  scale_y_continuous("proportion of reads passing filter") -> a



ggplot(all.reads %>% filter (process == "pass"),aes(x = reads.in, y = nreads)) +
  geom_point(aes(color = process)) +
  scale_x_log10() +
  scale_y_continuous("proportion of reads passing filter") -> b

return(list (output.dada2, all.reads,a, b ))
}

map (all.possibilities, trimmer.test) -> all.trimming.options

write_rds(all.trimming.options, "/home/rgallegosimon/data/noprimers/output/all.outputs.rds")


output.dada2 %>%
  pull(outFs) %>% map(~ as.data.frame(.x) %>%
                        rownames_to_column("sample")) %>%  bind_rows() %>%
  mutate(pass = reads.out/reads.in,
         sample = fct_reorder(sample, pass),
         fail = 1-pass) %>%
  pivot_longer(cols = c(fail,pass), names_to = "process", values_to = "nreads") -> all.reads



  ggplot(all.reads,aes(x = sample, y = nreads)) +
    geom_col(aes(fill = process), position="fill") +
    scale_y_continuous("proportion of reads passing filter") -> a

# 
# 
#   ggplot(all.reads %>% filter (process == "pass"),aes(x = reads.in, y = nreads)) +
#     geom_point(aes(color = process)) +
#     scale_x_log10() +
#     scale_y_continuous("proportion of reads passing filter") -> b
#   