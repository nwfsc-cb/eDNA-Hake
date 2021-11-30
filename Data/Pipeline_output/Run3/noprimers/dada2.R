params <-
list(folder = "/Users/ramon.gallegosimon/Projects/eDNA-Hake/Data/Pipeline_output/Run3/noprimers/", 
    hash = TRUE, trimming.length.Read1 = 220L, trimming.length.Read2 = 150L, 
    metadata = "output.metadata.csv", output.folder = "~/test_dad2", 
    keep.mid.files = TRUE)

## ----setup, include=FALSE----------------------------------------------------------------------

knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = params$folder)


## ----loading packages, echo=TRUE ,message=FALSE------------------------------------------------

library (tidyverse)
library (dada2)
library (digest)
library (seqinr)
library (knitr)
library (kableExtra)

sample.metadata <- read_csv(params$metadata)
filt_path <- file.path(params$folder, "filtered")
getwd()
filt_path
# Check if output directory exists - if not create as a subfolder of input dir

if(!dir.exists(file.path(params$output.folder))){
  dir.create(path = file.path(params$output.folder),recursive = T)
  output.dir <- file.path(params$output.folder)
}else{
  output.dir <- file.path(params$output.folder)
}
output.dir



## ----dadaing-----------------------------------------------------------------------------------

sample.metadata %>%
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
                       } )) %>% 
  filter (outFs$reads.out > 0)  %>% # Keep only cases in which there are sequences passing filter
  mutate(
          errF1 = map(filtF1, ~ learnErrors(.x, multithread=TRUE,verbose = 0)),     # Calculate errors
          errR1 = map(filtR1, ~ learnErrors(.x, multithread=TRUE,verbose = 0)),
          derepF1 = map(filtF1, derepFastq),                   # dereplicate seqs
          derepR1 = map(filtR1, derepFastq),
          dadaF1  = map2(derepF1,errF1, ~ dada(.x, err = .y, multithread = TRUE)),  # dada2
          dadaR1  = map2(derepR1,errR1, ~ dada(.x, err = .y, multithread = TRUE)),
          mergers = pmap(.l = list(dadaF1,derepF1, dadaR1,derepR1),                 # merge things
                         .f = mergePairs,
                         minOverlap = 10)) -> output.dada2
if ( params$keep.mid.files==TRUE){
write_rds(output.dada2, path = "output.halfway.rds")}


seqtabF <- makeSequenceTable(output.dada2$mergers)
dim(seqtabF)

table(nchar(getSequences(seqtabF)))




