library(tidyverse)
library(Biostrings)

sampleID <- commandArgs(trailingOnly = T)[1] %>% str_remove(".FinalBinSet.tar.gz") %>% str_remove(".*/")
print(sampleID)

untar(commandArgs(trailingOnly = T)[1])

list.files("Final_binset/", full.names = T, pattern = ".fa") %>% 
  map(~{
    dna <- readDNAStringSet(.x)
    bindID <- .x %>% str_extract("bin[0-9]+")
    newnames <- names(dna) %>% str_remove(".*.fa") %>% paste0(sampleID, "-", bindID ,.)
    names(dna) <- newnames  
    if(!dir.exists(sampleID)){ dir.create(sampleID) }
    writeXStringSet(dna, paste0(sampleID, "/", sampleID, "-" ,bindID, ".fa"))
  })

read_tsv("Final_binset/Best_binset_bin_stats_ext.tsv", col_names = c("genome", "info")) %>% 
  mutate(info = str_remove_all(info, "\\{|\\}")) %>% 
  separate(info, c("contamination", "lineage", "completeness", "glength"), sep = ", ") %>% 
  mutate(across(contamination:glength, ~{str_remove(.x, "\'.*: ")})) %>% 
  select(-lineage, -glength) %>% 
  mutate(genome = paste0(sampleID, "-", str_extract(genome, "bin[0-9]+"), ".fa")) %>% 
  write_csv(paste0(sampleID, "/bins_checkm.csv"))

tar(paste0(sampleID, ".ren.tar.gz"), files = sampleID, compression = "gzip")

unlink("Final_binset/", recursive = T)
unlink(sampleID, recursive = T)

