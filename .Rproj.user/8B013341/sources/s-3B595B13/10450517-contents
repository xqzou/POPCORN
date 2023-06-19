# Process genomic region coordinate files
library(tidyverse)

usethis::use_pipe()
telomere <- read.table("./telomere_hg38.tsv", sep = "\t", header = T, as.is = T)
usethis::use_data(telomere,overwrite = TRUE)

centromere <- read.table("./centromere_hg38_combined.txt", sep = "\t", header = T, as.is = T)
usethis::use_data(centromere,overwrite = TRUE)

gap <- read.table("./gap.txt", sep = "\t", header = T, as.is = T)
gap <- gap %>%
  mutate('ID'=paste0(Chrom, ":", Start, "-", End)) %>%
  dplyr::rename(CHROM=Chrom, POS=Start, END=End)
usethis::use_data(gap,overwrite = TRUE)

sg <- read.table("./20211015_SegDup.sorted.bed", sep = "\t", header = F, as.is = T)
names(sg) <- c("CHROM", "POS","END","ID")
sg <- sg %>%
  mutate('CHROM'=paste0("chr",CHROM)) %>%
  mutate('ID'=paste0(CHROM, ":", POS, "-", END)) %>%
  unique()
usethis::use_data(sg,overwrite = TRUE)

srep <- read.table("./simpleRepeat_hg38.txt", sep = "\t", header = F, as.is = T)
srep <- srep[,2:4]
names(srep) <- c("CHROM", "POS","END")
srep <- srep %>%
  mutate('ID'=paste0(CHROM, ":", POS, "-", END)) %>%
  unique()
usethis::use_data(srep,overwrite = TRUE)

