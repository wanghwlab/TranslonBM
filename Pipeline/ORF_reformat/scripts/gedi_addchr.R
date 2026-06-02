.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
options(scipen = 999)

prepare_gppy <- possibly(function(raw_bed,addchr_bed) {
    tmp_file <- read_tsv(raw_bed, col_names = F, col_types = cols(.default = "c")) |>
        mutate(X1 = str_c("chr", X1, sep = ""))

    readr::write_tsv(tmp_file, addchr_bed, col_names = F, quote = "none")
})

prepare_gppy(snakemake@input[[1]],snakemake@output[[1]])
