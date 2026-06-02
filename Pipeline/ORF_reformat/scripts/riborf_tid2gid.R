.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
options(scipen = 999)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 80 * 1024 * 1024^2)
library(tidyverse)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

riborf_extract_ORF <- possibly(function(withseq, tid_withseq, tid_noseq) {
    tmp <- read_tsv(withseq, col_names = T) |>
        select(transcript_id, ORF_tstart, ORF_tstop, ORF_sequence) |>
        mutate(
            ORF_tstart = ORF_tstart,
            ORF_tstop = ORF_tstop - 3 - 1, ## riborf gave wrong id, fix stop pos
            ORF_sequence = stringr::str_sub(ORF_sequence, end = -4L)
        ) |>
        distinct()

    tmp2 <- tmp |> select(-ORF_sequence)

    output_func(tmp, tid_withseq)
    output_func(tmp2, tid_noseq)
})


riborf_extract_ORF(
    withseq = snakemake@input[[1]],
    tid_withseq = snakemake@output[[1]],
    tid_noseq = snakemake@output[[2]]
)
