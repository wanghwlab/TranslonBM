.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
options(scipen = 999)
library(tidyverse)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

ribotricer_index <- read_tsv(
    "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_1core/ribotricer/ribotricer_chrN/index/Human.fasta",
    col_names = T
)

ribotricer_ORF_coordinate <- read_tsv(
    "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_1core/ribotricer/ribotricer_chrN/index/Human_candidate_orfs.tsv",
    col_names = T
)
## need improvement

ribotricer_ORF_metainfo <- ribotricer_index |>
    inner_join(ribotricer_ORF_coordinate, by = c("ORF_ID" = "ORF_ID")) |>
    dplyr::rename(ORF_sequence = sequence) |>
    select(ORF_ID, ORF_sequence, coordinate)

ribotricer_extract_ORF <- possibly(function(raw_file, withseq) {
    tmp <- read_tsv(
        raw_file,
        col_names = T,
        col_select = -c("profile")
    ) |>
        filter(phase_score > 0.440) |>
        inner_join(ribotricer_ORF_metainfo, by = c("ORF_ID" = "ORF_ID")) |>
        mutate(ORF_id_tmp = ORF_ID, ORF_ID_raw = ORF_ID) |>
        separate(ORF_id_tmp, c("transcript_id", "ORF_gstart", "ORF_gstop", "ORF_length"), sep = "_")

    output_func(tmp, withseq)
})

## stand_result:: ORF_ID_raw gene_id transcript_id chrom strand ORF_gstart ORF_gstop ORF_tstart ORF_tstop ORF_length start_codon ORF_sequence ORF_sequence_aa annotated_gstart annotated_gstop annotated_tstart annotated_tstop
## ribotricer lack ORF_tstart and ORF_tstop<U+FF0C> and orf_sequence do not contain stop codon, need fix

ribotricer_extract_ORF(
    raw_file = snakemake@input[[1]],
    withseq = snakemake@output[[1]]
)
