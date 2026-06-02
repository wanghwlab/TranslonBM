options(scipen = 999)
library(tidyverse)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

gencode_transcript_seq <- read_tsv(
    "/home/chengennong/Ribo_benchmark/ref/Human/gencode.v43.transcripts.seq",
    col_names = F,
    col_select = c(
        "transcript_id" = "X1",
        "transcript_seq" = "X2"
    )
)

ribocode_extract_ORF <- possibly(function(raw_file,withseq) {
    tmp <- read_tsv(raw_file, col_names = T) |>
        inner_join(
            gencode_transcript_seq,
            by = c("transcript_id" = "transcript_id")
        ) |>
        mutate(ORF_sequence = stringr::str_sub(transcript_seq, ORF_tstart, ORF_tstop)) |>
        select(-transcript_seq) |>
        dplyr::rename(ORF_ID_raw = ORF_ID, ORF_sequence_aa = AAseq)

    output_func(tmp, withseq)
})

## stand_result:: ORF_ID_raw gene_id transcript_id chrom strand ORF_gstart ORF_gstop ORF_tstart ORF_tstop ORF_length start_codon ORF_sequence ORF_sequence_aa annotated_gstart annotated_gstop annotated_tstart annotated_tstop

ribocode_extract_ORF(
    raw_file = snakemake@input[[1]],
    withseq = snakemake@output[[1]]
)