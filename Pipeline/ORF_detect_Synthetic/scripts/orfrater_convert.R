library(tidyverse)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 60 * 1024 * 1024^2)
options(scipen = 999)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = F, quote = "none")
}

gencode_transcript_seq <- read_tsv(
    "/home/chengennong/Ribo_benchmark/ref/Human/gencode.v43.transcripts.seq",
    col_names = F,
    col_select = c(
        "transcript_id" = "X1",
        "transcript_sequence" = "X2"
    )
)

prepare_gppy <- possibly(function(raw_file,gppy) {

    tmp_file <- read_csv(raw_file, col_names = T, col_types = cols(.default = "c"))
    tmp_file <- tmp_file |>
        dplyr::filter(as.numeric(orfrating) > 0.8, as.numeric(AAlen) > 8) |>
        dplyr::inner_join(gencode_transcript_seq, by = c("tid" = "transcript_id")) |>
        dplyr::rename(
            transcript_id = tid,
            start_codon = codon,
            ORF_tstart = tcoord,
            ORF_tstop = tstop,
        ) |>
        dplyr::mutate(
            ORF_tstart_1base = as.numeric(ORF_tstart) + 1,
            ORF_tstop = as.numeric(ORF_tstop),
            ORF_sequence = stringr::str_sub(transcript_sequence, ORF_tstart_1base, ORF_tstop),
            ORF_tstop_1base = if_else(
                (stringr::str_length(ORF_sequence) %% 3 == 0),
                ORF_tstop - 3,
                ORF_tstop
            ),
            ORF_sequence_correct = if_else(
                (stringr::str_length(ORF_sequence) %% 3 == 0),
                stringr::str_sub(ORF_sequence, end = -4L),
                ORF_sequence
            )
        ) |>
        dplyr::select(transcript_id, ORF_tstart_1base, ORF_tstop_1base, ORF_sequence_correct, start_codon) |>
        dplyr::distinct()

    ## if orf is completed ,need remove stop codon .
    output_func(tmp_file, gppy)
})

prepare_gppy(
    raw_file = snakemake@input[[1]],
    gppy = snakemake@output[[1]]
)
