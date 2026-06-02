library(tidyverse)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = F, quote = "none")
}

prepare_gppy <- possibly(function(withseq, gppy) {
    tmp_file <- read_tsv(withseq, col_names = T)
    tmp_file <- tmp_file |>
        dplyr::mutate(
            ORF_tstart_1base = as.numeric(ORF_tstart),
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
    withseq = snakemake@input[[1]],
    gppy = snakemake@output[[1]]
)

## gppy need 1-based coord
## ribocode result is 1-based coord
