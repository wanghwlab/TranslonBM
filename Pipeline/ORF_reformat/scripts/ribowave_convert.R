.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 60 * 1024 * 1024^2)
options(scipen = 999)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = F, quote = "none")
}

gencode_transcript_seq <- read_tsv(
    "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.transcripts.seq",
    col_names = F,
    col_select = c(
        "transcript_id" = "X1",
        "transcript_sequence" = "X2"
    )
)

## contains no stop codon ,do not need to modify ORF_stop

#prepare_gppy <- possibly(function(raw_file, gppy) {
#    tmp_file <- read_tsv(raw_file, col_names = T)
#    tmp_file <- tmp_file |>
##        dplyr::inner_join(
#            gencode_transcript_seq,
#            by = c("transcript_id" = "transcript_id")
#        ) |>
#        dplyr::mutate(
#            ORF_sequence = stringr::str_sub(transcript_sequence, ORF_tstart, ORF_tstop),
#            ORF_tstart_1base = ORF_tstart,
#            ORF_tstop_1base = ORF_tstop,
#            ORF_sequence_correct = stringr::str_sub(ORF_sequence),
#            start_codon = stringr::str_sub(ORF_sequence_correct, 1, 3)
#        ) |>
#        dplyr::select(transcript_id, ORF_tstart_1base, ORF_tstop_1base, ORF_sequence_correct, start_codon) |>
#        dplyr::distinct()

#    output_func(tmp_file, gppy)
#})

# Remove possibly during debugging, or keep it while ensuring that the logic remains correct
prepare_gppy <- possibly(function(raw_file, gppy) {
    
    # 1. Read the file with automatic header detection disabled (col_names = F)
    tmp_file <- read_tsv(raw_file, col_names = FALSE, col_types = cols(.default = "c"))
    
    # 2. Key fix: split columns
    # Assumed format: TranscriptID_Frame_Start_Stop
    # Use separate to split X1 into the required columns
    tmp_file <- tmp_file |>
        tidyr::separate(
            col = X1, 
            into = c("transcript_id", "frame", "ORF_tstart", "ORF_tstop"), 
            sep = "_", 
            remove = FALSE, # Retain the original column as a safeguard
            convert = TRUE  # Convert numeric strings to numeric values automatically
        ) |>
        # Ensure Start and Stop are numeric to prevent downstream errors
        dplyr::mutate(
            ORF_tstart = as.numeric(ORF_tstart),
            ORF_tstop = as.numeric(ORF_tstop)
        )

    # 3. The subsequent logic is unchanged
    tmp_file <- tmp_file |>
        dplyr::inner_join(
            gencode_transcript_seq,
            by = c("transcript_id" = "transcript_id")
        ) |>
        dplyr::mutate(
            # Note: substr indices in R start at 1. If the mx file uses 0-based coordinates, this may require adjustment (+1)
            # Ribowave output may be 0-based or 1-based; confirm this in the tool documentation.
            # Assume the input is 1-based and use it directly:
            ORF_sequence = stringr::str_sub(transcript_sequence, ORF_tstart, ORF_tstop),
            ORF_tstart_1base = ORF_tstart,
            ORF_tstop_1base = ORF_tstop,
            ORF_sequence_correct = stringr::str_sub(ORF_sequence),
            start_codon = stringr::str_sub(ORF_sequence_correct, 1, 3)
        ) |>
        dplyr::select(transcript_id, ORF_tstart_1base, ORF_tstop_1base, ORF_sequence_correct, start_codon) |>
        dplyr::distinct()

    output_func(tmp_file, gppy)
}, otherwise = NULL)

prepare_gppy(
    raw_file = snakemake@input[[1]],
    gppy = snakemake@output[[1]]
)
