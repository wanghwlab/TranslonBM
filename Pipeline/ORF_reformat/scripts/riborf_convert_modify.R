.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
options(scipen = 999)

# --- Load required packages ---
library(tidyverse)
# Added: load Biostrings to read FASTA files
if (!require("Biostrings", quietly = TRUE)) {
    stop("Ensure that Biostrings is installed in the environment: conda install bioconductor-biostrings")
}
library(Biostrings)

# --- Modification begins: read the FASTA file correctly ---
# The original read_tsv call fails for .fa files; use readDNAStringSet instead
fa_path <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_1core/riborf/riborf_chrN/index/Human/candidateORF.fa"

# Read FASTA
fa_data <- readDNAStringSet(fa_path)

# Convert the FASTA object to a Tibble (DataFrame)
# names(fa_data)  is the ID, and as.character(fa_data)  is the sequence
riborf_ORF_sequence <- tibble(
    ORF_ID = names(fa_data) |> stringr::str_split_i(" ", 1), # Extract the ID before the first space
    ORF_sequence = as.character(fa_data)
)
# --- Modification ends ---


# --- The remaining code is unchanged (the index-reading step is adjusted slightly for safety) ---

riborf_index <- read_tsv(
    "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_1core/riborf/riborf_chrN/index/Human/candidateORF.genepred.txt",
    col_names = F,
    col_select = c(
        "ORF_ID" = "X1",
        "ORF_gstart" = "X6",
        "ORF_gstop" = "X7"
    )
)

riborf_metainfo <- riborf_ORF_sequence |>
    inner_join(riborf_index, by = c("ORF_ID" = "ORF_ID"))

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

tid_gid_looktable <- read_tsv(
    "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.looktable",
    col_names = F
) |>
    dplyr::rename(
        "transcript_id" = "X1",
        "gene_id" = "X2"
    )

# Remove possibly while debugging, or keep it as written; it is retained here to match the existing style
riborf_extract_ORF <- possibly(function(raw_file, withseq) {
    tmp <- read_tsv(raw_file, col_names = T) |>
        filter(pred.pvalue > 0.7 & length >= 27) |>
        mutate(ORF_ID_tmp = orfID, ORF_ID_raw = orfID) |>
        separate(ORF_ID_tmp, c("tid_tmp", "RankNumber", "transcript_tmp", "ORF_type", "start_codon"), sep = "\\|") |>
        separate(tid_tmp, c("transcript_id", "chrom_tmp", "strand_tmp"), sep = ":") |>
        separate(transcript_tmp, c("transcript_length", "ORF_tstart", "ORF_tstop")) |>
        inner_join(riborf_metainfo, by = c("orfID" = "ORF_ID")) |>
        dplyr::rename(
            "ORF_ID" = "orfID",
            "ORF_length" = "length"
        ) |>
        inner_join(tid_gid_looktable, by = c("transcript_id" = "transcript_id"))

    output_func(tmp, withseq)
}, otherwise = NULL)

## stand_result:: ORF_ID_raw gene_id transcript_id chrom strand ORF_gstart ORF_gstop ORF_tstart ORF_tstop ORF_length start_codon ORF_sequence ORF_sequence_aa annotated_gstart annotated_gstop annotated_tstart annotated_tstop

riborf_extract_ORF(
    raw_file = snakemake@input[[1]],
    withseq = snakemake@output[[1]]
)
