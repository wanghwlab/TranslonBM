options(scipen = 999)
library(tidyverse)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

ribotish_extract_ORF <- possibly(function(raw_file,withseq) {

    tmp <- read_tsv(raw_file, col_names = T) |>
        dplyr::rename(
            "gene_id" = "Gid",
            "transcript_id" = "Tid",
            "start_codon" = "StartCodon",
            "ORF_sequence" = "Seq",
            "ORF_sequence_aa" = "AASeq",
            "ORF_tstart" = "Start",
            "ORF_tstop" = "Stop"
        ) |>
        mutate(ORF_ID_raw = "not_set", ORF_length = str_length(ORF_sequence)) |>
        separate(GenomePos, c("chrom", "gene_pos", "strand"), sep = ":") |>
        separate(gene_pos, c("ORF_gstart", "ORF_gstop"), sep = "-")

    output_func(tmp, withseq)
})

## stand_result:: ORF_ID_raw gene_id transcript_id chrom strand ORF_gstart ORF_gstop ORF_tstart ORF_tstop ORF_length start_codon ORF_sequence ORF_sequence_aa annotated_gstart annotated_gstop annotated_tstart annotated_tstop
## ribotricer lack ORF_tstart and ORF_tstop<U+FF0C> and orf_sequence do not contain stop codon, need fix

ribotish_extract_ORF(
    raw_file = snakemake@input[[1]],
    withseq = snakemake@output[[1]]
)
