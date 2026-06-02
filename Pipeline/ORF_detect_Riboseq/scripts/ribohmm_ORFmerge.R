library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 80 * 1024 * 1024^2)
options(scipen = 999)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

get_coordinate <- function(tid_block_length, tid_block_start, ORF_gstart_0base, ORF_gstop_0base) {
    ORF_gstart_0base <- as.numeric(ORF_gstart_0base)
    ORF_gstop_0base <- as.numeric(ORF_gstop_0base)

    tid_block_length <- tid_block_length |>
        str_split_1(",") |>
        map_int(\(x) as.numeric(x))

    tid_exon_start <- tid_block_start |>
        str_split_1(",") |>
        map_int(\(x) as.numeric(x))

    tid_exon_stop <- tid_exon_start + tid_block_length
    ## need remove block which do not belong to orf

    tiv_start <- purrr::map2_lgl(tid_exon_start, tid_exon_stop, \(x, y) x <= ORF_gstart_0base & ORF_gstart_0base <= y)
    tiv_stop <- purrr::map2_lgl(tid_exon_start, tid_exon_stop, \(x, y) x <= ORF_gstop_0base & ORF_gstop_0base <= y)

    tid_exon_start <- tid_exon_start[which(tiv_start):which(tiv_stop)] |> map_chr(\(x) as.character(x))
    tid_exon_stop <- tid_exon_stop[which(tiv_start):which(tiv_stop)] |> map_chr(\(x) as.character(x))
    tid_exon_start[1] <- ORF_gstart_0base
    tid_exon_stop[length(tid_exon_stop)] <- ORF_gstop_0base
    tid_exon_start <- tid_exon_start |> map_chr(\(x) as.character(x))
    tid_exon_stop <- tid_exon_stop |> map_chr(\(x) as.character(x))

    stringr::str_c(tid_exon_start, tid_exon_stop, sep = "-", collapse = ",")
}

## unique func for ribohmm,remove duplicates and get_coor

remove_duplicates <- function(df) {
    tmp <- df |>
        dplyr::mutate(
            start_codon = case_when(
                stringr::str_sub(ORF_sequence_aa, 1, 1) == "M" ~ "ATG",
                stringr::str_sub(ORF_sequence_aa, 1, 1) == "L" ~ "CTG/TTG",
                stringr::str_sub(ORF_sequence_aa, 1, 1) == "V" ~ "GTG",
                .default = "others"
            )
        ) |>
        dplyr::select(chrom, coordinate_0base, strand, ORF_sequence_aa, start_codon, transcript_id) |>
        dplyr::distinct() |>
        tidyr::unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))

    grange_list <- tmp |>
        dplyr::select(coordinate_0base, chrom, strand, coordinate_id) |>
        tidyr::separate_longer_delim(coordinate_0base, delim = ",") |>
        tidyr::separate(coordinate_0base, c("start", "stop"), sep = "-") |>
        GenomicRanges::makeGRangesListFromDataFrame(
            split.field = "coordinate_id",
            keep.extra.columns = TRUE,
            starts.in.df.are.0based = TRUE
        )

    tmp_coor_aa <- tmp |>
        dplyr::select(coordinate_id, ORF_sequence_aa) |>
        distinct()

    hits <- findOverlaps(grange_list, grange_list, type = "within")

    if (length(grange_list) == length(hits)) {
        return(tmp |> dplyr::select(-coordinate_id))
    }

    total <- grange_list[-hits@to[duplicated(hits@from)]]
    ## check if aa_seq contains stop *
    giv_candidate_orf <- names(total)

    hits <- hits |>
        as_tibble() |>
        filter(!(queryHits == subjectHits))
    ## grange_list order do not match tmp order,using names to match
    check_outframe <- hits |>
        dplyr::mutate(
            queryHits_orf = names(grange_list[queryHits]),
            subjectHits_orf = names(grange_list[subjectHits])
        ) |>
        dplyr::inner_join(tmp_coor_aa, by = c("queryHits_orf" = "coordinate_id")) |>
        dplyr::rename("queryHits_orf_aa" = "ORF_sequence_aa") |>
        dplyr::inner_join(tmp_coor_aa, by = c("subjectHits_orf" = "coordinate_id")) |>
        dplyr::rename("subjectHits_orf_aa" = "ORF_sequence_aa") |>
        dplyr::mutate(included_info = str_detect(subjectHits_orf_aa, queryHits_orf_aa)) |>
        dplyr::group_by(queryHits_orf) |>
        dplyr::filter(base::all(included_info == FALSE)) |>
        dplyr::ungroup() |>
        dplyr::select(queryHits_orf) |>
        dplyr::distinct()

    giv_candidate_orf <- c(giv_candidate_orf, check_outframe$queryHits_orf)

    giv_candidate_orf_df <- tmp |>
        dplyr::filter(coordinate_id %in% giv_candidate_orf) |>
        dplyr::select(-coordinate_id)

    return(giv_candidate_orf_df)
}

tid2gid <- read_tsv(
    "/home/chengennong/Ribo_benchmark/ref/Human/gencode.v43.looktable",
    col_names = F,
    col_select = c(
        "transcript_id" = "X1",
        "gene_id" = "X2"
    )
)

preprocessing <- possibly(function(block) {
    ORF_blocks <- read_delim(block, col_names = T, col_types = cols(.default = "c")) |>
        dplyr::filter(as.numeric(posterior) > 8000) |>
        dplyr::rename(
            chrom = chromosome,
            ORF_gstart_0base = cdstart,
            ORF_gstop_0base = cdstop,
            tid_block_length = exon_sizes,
            tid_block_start = exon_starts,
            ORF_sequence_aa = protein_seq
        ) |>
        dplyr::mutate(
            coordinate_0base = furrr::future_pmap_chr(list(tid_block_length, tid_block_start, ORF_gstart_0base, ORF_gstop_0base), get_coordinate)
        ) |>
        dplyr::inner_join(tid2gid, by = c("transcript_id" = "transcript_id"))
})

format_ORF <- function(ORF_blocks, formatted_ORF) {
    ORF_blocks <- ORF_blocks |>
        dplyr::mutate(
            start_codon = case_when(
                stringr::str_sub(ORF_sequence_aa, 1, 1) == "M" ~ "ATG",
                stringr::str_sub(ORF_sequence_aa, 1, 1) == "L" ~ "CTG/TTG",
                stringr::str_sub(ORF_sequence_aa, 1, 1) == "V" ~ "GTG",
                .default = "others"
            )
        ) |>
        dplyr::select(gene_id, transcript_id, chrom, coordinate_0base, strand, start_codon, ORF_sequence_aa) |>
        dplyr::distinct() |>
        tidyr::unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::mutate(ORF_sequence_correct = "NA") |>
        dplyr::filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))
    output_func(ORF_blocks, formatted_ORF)
}

merge_ORF <- function(ORF_blocks, merged_ORF) {
    ORF_blocks <- ORF_blocks |>
        dplyr::group_nest(gene_id) |>
        dplyr::mutate(duplicated_data = data |> furrr::future_map(remove_duplicates)) |>
        dplyr::select(-data) |>
        tidyr::unnest(duplicated_data) |>
        dplyr::mutate(ORF_sequence_correct = "NA")
    output_func(ORF_blocks, merged_ORF)
}

ORF_processed <- preprocessing(
    block = snakemake@input[[1]]
)

format_ORF(
    ORF_blocks = ORF_processed,
    formatted_ORF = snakemake@output[[1]]
)

merge_ORF(
    ORF_blocks = ORF_processed,
    merged_ORF = snakemake@output[[2]]
)
