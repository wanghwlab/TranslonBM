.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 80 * 1024 * 1024^2)
options(scipen = 999)

get_coordinate <- function(ORF_gstart_0base, block_length, block_start, strand) {
    ORF_gstart_0base <- as.numeric(ORF_gstart_0base)
    block_length <- block_length |>
        str_split_1(",") |>
        map_int(\(x) as.numeric(x))
    ORF_exon_start <- block_start |>
        str_split_1(",") |>
        map_int(\(x) as.numeric(x) + ORF_gstart_0base)
    ORF_exon_stop <- ORF_exon_start + block_length

    ## need remove stop codon
    if (strand == "-") {
        if (block_length[1] < 4) {
            ORF_exon_start[2] <- ORF_exon_start[2] + (3 - block_length[1])
            ORF_exon_start <- ORF_exon_start[-1]
            ORF_exon_stop <- ORF_exon_stop[-1]
        } else {
            ORF_exon_start[1] <- ORF_exon_start[1] + 3
        }
    } else {
        block_nums <- length(block_length)
        if (block_length[block_nums] < 4) {
            ORF_exon_stop[block_nums - 1] <- ORF_exon_stop[block_nums - 1] - (3 - block_length[block_nums])
            ORF_exon_start <- ORF_exon_start[-block_nums]
            ORF_exon_stop <- ORF_exon_stop[-block_nums]
        } else {
            ORF_exon_stop[block_nums] <- ORF_exon_stop[block_nums] - 3
        }
    }

    stringr::str_c(ORF_exon_start, ORF_exon_stop, sep = "-", collapse = ",")
}

remove_duplicates <- function(df) {
    tmp <- df |>
        dplyr::select(chrom, coordinate_0base, strand, ORF_sequence_correct, start_codon) |>
        distinct() |>
        tidyr::unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence_correct |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        dplyr::filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))
    ## check if aa_seq contains stop *

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
    "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.looktable",
    col_names = F,
    col_select = c(
        "transcript_id" = "X1",
        "gene_id" = "X2"
    )
)

## need improve

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

reformat_ORF <- possibly(function(addchr_bed, ORF_sequence_tsv, ORF_gcoor_reformat) {
    rawfile_seq <- read_tsv(
        ORF_sequence_tsv,
        col_names = F,
        col_types = cols(.default = "c"),
        col_select = c(
            "ORF_id" = "X1",
            "ORF_sequence" = "X2"
        )
    ) |>
        dplyr::mutate(
            ORF_id = stringr::str_sub(ORF_id, end = -4L)
        )

    ORF_blocks <- read_tsv(
        addchr_bed,
        col_names = F,
        col_types = cols(.default = "c")
    ) |>
        dplyr::rename(
            chrom = X1,
            ORF_gstart_0base = X2,
            ORF_gstop_0base = X3,
            ORF_id = X4,
            strand = X6,
            block_length = X11,
            block_start = X12
        ) |>
        dplyr::inner_join(rawfile_seq, by = c("ORF_id" = "ORF_id")) |>
        dplyr::mutate(
            block_length = stringr::str_sub(block_length, end = -2L),
            block_start = stringr::str_sub(block_start, end = -2L),
            coordinate_0base = furrr::future_pmap_chr(list(ORF_gstart_0base, block_length, block_start, strand), get_coordinate)
        ) |>
        tidyr::separate(ORF_id, c("transcript_id", "ORF_type", "tmp"), sep = "_", remove = F) |>
        dplyr::inner_join(tid2gid, by = c("transcript_id" = "transcript_id")) |>
        dplyr::mutate(
            ORF_sequence_correct = stringr::str_sub(ORF_sequence, end = -4L),
            start_codon = stringr::str_sub(ORF_sequence, 1, 3)
        )

    ORF_blocks <- ORF_blocks |>
        dplyr::select(gene_id, chrom, coordinate_0base, strand, ORF_sequence_correct, start_codon) |>
        distinct() |>
        tidyr::unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence_correct |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        dplyr::filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))

    ## here we do not need to combine ORF in block format
    output_func(ORF_blocks, ORF_gcoor_reformat)
})

reformat_ORF(
    addchr_bed = snakemake@input[[1]],
    ORF_sequence_tsv = snakemake@input[[2]],
    ORF_gcoor_reformat = snakemake@output[[1]]
)
