.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 100 * 1024 * 1024^2)
options(scipen = 999)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

base1_to_base0 <- function(coordinate_1base) {
    tmp <- str_split(coordinate_1base, ",") |> unlist()
    start <- tmp |>
        str_split_i("-", 1) |>
        map_int(\(x) as.numeric(x) - 1) |>
        map_chr(as.character)
    stop <- tmp |> str_split_i("-", 2)
    stringr::str_c(start, stop, sep = "-", collapse = ",")
}

## ribotricer use 1-based coordinate; need convert
remove_duplicates <- function(df) {
    tmp <- df |>
        dplyr::select(chrom, coordinate_0base, strand, ORF_sequence, start_codon, transcript_id) |>
        dplyr::distinct() |>
        tidyr::unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE)

    tmp_coor_aa <- tmp |>
        dplyr::select(coordinate_id, ORF_sequence) |>
        dplyr::distinct() |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        filter(!stringr::str_detect(ORF_sequence_aa, "\\*")) |>
        dplyr::select(-ORF_sequence) |>
        dplyr::distinct()

    tmp_coordinate_id <- tmp_coor_aa |> dplyr::pull(coordinate_id)
    tmp <- tmp |>
        filter(coordinate_id %in% tmp_coordinate_id)

    grange_list <- tmp |>
        dplyr::select(coordinate_0base, chrom, strand, coordinate_id) |>
        dplyr::distinct() |>
        tidyr::separate_longer_delim(coordinate_0base, delim = ",") |>
        tidyr::separate(coordinate_0base, c("start", "stop"), sep = "-") |>
        GenomicRanges::makeGRangesListFromDataFrame(
            split.field = "coordinate_id",
            keep.extra.columns = TRUE,
            starts.in.df.are.0based = TRUE
        )

    hits <- findOverlaps(grange_list, grange_list, type = "within")

    if (length(grange_list) == length(hits)) {
        tmp <- tmp |>
            dplyr::mutate(
                ORF_sequence_aa = ORF_sequence |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
            ) |>
            filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))
        return(tmp)
    }

    ## <U+4F7F><U+7528>findOverlaps <U+5E76><U+6307><U+5B9A>within, <U+5F97><U+5230><U+5404><U+4E2A>block<U+4E4B><U+95F4><U+7684><U+60C5><U+51B5>
    ## <U+4F46><U+662F><U+6CA1><U+529E><U+6CD5><U+5224><U+5B9A><U+7A76><U+7ADF><U+662F><U+54EA><U+4E2A><U+8D77><U+59CB><U+4E86>
    total <- grange_list[-hits@to[duplicated(hits@from)]]
    ## check if aa_seq contains stop *
    giv_candidate_orf <- names(total)

    hits <- hits |>
        as_tibble() |>
        filter(!(queryHits == subjectHits))
    ## grange_list order do not match tmp order,using names to match
    check_outframe <- hits |>
        mutate(
            queryHits_orf = names(grange_list[queryHits]),
            subjectHits_orf = names(grange_list[subjectHits])
        ) |>
        inner_join(tmp_coor_aa, by = c("queryHits_orf" = "coordinate_id")) |>
        dplyr::rename("queryHits_orf_aa" = "ORF_sequence_aa") |>
        inner_join(tmp_coor_aa, by = c("subjectHits_orf" = "coordinate_id")) |>
        dplyr::rename("subjectHits_orf_aa" = "ORF_sequence_aa") |>
        mutate(included_info = str_detect(subjectHits_orf_aa, queryHits_orf_aa)) |>
        group_by(queryHits_orf) |>
        filter(base::all(included_info == FALSE)) |>
        ungroup() |>
        select(queryHits_orf) |>
        distinct()

    ## <U+5BF9><U+627E><U+5230><U+7684><U+5305><U+542B><U+5173><U+7CFB><U+90FD><U+505A><U+4E00><U+6B21><U+5224><U+65AD>,<U+53EA><U+8981><U+6709><U+4E00><U+6B21><U+88AB><U+5305><U+542B><U+5C31><U+4E22><U+6389>

    giv_candidate_orf <- c(giv_candidate_orf, check_outframe$queryHits_orf)

    giv_candidate_orf_df <- tmp |>
        filter(coordinate_id %in% giv_candidate_orf)

    rm(tmp)
    gc()

    giv_candidate_orf_df <- giv_candidate_orf_df |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))
    return(giv_candidate_orf_df)
}

preprocessing <- possibly(function(withseq) {
    tmp_file <- read_tsv(withseq, col_names = T)
    tmp_file <- tmp_file |>
        dplyr::select(gene_id, transcript_id, chrom, coordinate, strand, ORF_sequence, start_codon) |>
        dplyr::distinct() |>
        dplyr::mutate(
            coordinate_0base = furrr::future_map_chr(coordinate, base1_to_base0),
        ) |>
        dplyr::select(-coordinate)

    gc()
    return(tmp_file)
})


format_ORF <- function(ORF_blocks, formatted_ORF) {
    ORF_blocks <- ORF_blocks |>
        select(gene_id, transcript_id, chrom, coordinate_0base, strand, ORF_sequence, start_codon) |>
        dplyr::distinct() |>
        unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        filter(!stringr::str_detect(ORF_sequence_aa, "\\*")) |>
        dplyr::rename(ORF_sequence_correct = ORF_sequence)
    output_func(ORF_blocks, formatted_ORF)
}

merge_ORF <- function(ORF_blocks, merged_ORF) {
    ORF_blocks <- ORF_blocks |>
        group_nest(gene_id) |>
        mutate(
            duplicated_data = data |> furrr::future_map(remove_duplicates)
        ) |>
        select(-data) |>
        unnest(duplicated_data) |>
        dplyr::rename(ORF_sequence_correct = ORF_sequence)
    output_func(ORF_blocks, merged_ORF)
}

ORF_processed <- preprocessing(withseq = snakemake@input[[1]])

format_ORF(
    ORF_blocks = ORF_processed,
    formatted_ORF = snakemake@output[[1]]
)

merge_ORF(
    ORF_blocks = ORF_processed,
    merged_ORF = snakemake@output[[2]]
)
