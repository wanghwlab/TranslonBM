.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 80 * 1024 * 1024^2)
options(scipen = 999)

mapping_software <- c("STAR", "hisat2", "tophat2")
bench <- c("rpbp")
bench <- str_c(bench, "_chrN")
pred <- c("orf_pred_default", "orf_pred_undisputed")
path <- "/home/chengennong/Ribo_benchmark/unique_mapping/Human/preprocessing"
RPF_samples <- tibble::tribble(
    ~samples,
    "SRX11812007_SRX11812008_SRX11812009",
    "SRX1254413",
    "SRX1447296",
    "SRX5256543_SRX5256555",
    "SRX5256553_SRX5256554",
    "SRX5887328_SRX5887329_SRX5887330",
    "SRX740748",
    "SRX7666669-73",
    "SRX7666674-78",
    "SRX7666679-83",
    "SRX7666684-88",
    "SRX7666689-93",
    "SRX7666694-98",
    "SRX876063_SRX876069"
)
sample_info <- tidyr::expand_grid(RPF_samples, bench, mapping_software, pred, path)

sample_info_default <- sample_info |> filter(pred == "orf_pred_default")
sample_info_undisputed <- sample_info |> filter(pred == "orf_pred_undisputed")

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

get_coordinate <- function(ORF_gstart_0base, block_length, block_start) {
    ORF_gstart_0base <- as.numeric(ORF_gstart_0base)
    block_length <- block_length |>
        str_split_1(",") |>
        map_int(\(x) as.numeric(x))
    ORF_exon_start <- block_start |>
        str_split_1(",") |>
        map_int(\(x) as.numeric(x) + ORF_gstart_0base)
    ORF_exon_stop <- ORF_exon_start + block_length

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
    "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.looktable",
    col_names = F,
    col_select = c(
        "transcript_id" = "X1",
        "gene_id" = "X2"
    )
)

duplicate_ORF <- possibly(function(raw_file, ORF_fa, ORF_gcoor_reformat) {
    ORF_sequence <- read_tsv(
        ORF_fa,
        col_names = F,
        col_select = c(
            "id" = "X1",
            "ORF_sequence" = "X2"
        )
    )

    ORF_blocks <- read_tsv(
        raw_file,
        col_names = T,
        col_types = cols(.default = "c")
    ) |>
        dplyr::rename_with(~ stringr::str_sub(.x, 2, -1), everything()) |>
        dplyr::inner_join(ORF_sequence, by = c("id" = "id")) |>
        dplyr::rename(
            chrom = seqname,
            ORF_gstart_0base = start,
            ORF_gstop_0base = end,
            strand = strand,
            block_length = exon_lengths,
            block_start = exon_genomic_relative_starts
        ) |>
        dplyr::mutate(
            coordinate_0base = furrr::future_pmap_chr(list(ORF_gstart_0base, block_length, block_start), get_coordinate),
            ORF_sequence_correct = ORF_sequence,
            start_codon = stringr::str_sub(ORF_sequence_correct, 1, 3)
        )

    ORF_blocks <- ORF_blocks |>
        tidyr::separate(id, c("transcript_id", "pos"), sep = "_") |>
        dplyr::inner_join(tid2gid, by = c("transcript_id" = "transcript_id")) |>
        dplyr::select(gene_id, chrom, coordinate_0base, strand, ORF_sequence_correct, start_codon) |>
        distinct() |>
        tidyr::unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence_correct |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        dplyr::filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))

    output_func(ORF_blocks, ORF_gcoor_reformat)
})

duplicate_ORF(
    raw_file = snakemake@input[[1]],
    ORF_fa = snakemake@input[[2]],
    ORF_gcoor_reformat = snakemake@output[[1]]
)
