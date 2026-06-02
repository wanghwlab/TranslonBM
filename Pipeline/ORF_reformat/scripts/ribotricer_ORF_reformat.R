.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 60 * 1024 * 1024^2)
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
        select(chrom, coordinate_0base, strand, ORF_sequence, start_codon) |>
        distinct() |>
        unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))

    ## translating<U+91CC><U+9700><U+8981><U+6307><U+5B9A>no.init.codon=T
    ## <U+5426><U+5219>TTG<U+8D77><U+59CB><U+7FFB><U+8BD1><U+51FA><U+6765><U+4E0D><U+662F>L,<U+800C><U+662F>M
    ## <U+8FC7><U+6EE4><U+6389><U+542B><U+6709><U+7EC8><U+6B62><U+5BC6><U+7801><U+5B50><U+7684><U+5E8F><U+5217>

    grange_list <- tmp |>
        select(-ORF_sequence_aa, -start_codon) |>
        tidyr::separate_longer_delim(coordinate_0base, delim = ",") |>
        tidyr::separate(coordinate_0base, c("start", "stop"), sep = "-") |>
        GenomicRanges::makeGRangesListFromDataFrame(
            split.field = "coordinate_id",
            keep.extra.columns = TRUE,
            starts.in.df.are.0based = TRUE
        )

    tmp_coor_aa <- tmp |> select(coordinate_id, ORF_sequence_aa)

    hits <- findOverlaps(grange_list, grange_list, type = "within")

    if (length(grange_list) == length(hits)) {
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
        filter(coordinate_id %in% giv_candidate_orf) |>
        select(-coordinate_id)

    return(giv_candidate_orf_df)
}

duplicate_ORF <- possibly(function(withseq, ORF_gcoor_reformat) {
    tmp_file <- read_tsv(withseq, col_names = T)
    tmp_file <- tmp_file |>
        select(gene_id, chrom, coordinate, strand, ORF_sequence, start_codon) |>
        distinct() |>
        dplyr::rename(coordinate_1base = coordinate) |>
        dplyr::mutate(
            coordinate_0base = furrr::future_map_chr(coordinate_1base, base1_to_base0),
        ) |>
        dplyr::select(-coordinate_1base)

    gc()

    tmp_file <- tmp_file |>
        select(gene_id, chrom, coordinate_0base, strand, ORF_sequence, start_codon) |>
        distinct() |>
        unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))

    output_func(tmp_file, ORF_gcoor_reformat)
})

duplicate_ORF(
    withseq = snakemake@input[[1]],
    ORF_gcoor_reformat = snakemake@output[[1]]
)
