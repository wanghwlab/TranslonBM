.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(furrr)
plan(multisession, workers = 24)
options("future.globals.maxSize" = 80 * 1024 * 1024^2)
options(scipen = 999)

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
        dplyr::select(chrom, coordinate_0base, strand, ORF_sequence_correct, start_codon, transcript_id) |>
        dplyr::distinct() |>
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

# Load Biostrings explicitly at the start of the script (GenomicRanges usually depends on it)
library(Biostrings)

preprocessing <- possibly(function(raw_file, ORF_fa) {
    
    # --- Modification begins: read the FASTA file correctly ---
    # Read FASTA as a DNAStringSet object
    fa_obj <- readDNAStringSet(ORF_fa)
    
    # Convert to a Tibble for the subsequent join
    # Note: names(fa_obj) usually contains header information; ensure it matches the id in raw_file
    # Assume the id in raw_file matches the first space-delimited field of the FASTA header
    ORF_sequence <- tibble(
        id = names(fa_obj) |> str_split_i(" ", 1), # Extract the ID before the first space; adjust for the actual header format
        ORF_sequence = as.character(fa_obj)
    )
    # --- Modification ends ---

    ORF_blocks <- read_tsv(
        raw_file,
        col_names = T,
        col_types = cols(.default = "c")
    ) |>
        dplyr::rename_with(~ stringr::str_sub(.x, 2, -1), everything()) |>
        # The original inner_join logic is unchanged
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
        dplyr::inner_join(tid2gid, by = c("transcript_id" = "transcript_id"))
    
    return(ORF_blocks) # Return explicitly
}, otherwise = NULL) # Specify otherwise explicitly

format_ORF <- function(ORF_blocks, formatted_ORF) {
    ORF_blocks <- ORF_blocks |>
        dplyr::select(gene_id, transcript_id, chrom, coordinate_0base, strand, ORF_sequence_correct, start_codon) |>
        dplyr::distinct() |>
        tidyr::unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE) |>
        dplyr::mutate(
            ORF_sequence_aa = ORF_sequence_correct |> Biostrings::DNAStringSet() |> Biostrings::translate(no.init.codon = T) |> map_chr(\(x) toString(x))
        ) |>
        dplyr::filter(!stringr::str_detect(ORF_sequence_aa, "\\*"))
    output_func(ORF_blocks, formatted_ORF)
}

merge_ORF <- function(ORF_blocks, merged_ORF) {
    ORF_blocks <- ORF_blocks |>
        dplyr::group_nest(gene_id) |>
        dplyr::mutate(duplicated_data = data |> furrr::future_map(remove_duplicates)) |>
        dplyr::select(-data) |>
        tidyr::unnest(duplicated_data)
    output_func(ORF_blocks, merged_ORF)
}


#target_dir <- snakemake@params[["fa_dir"]]

# Modification: stop matching the sample prefix and find any .fa file in this directory
# ignore.case = TRUE matches both .fa and .FA
#found_fa <- list.files(target_dir, pattern = "\\.fa$", full.names = TRUE, ignore.case = TRUE)[1]

# Print the path for debugging so the log shows where the file search occurs
#message("Searching the directory for an FA file: ", target_dir)
#message("Found file: ", found_fa)

#if (is.na(found_fa)) {
##    # If no file is found, list the directory contents for troubleshooting
#    existing_files <- list.files(target_dir)
#    stop(paste0("Error: no .fa file was found in directory ", target_dir, ".\nFiles in the directory: ", paste(existing_files, collapse = ", ")))
#}

ORF_processed <- preprocessing(
    raw_file = '/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_1core/rpbp/rpbp_chrN/orf_pred_default/SRX876063_SRX876069_STAR/orf-predictions/SRX876063_SRX876069_STAR_rpbp_filtered_raw.txt',
    #ORF_fa = snakemake@input[[2]]
    ORF_fa = '/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_1core/rpbp/rpbp_chrN/orf_pred_default/SRX876063_SRX876069_STAR/orf-predictions/SRX876063_SRX876069_STAR.length-23-25-28-29.offset-9-9-12-12.filtered.predicted-orfs.dna.fa'
)

format_ORF(
    ORF_blocks = ORF_processed,
    formatted_ORF = 'SRX876063_SRX876069_STAR_rpbp_gcoor.tsv.gz'
)

merge_ORF(
    ORF_blocks = ORF_processed,
    merged_ORF = 'SRX876063_SRX876069_STAR_rpbp_merged_gcoor.tsv.gz'
)
