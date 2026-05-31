# ==============================================================================
# Simulation.R — A Flexible Ribo-seq Data Simulator
# Version : 1.0.0
# Requires: R >= 4.1, argparse, furrr, tidyverse; Python >= 3.6 + gppy
# Usage：
#   Rscript Simulation.R --help
#   Rscript Simulation.R \
#       --counts      6000000 \
#       --replicate   T1 \
#       --Path_dir    $NMHOME/SimulationInput \
#       --threads     8
#       --out_dir     $NMHOME/SimulationOutput
# ==============================================================================

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(
    description = "A Flexible Ribo-seq Data Simulator"
)

# Required Parameters
parser$add_argument(
    "--counts",
    type     = "integer",
    required = TRUE,
    metavar  = "INT",
    help     = "Simulate number of reads, such as 6,000,000 or 60,000,000"
)
parser$add_argument(
    "--replicate",
    type     = "character",
    required = TRUE,
    metavar  = "Label",
    help     = "Replicate label, e.g. T1, T2, T3"
)
parser$add_argument(
    "--Path_dir",
    type     = "character",
    required = TRUE,
    metavar  = "PATH",
    help     = "Project root directory, containing all required input files"
)
parser$add_argument(
    "--threads",
    type     = "integer",
    required = TRUE,
    metavar  = "INT",
    help     = "Number of threads"
)

# Optional Parameters
parser$add_argument(
    "--gppy_script",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "gppy gtf.py Script path"
)
parser$add_argument(
    "--out_dir",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Output directory (specified separately; if not specified, it defaults to Path_dir)"
)

# Simulated quality parameters
parser$add_argument(
    "--seed",
    type    = "integer",
    default = 42,
    metavar = "INT",
    help    = "Global random seed to ensure reproducible results [Default: 42]"
)
parser$add_argument(
    "--n_transcripts",
    type    = "integer",
    default = 12000,
    metavar = "INT",
    help    = "Number of randomly selected transcripts [Default: 12,000]"
)
parser$add_argument(
    "--orfs_num",
    type    = "integer",
    default = 12000,
    metavar = "INT",
    help    = "Simulated number of ORFs, typically matching n_transcripts [default: 12,000]"
)
parser$add_argument(
    "--perturb_cutoff",
    type    = "double",
    default = 0.02,
    metavar = "FLOAT",
    help    = "Upper limit on the percentage of non-periodic reads [Default: 0.02]"
)
parser$add_argument(
    "--snp_ratio",
    type    = "double",
    default = 0.001,
    metavar = "FLOAT",
    help    = "Probability of an SNP occurring [Default: 0.001]"
)
parser$add_argument(
    "--indel_ratio",
    type    = "double",
    default = 0.0005,
    metavar = "FLOAT",
    help    = "Probability of indels [Default: 0.0005]"
)
parser$add_argument(
    "--insertion_ratio",
    type    = "double",
    default = 0.00025,
    metavar = "FLOAT",
    help    = "Probability of an insertion event in the indel [default: 0.00025]"
)
parser$add_argument(
    "--seq_error_ratio",
    type    = "double",
    default = 0.005,
    metavar = "FLOAT",
    help    = "Sequencing error probability [Default: 0.005]"
)
parser$add_argument(
    "--max_iter",
    type    = "integer",
    default = 1000,
    metavar = "INT",
    help    = "Maximum number of iterations for `generate_ORF_counts` to prevent an infinite loop [Default: 1,000]"
)

args <- parser$parse_args()

options(scipen = 999)
options("future.globals.maxSize" = 80 * 1024 * 1024^2)

`%||%` <- function(a, b) if (!is.null(a)) a else b

Path_dir <- args$Path_dir

# Path
out_dir   <- args$out_dir   %||% Path_dir
gppy_script <- args$gppy_script %||% file.path(args$Path_dir, "gtf.py")

# Specify the input file path
path_coding_transcripts <- file.path(Path_dir, "coding_transcripts.tsv")
path_transcript_seq     <- file.path(Path_dir,   "gencode.v43.transcripts.seq")
path_transcripts_start  <- file.path(Path_dir, "coding_transcripts_start.tsv")
path_transcripts_stop   <- file.path(Path_dir, "coding_transcripts_stop.tsv")
path_codons_freq        <- file.path(Path_dir,  "codons_freq.tsv")
path_gtf                <- file.path(Path_dir,   "gencode.v43.annotation.gtf")

# Lengths corresponding to the 5'-end offset positions: a1=9, a2=10, a3=11, a4=12, a5=13, a6=14 nt
vec_d5 <- c(0.06018293, 0.05040521, 0.09572164, 0.74742327, 0.02618989, 0.02007706)
# Lengths corresponding to the 3'-end offset positions: b1=13, b2=14, b3=15, b4=16, b5=17, b6=18 nt
vec_d3 <- c(0.19138088, 0.53849728, 0.21173048, 0.03839137, 0.01000000, 0.01000000)

output_prefix <- paste0("simulation_", args$counts / 1e6, "M_", args$replicate)


message("========================================")
message(" counts_size  : ", args$counts)
message(" replicate    : ", args$replicate)
message(" threads      : ", args$threads)
message(" seed         : ", args$seed)
message(" Path_dir     : ", args$Path_dir)
message(" out_dir      : ", out_dir)
message("========================================\n")

suppressPackageStartupMessages({
    library(furrr)
    library(tidyverse)
})

# Detect the number of available cores to prevent the number of requests from exceeding actual resources
available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores)) available_cores <- 1L

safe_threads <- min(args$threads, max(1L, available_cores - 1L))

if (safe_threads < args$threads) {
    warning(
        "Requested threads = ", args$threads,
        ",Machine has only ", available_cores, " available_cores ",
        "\nReduced to ", safe_threads, " threads"
    )
}

if (safe_threads <= 1L) {
    plan(sequential)
    message("Sequential")
} else {
    plan(multisession, workers = safe_threads)
    message("Workers: ", safe_threads, " / ", available_cores, " cores available")
}

set.seed(args$seed)
message("Global seed: ", args$seed)

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message("Output directory: ", out_dir)
}

setwd(Path_dir)
message("Processing directory: ", Path_dir)

required_files <- c(
    path_coding_transcripts,
    path_transcript_seq,
    path_transcripts_start,
    path_transcripts_stop,
    path_codons_freq,
    path_gtf
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
    stop(
        "The following required input files are missing:\n",
        paste(" -", missing_files, collapse = "\n"),
        "\nPlease check --Path_dir. to override individual paths.",
        call. = FALSE
)
}


message("Input files are present; Continuing")


nucleotides <- c("A", "C", "G", "T")
codons      <- expand.grid(nucleotides, nucleotides, nucleotides)
codons      <- apply(codons, 1, paste, collapse = "")
codons      <- codons[!codons %in% c("TAA", "TAG", "TGA")]

# codons_freq is used to weight the sampling of P sites based on codon frequency
codons_freq <- read_tsv(path_codons_freq, col_names = TRUE, show_col_types = FALSE)


message("[Step 1 of 7] Loading transcript data...")

selected_transcripts <- read_tsv(
    path_coding_transcripts,
    col_names      = FALSE,
    col_select     = c(transcript_id = X1),
    show_col_types = FALSE
)

transcripts_seq <- read_tsv(
    path_transcript_seq,
    col_names      = FALSE,
    col_select     = c(transcript_id = X1, transcript_seq = X2),
    show_col_types = FALSE
)

transcripts_start <- read_tsv(
    path_transcripts_start,
    col_names      = TRUE,
    col_select     = c(transcript_id = Key, ORF_tstart = Value),
    show_col_types = FALSE
)

transcripts_stop <- read_tsv(
    path_transcripts_stop,
    col_names      = TRUE,
    col_select     = c(transcript_id = Key, ORF_tstop = Value),
    show_col_types = FALSE
)

selected_transcripts <- selected_transcripts |>
    inner_join(transcripts_seq,   by = "transcript_id") |>
    inner_join(transcripts_start, by = "transcript_id") |>
    inner_join(transcripts_stop,  by = "transcript_id") |>
    filter(
        as.numeric(ORF_tstart) > 14,
        (str_length(transcript_seq) - as.numeric(ORF_tstop)) > 25
    )

message("[Step 1 of 7] Number of transcripts after filtering: ", nrow(selected_transcripts))

selected_transcripts <- selected_transcripts[
    sample(nrow(selected_transcripts), args$n_transcripts, replace = FALSE),
]
message("[Step 1 of 7] Number of sampled transcripts: ", nrow(selected_transcripts))

rm(transcripts_seq, transcripts_start, transcripts_stop)
gc()


message("[Step 2 of 7] Building a model of read length distribution...")

names(vec_d5) <- str_c("a", 1:6)
names(vec_d3) <- str_c("b", 1:6)

prop_mat <- outer(vec_d5, vec_d3)

vec_detail <- expand_grid(vec_d5 = 1:6, vec_d3 = 1:6) |>
    mutate(index = vec_d5 + vec_d3) |>
    filter(between(index, 2, 12)) |>
    mutate(
        vec_d5 = str_c("a", vec_d5),
        vec_d3 = str_c("b", vec_d3)
    ) |>
    arrange(index)

extract_element <- function(a, b) prop_mat[a, b]

vec_detail <- vec_detail |>
    mutate(
        prop         = map2_dbl(vec_d5, vec_d3, extract_element),
        reads_length = index + 23
    ) |>
    group_by(reads_length) |>
    mutate(reads_length_ratio = sum(prop)) |>
    ungroup() |>
    mutate(
        reads_length_d5 = case_match(
            vec_d5,
            "a1" ~ 9,  "a2" ~ 10, "a3" ~ 11,
            "a4" ~ 12, "a5" ~ 13, "a6" ~ 14
        ),
        reads_length_d3 = case_match(
            vec_d3,
            "b1" ~ 13, "b2" ~ 14, "b3" ~ 15,
            "b4" ~ 16, "b5" ~ 17, "b6" ~ 18
        )
    )

select_offset <- function(df) {
    offset_temp <- df |> arrange(desc(prop)) |> slice(1) |> pull(reads_length_d5)
    df |> mutate(offset = offset_temp)
}

vec_detail_qc_frame <- vec_detail |>
    group_nest(reads_length) |>
    mutate(data = map(data, select_offset)) |>
    unnest(data) |>
    mutate(
        frame_shift = offset - reads_length_d5,
        frame_qc    = (offset - reads_length_d5) %% 3
    ) |>
    group_by(reads_length, frame_qc) |>
    summarise(frame_count = sum(prop), .groups = "drop") |>
    pivot_wider(
        names_from   = frame_qc,
        names_expand = TRUE,
        names_prefix = "frame",
        values_from  = frame_count,
        values_fill  = 0
    ) |>
    rowwise() |>
    mutate(reads_length_ratio = sum(frame0, frame1, frame2)) |>
    ungroup() |>
    mutate(across(c(frame0, frame1, frame2), ~ .x / reads_length_ratio)) |>
    rowwise() |>
    mutate(max_frame = max(c_across(c(frame0, frame1, frame2))))


message("[Step 3 of 7] Allocate periodic / non-periodic reads...")

vec_detail_for_simulation <- vec_detail |>
    mutate(
        perturb = ifelse(
            args$perturb_cutoff * (prop / reads_length_ratio) > prop,
            prop,
            args$perturb_cutoff * (prop / reads_length_ratio)
        ),
        periodic              = prop - perturb,
        non_periodic_insitu   = perturb * 1 / 3,
        non_periodic_forward  = perturb * 1 / 3,
        non_periodic_backward = perturb * 1 / 3
    ) |>
    mutate(across(
        c(periodic, non_periodic_insitu, non_periodic_forward, non_periodic_backward),
        ~ round(.x * args$counts)
    ))

frame_expand <- function(periodic, non_periodic_insitu, non_periodic_forward, non_periodic_backward) {
    bind_rows(
        tibble(frame_shift = rep( 0, periodic),              periodic_tag = "periodic"),
        tibble(frame_shift = rep( 0, non_periodic_insitu),   periodic_tag = "non_periodic"),
        tibble(frame_shift = rep( 1, non_periodic_forward),  periodic_tag = "non_periodic"),
        tibble(frame_shift = rep(-1, non_periodic_backward), periodic_tag = "non_periodic")
    ) |> slice_sample(prop = 1)
}

vec_df <- vec_detail_for_simulation |>
    mutate(
        final = pmap(
            list(periodic, non_periodic_insitu, non_periodic_forward, non_periodic_backward),
            frame_expand
        )
    ) |>
    unnest(final)

vec_df <- vec_df[sample(nrow(vec_df)), ]

vec_df_non_periodic <- vec_df |> filter(periodic_tag == "non_periodic")
vec_df_periodic     <- vec_df |> filter(periodic_tag == "periodic")

message("[Step 3 of 7] periodic reads: ", nrow(vec_df_periodic),
        " | non-periodic reads: ", nrow(vec_df_non_periodic))


generate_ORF_counts <- function(ORFs_num, target_df, max_iter = 1000) {
    target_n <- nrow(target_df)
    lambda   <- target_n / ORFs_num

    for (i in seq_len(max_iter)) {
        tmp <- ceiling(rexp(ORFs_num, rate = 1 / lambda))
        if (sum(tmp) < target_n && (sum(tmp) / target_n) > 0.99) return(tmp)
    }

    warning(
        "generate_ORF_counts: Reached the maximum number of iterations (", max_iter, ")，"
    )
    base_count <- floor(target_n / ORFs_num)
    remainder  <- target_n - base_count * ORFs_num
    counts     <- rep(base_count, ORFs_num)
    counts[seq_len(remainder)] <- counts[seq_len(remainder)] + 1
    return(counts)
}

periodic_count_distr     <- generate_ORF_counts(args$orfs_num, vec_df_periodic,     args$max_iter)
non_periodic_count_distr <- generate_ORF_counts(args$orfs_num, vec_df_non_periodic, args$max_iter)

# Precisely align the total count (add/subtract the excess/shortfall to/from the last ORF).
periodic_count_distr[args$orfs_num]     <- periodic_count_distr[args$orfs_num]     - (sum(periodic_count_distr)     - nrow(vec_df_periodic))
non_periodic_count_distr[args$orfs_num] <- non_periodic_count_distr[args$orfs_num] - (sum(non_periodic_count_distr) - nrow(vec_df_non_periodic))

if (sum(periodic_count_distr) != nrow(vec_df_periodic)) {
    stop("Internal error: periodic read count mismatch after ORF assignment. ",
         "Got ", sum(periodic_count_distr), ", expected ", nrow(vec_df_periodic),
         call. = FALSE)
}
if (sum(non_periodic_count_distr) != nrow(vec_df_non_periodic)) {
    stop("Internal error: non-periodic read count mismatch after ORF assignment. ",
         "Got ", sum(non_periodic_count_distr), ", expected ", nrow(vec_df_non_periodic),
         call. = FALSE)
}

vec_df_periodic     <- vec_df_periodic[sample(nrow(vec_df_periodic)), ]
vec_df_non_periodic <- vec_df_non_periodic[sample(nrow(vec_df_non_periodic)), ]

periodic_indices     <- rep(seq_along(periodic_count_distr),     times = periodic_count_distr)
non_periodic_indices <- rep(seq_along(non_periodic_count_distr), times = non_periodic_count_distr)

vec_df_periodic_indexed     <- base::split(vec_df_periodic,     factor(periodic_indices,     levels = seq_along(periodic_count_distr)))
vec_df_non_periodic_indexed <- base::split(vec_df_non_periodic, factor(non_periodic_indices, levels = seq_along(non_periodic_count_distr)))

vec_df_indexed <- map2(
    vec_df_non_periodic_indexed,
    vec_df_periodic_indexed,
    ~ bind_rows(.x, .y)
)

message("[Step 4 of 7] All ORFs have been assigned; total reads: ",
        sum(map_int(vec_df_indexed, nrow)))

selected_transcripts <- selected_transcripts |>
    mutate(RPF_count_info = vec_df_indexed)

detail_path <- file.path(out_dir, paste0(output_prefix, "_detail.tsv"))
selected_transcripts |>
    mutate(
        ORF_tstart = as.numeric(ORF_tstart) + 1,
        ORF_tstop  = as.numeric(ORF_tstop)
    ) |>
    select(-transcript_seq) |>
    unnest(cols = RPF_count_info) |>
    write_tsv(detail_path, col_names = TRUE, quote = "none")
message("[Step 4 of 7] Output intermediate detail file: ", detail_path)

gc()


get_codon <- function(ORF) {
    ORF_split   <- str_split_1(ORF, "")
    codon_index <- seq(1, str_length(ORF), by = 3)
    map_chr(codon_index, ~ str_c(ORF_split[.x:min(.x + 2, length(ORF_split))], collapse = ""))
}

add_snp <- function(reads_sequence) {
    nucs    <- c("A", "T", "C", "G")
    mut_pos <- which(rbinom(str_length(reads_sequence), 1, prob = args$snp_ratio) == 1)

    if (length(mut_pos) == 0) {
        return(tibble(reads_sequence = reads_sequence, snp_info = ""))
    }

    mut_raw  <- map_chr(mut_pos, ~ str_sub(reads_sequence, .x, .x))
    mut      <- map_chr(mut_raw, ~ sample(nucs[nucs != .x], 1))
    mut_info <- str_c(mut_raw, mut_pos, mut, sep = ":", collapse = ",")

    seq_vec <- str_split_1(reads_sequence, "")
    seq_vec[mut_pos] <- mut
    tibble(reads_sequence = str_c(seq_vec, collapse = ""), snp_info = mut_info)
}

add_indel <- function(reads_sequence) {
    nucs    <- c("A", "T", "C", "G")
    mut_pos <- which(rbinom(str_length(reads_sequence), 1, prob = args$indel_ratio) == 1)

    if (length(mut_pos) == 0) {
        return(tibble(reads_sequence = reads_sequence, indel_info = ""))
    }

    mut_raw  <- map_chr(mut_pos, ~ str_sub(reads_sequence, .x, .x))
    mut      <- map_chr(mut_raw, ~ ifelse(
        runif(1) < (args$insertion_ratio / args$indel_ratio),
        paste0(.x, sample(nucs, 1)),
        ""
    ))
    mut_info <- str_c(mut_raw, mut_pos, mut, sep = ":", collapse = ",")

    seq_vec <- str_split_1(reads_sequence, "")
    seq_vec[mut_pos] <- mut
    tibble(reads_sequence = str_c(seq_vec, collapse = ""), indel_info = mut_info)
}

add_sequencing_error <- function(reads_sequence) {
    nucs    <- c("A", "T", "C", "G", "N")
    mut_pos <- which(rbinom(str_length(reads_sequence), 1, prob = args$seq_error_ratio) == 1)

    if (length(mut_pos) == 0) {
        return(tibble(reads_sequence_final = reads_sequence, sequencing_error_info = ""))
    }

    mut_raw  <- map_chr(mut_pos, ~ str_sub(reads_sequence, .x, .x))
    mut      <- map_chr(mut_raw, ~ sample(nucs[nucs != .x], 1))
    mut_info <- str_c(mut_raw, mut_pos, mut, sep = ":", collapse = ",")

    seq_vec <- str_split_1(reads_sequence, "")
    seq_vec[mut_pos] <- mut
    tibble(reads_sequence_final = str_c(seq_vec, collapse = ""), sequencing_error_info = mut_info)
}

# Generate the sequences and coordinates of all reads mapping to this ORF
get_raw_reads <- function(ORF_tstart, ORF_tstop, transcript_seq, RPF_count_info, mutation = NULL) {
    cds <- str_sub(transcript_seq, ORF_tstart, ORF_tstop)

    cds_codon_prob <- get_codon(cds) |>
        as_tibble() |>
        rename(codon = value) |>
        inner_join(codons_freq, by = "codon") |>
        mutate(value = value / sum(value)) |>
        group_by(codon) |>
        summarise(value = sum(value))

    RPF_count_size <- nrow(RPF_count_info)

    choosen_codon <- sample(
        cds_codon_prob$codon,
        size    = RPF_count_size,
        replace = TRUE,
        prob    = cds_codon_prob$value
    )

    codon_index <- choosen_codon |>
        unique() |>
        map(~ c(
            which(get_codon(cds) == .x)[1],
            sample(
                which(get_codon(cds) == .x),
                size    = max(sum(choosen_codon == .x) - 1, 0),
                replace = TRUE
            )
        ))

    RPF_count_info |>
        mutate(
            codon_index      = reduce(codon_index, c),
            codon_frame0_pos = (codon_index - 1) * 3 + 1,
            codon_pos        = codon_frame0_pos + frame_shift,
            reads_start_pos  = (ORF_tstart - 1) + codon_pos - reads_length_d5,
            reads_stop_pos   = (ORF_tstart - 1) + codon_pos + 2 + reads_length_d3,
            reads_sequence   = str_sub(transcript_seq, reads_start_pos, reads_stop_pos)
        )
}

get_alt_reads <- function(RPF_count_info, mutation_data) {
    single_read_alt <- function(reads_start_pos, reads_stop_pos, reads_sequence, mutation_data) {
        mutation_data <- filter(mutation_data, between(pos, reads_start_pos, reads_stop_pos))
        if (nrow(mutation_data) == 0) {
            return(tibble(reads_sequence_alt = reads_sequence, mutation_info = ""))
        }
        mut_pos  <- pull(mutation_data, pos)
        alt      <- pull(mutation_data, alt)
        mut_info <- str_c(pull(mutation_data, mutation_meta), collapse = ",")
        seq_vec  <- str_split_1(reads_sequence, "")
        seq_vec[mut_pos - reads_start_pos + 1] <- alt
        tibble(reads_sequence_alt = str_c(seq_vec, collapse = ""), mutation_info = mut_info)
    }

    if (is.null(mutation_data)) {
        return(RPF_count_info |> mutate(reads_sequence_alt = reads_sequence, mutation_info = ""))
    }

    mutation_data <- mutation_data |>
        mutate(mutation_meta = mutation) |>
        separate_wider_delim(mutation, delim = ":", names = c("ref", "pos", "alt")) |>
        mutate(pos = as.numeric(pos))

    RPF_count_info |>
        mutate(
            alt_reads = pmap(
                list(reads_start_pos, reads_stop_pos, reads_sequence),
                single_read_alt,
                mutation_data = mutation_data
            )
        ) |>
        unnest(cols = alt_reads)
}


selected_transcripts <- selected_transcripts |>
    mutate(
        ORF_tstart = as.numeric(ORF_tstart) + 1,
        ORF_tstop  = as.numeric(ORF_tstop)
    )

mutation_info <- selected_transcripts |>
    mutate(
        indel = map(transcript_seq, add_indel),
        snp   = map(transcript_seq, add_snp)
    ) |>
    unnest(cols = indel) |>
    select(-reads_sequence) |>
    unnest(cols = snp) |>
    unite("mutation", c("indel_info", "snp_info"), sep = ",") |>
    separate_rows(mutation, sep = ",") |>
    select(transcript_id, mutation) |>
    distinct() |>
    filter(mutation != "") |>
    mutate(mutation_pos = str_extract(mutation, "[0-9]+")) |>
    select(transcript_id, mutation_pos, mutation)

mut_tsv_path   <- file.path(out_dir, paste0(output_prefix, "_mut.tsv"))
mut_block_path <- file.path(out_dir, paste0(output_prefix, "_mut_block.tsv"))

write_tsv(mutation_info, mut_tsv_path, col_names = FALSE)

# Run python script: Transcript coordinates to genomic coordinates
gppy_command <- paste(
    "python", gppy_script, "t2g",
    "-g", path_gtf,
    "-i", mut_tsv_path,
    ">",  mut_block_path
)
result <- try(system(gppy_command, wait = TRUE, intern = TRUE), silent = TRUE)
if (inherits(result, "try-error") || !file.exists(mut_block_path) || file.size(mut_block_path) == 0) {
    stop("gppy execution failed. Please verify that the installation package is intact and that the output is not empty.\n",
         "Command: ", gppy_command, call. = FALSE)
}

message("[Step 5 of 7] Coordinate transformation complete: ", mut_block_path)

mutation_info_gpos <- read_tsv(
    mut_block_path,
    col_names      = c("transcript_id", "mutation_pos", "mutation", "chrom", "strand", "genomic_position"),
    show_col_types = FALSE
)

mutation_intergrate <- function(df) {
    if (nrow(df) == 1) return(select(df, transcript_id, mutation))

    selected_mutation <- sample(df$mutation, 1)
    alt_para          <- str_split_1(selected_mutation, ":")[3]

    df |>
        separate_wider_delim(mutation, delim = ":", names = c("ref", "mutation_pos_2", "alt")) |>
        select(-mutation_pos_2) |>
        mutate(alt = alt_para) |>
        unite("mutation", c(ref, mutation_pos, alt), sep = ":")
}

mutation_info_gpos_nest <- mutation_info_gpos |>
    group_nest(chrom, strand, genomic_position) |>
    mutate(
        data = future_map(data, mutation_intergrate, .options = furrr_options(seed = TRUE))
    ) |>
    unnest(cols = data)


message("[Step 6 of 7] Generate read sequences (mutations + sequencing errors)...")

mutation_info_gpos_nest <- mutation_info_gpos_nest |>
    select(transcript_id, mutation) |>
    distinct() |>
    group_nest(transcript_id) |>
    right_join(selected_transcripts, by = "transcript_id") |>
    mutate(
        RPF_count_info = future_pmap(
            list(ORF_tstart, ORF_tstop, transcript_seq, RPF_count_info),
            get_raw_reads,
            .options = furrr_options(seed = TRUE)
        )
    ) |>
    select(-transcript_seq)

final <- mutation_info_gpos_nest |>
    mutate(
        RPF_count_info = future_pmap(
            list(RPF_count_info, data),
            get_alt_reads,
            .options = furrr_options(seed = TRUE)
        )
    ) |>
    select(-data) |>
    unnest(cols = RPF_count_info) |>
    mutate(
        reads_sequence_alt    = future_map(reads_sequence_alt, add_sequencing_error, .options = furrr_options(seed = TRUE)),
        reads_sequence_final  = map_chr(reads_sequence_alt, ~ .x[[1]]),
        sequencing_error_info = map_chr(reads_sequence_alt, ~ .x[[2]])
    ) |>
    select(-reads_sequence_alt)


tsv_path <- file.path(out_dir, paste0(output_prefix, ".tsv"))
write_tsv(final, tsv_path, col_names = TRUE, quote = "none")
message("[Step 7 of 7] Output TSV: ", tsv_path)

final_fa <- final |>
    select(
        transcript_id, ORF_tstart, ORF_tstop,
        vec_d5, vec_d3, periodic_tag, frame_shift,
        codon_frame0_pos, reads_start_pos, reads_stop_pos,
        reads_length, reads_length_d5, reads_length_d3,
        reads_sequence_final
    ) |>
    mutate(row_id = paste0("num", row_number())) |>
    unite(
        "header",
        c(transcript_id, ORF_tstart, ORF_tstop, vec_d5, vec_d3,
          periodic_tag, frame_shift, codon_frame0_pos,
          reads_start_pos, reads_stop_pos, reads_length,
          reads_length_d5, reads_length_d3, row_id),
        sep = ":"
    ) |>
    mutate(
        phred  = str_length(reads_sequence_final),
        tag    = "+",
        header = map_chr(header, ~ str_c("@", .x))
    ) |>
    group_nest(phred) |>
    mutate(phred = map_chr(phred, ~ str_c(rep("?", .x), collapse = ""))) |>
    unnest(data)

write_func <- function(header, reads_sequence_final, tag, phred) {
    list(header, reads_sequence_final, tag, phred)
}

fq_file <- future_pmap(final_fa, write_func, .options = furrr_options(seed = TRUE)) |> unlist()

fq_path <- file.path(out_dir, paste0(output_prefix, ".fq"))
write_lines(fq_file, fq_path)
message("[Step 7 of 7] Output FASTQ: ", fq_path)

message(" Simulation complete! ")

