.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(colorspace)

options(readr.show_col_types = FALSE)

base_font_family <- "sans"

setwd("/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/Other_plots/")
workdir <- "/home/tangyuewen/ORF_benchmark/simulation_preprocessing"
gcoor_dir <- "/home/tangyuewen/ORF_benchmark/final_ORFs_2025.7/final_ORFs/merged/orf_pred_default_trim"

out_dir <- "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/Other_plots/Sensitivity_Analysis_Normalized_lineplot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

use_cache <- FALSE
cache_file <- file.path(out_dir, "Sensitivity_Analysis_Normalized_processed_data.rds")

samples_vec <- c("simulation_6M_T1", "simulation_6M_T3", "simulation_60M_T1", "simulation_60M_T3")
aligner_vec <- c("tophat2", "hisat2", "STAR")
orf_detector_vec <- c("ribocode", "ribotish", "riborf", "ribowave", "rpbp", "ribotricer", "ORFquant", "orfrater", "gedi")
codon_group_vec <- c("ATG", "NTG")

gene_group_file <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/featureCount/gene_expression_groups.csv"
gene_groups_raw <- read_csv(gene_group_file, show_col_types = FALSE)

get_coordinate <- function(ORF_gstart_0base, block_length, block_start) {
    ORF_gstart_0base <- as.numeric(ORF_gstart_0base)
    block_length <- block_length |> str_split_1(",") |> map_int(\(x) as.numeric(x))
    ORF_exon_start <- block_start |> str_split_1(",") |> map_int(\(x) as.numeric(x) + ORF_gstart_0base)
    ORF_exon_stop <- ORF_exon_start + block_length
    stringr::str_c(ORF_exon_start, ORF_exon_stop, sep = "-", collapse = ",")
}

parse_ref <- function(sample_id) {
    file_path <- paste0(workdir, "/", sample_id, "_block.tsv")
    read_tsv(
        file_path,
        col_names = FALSE,
        col_types = stringr::str_c(rep("c", 15), collapse = ""),
        col_select = c(
            chrom = "X1", ORF_gstart_0base = "X2", ORF_gstop_0base = "X3", gene_id = "X5",
            strand = "X6", block_length = "X11", block_start = "X12", transcript_id = "X13",
            ORF_tstart_1base = "X14", ORF_tstop_1base = "X15"
        )
    ) |>
        mutate(
            block_length = stringr::str_sub(block_length, end = -2L),
            block_start = stringr::str_sub(block_start, end = -2L),
            coordinate_0base = purrr::pmap_chr(list(ORF_gstart_0base, block_length, block_start), get_coordinate)
        ) |>
        unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE)
}

parse_unique_alignments <- function(sample_id, aligner_id) {
    file_path <- paste0(workdir, "/", sample_id, "_", aligner_id, ".header")
    read_tsv(file_path, col_names = "reads_id", show_col_types = FALSE) |>
        separate_wider_delim(
            reads_id,
            names = c(
                "transcript_id", "ORF_tstart", "ORF_tstop", "vec_d5", "vec_d3", "periodic_tag", "frame_shift",
                "codon_frame0_pos", "reads_start_pos", "reads_stop_pos", "reads_length", "reads_length_d5", "reads_length_d3", "id"
            ),
            delim = ":"
        )
}

get_ORF_inframe_ratio <- function(df) {
    df |>
        count(transcript_id, ORF_tstart, ORF_tstop, periodic_tag, frame_shift, name = "count") |>
        group_by(transcript_id, ORF_tstart, ORF_tstop) |>
        mutate(ORF_reads_count = sum(count)) |>
        ungroup() |>
        filter(frame_shift == 0) |>
        group_by(transcript_id, ORF_tstart, ORF_tstop) |>
        summarise(
            ORF_reads_count = first(ORF_reads_count),
            ORF_reads_count_inframe = sum(count),
            .groups = "drop"
        ) |>
        mutate(inframe_ratio = ORF_reads_count_inframe / ORF_reads_count)
}

parse_gcoor <- function(sample_id, aligner_id, detector_id) {
    file_path <- paste0(gcoor_dir, "/", sample_id, "_", aligner_id, "_", detector_id, "_gcoor.tsv.gz")
    if (!file.exists(file_path)) {
        return(tibble(gene_id = character(), coordinate_id = character(), start_codon = character()))
    }
    read_tsv(file_path, col_names = TRUE, col_select = c(gene_id, coordinate_id, start_codon), show_col_types = FALSE) |>
        distinct()
}

get_final_ref_block <- function(ref, ORF_inframe_ratio, current_gene_groups) {
    ORF_inframe_ratio |>
        inner_join(
            ref,
            by = c(
                "transcript_id" = "transcript_id",
                "ORF_tstart" = "ORF_tstart_1base",
                "ORF_tstop" = "ORF_tstop_1base"
            )
        ) |>
        left_join(current_gene_groups, by = "gene_id") |>
        filter(!is.na(gene_group)) |>
        mutate(
            frame_periodicity = case_when(
                inframe_ratio <= 0.35 ~ 0.35,
                inframe_ratio > 0.35 & inframe_ratio <= 0.40 ~ 0.40,
                inframe_ratio > 0.40 & inframe_ratio <= 0.45 ~ 0.45,
                inframe_ratio > 0.45 & inframe_ratio <= 0.50 ~ 0.50,
                inframe_ratio > 0.50 & inframe_ratio <= 0.55 ~ 0.55,
                inframe_ratio > 0.55 & inframe_ratio <= 0.60 ~ 0.60,
                inframe_ratio > 0.60 & inframe_ratio <= 0.65 ~ 0.65,
                inframe_ratio > 0.65 & inframe_ratio <= 0.70 ~ 0.70,
                inframe_ratio > 0.70 & inframe_ratio <= 0.75 ~ 0.75,
                inframe_ratio > 0.75 & inframe_ratio <= 0.80 ~ 0.80,
                inframe_ratio > 0.80 & inframe_ratio <= 0.85 ~ 0.85,
                inframe_ratio > 0.85 & inframe_ratio <= 0.90 ~ 0.90,
                inframe_ratio > 0.90 & inframe_ratio <= 0.95 ~ 0.95,
                inframe_ratio > 0.95 & inframe_ratio <= 1.00 ~ 1.00,
                TRUE ~ NA_real_
            ),
            inframe_ratio_bin = case_when(
                frame_periodicity == 0.35 ~ "<=35%",
                frame_periodicity == 0.40 ~ "35~40%",
                frame_periodicity == 0.45 ~ "40~45%",
                frame_periodicity == 0.50 ~ "45~50%",
                frame_periodicity == 0.55 ~ "50~55%",
                frame_periodicity == 0.60 ~ "55~60%",
                frame_periodicity == 0.65 ~ "60~65%",
                frame_periodicity == 0.70 ~ "65~70%",
                frame_periodicity == 0.75 ~ "70~75%",
                frame_periodicity == 0.80 ~ "75~80%",
                frame_periodicity == 0.85 ~ "80~85%",
                frame_periodicity == 0.90 ~ "85~90%",
                frame_periodicity == 0.95 ~ "90~95%",
                frame_periodicity == 1.00 ~ "95~100%",
                TRUE ~ NA_character_
            ),
            gene_group = factor(gene_group, levels = c("Low", "Mid", "High"))
        ) |>
        filter(!is.na(frame_periodicity), !is.na(inframe_ratio_bin))
}

valid_ratio_fast <- function(ref_coords, pred_coords) {
    if (length(ref_coords) == 0) return(0)
    length(intersect(ref_coords, pred_coords)) / length(ref_coords)
}

unvalid_ratio_fast <- function(ref_coords_all, pred_coords) {
    if (length(pred_coords) == 0) return(0)
    length(setdiff(pred_coords, ref_coords_all)) / length(pred_coords)
}

calculate_discrete_auc <- function(x, y) {
    ok <- !is.na(x) & !is.na(y)
    x <- x[ok]
    y <- y[ok]
    if (length(unique(x)) < 2) return(NA_real_)
    ord <- order(x)
    x_sorted <- x[ord]
    y_sorted <- y[ord]
    sum(diff(x_sorted) * (y_sorted[-1] + y_sorted[-length(y_sorted)]) / 2)
}

format_factor_levels <- function(df) {
    df |>
        mutate(
            orf_detector = case_match(
                orf_detector,
                "plastid" ~ "Plastid",
                "riboseqc" ~ "Ribo-seQC",
                "ribowaltz" ~ "riboWaltz",
                "shoelaces" ~ "Shoelaces",
                "gedi" ~ "PRICE",
                "ribocode" ~ "RiboCode",
                "ribotish" ~ "Ribo-TISH",
                "ribowave" ~ "RiboWave",
                "rpbp" ~ "RP-BP",
                "ribotricer" ~ "Ribotricer",
                "ribohmm" ~ "riboHMM",
                "riborf" ~ "RibORF",
                "ribotaper" ~ "RiboTaper",
                "ORFquant" ~ "ORFquant",
                "orfrater" ~ "ORF-RATER",
                .default = orf_detector
            ),
            orf_detector = factor(
                orf_detector,
                levels = c("Plastid", "Shoelaces", "riboWaltz", "Ribo-seQC", "PRICE", "RiboCode", "Ribo-TISH", "Ribotricer", "RiboTaper", "RiboWave", "RibORF", "RP-BP", "ORFquant", "ORF-RATER", "riboHMM")
            ),
            aligner = case_match(aligner, "STAR" ~ "STAR", "hisat2" ~ "HISAT2", "tophat2" ~ "TopHat2", .default = aligner),
            aligner = factor(aligner, levels = c("STAR", "HISAT2", "TopHat2")),
            samples = factor(samples, levels = samples_vec),
            codon_group = factor(codon_group, levels = codon_group_vec)
        )
}

process_one_condition <- function(sample_id, aligner_id) {
    message("Processing: ", sample_id, " | ", aligner_id)

    target_col <- paste0(sample_id, "_", aligner_id)
    if (target_col %in% colnames(gene_groups_raw)) {
        current_gene_groups <- gene_groups_raw |>
            select(gene_id = Geneid, gene_group = all_of(target_col))
    } else {
        warning("在 gene_expression_groups.csv 中未找到列名: ", target_col)
        current_gene_groups <- tibble(gene_id = character(), gene_group = character())
    }

    ref_block <- parse_ref(sample_id)
    mapped_ORF <- parse_unique_alignments(sample_id, aligner_id) |> get_ORF_inframe_ratio()
    mapped_ORF_ref <- get_final_ref_block(ref_block, mapped_ORF, current_gene_groups) |>
        mutate(samples = sample_id, aligner = aligner_id)

    gcoor_raw <- set_names(
        map(orf_detector_vec, \(detector_id) parse_gcoor(sample_id, aligner_id, detector_id)),
        orf_detector_vec
    )

    gcoor_sets <- crossing(orf_detector = orf_detector_vec, codon_group = codon_group_vec) |>
        mutate(
            pred_coords = map2(orf_detector, codon_group, \(detector_id, codon_id) {
                df <- gcoor_raw[[detector_id]]
                if (nrow(df) == 0) return(character(0))
                df |>
                    mutate(codon_group_eval = if_else(start_codon == "ATG", "ATG", "NTG")) |>
                    filter(codon_group_eval == codon_id) |>
                    pull(coordinate_id) |>
                    unique()
            })
        )

    all_ref_coords <- unique(mapped_ORF_ref$coordinate_id)

    ref_1D <- mapped_ORF_ref |>
        group_by(samples, aligner, inframe_ratio_bin, frame_periodicity) |>
        summarise(
            ref_coords = list(unique(coordinate_id)),
            ref_ORF_nums = n(),
            .groups = "drop"
        )

    fig_1D_valid <- tidyr::expand_grid(ref_1D, gcoor_sets) |>
        mutate(block_valid = map2_dbl(ref_coords, pred_coords, valid_ratio_fast)) |>
        select(samples, aligner, inframe_ratio_bin, frame_periodicity, ref_ORF_nums, orf_detector, codon_group, block_valid)

    unvalid_df <- gcoor_sets |>
        mutate(block_unvalid = map_dbl(pred_coords, \(x) unvalid_ratio_fast(all_ref_coords, x))) |>
        select(orf_detector, codon_group, block_unvalid)

    fig_1D <- fig_1D_valid |>
        left_join(unvalid_df, by = c("orf_detector", "codon_group"))

    ref_2D <- mapped_ORF_ref |>
        group_by(samples, aligner, frame_periodicity, gene_group) |>
        summarise(
            ref_coords = list(unique(coordinate_id)),
            .groups = "drop"
        )

    fig_2D <- tidyr::expand_grid(ref_2D, gcoor_sets) |>
        mutate(block_valid = map2_dbl(ref_coords, pred_coords, valid_ratio_fast)) |>
        select(samples, aligner, frame_periodicity, gene_group, orf_detector, codon_group, block_valid)

    list(fig_1D = fig_1D, fig_2D = fig_2D)
}

if (use_cache && file.exists(cache_file)) {
    message("Reading cached data: ", cache_file)
    cached <- readRDS(cache_file)
    fig_data <- cached$fig_data
    fig_data_2D <- cached$fig_data_2D
    auc_results <- cached$auc_results
    auc_results_2D <- cached$auc_results_2D
} else {
    condition_grid <- crossing(samples = samples_vec, aligner = aligner_vec)
    res_all <- pmap(condition_grid, \(samples, aligner) process_one_condition(samples, aligner))

    fig_data <- map_dfr(res_all, "fig_1D") |> format_factor_levels()
    fig_data_2D <- map_dfr(res_all, "fig_2D") |> format_factor_levels()

    auc_results <- fig_data |>
        group_by(samples, aligner, codon_group, orf_detector) |>
        summarise(
            Absolute_AUC = calculate_discrete_auc(frame_periodicity, block_valid),
            X_Range = max(frame_periodicity, na.rm = TRUE) - min(frame_periodicity, na.rm = TRUE),
            AUC_Sensitivity = if_else(X_Range > 0, Absolute_AUC / X_Range, NA_real_),
            .groups = "drop"
        )

    auc_results_2D <- fig_data_2D |>
        group_by(samples, aligner, codon_group, orf_detector, gene_group) |>
        summarise(
            Absolute_AUC = calculate_discrete_auc(frame_periodicity, block_valid),
            X_Range = max(frame_periodicity, na.rm = TRUE) - min(frame_periodicity, na.rm = TRUE),
            AUC_Sensitivity = if_else(X_Range > 0, Absolute_AUC / X_Range, NA_real_),
            .groups = "drop"
        ) |>
        mutate(
            gene_group = factor(gene_group, levels = c("Low", "Mid", "High")),
            AUC_label = paste0(as.character(gene_group), ": ", sprintf("%.2f", AUC_Sensitivity)),
            auc_y = case_when(
                gene_group == "High" ~ 0.96,
                gene_group == "Mid" ~ 0.88,
                gene_group == "Low" ~ 0.80,
                TRUE ~ 0.72
            )
        )

    saveRDS(
        list(
            fig_data = fig_data,
            fig_data_2D = fig_data_2D,
            auc_results = auc_results,
            auc_results_2D = auc_results_2D
        ),
        cache_file
    )
}


gene_group_color <- c(
    "High" = "#D55E00",
    "Mid"  = "#0072B2",
    "Low"  = "#009E73"
)

fig_data_2D <- fig_data_2D |>
    mutate(
        samples = as.character(samples),
        aligner = as.character(aligner),
        orf_detector = as.character(orf_detector),
        codon_group = as.character(codon_group),
        gene_group = factor(as.character(gene_group), levels = c("Low", "Mid", "High"))
    )

auc_gene_group_results <- fig_data_2D |>
    group_by(samples, aligner, codon_group, orf_detector, gene_group) |>
    summarise(
        Absolute_AUC = calculate_discrete_auc(frame_periodicity, block_valid),
        X_Range = max(frame_periodicity, na.rm = TRUE) - min(frame_periodicity, na.rm = TRUE),
        AUC_Sensitivity = Absolute_AUC / X_Range,
        .groups = "drop"
    ) |>
    mutate(
        AUC_label = paste0(as.character(gene_group), ": ", sprintf("%.2f", AUC_Sensitivity))
    )

get_auc_label_data <- function(plot_data, auc_data) {
    label_data <- auc_data |>
        group_by(samples, aligner, codon_group, orf_detector) |>
        summarise(
            AUC_text = paste(AUC_label, collapse = "\n"),
            .groups = "drop"
        ) |>
        left_join(
            plot_data |>
                group_by(samples, aligner, codon_group, orf_detector) |>
                summarise(
                    x_pos = min(frame_periodicity, na.rm = TRUE),
                    y_pos = 0.98,
                    .groups = "drop"
                ),
            by = c("samples", "aligner", "codon_group", "orf_detector")
        )

    return(label_data)
}

plot_gene_group_line <- function(sample_id, codon_id) {

    plot_data <- fig_data_2D |>
        filter(samples == sample_id, codon_group == codon_id)

    auc_label_data <- auc_gene_group_results |>
        filter(samples == sample_id, codon_group == codon_id) |>
        get_auc_label_data(plot_data = plot_data)

    p <- ggplot(
        plot_data,
        aes(
            x = frame_periodicity,
            y = block_valid,
            color = gene_group,
            group = gene_group
        )
    ) +
        geom_line(linewidth = 0.7) +
        geom_point(size = 1.5) +
        geom_text(
            data = auc_label_data,
            aes(x = x_pos, y = y_pos, label = AUC_text),
            inherit.aes = FALSE,
            hjust = 0,
            vjust = 1,
            size = 2.4,
            color = "black",
            lineheight = 0.9
        ) +
        facet_grid(aligner ~ orf_detector) +
        scale_color_manual(values = gene_group_color, name = "Host Gene Group") +
        scale_y_continuous(
            labels = scales::percent,
            limits = c(0, 1),
            expand = expansion(mult = c(0.02, 0.05))
        ) +
        scale_x_continuous(
            breaks = sort(unique(plot_data$frame_periodicity))
        ) +
        labs(
            title = paste("Sensitivity by Gene Expression:", sample_id, "(", codon_id, ")"),
            x = "Frame periodicity strength",
            y = "Sensitivity = TP / (TP + FN)"
        ) +
        theme_bw(base_family = base_font_family) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            axis.text.y = element_text(size = 7),
            axis.title = element_text(size = 10),
            strip.text.x = element_text(size = 7),
            strip.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            panel.grid = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 10),
            aspect.ratio = 0.4
        )

    return(p)
}

sample_levels <- c(
    "simulation_6M_T1",
    "simulation_6M_T3",
    "simulation_60M_T1",
    "simulation_60M_T3"
)

for (sample_id in sample_levels) {

    message("Plotting GeneGroup all codon groups: ", sample_id)

    pdf_file <- file.path(
        out_dir,
        paste0("Sensitivity_Heatmap_GeneGroup_", sample_id, ".pdf")
    )

    pdf(pdf_file, width = 18, height = 6, onefile = TRUE)

    for (codon_id in c("ATG", "NTG")) {
        p <- plot_gene_group_line(sample_id = sample_id, codon_id = codon_id)
        print(p)
    }

    dev.off()
}

for (sample_id in sample_levels) {

    message("Plotting GeneGroup ATG only: ", sample_id)

    pdf_file <- file.path(
        out_dir,
        paste0("Sensitivity_Heatmap_GeneGroup_ATG_", sample_id, ".pdf")
    )

    p <- plot_gene_group_line(sample_id = sample_id, codon_id = "ATG")

    ggsave(
        filename = pdf_file,
        plot = p,
        width = 18,
        height = 6,
        limitsize = FALSE,
        device = cairo_pdf
    )
}

write_csv(fig_data, file.path(out_dir, "Sensitivity_Curve_1D_Data.csv"))
write_csv(fig_data_2D, file.path(out_dir, "Sensitivity_Heatmap_GeneGroup_Data.csv"))
write_csv(auc_results, file.path(out_dir, "Sensitivity_AUC_Results.csv"))
write_csv(auc_results_2D, file.path(out_dir, "Sensitivity_AUC_GeneGroup_Results.csv"))

message("=== 完成：已输出 CSV、Sensitivity_Heatmap_GeneGroup.pdf 和 Sensitivity_Heatmap_GeneGroup_ATG.pdf。===")
