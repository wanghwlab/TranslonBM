.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(cowplot)


# --- 数据读取与预处理 ---
setwd('/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/time_benchmark_cores')

#csv_file_path <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/time_benchmark/version2026.1/ORF_time_benchmark_summary_10.23.csv"
#csv_file_path <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/time_benchmark_cores/ORF_time_benchmark_summary_12.29_merge.csv"
csv_file_path <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/time_benchmark_cores/ORF_time_benchmark_summary_12.29_11tools.csv"

raw_data <- read_csv(csv_file_path)

#  --- 1. 在原始数据上创建 basename 列 --- 
data_with_basename <- raw_data %>%
    mutate(
        basename = case_when(
            TRUE ~ basename(doc_path)
        )
    )

time_cols <- c("s", "cpu_time")
space_cols <- c("max_rss", "max_vms", "max_uss", "max_pss", "io_in", "io_out", "mean_load")


main_aligners <- c("STAR", "hisat2", "tophat2")

indexing_steps <- data_with_basename %>%
    filter(!aligner %in% main_aligners)

analysis_steps <- data_with_basename %>%
    filter(aligner %in% main_aligners)

expansion_grid <- analysis_steps %>%
    distinct(sample, aligner)

expanded_indexing_steps <- indexing_steps %>%
    select(-sample, -aligner) %>%
    crossing(expansion_grid)

combined_for_final_agg <- bind_rows(analysis_steps, expanded_indexing_steps)


# --- 2. 软件聚合 ---
final_data <- combined_for_final_agg %>%
  group_by(sample, aligner, software, chr, parameter) %>%
  summarise(
    across(all_of(time_cols), sum, na.rm = TRUE),
    across(all_of(space_cols), max, na.rm = TRUE),
    .groups = 'drop'
  )

# --- 3. 添加绘图所需的新分类和标签列 ---
study_name_map <- c(
    "SRX876063_SRX876069" = "Ji et al. (2015)", "SRX740748" = "Gao et al. (2015)",
    "SRX1254413" = "Calviello et al. (2016)", "SRX1447296" = "Raj et al. (2016)",
    "SRX5256543_SRX5256555" = "Martinez et al. (2020)", "SRX5887328_SRX5887329_SRX5887330" = "Chen et al. (2020)",
    "SRX11812007_SRX11812008_SRX11812009" = "Chothani et al. (2022)",
    "SRX7666669-73"="SRX7666669-73", "SRX7666674-78"="SRX7666674-78",
    "SRX7666679-83"="SRX7666679-83", "SRX7666684-88"="SRX7666684-88",
    "SRX7666689-93"="SRX7666689-93", "SRX7666694-98"="SRX7666694-98"
)
sample_order <- c(
    "Ji et al. (2015)", "Gao et al. (2015)", "Calviello et al. (2016)", "Raj et al. (2016)",
    "Martinez et al. (2020)", "Chen et al. (2020)", "Chothani et al. (2022)",
    "SRX7666669-73", "SRX7666674-78", "SRX7666679-83", "SRX7666684-88",
    "SRX7666689-93", "SRX7666694-98"
)

plot_ready_data <- final_data %>%
    mutate(
        dataset_type = case_when(
            startsWith(sample, "simulation") ~ "Simulation",
            startsWith(sample, "SRX7666") ~ "Replicate",
            startsWith(sample, "SRX") ~ "Real",
            TRUE ~ "Other"
        ),
        elapsed_time = s / 60,
        max_rss_GiB = max_rss / 1024,
        `uss/rss` = if_else(max_rss > 0, (max_uss / max_rss) * 100, 0),
        cpu_load_ave = mean_load,
        study_name = coalesce(study_name_map[sample], sample)
    ) %>%
    mutate(study_name = factor(study_name, levels = unique(c(sample_order, sample)))) %>%
    mutate(aligner = factor(aligner, levels = c("STAR", "hisat2", "tophat2")))

# --- 4. 定义自定义颜色 ---
software_colors <- c(
  "riboseqc" = "#FAE17D", "gedi" = "#93C681", "ribocode" = "#58AC52", "ribotish" = "#E28792",
  "ribotricer" = "#4E8BC6", "ribotaper" = "#84CFAB", "ribowave" = "#836AB8", "riborf" = "#C2A2B8",
  "rpbp" = "#79A690", "ORFquant" = "#6E8486", "orfrater" = "#6F5E50", "ribohmm" = "#EB8730"
)

# --- 5. 循环绘图与保存 ---
print("Step 2: Generating plots...")

dataset_types <- unique(plot_ready_data$dataset_type)

for (current_type in dataset_types) {

    type_data <- plot_ready_data %>% filter(dataset_type == current_type)
    if (nrow(type_data) == 0) next
    
    csv_filename <- paste0("ORFdetector_resources_", current_type, ".csv")
    data_to_save <- type_data %>%
        select(study_name, sample, aligner, software, max_rss_GiB, elapsed_time, `uss/rss`, cpu_load_ave, max_rss, s, cpu_time)
    write_csv(data_to_save, csv_filename)
    print(paste("Saved data table to:", csv_filename))
    
    fig_memory_usage <- type_data %>%
        ggplot(aes(x = study_name, y = max_rss_GiB, color = software, shape = aligner)) +
        geom_jitter(width = 0.2, alpha = 0.8, stroke = 0, size = 2) +
        scale_color_manual(values = software_colors) +
        scale_shape_manual(values = c("STAR" = 15, "hisat2" = 16, "tophat2" = 17)) +
        scale_y_continuous(breaks = function(x) seq(0, max(x, na.rm = TRUE) + 20, by = 200)) + 
        labs(
            x = NULL, 
            y = "Maximum resident set size (GiB)",
            color = "ORF Detectors",
            shape = "Aligner"
        ) +
        theme_classic() +
        theme(
            text = element_text(color = "black"), axis.text = element_text(color = "black"),
            axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
            legend.position = "right", axis.title = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    fig_cpu_load <- type_data %>%
        ggplot(aes(x = study_name, y = elapsed_time, color = software, shape = aligner)) +
        geom_jitter(width = 0.2, alpha = 0.8, stroke = 0, size = 2) +
        scale_color_manual(values = software_colors) +
        scale_shape_manual(values = c("STAR" = 15, "hisat2" = 16, "tophat2" = 17)) +
        scale_y_continuous(breaks = function(x) seq(0, max(x, na.rm = TRUE) + 1500, by = 1500)) +
        labs(
            x = NULL, 
            y = "Elapsed real time (mins)",
            color = "ORF Detectors",
            shape = "Aligner"
        ) +
        theme_classic() +
        theme(
            text = element_text(color = "black"), axis.text = element_text(color = "black"),
            axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
            legend.position = "right", axis.title = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    final_plot <- cowplot::plot_grid(
        fig_memory_usage, fig_cpu_load,
        ncol = 2, 
        align = "hv", 
        labels = "AUTO" 
    )

    output_filename <- paste0("ORFdetector_resources_real_merged_", current_type, ".pdf")

    ggsave(
        plot = final_plot,
        filename = output_filename,
        width = 12, 
        height = 3, 
        units = "in",
        limitsize = FALSE
    )

    print(paste("Saved plot to:", output_filename))
}
