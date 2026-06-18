.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(RColorBrewer)
library(cowplot)

# ==============================================================================
# 1. 设置文件路径
# ==============================================================================
rld_file <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_simu/merged_RLD_stats_final.csv"
offset_file <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_simu/merged_offsets_all.csv"
output_pdf <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_simu/rld_offsets_chrN_simu_untrim.pdf"

# ==============================================================================
# 2. 读取并处理 RLD 数据
# ==============================================================================
rld_data <- read_csv(rld_file, show_col_types = FALSE) %>%
  rename(mapping_software = ORFtools, samples = sample) %>% 
  group_by(samples, mapping_software) %>%
  mutate(
    fraction = read_count / sum(read_count)
  ) %>%
  ungroup() %>%
  filter(between(read_length, 25, 35))

# ==============================================================================
# 3. 读取并处理 Offset 数据 
# ==============================================================================
offset_data_raw <- read_csv(offset_file, show_col_types = FALSE)

colnames(offset_data_raw) <- c("sample_combined", "psite_software", "read_length", "offset")

offset_data <- offset_data_raw %>%
  extract(sample_combined, into = c("samples", "mapping_software"), regex = "^(.*)_([^_]+)$") %>%
  filter(between(read_length, 25, 35))

print("Data Split Check:")
print(head(offset_data))

# ==============================================================================
# 4. 合并数据
# ==============================================================================
full_data <- full_join(
  rld_data,
  offset_data,
  by = c("samples", "mapping_software", "read_length")
) %>%
  mutate(
    fraction = replace_na(fraction, 0),
    pipeline = paste(psite_software, mapping_software, sep = ":"),
    paper_info = gsub("simulation_", "", samples)
  ) %>%
  filter(!is.na(psite_software)) 

# ==============================================================================
# 5. 设置因子水平 
# ==============================================================================
# 1. Mapping 顺序
full_data <- full_data %>%
  mutate(mapping_software = factor(mapping_software, levels = c("STAR", "hisat2", "tophat2")))

# 2. P-site Software 顺序
psite_levels <- c(
  "plastid", "shoelaces", "ribowaltz", "riboseqc", 
  "price", "ribocode", "ribotish", "ribotricer", 
  "ribotaper", "ribowave", "riborf", "rpbp", 
  "orfquant", "orfrater", "ribohmm"
)

full_data <- full_data %>%
  mutate(psite_software = factor(psite_software, levels = psite_levels)) %>%
  arrange(psite_software, mapping_software) %>%
  mutate(pipeline = factor(pipeline, levels = unique(pipeline) %>% rev()))

# ==============================================================================
# 6. 绘图函数 
# ==============================================================================
plot_rld_offset <- function(df, aes_x, aes_y, aes_fill, text_para, label_para) {
  
  df_complete <- df %>%
    group_by(!!sym(label_para), pipeline) %>%
    complete(read_length = 25:35, fill = list(fraction = 0)) %>% 
    ungroup()
  
  ggplot(df_complete, aes(x = {{ aes_x }}, y = {{ aes_y }}, fill = {{ aes_fill }})) +
    geom_tile(color = NA, linewidth = 0.2) + 
    geom_text(aes(label = {{ text_para }}), color = "black", size = 3, na.rm = TRUE) +
    facet_wrap(as.formula(paste0("~", label_para)), ncol = 4, scale = "free") +
    labs(
      x = "Length distribution of uniquely mapped reads",
      y = NULL, 
      fill = "Fraction of reads"
    ) +

    scale_fill_gradient(low = "#d9d9d9", high = "#bfbbd9", labels = scales::label_percent()) +
    scale_x_continuous(breaks = seq(25, 35, 1), expand = c(0, 0)) +
    theme_bw() + 
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      axis.title.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.5, colour = "black"),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.key.width = unit(30, "pt"),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      strip.background = element_blank(),
      strip.text = element_text(size = 11, color = "black", face = "bold") # 移除 family
    )
}

# ==============================================================================
# 7. 生成并保存
# ==============================================================================

p <- plot_rld_offset(
  df = full_data,
  aes_x = read_length,
  aes_y = pipeline,
  aes_fill = fraction,
  text_para = offset,
  label_para = "paper_info"
)

n_pipelines <- length(unique(full_data$pipeline))
calc_height <- max(10, n_pipelines * 0.1) 

ggsave(
  filename = output_pdf,
  plot = p,
  width = 12, 
  height = calc_height,
  limitsize = FALSE
)

print(paste("Done! Saved to:", output_pdf))
