.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(ggpubr)
library(extrafont)
library(colorspace)

# install.packages("extrafont")
library(extrafont)
font_import()
loadfonts(device = "pdf")

cowplot::set_null_device("pdf")


loadfonts(device = "pdf")
cowplot::set_null_device("pdf")

# --- 1. 设置文件路径 ---
setwd("/home/tangyuewen/ORF_benchmark/final_ORFs_2025.7/plots/Other_plots/")
workdir <- "/home/tangyuewen/ORF_benchmark/simulation_preprocessing"
#gcoor_dir <- "/home/tangyuewen/ORF_benchmark/final_ORFs_2025.7/final_ORFs/merged/orf_pred_default_trim"
gcoor_dir <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/orf_pred_default_untrim"


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

samples <- c("simulation_6M_T1","simulation_6M_T3","simulation_60M_T1","simulation_60M_T3")
aligner <- c("tophat2","hisat2","STAR")
orf_detector <- c("ribocode", "ribotish", "riborf", "ribowave", "rpbp", "ribotricer", "ORFquant", "orfrater", "gedi")


get_gcoor_data <- function(df,samples,aligner,orf_detector){
    tmp <- df |> 
        filter(
            samples == {{samples}},
            aligner == {{aligner}},
            orf_detector == {{orf_detector}}
        ) |> 
        pull(gcoor_data) |> 
        pluck(1)
    return(tmp)
}

parse_gcoor <- function(samples,aligner,orf_detector){
    file_path <- paste0(gcoor_dir,"/",samples,"_",aligner,"_",orf_detector,"_gcoor.tsv.gz")
    if (!file.exists(file_path)) {
        return(tibble(gene_id = character(), coordinate_id = character(), start_codon = character()))
    }
    tmp <- read_tsv(file_path,
        col_names = T,
        col_select = c(
            gene_id, coordinate_id, start_codon
        )) |> 
        distinct()
    return(tmp)
}

parse_ref <- function(samples){
    file_path <- paste0(workdir,"/",samples,"_block.tsv")
    tmp <- read_tsv(
        file_path, col_names = F, col_types = stringr::str_c(rep("c", 15), collapse = ""),
        col_select = c(
            chrom = "X1", ORF_gstart_0base = "X2", ORF_gstop_0base = "X3", gene_id = "X5",
            strand = "X6", block_length = "X11", block_start = "X12", transcript_id = "X13",
            ORF_tstart_1base = "X14", ORF_tstop_1base = "X15"
        )) |> 
        dplyr::mutate(
            block_length = stringr::str_sub(block_length, end = -2L),
            block_start = stringr::str_sub(block_start, end = -2L),
            coordinate_0base = purrr::pmap_chr(list(ORF_gstart_0base, block_length, block_start), get_coordinate)
            ) |>
        tidyr::unite("coordinate_id", c(chrom, strand, coordinate_0base), sep = ":", remove = FALSE)
    return(tmp)
}

parse_unique_alignments <- function(samples,aligner){
    file_path <- paste0(workdir,"/",samples,"_",aligner,".header")
    tmp <- read_tsv(file_path,col_names = "reads_id") |> 
        separate_wider_delim(
            "reads_id",
            names = c("transcript_id", "ORF_tstart", "ORF_tstop", "vec_d5", "vec_d3", "periodic_tag", "frame_shift", 
                      "codon_frame0_pos", "reads_start_pos", "reads_stop_pos", "reads_length", "reads_length_d5", 
                      "reads_length_d3", "id"),
            delim = ":"
        )
    return(tmp)
}

get_ORF_inframe_ratio <- function(df){
    tmp <- df |> 
        group_by(transcript_id, ORF_tstart, ORF_tstop, periodic_tag, frame_shift) |>
        summarise(count = n(), .groups = 'drop') |>
        group_by(transcript_id, ORF_tstart, ORF_tstop) |>
        mutate(ORF_reads_count = sum(count)) |>
        ungroup() |>
        filter(frame_shift == 0) |>
        group_by(transcript_id, ORF_tstart, ORF_tstop) |>
        mutate(ORF_reads_count_inframe = sum(count)) |>
        distinct() |>
        ungroup() |>
        mutate(inframe_ratio = ORF_reads_count_inframe / ORF_reads_count) |>
        select(-c(periodic_tag,count)) |>
        distinct()
    return(tmp)
}

get_final_ref_block <- function(ref,ORF_inframe_ratio){
    tmp <- ORF_inframe_ratio |> 
        inner_join(
            ref,
            by = c("transcript_id" = "transcript_id", "ORF_tstart" = "ORF_tstart_1base", "ORF_tstop" = "ORF_tstop_1base")
        ) |> 
        arrange(desc(inframe_ratio)) |> 
        mutate(
            inframe_ratio_bin = case_when(
            inframe_ratio <= 0.35 ~ "<=35%", inframe_ratio > 0.35 & inframe_ratio <= 0.40 ~ "35~40%",
            inframe_ratio > 0.4 & inframe_ratio <= 0.45 ~ "40~45%", inframe_ratio > 0.45 & inframe_ratio <= 0.50 ~ "45~50%",
            inframe_ratio > 0.5 & inframe_ratio <= 0.55 ~ "50~55%", inframe_ratio > 0.55 & inframe_ratio <= 0.60 ~ "55~60%",
            inframe_ratio > 0.6 & inframe_ratio <= 0.65 ~ "60~65%", inframe_ratio > 0.65 & inframe_ratio <= 0.70 ~ "65~70%",
            inframe_ratio > 0.7 & inframe_ratio <= 0.75 ~ "70~75%", inframe_ratio > 0.75 & inframe_ratio <= 0.80 ~ "75~80%",
            inframe_ratio > 0.8 & inframe_ratio <= 0.85 ~ "80~85%", inframe_ratio > 0.85 & inframe_ratio <= 0.90 ~ "85~90%",
            inframe_ratio > 0.9 & inframe_ratio <= 0.95 ~ "90~95%", inframe_ratio > 0.95 & inframe_ratio <= 1 ~ "95~100%",
    )) 
    return(tmp)
}

get_valid_ratio <- function(df1,df2){
    tmp1 <- df1 |> pull(coordinate_id) |> unique()
    tmp2 <- df2 |> pull(coordinate_id) |> unique()
    if (length(tmp1) == 0) return(0)
    valid_tmp <- intersect(tmp1,tmp2) |> unique() |> length()
    valid_ratio <- valid_tmp / length(tmp1)
    return(valid_ratio)
}

get_unvalid_ratio <- function(df1,df2){
    tmp1 <- df1 |> pull(coordinate_id) |> unique()
    tmp2 <- df2 |> pull(coordinate_id) |> unique()
    if (length(tmp2) == 0) return(0)
    unvalid_tmp <- setdiff(tmp2,tmp1) |> unique() |> length()
    unvalid_ratio <- unvalid_tmp / length(tmp2)
    return(unvalid_ratio)
}

main_func <- function(samples){
    ref_block <- parse_ref(samples)
    
    aligner_func <- function(aligner){
        mapped_ORF <- parse_unique_alignments(samples,aligner) |> 
            get_ORF_inframe_ratio()
        mapped_ORF_ref <- get_final_ref_block(ref_block,mapped_ORF) |> 
            mutate(samples = {{samples}}, aligner = {{aligner}})
        
        gcoor_df <- crossing(samples,aligner,orf_detector) |> 
            mutate(gcoor_data = pmap(pick(everything()), parse_gcoor))

        get_valid_ratio_wapper <- function(samples,aligner,inframe_ratio_bin,ref_data,orf_detector, codon_group){
            df1 <- ref_data
            df2_all <- get_gcoor_data(df = gcoor_df, samples = {{samples}}, aligner = {{aligner}}, orf_detector = {{orf_detector}})
            
            df2 <- df2_all |> 
                mutate(codon_group_eval = ifelse(start_codon == "ATG", "ATG", "NTG")) |>
                filter(codon_group_eval == {{codon_group}})

            tmp <- get_valid_ratio(df1,df2)
            return(tmp)
        }

        get_unvalid_ratio_wapper <- function(samples,aligner,ref_data,orf_detector, codon_group){
            df1 <- ref_data
            df2_all <- get_gcoor_data(df = gcoor_df, samples = {{samples}}, aligner = {{aligner}}, orf_detector = {{orf_detector}})

            df2 <- df2_all |> 
                mutate(codon_group_eval = ifelse(start_codon == "ATG", "ATG", "NTG")) |>
                filter(codon_group_eval == {{codon_group}})
            
            tmp <- get_unvalid_ratio(df1,df2)
            return(tmp)
        }
        
        mapped_ORF_valid <- mapped_ORF_ref |> 
            group_nest(samples,aligner,inframe_ratio_bin,.key = "ref_data") |> 
            crossing(orf_detector) |>
            crossing(codon_group = c("ATG", "NTG")) |> 
            mutate(block_valid = pmap_dbl(pick(everything()), get_valid_ratio_wapper)) |> 
            mutate(ref_ORF_nums = ref_data |> map_dbl(\(x) nrow(x))) |> 
            select(-"ref_data") |> 
            distinct()

        mapped_ORF_unvalid <- mapped_ORF_ref |> 
            group_nest(samples,aligner,.key = "ref_data") |> 
            crossing(orf_detector) |>
            crossing(codon_group = c("ATG", "NTG")) |> 
            mutate(block_unvalid = pmap_dbl(pick(everything()), get_unvalid_ratio_wapper)) |>
            select(-"ref_data") |> 
            distinct()

        fig_data <- inner_join(mapped_ORF_valid, mapped_ORF_unvalid, 
                               by=c("samples", "aligner", "orf_detector", "codon_group"))
        return(fig_data)
    }
    aligner |> map(aligner_func) |> list_rbind()
}

fig_data <- samples |> map(main_func) |> list_rbind()

orf_detector_color <- c(
  "Consensus offset" = "#dc575a", "Plastid" = "#9fc7da", "Shoelaces" = "#f1ba7e", "riboWaltz" = "#c29e86",
  "Ribo-seQC" = "#73aaa4", "PRICE" = "#add698", "RiboCode" = "#51a54f", "Ribo-TISH" = "#ea8a87",
  "Ribotricer" = "#3c80af", "RiboTaper" = "#73aaa4", "RiboWave" = "#73539e", "RibORF" = "#d1a3a3",
  "RP-BP" = "#5f9b8c", "ORFquant" = "#73aaa4", "ORF-RATER" = "#5a4b50", "riboHMM" = "#f78a41"
)

get_valid_fig <- function(samples, aligner, simulation_data){
    get_ORF_num_plot <- function(simulation_data){
        tmp <- ggplot() + 
        geom_col(data = simulation_data |> select(inframe_ratio_bin,ref_ORF_nums) |> distinct(),
            mapping = aes(x = inframe_ratio_bin,y = ref_ORF_nums), fill = "#595959") + 
        scale_y_continuous(labels = scales::label_comma(), position = "right") +
        labs(title = paste0(samples, ":", aligner), x = "Frame percentage") +
        coord_cartesian(ylim = c(0, 5000), clip = "off") + 
        theme_tufte(base_family = "Arial") +
        theme(
            axis.text = element_text(size = 11,color = "black"), 
            axis.text.x = element_text(angle = 45,hjust = 1),
            axis.title.y = element_blank(), 
            plot.title = element_text(hjust = 0.5),
            axis.line.y.right = element_line(color = "black", size = 0.5) 
        )
        return(tmp)
    }
    
    ORF_num_plot <- get_ORF_num_plot(simulation_data)

    valid_plot <- ggplot( 
        data = simulation_data |> select(inframe_ratio_bin,block_valid,orf_detector) |> distinct(),
        aes(x = inframe_ratio_bin,y = block_valid,group = orf_detector)) + 
        geom_point(aes(color = orf_detector)) +
        geom_line(aes(group = orf_detector,color = orf_detector),linetype = 2) + 
        coord_cartesian(clip = "off") + 
        scale_y_continuous(labels = scales::label_percent(), limits=c(0,1)) + 
        scale_color_manual(values = orf_detector_color)+
        theme_tufte(base_family = "Arial")+
        labs(x = "Frame percentage", y = "Block level", color = "ORF Detectors") + 
        theme(
            axis.text = element_text(size = 11,color = "black"), 
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.text = element_text(family = "Arial",size = 11.5),
            axis.line = element_line(color = "black", size = 0.5) 
        )
    
    tmp <- align_plots(ORF_num_plot, valid_plot + theme(legend.position = "none"), align = "hv", axis = "tblr")
    ggdraw(tmp[[1]]) + draw_plot(tmp[[2]])
}


fig_data <- fig_data |>
    mutate(
        orf_detector = case_match(
            orf_detector, "plastid" ~ "Plastid", "riboseqc" ~ "Ribo-seQC", "ribowaltz" ~ "riboWaltz",
            "shoelaces" ~ "Shoelaces", "gedi" ~ "PRICE", "ribocode" ~ "RiboCode", "ribotish" ~ "Ribo-TISH",
            "ribowave" ~ "RiboWave", "rpbp" ~ "RP-BP", "ribotricer" ~ "Ribotricer", "ribohmm" ~ "riboHMM",
            "riborf" ~ "RibORF", "ribotaper" ~ "RiboTaper", "ORFquant" ~ "ORFquant", "orfrater" ~ "ORF-RATER"
        ),
        orf_detector = factor(orf_detector, levels = c("Plastid", "Shoelaces", "riboWaltz", "Ribo-seQC", "PRICE", "RiboCode", "Ribo-TISH", "Ribotricer", "RiboTaper", "RiboWave", "RibORF", "RP-BP", "ORFquant", "ORF-RATER", "riboHMM")),
        aligner = case_match(aligner, "STAR" ~ "STAR", "hisat2" ~ "HISAT2", "tophat2" ~ "TopHat2"),
        aligner = factor(aligner, levels = c("STAR", "HISAT2", "TopHat2")),
        samples = factor(samples, levels = c("simulation_6M_T1","simulation_6M_T3","simulation_60M_T1","simulation_60M_T3"))
    )

get_gcoor_num <- function(samples,aligner,orf_detector, codon_group){
    file_path <- paste0(gcoor_dir,"/",samples,"_",aligner,"_",orf_detector,"_gcoor.tsv.gz")
    if (!file.exists(file_path)) {
        return(tibble(coordinate_id_n = 0))
    }
    tmp <- read_tsv(file_path, col_names = T, col_select = c(coordinate_id, start_codon)) |>
        mutate(codon_group_eval = ifelse(start_codon == "ATG", "ATG", "NTG")) |>
        filter(codon_group_eval == {{codon_group}}) |>
        summarise(coordinate_id_n = n_distinct(coordinate_id))
    return(tmp)
}


df <- crossing(samples, aligner, orf_detector, codon_group = c("ATG", "NTG"))
df_nums <- df |> 
    mutate(counts = pmap(df, get_gcoor_num)) |>
    unnest(cols = counts)


ORF_unvalid_data <- fig_data |>
    select(samples, aligner, orf_detector, codon_group, block_unvalid) |>
    distinct() |>
    inner_join(df_nums, by=c("samples", "aligner", "orf_detector", "codon_group"))

legend_fig_data <- tibble(tmp1 = names(orf_detector_color), tmp2 = orf_detector_color)
legend_fig <- legend_fig_data |> 
    filter(tmp1 %in% c("PRICE","ORFquant","ORF-RATER","RiboCode","RibORF","Ribo-TISH","Ribotricer","RiboWave","RP-BP")) |> 
    mutate(tmp1 = factor(tmp1, levels = c("Plastid", "Shoelaces", "riboWaltz", "Ribo-seQC", "PRICE", "RiboCode", "Ribo-TISH", "Ribotricer", "RiboTaper", "RiboWave", "RibORF", "RP-BP", "ORFquant", "ORF-RATER", "riboHMM"))) |> 
    ggplot() + geom_point(aes(x = tmp1,y = tmp2,color = tmp1)) + 
    scale_color_manual(values = orf_detector_color) + labs(color = "ORF Detectors") + theme_bw() +
    guides(color = guide_legend(nrow = 2)) + theme(legend.position = "bottom",legend.text = element_text(size = 10))
legend_plot <- cowplot::get_legend(legend_fig)


list(codon_group = c("ATG", "NTG")) |>
  pwalk(\(codon_group) {
    
    # -- 1. 验证率图 --
    fig_data_subset <- fig_data |> filter(codon_group == !!codon_group)
    
    fig_data_df <- fig_data_subset |> 
        group_nest(samples,aligner,.key = "simulation_data") |>
        arrange(samples, aligner) |>
        mutate(fig = pmap(list(samples, aligner, simulation_data), get_valid_fig))
        
    block_valid_plots <- fig_data_df |> pull(fig)
    tmp <- cowplot::plot_grid(plotlist = block_valid_plots, ncol = 3, nrow = 4)
    block_valid_fig <- cowplot::plot_grid(tmp, legend_plot, ncol = 1, rel_heights = c(1, 0.1))
    
    validation_filename <- file.path(workdir, paste0("block_level_validation_", codon_group, ".pdf"))
    ggsave(filename = validation_filename, plot = block_valid_fig, width = 12, height = 16, limitsize = FALSE, device = cairo_pdf)
    print(paste("Saved:", validation_filename))
    
    # -- 2. 非验证率图 --
    unvalid_data_subset <- ORF_unvalid_data |> filter(codon_group == !!codon_group)
    
    block_unvalid_fig <- unvalid_data_subset |> 
        mutate(orf_detector = forcats::fct_rev(orf_detector)) |> 
        ggplot() +
        geom_point(aes(x = aligner, y = orf_detector, color = block_unvalid, size = coordinate_id_n)) + 
        facet_wrap(~ samples) + 
        scale_size(name = "Detected ORFs (Block)", labels = scales::label_number(scale = 1, big.mark = ",", suffix = "")) +
        scale_color_continuous_sequential(
            name = "Unvalid Ratio (Block)", limit = c(0,1), labels = scales::label_percent(), palette = "Viridis", rev = FALSE) +
        labs(title = paste("Block Level Unvalidation -", codon_group)) + 
        theme_bw() +
        theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5))
        
    unvalidation_filename <- file.path(workdir, paste0("block_level_unvalidation_", codon_group, ".pdf"))
    ggsave(filename = unvalidation_filename, plot = block_unvalid_fig, width = 12, height = 10, limitsize = FALSE)
    print(paste("Saved:", unvalidation_filename))
  })
  
  
  
fig_data_ATG <- fig_data |> filter(codon_group == "ATG")
fig_data_NTG <- fig_data |> filter(codon_group == "NTG")

df_ATG <- crossing(samples, aligner, orf_detector) |> mutate(codon_group = "ATG")
df_nums_ATG <- df_ATG |> 
    mutate(counts = pmap(df_ATG, get_gcoor_num)) |>
    unnest(cols = counts)
ORF_unvalid_data_ATG <- fig_data_ATG |>
    select(samples, aligner, orf_detector, codon_group, block_unvalid) |>
    distinct() |>
    inner_join(df_nums_ATG, by=c("samples", "aligner", "orf_detector", "codon_group"))

df_NTG <- crossing(samples, aligner, orf_detector) |> mutate(codon_group = "NTG")
df_nums_NTG <- df_NTG |> 
    mutate(counts = pmap(df_NTG, get_gcoor_num)) |>
    unnest(cols = counts)
ORF_unvalid_data_NTG <- fig_data_NTG |>
    select(samples, aligner, orf_detector, codon_group, block_unvalid) |>
    distinct() |>
    inner_join(df_nums_NTG, by=c("samples", "aligner", "orf_detector", "codon_group"))
	
	
# ===================================================
# ===========   为 ATG 组绘图并保存   ==============
# ===================================================

# -- 1. 验证率图 (ATG) --
fig_data_df_ATG <- fig_data_ATG |> 
    group_nest(samples, aligner, .key = "simulation_data") |>
    arrange(samples, aligner) |>
    mutate(fig = pmap(list(samples, aligner, simulation_data), get_valid_fig))
    
block_valid_plots_ATG <- fig_data_df_ATG |> pull(fig)
tmp_ATG <- cowplot::plot_grid(plotlist = block_valid_plots_ATG, ncol = 3, nrow = 4)
block_valid_fig_ATG <- cowplot::plot_grid(tmp_ATG, legend_plot, ncol = 1, rel_heights = c(1, 0.1))

validation_filename_ATG <- file.path(workdir, "block_level_validation_ATG.pdf")
ggsave(filename = validation_filename_ATG, plot = block_valid_fig_ATG, width = 12, height = 16, limitsize = FALSE, device = cairo_pdf)
print(paste("Saved:", validation_filename_ATG))

# -- 2. 非验证率图 (ATG) --
block_unvalid_fig_ATG <- ORF_unvalid_data_ATG |> 
    mutate(orf_detector = forcats::fct_rev(orf_detector)) |> 
    ggplot() +
    geom_point(aes(x = aligner, y = orf_detector, color = block_unvalid, size = coordinate_id_n)) + 
    facet_wrap(~ samples) + 
    scale_size(name = "Detected ORFs (Block)", labels = scales::label_number(scale = 1, big.mark = ",", suffix = "")) +
    scale_color_continuous_sequential(
        name = "Unvalid Ratio (Block)", limit = c(0,1), labels = scales::label_percent(), palette = "Viridis", rev = FALSE) +
    labs(title = "Block Level Unvalidation - ATG") + 
    theme_bw() +
    theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5))
    
unvalidation_filename_ATG <- file.path(workdir, "block_level_unvalidation_ATG.pdf")
ggsave(filename = unvalidation_filename_ATG, plot = block_unvalid_fig_ATG, width = 12, height = 10, limitsize = FALSE)
print(paste("Saved:", unvalidation_filename_ATG))


# ===================================================
# ===========   为 NTG 组绘图并保存   ==============
# ===================================================

# -- 1. 验证率图 (NTG) --
fig_data_df_NTG <- fig_data_NTG |> 
    group_nest(samples, aligner, .key = "simulation_data") |>
    arrange(samples, aligner) |>
    mutate(fig = pmap(list(samples, aligner, simulation_data), get_valid_fig))
    
block_valid_plots_NTG <- fig_data_df_NTG |> pull(fig)
tmp_NTG <- cowplot::plot_grid(plotlist = block_valid_plots_NTG, ncol = 3, nrow = 4)
block_valid_fig_NTG <- cowplot::plot_grid(tmp_NTG, legend_plot, ncol = 1, rel_heights = c(1, 0.1))

validation_filename_NTG <- file.path(workdir, "block_level_validation_NTG.pdf")
ggsave(filename = validation_filename_NTG, plot = block_valid_fig_NTG, width = 12, height = 16, limitsize = FALSE, device = cairo_pdf)
print(paste("Saved:", validation_filename_NTG))

# -- 2. 非验证率图 (NTG) --
block_unvalid_fig_NTG <- ORF_unvalid_data_NTG |> 
    mutate(orf_detector = forcats::fct_rev(orf_detector)) |> 
    ggplot() +
    geom_point(aes(x = aligner, y = orf_detector, color = block_unvalid, size = coordinate_id_n)) + 
    facet_wrap(~ samples) + 
    scale_size(name = "Detected ORFs (Block)", labels = scales::label_number(scale = 1, big.mark = ",", suffix = "")) +
    scale_color_continuous_sequential(
        name = "Unvalid Ratio (Block)", limit = c(0,1), labels = scales::label_percent(), palette = "Viridis", rev = FALSE) +
    labs(title = "Block Level Unvalidation - NTG") + 
    theme_bw() +
    theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5))
    
unvalidation_filename_NTG <- file.path(workdir, "block_level_unvalidation_NTG.pdf")
ggsave(filename = unvalidation_filename_NTG, plot = block_unvalid_fig_NTG, width = 12, height = 10, limitsize = FALSE)
print(paste("Saved:", unvalidation_filename_NTG))
