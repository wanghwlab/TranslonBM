.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
options(scipen = 999)

# --- 加载必要的包 ---
library(tidyverse)
# [新增] 加载 Biostrings 用于读取 FASTA 文件
if (!require("Biostrings", quietly = TRUE)) {
    stop("请确保环境中已安装 Biostrings 包: conda install bioconductor-biostrings")
}
library(Biostrings)

# --- 修改开始: 正确读取 FASTA 文件 ---
# 原代码使用 read_tsv 读取 .fa 会报错，改为使用 readDNAStringSet
fa_path <- "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_1core/riborf/riborf_chrN/index/Human/candidateORF.fa"

# 读取 FASTA
fa_data <- readDNAStringSet(fa_path)

# 将 FASTA 对象转换为 Tibble (DataFrame)
# names(fa_data) 是 ID，as.character(fa_data) 是序列
riborf_ORF_sequence <- tibble(
    ORF_ID = names(fa_data) |> stringr::str_split_i(" ", 1), # 提取第一个空格前的 ID
    ORF_sequence = as.character(fa_data)
)
# --- 修改结束 ---


# --- 以下代码保持原样 (仅 index 读取部分微调以确保安全) ---

riborf_index <- read_tsv(
    "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_1core/riborf/riborf_chrN/index/Human/candidateORF.genepred.txt",
    col_names = F,
    col_select = c(
        "ORF_ID" = "X1",
        "ORF_gstart" = "X6",
        "ORF_gstop" = "X7"
    )
)

riborf_metainfo <- riborf_ORF_sequence |>
    inner_join(riborf_index, by = c("ORF_ID" = "ORF_ID"))

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = T, quote = "none")
}

tid_gid_looktable <- read_tsv(
    "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.looktable",
    col_names = F
) |>
    dplyr::rename(
        "transcript_id" = "X1",
        "gene_id" = "X2"
    )

# 建议去掉 possibly 以便调试，或者保持原样。这里为了配合你的习惯保留 possibly
riborf_extract_ORF <- possibly(function(raw_file, withseq) {
    tmp <- read_tsv(raw_file, col_names = T) |>
        filter(pred.pvalue > 0.7 & length >= 27) |>
        mutate(ORF_ID_tmp = orfID, ORF_ID_raw = orfID) |>
        separate(ORF_ID_tmp, c("tid_tmp", "RankNumber", "transcript_tmp", "ORF_type", "start_codon"), sep = "\\|") |>
        separate(tid_tmp, c("transcript_id", "chrom_tmp", "strand_tmp"), sep = ":") |>
        separate(transcript_tmp, c("transcript_length", "ORF_tstart", "ORF_tstop")) |>
        inner_join(riborf_metainfo, by = c("orfID" = "ORF_ID")) |>
        dplyr::rename(
            "ORF_ID" = "orfID",
            "ORF_length" = "length"
        ) |>
        inner_join(tid_gid_looktable, by = c("transcript_id" = "transcript_id"))

    output_func(tmp, withseq)
}, otherwise = NULL)

## stand_result:: ORF_ID_raw gene_id transcript_id chrom strand ORF_gstart ORF_gstop ORF_tstart ORF_tstop ORF_length start_codon ORF_sequence ORF_sequence_aa annotated_gstart annotated_gstop annotated_tstart annotated_tstop

riborf_extract_ORF(
    raw_file = snakemake@input[[1]],
    withseq = snakemake@output[[1]]
)
