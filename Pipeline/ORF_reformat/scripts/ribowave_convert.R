.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(furrr)
plan(multisession, workers = snakemake@threads)
options("future.globals.maxSize" = 60 * 1024 * 1024^2)
options(scipen = 999)

output_func <- function(df, output_path) {
    readr::write_tsv(df, output_path, col_names = F, quote = "none")
}

gencode_transcript_seq <- read_tsv(
    "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.transcripts.seq",
    col_names = F,
    col_select = c(
        "transcript_id" = "X1",
        "transcript_sequence" = "X2"
    )
)

## contains no stop codon ,do not need to modify ORF_stop

#prepare_gppy <- possibly(function(raw_file, gppy) {
#    tmp_file <- read_tsv(raw_file, col_names = T)
#    tmp_file <- tmp_file |>
##        dplyr::inner_join(
#            gencode_transcript_seq,
#            by = c("transcript_id" = "transcript_id")
#        ) |>
#        dplyr::mutate(
#            ORF_sequence = stringr::str_sub(transcript_sequence, ORF_tstart, ORF_tstop),
#            ORF_tstart_1base = ORF_tstart,
#            ORF_tstop_1base = ORF_tstop,
#            ORF_sequence_correct = stringr::str_sub(ORF_sequence),
#            start_codon = stringr::str_sub(ORF_sequence_correct, 1, 3)
#        ) |>
#        dplyr::select(transcript_id, ORF_tstart_1base, ORF_tstop_1base, ORF_sequence_correct, start_codon) |>
#        dplyr::distinct()

#    output_func(tmp_file, gppy)
#})

# 建议在调试阶段去掉 possibly，或者像这样保留但确保逻辑正确
prepare_gppy <- possibly(function(raw_file, gppy) {
    
    # 1. 读取文件：关闭表头自动识别 (col_names = F)
    tmp_file <- read_tsv(raw_file, col_names = FALSE, col_types = cols(.default = "c"))
    
    # 2. 关键修复：拆分列
    # 假设格式为: TranscriptID_Frame_Start_Stop
    # 使用 separate 将 X1 拆分为需要的列
    tmp_file <- tmp_file |>
        tidyr::separate(
            col = X1, 
            into = c("transcript_id", "frame", "ORF_tstart", "ORF_tstop"), 
            sep = "_", 
            remove = FALSE, # 保留原始列以防万一
            convert = TRUE  # 自动将数字字符串转换为数值型
        ) |>
        # 确保 Start/Stop 是数字，防止后续计算报错
        dplyr::mutate(
            ORF_tstart = as.numeric(ORF_tstart),
            ORF_tstop = as.numeric(ORF_tstop)
        )

    # 3. 后续逻辑保持不变
    tmp_file <- tmp_file |>
        dplyr::inner_join(
            gencode_transcript_seq,
            by = c("transcript_id" = "transcript_id")
        ) |>
        dplyr::mutate(
            # 注意：R中substr索引从1开始。如果你的mx文件是0-based坐标，这里可能需要调整 (+1)
            # 通常 Ribowave 输出可能是 0-based 或 1-based，请根据具体工具文档确认。
            # 假设输入是 1-based，直接使用：
            ORF_sequence = stringr::str_sub(transcript_sequence, ORF_tstart, ORF_tstop),
            ORF_tstart_1base = ORF_tstart,
            ORF_tstop_1base = ORF_tstop,
            ORF_sequence_correct = stringr::str_sub(ORF_sequence),
            start_codon = stringr::str_sub(ORF_sequence_correct, 1, 3)
        ) |>
        dplyr::select(transcript_id, ORF_tstart_1base, ORF_tstop_1base, ORF_sequence_correct, start_codon) |>
        dplyr::distinct()

    output_func(tmp_file, gppy)
}, otherwise = NULL)

prepare_gppy(
    raw_file = snakemake@input[[1]],
    gppy = snakemake@output[[1]]
)
