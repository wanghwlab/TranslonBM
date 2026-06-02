#! /usr/bin/env R
.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(furrr)
options(scipen = 999)
options("future.globals.maxSize" = 20 * 1024 * 1024^2)
plan(multisession, workers = 16)



trim_mismatch_minus <- function(fq_minus, trimmed_fq_minus) {
    ## trim mismatch and filter reads length below 25nt
    strand_minus <- read_lines(fq_minus) |>
        matrix(byrow = T, ncol = 4) |>
        as_tibble(.name_repair) |>
        rename(
            header_raw = V1,
            sequence = V2,
            info = V3,
            quality = V4
        ) |>
        separate_wider_delim(header_raw, names = c("header", "MD_tag"), delim = "\t") |>
        mutate(
            five_trim_num = str_count(MD_tag, "[ATCG]0$") + str_count(MD_tag, "[ATCG]0[ATCG]0$"), ## minus strand
            five_secondary_flag = ifelse(str_sub(sequence, start = 2,end = 2) == "N",-1,0)
        ) |>
        mutate(
            five_trim_num = five_trim_num + five_secondary_flag
        )

    ## max mismatch is set to 2
    ## trim mismatch 
    ## note: do not to trim N

    strand_minus_trimmed <- strand_minus |>
        group_by(five_trim_num) |>
        mutate(
            trimmed_sequence = str_sub(sequence, start = five_trim_num + 1),
            trimmed_quality = str_sub(quality, start = five_trim_num + 1)
        ) |>
        ungroup() |>
        select(
            header, trimmed_sequence, MD_tag, trimmed_quality
        ) |>
        mutate(
            MD_tag = paste0("+", MD_tag)
        )

    write_func <- function(header, trimmed_sequence, MD_tag, trimmed_quality) {
        tmp <- list(header, trimmed_sequence, MD_tag, trimmed_quality)
        return(tmp)
    }

    final_minus <- strand_minus_trimmed |>
        furrr::future_pmap(write_func) |>
        unlist()

    write_lines(
        final_minus,
        file = trimmed_fq_minus
    )
}

#trim_mismatch_minus(
#    fq_minus = snakemake@input[[1]],
#    trimmed_fq_minus = snakemake@output[[1]]
#)


# 获取所有minus.fq文件
input_files <- list.files("trim_five_prime_mismatch", pattern = "*minus.fq", full.names = TRUE)
for(input_file in input_files) {
  # 从文件名提取merge和mpsf
  filename <- basename(input_file)
  parts <- strsplit(sub("_minus.fq", "", filename), "_")[[1]]
  merge <- paste(parts[1:(length(parts)-1)], collapse = "_")
  mpsf <- parts[length(parts)]

  # 构建输出文件名
  output_file <- file.path("trim_five_prime_mismatch", paste0(merge, "_", mpsf, "_minus_trimmed.fq"))

  # 调用处理函数
  trim_mismatch_minus(
    fq_minus = input_file,
    trimmed_fq_minus = output_file
  )
}

