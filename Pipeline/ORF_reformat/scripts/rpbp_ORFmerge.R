.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(furrr)
plan(multisession, workers = snakemake@threads)
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

library(Biostrings)

preprocessing <- possibly(function(raw_file, ORF_fa) {
    
    # --- 修改开始: 正确读取 FASTA 文件 ---
    # 读取 FASTA 为 DNAStringSet 对象
    fa_obj <- readDNAStringSet(ORF_fa)
    
    # 转换为 Tibble 用于后续 join
    # 注意：names(fa_obj) 通常包含 header 信息，需要确保它能和 raw_file 中的 id 对应
    # 这里假设 raw_file 中的 id 与 fasta header 的第一部分(空格前)一致
    ORF_sequence <- tibble(
        id = names(fa_obj) |> str_split_i(" ", 1), # 提取空格前的 ID，视具体 header 格式调整
        ORF_sequence = as.character(fa_obj)
    )
    # --- 修改结束 ---

    ORF_blocks <- read_tsv(
        raw_file,
        col_names = T,
        col_types = cols(.default = "c")
    ) |>
        dplyr::rename_with(~ stringr::str_sub(.x, 2, -1), everything()) |>
        # 这里原来的 inner_join 逻辑保持不变
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
    
    return(ORF_blocks) # 显式返回，是个好习惯
}, otherwise = NULL) # 显式写出 otherwise

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


target_dir <- snakemake@params[["fa_dir"]]

# 【修改点】不再匹配样本名前缀，直接找该目录下任何以 .fa 结尾的文件
# ignore.case = TRUE 可以同时匹配 .fa 和 .FA
found_fa <- list.files(target_dir, pattern = "\\.fa$", full.names = TRUE, ignore.case = TRUE)[1]

# 打印一下路径方便调试 (如果还报错，查看日志就能知道它到底在哪找文件)
message("正在目录查找 FA 文件: ", target_dir)
message("找到文件: ", found_fa)

if (is.na(found_fa)) {
    # 如果找不到，列出该目录下有哪些文件，方便排查
    existing_files <- list.files(target_dir)
    stop(paste0("Error: 在目录 ", target_dir, " 中未找到 .fa 文件。\n目录下的文件有: ", paste(existing_files, collapse = ", ")))
}

ORF_processed <- preprocessing(
    raw_file = snakemake@input[[1]],
    #ORF_fa = snakemake@input[[2]]
    ORF_fa = found_fa
)

format_ORF(
    ORF_blocks = ORF_processed,
    formatted_ORF = snakemake@output[[1]]
)

merge_ORF(
    ORF_blocks = ORF_processed,
    merged_ORF = snakemake@output[[2]]
)
