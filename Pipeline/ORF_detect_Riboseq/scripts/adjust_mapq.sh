#!/usr/bin/env bash

# 确保脚本在遇到任何错误时立即退出，这是一种良好的编程习惯
set -euo pipefail

# 将BAM转换为BED格式，输出到Snakemake现在已知的临时bed文件
~/miniconda3/bin/bedtools bamtobed -i "${snakemake_input_raw_bam}" \
    -bed12 -split -cigar >"${snakemake_output_bed}"

# 将临时BED文件转换回BAM格式，并设置MAPQ值
~/miniconda3/bin/bedtools bedtobam -i "${snakemake_output_bed}" \
    -mapq 51 -g "${snakemake_input_genome_file}" -bed12 -ubam >"${snakemake_output_adjusted_bam}"
