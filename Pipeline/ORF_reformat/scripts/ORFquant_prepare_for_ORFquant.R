.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
#! /usr/bin/env R
library(ORFquant)
# ============================================================================ #
# prepare input file for ORFquant

prepare_for_ORFquant(
    annotation_file = snakemake@input[[1]],
    bam_file = snakemake@input[[2]],
    path_to_rl_cutoff_file = snakemake@input[[3]],
    chunk_size = 5e+06,
    dest_name = snakemake@params[[1]]
)