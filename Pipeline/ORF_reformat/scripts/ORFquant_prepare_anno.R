.libPaths(c('/home/tangyuewen/R/4.3_user/','/home/tangyuewen/R/4.3','/home/tangyuewen/miniconda3/envs/r_deseq2/lib/R/library'))
#! /usr/bin/env R
library(ORFquant)
# ============================================================================ #
# prepare annotation files

prepare_annotation_files(
    annotation_directory = snakemake@params[[1]],
    twobit_file = snakemake@input[[2]],
    gtf_file = snakemake@input[[1]],
    scientific_name = snakemake@params[[2]],
    annotation_name = snakemake@params[[3]],
    export_bed_tables_TxDb = TRUE,
    forge_BSgenome = FALSE,
    genome_seq = snakemake@input[[3]],
    circ_chroms = DEFAULT_CIRC_SEQS,
    create_TxDb = TRUE
)