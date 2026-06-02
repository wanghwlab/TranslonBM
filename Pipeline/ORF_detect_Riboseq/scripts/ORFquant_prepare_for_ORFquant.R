#! /usr/bin/env R
library(ORFquant)
source('/home/tangyuewen/software/disjointExons.R')
elementLengths <- S4Vectors::elementNROWS
library('magrittr')
# ============================================================================ #
# prepare input file for ORFquant

prepare_for_ORFquant(
    annotation_file = snakemake@input[[1]],
    bam_file = snakemake@input[[2]],
    path_to_rl_cutoff_file = snakemake@input[[3]],
    chunk_size = 5e+06,
    dest_name = snakemake@params[[1]]
)
