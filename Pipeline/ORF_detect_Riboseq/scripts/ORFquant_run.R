#! /usr/bin/env R
library(ORFquant)
source('/home/tangyuewen/software/disjointExons.R')
elementLengths <- S4Vectors::elementNROWS
library('magrittr')

run_ORFquant(
    for_ORFquant_file = snakemake@input[[1]],
    annotation_file = snakemake@input[[2]],
    n_cores = snakemake@threads,
    prefix = snakemake@params[[1]],
    gene_name = NA,
    gene_id = NA,
    genomic_region = NA,
    write_temp_files = T,
    write_GTF_file = T,
    write_protein_fasta = T,
    interactive = T,
    stn.orf_find.all_starts = T,
    stn.orf_find.nostarts = F,
    stn.orf_find.start_sel_cutoff = NA,
    stn.orf_find.start_sel_cutoff_ave = 0.5,
    stn.orf_find.cutoff_fr_ave = 0.5,
    stn.orf_quant.cutoff_cums = NA,
    stn.orf_quant.cutoff_pct = 2,
    stn.orf_quant.cutoff_P_sites = NA,
    unique_reads_only = FALSE,
    canonical_start_only = FALSE
    #stn.orf_quant.scaling = "total_Psites"
)
