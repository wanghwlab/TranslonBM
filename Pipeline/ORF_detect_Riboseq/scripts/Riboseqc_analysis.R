#! /usr/bin/env R
library("RiboseQC")

RiboseQC_analysis(
    annotation_file = snakemake@input[[2]],
    bam_files = snakemake@input[[1]],
    # genome_seq = snakemake@input[[3]], ## for env riboseqc we do not need it
    read_subset = T,
    readlength_choice_method = "max_coverage",
    chunk_size = 5000000L,
    write_tmp_files = TRUE,
    dest_names = snakemake@params[[1]],
    rescue_all_rls = FALSE,
    fast_mode = TRUE,
    create_report = TRUE,
    sample_names = NA,
    report_file = NA,
    extended_report = FALSE,
    pdf_plots = FALSE
)
