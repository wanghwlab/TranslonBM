#!/usr/bin/env bash
reads=$(awk -F '\t' '{print $1}' "${snakemake_input[offsets]}" | tr '\n' ',' | sed 's/,$//')
offsets=$(awk -F '\t' '{print $2}' "${snakemake_input[offsets]}" | tr '\n' ',' | sed 's/,$//')

(cd ${snakemake_params[out_prefix]} && \
    Ribotaper.sh "${snakemake_input[ribo_bam]}" \
    "${snakemake_input[rna_bam]}" \
    ${snakemake_params[ribotaper_index]} \
    $reads \
    $offsets \
    ${snakemake[threads]})