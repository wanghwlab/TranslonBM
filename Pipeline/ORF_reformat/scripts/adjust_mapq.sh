#!/usr/bin/env bash

bedtools bamtobed -i "${snakemake_input[raw_bam]}" \
    -bed12 -split -cigar >"${snakemake_output[bed]}" && \
bedtools bedtobam -i "${snakemake_output[bed]}" \
    -mapq 51 -g "${snakemake_input[genome_file]}" -bed12 -ubam >"${snakemake_output[adjusted_bam]}"
