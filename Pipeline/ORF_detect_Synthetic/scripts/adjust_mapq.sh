#!/usr/bin/env bash

# Exit immediately if any command fails
set -euo pipefail

# Convert BAM to BED and write to the temporary BED file declared by Snakemake
bedtools bamtobed -i "${snakemake_input_raw_bam}" \
    -bed12 -split -cigar >"${snakemake_output_bed}"

# Convert the temporary BED file back to BAM and set the MAPQ value
bedtools bedtobam -i "${snakemake_output_bed}" \
    -mapq 51 -g "${snakemake_input_genome_file}" -bed12 -ubam >"${snakemake_output_adjusted_bam}"
