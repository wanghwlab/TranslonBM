#!/usr/bin/env bash

samtools view -@ ${snakemake[threads]} -bh \
    -e '[NH]==1 && rname =~ "chr[^M]+"' \
    -o "${snakemake_output[unique_bam_unsort]}" "${snakemake_input[raw_bam]}" && \
samtools sort -@ ${snakemake[threads]} -o "${snakemake_output[unique_bam]}" "${snakemake_output[unique_bam_unsort]}" && \
samtools index -@ ${snakemake[threads]} -o "${snakemake_output[unique_bam_index]}" "${snakemake_output[unique_bam]}"

## require alignments which are uniquely mapped and chrN-sourced
## there is still improvement
