#!/usr/bin/env bash

samtools view -@ ${snakemake[threads]} -bh -F 16 \
    -e '[MD] =~ "^0[ATCG].*"' \
    -o "${snakemake_output[mismatch_bam]}" "${snakemake_input[raw_bam]}" && \
samtools view -@ ${snakemake[threads]} -bh -F 16 \
    -e '[MD] !~ "^0[ATCG].*"' \
    -o "${snakemake_output[perfect_bam_unsort]}" "${snakemake_input[raw_bam]}" && \
samtools sort -@ ${snakemake[threads]} \
    -o "${snakemake_output[perfect_bam]}" "${snakemake_output[perfect_bam_unsort]}" && \
samtools index -o "${snakemake_output[perfect_bam_index]}" "${snakemake_output[perfect_bam]}"
