#!/usr/bin/env bash

grep ^RL "${snakemake_input[0]}" | cut -f 2 | tr ',' '\n' | sort -n | uniq | awk 'BEGIN{ORS=","} {print}' > "${snakemake_output[0]}"
