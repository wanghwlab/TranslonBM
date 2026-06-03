# TranslonBM

**A unified framework for benchmarking translon detection**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3?logo=r)
![Python](https://img.shields.io/badge/Python-%3E%3D3.6-3776AB?logo=python)

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Modules](#modules)
  - [RiboSim](#ribosim)
  - [Pipeline](#pipeline)
  - [Evaluation](#evaluation)
- [Citation](#citation)

---

## Overview

TranslonBM is a unified benchmarking framework for evaluating computational tools that detect translons — actively translated open reading frames (ORFs) identified from Ribo-seq data. The framework covers three stages of a complete benchmarking workflow:

1. **RiboSim** — A Flexible Ribo-seq Data Simulator
2. **Pipeline** — Ribo-seq data preprocessing, alignment, ORF prediction, and output standardization
3. **Evaluation** — Statistical analysis, benchmark evaluation, and figure plotting

---

## Repository Structure

```
TranslonBM/
├── RiboSim/                              # Ribo-seq RPF data simulator
│   ├── Simulation.R                      # Main simulation script
│   └── Require Input Data/               # Required input files
│
├── Pipeline/                             # Ribo-seq data preprocessing and ORF prediction
│   ├── ORF_detect_Riboseq/               # ORF prediction on real Ribo-seq datasets
│   │   └── scripts/
│   ├── ORF_detect_Synthetic/             # ORF prediction on simulation datasets
│   │   └── scripts/
│   ├── ORF_reformat/                     # Standardize ORF prediction output format
│   │   └── scripts/
│   └── trim5prime_scripts/               # 5' end mismatch trimming and remapping
│       └── scripts/
│
├── Evaluation/                           # Benchmark analysis and figure plotting scripts
│   ├── Fig1/
│   ├── Fig2/
│   ├── Fig3/
│   └── Fig4/
│
└── README.md
```

---

## Modules

### RiboSim

RiboSim simulates ribosome-protected fragment (RPF) sequencing data from Ribo-seq experiments. Starting from real transcript sequences and ORF annotations, it generates synthetic FASTQ files that closely mimic real data, including:

- Empirical read length distribution (5'/3' offset model fitted from real data)
- 3-nucleotide periodicity signal with controllable noise level
- SNP, indel, and sequencing error simulation
- Full ground truth annotation in the output TSV

The output FASTQ encodes all ground truth information in the read header, enabling precise benchmarking of any downstream Ribo-seq analysis tool.

#### Requirements

Install RiboSim (recommended via conda):

```bash
conda create -n RiboSim python=3.8
conda activate RiboSim
pip install gppy
```

Install R packages:

```r
install.packages(c("argparse", "furrr", "tidyverse"))
```

#### Required Input Files

All input files should be placed in a single directory and passed via `--Path_dir`. The following files are required:

| File | Description |
|---|---|
| `coding_transcripts.tsv` | List of coding transcript IDs |
| `coding_transcripts_start.tsv` | ORF start positions on transcripts |
| `coding_transcripts_stop.tsv` | ORF stop positions on transcripts |
| `gencode.v43.transcripts.seq` | Transcript sequences |
| `gencode.v43.annotation.gtf` | Genome annotation GTF file |
| `codons_freq.tsv` | Codon usage frequency table |
| `gtf.py` | gppy coordinate conversion script |

Example files are provided in [`RiboSim/Require Input Data/`](RiboSim/Require%20Input%20Data/).

#### Usage

```bash
Rscript Simulation.R --help
```

**Minimal example:**

```bash
Rscript Simulation.R \
    --counts    6000000 \
    --replicate T1 \
    --Path_dir  /path/to/input/data \
    --out_dir   /path/to/output/ \
    --threads   8
```

**Full parameter example:**

```bash
Rscript Simulation.R \
    --counts          6000000 \
    --replicate       T1 \
    --Path_dir        /path/to/input/data \
    --threads         8 \
    --seed            42 \
    --out_dir         /path/to/output \
    --n_transcripts   12000 \
    --perturb_cutoff  0.02 \
    --snp_ratio       0.001 \
    --indel_ratio     0.0005 \
    --seq_error_ratio 0.005
```

#### Parameters

**Required:**

| Parameter | Type | Description |
|---|---|---|
| `--counts` | INT | Total number of simulated reads (e.g. `6000000`, `60000000`) |
| `--replicate` | Label | Replicate label, e.g. `T1`, `T2`, `T3` |
| `--Path_dir` | PATH | Directory containing all required input files |
| `--threads` | INT | Number of parallel workers |

**Optional:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--seed` | INT | `42` | Global random seed for reproducibility |
| `--n_transcripts` | INT | `12000` | Number of transcripts to sample |
| `--orfs_num` | INT | `12000` | Number of ORFs to simulate |
| `--perturb_cutoff` | FLOAT | `0.02` | Maximum proportion of non-periodic reads |
| `--snp_ratio` | FLOAT | `0.001` | Per-base SNP probability |
| `--indel_ratio` | FLOAT | `0.0005` | Per-base indel probability |
| `--insertion_ratio` | FLOAT | `0.00025` | Insertion probability within indels |
| `--seq_error_ratio` | FLOAT | `0.005` | Sequencing error probability |
| `--max_iter` | INT | `1000` | Max iterations in ORF count sampling |
| `--gppy_script` | PATH | `<Path_dir>/gtf.py` | Path to gppy gtf.py script |
| `--out_dir` | PATH | `<Path_dir>` | Output directory |

#### Output Files

All output files are prefixed with `simulation_<counts>M_<replicate>` (e.g. `simulation_6M_T1`):

| File | Description |
|---|---|
| `simulation_*_detail.tsv` | Per-ORF read assignment details |
| `simulation_*_mut.tsv` | Simulated mutations in transcript coordinates |
| `simulation_*_mut_block.tsv` | Simulated mutations in genomic coordinates |
| `simulation_*.tsv` | Full ground truth annotation table for all simulated reads |
| `simulation_*.fq` | Final simulated FASTQ file |


---

### Pipeline

The Pipeline module contains Snakemake workflows for Ribo-seq data preprocessing, alignment, ORF prediction, and output standardization.

#### ORF_detect_Riboseq

Snakemake workflows for running 11 ORF prediction tools on 6 real Ribo-seq datasets. All tools are run with default parameters and unified to use NTG start codons (ATG or tool default if NTG is not configurable).

- `scripts/` — Helper scripts called by the Snakemake workflows

#### ORF_detect_Synthetic

Snakemake workflows for running 9 ORF prediction tools on 4 simulation datasets. All tools are run with default parameters and unified to use NTG start codons (ATG or tool default if NTG is not configurable).

- `scripts/` — Helper scripts called by the Snakemake workflows

#### ORF_reformat

Standardizes the raw output of 11 ORF prediction tools into a unified genome-based exon block format. The reformatted results include gene information, genomic coordinates, start codon, nucleotide sequence, and amino acid sequence, and are used for downstream benchmark analysis.

- `scripts/` — Helper scripts called by the Snakemake workflows

#### trim5prime_scripts

Removes mismatched bases at the first 2 nt of the 5′ end of Ribo-seq reads, then remaps the trimmed reads using three alignment tools (STAR, HISAT2, TopHat2).

- `scripts/` — Helper scripts called by the Snakemake workflows

#### Environment Setup

```bash
conda activate snakemake_env

snakemake -s ./Pipeline/ORF_detect_Riboseq/snakemake_ribocode.smk -j 6 --keep-going -p --use-conda
```

#### Docker Images

Pre-built Docker images are provided for reproducibility. To use them:

```bash
# Import image from local file
gunzip -c translon-detector-aligner_v1.tar.gz | docker load

# Run image
docker run -it --rm translon-detector-aligner:v1 /bin/bash
```

| Image | Contents |
|:---|:---|
| `translon-detector-aligner:v1` | STAR 2.7.10b, HISAT2, TopHat2 |
| `translon-detector-ribohmm:v1` | RiboHMM |
| `translon-detector-orfquant:v1` | RiboseQC, ORFquant |
| `translon-detector-suite:v1` | GEDI, ORFrater, RiboCode, RiboTaper, ribotricer, Rp-Bp, RibORF, Ribo-TISH, RiboWave |

All images include `snakemake_env` for running Snakemake workflows.

---

### Evaluation

The Evaluation module contains statistical analysis and plotting scripts for the main benchmark results (Figures 1–4 in the manuscript).

| Directory | Contents |
|---|---|
| `Fig1/` | Data quality and preprocessing evaluation scripts |
| `Fig2/` | ORF prediction overlap analysis scripts |
| `Fig3/` | Benchmark scoring scripts (TIS, TISeq, MS) |
| `Fig4/` | Tool combination analysis and runtime benchmark scripts |

---

## Citation

If you use TranslonBM in your research, please cite:

> *Citation details will be provided upon publication.*

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
