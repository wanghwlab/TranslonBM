# TranslonBM
 A unified framework for benchmarking translon detection

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
2. **Pipeline** — 
3. **Evaluation** — 

---

## Repository Structure

```
TranslonBM/
├── RiboSim/                        # Ribo-seq RPF data simulator
│   ├── Simulation.R                # Main simulation script
│   └── Require Input Data/         # Required input files (see details below)
│
├── Pipeline/                       # 
│
├── Evaluation/                     # 
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

Install R packages:

```r
install.packages(c("argparse", "furrr", "tidyverse"))
```

Install RiboSim (recommended via conda):

```bash
conda create -n RiboSim python=3.8
conda activate RiboSim
pip install gppy
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

> 🚧 Documentation in progress.

---

### Evaluation

> 🚧 Documentation in progress.

---

## Citation

If you use TranslonBM in your research, please cite:

> *Citation details will be provided upon publication.*

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
