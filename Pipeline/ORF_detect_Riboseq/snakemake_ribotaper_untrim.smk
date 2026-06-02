# snakefile_ribotaper.smk 
import os
import pandas as pd

#################################
# 1. 设置全局参数和样本列表
#################################

SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_match_bak1.tsv"

try:
    SAMPLES_DF = pd.read_csv(SAMPLE_SHEET, sep="\t").set_index("ribo_sample", drop=False)
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含 'ribo_sample' 和 'rna_sample' 列头。错误: {e}")

SAMPLES = SAMPLES_DF.index.tolist()

SPE = 'Human'
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']

RIBOTAPER_ENV = "ribotaper_env"

#################################
# 2. 路径定义
#################################
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_untrim/ribotaper"
RIBO_BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping/merge_chrN"
RNA_BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.5/Mapping/rnamerge_chrN"
RIBOSEQ_RESULTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_untrim/orfquant"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"

#################################
# --- 脚本主体 ---
#################################

workdir: OUT_DIR

rule all:
  input:
    'ribotaper_chrN/index/' + SPE + '/all_exons.bed',
    expand('ribotaper_chrN/orf_pred_default/{sample}_{mpsf}/ORFs_max', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_RiboTaper_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule genome_fasta_index:
  input:
    fa = FA
  output:
    fa_index = FA + '.fai'
  conda: RIBOTAPER_ENV
  shell:
    r'samtools faidx {input.fa}'

rule ribotaper_index:
  input:
    gtf = GTF,
    fa = FA,
    fa_index = FA + '.fai'
  output: 'ribotaper_chrN/index/' + SPE + '/all_exons.bed'
  conda: RIBOTAPER_ENV
  params:
    ribotaper_index = 'ribotaper_chrN/index/' + SPE
  benchmark: 'time_benchmarks/ribotaper_chrN/ribotaper_index/' + SPE + '.txt'
  shell:
    r'create_annotations_files.bash {input.gtf} {input.fa} false false {params.ribotaper_index}'

rule ribotaper_pred_default:
  input:
    ribo_bam = os.path.join(RIBO_BAM_DIR, '{sample}_{mpsf}.bam'),
    rna_bam = lambda wildcards: os.path.join(RNA_BAM_DIR, SAMPLES_DF.loc[wildcards.sample, "rna_sample"] + f'_{wildcards.mpsf}.bam'),
    all_exons = 'ribotaper_chrN/index/' + SPE + '/all_exons.bed',
    offsets = os.path.join(RIBOSEQ_RESULTS_DIR, 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt')
  output: 'ribotaper_chrN/orf_pred_default/{sample}_{mpsf}/ORFs_max'
  params:
    ribotaper_index = os.path.abspath('ribotaper_chrN/index/' + SPE),
    out_prefix = 'ribotaper_chrN/orf_pred_default/{sample}_{mpsf}',
    rnamerge = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, "rna_sample"]
  threads: 10
  conda: RIBOTAPER_ENV
  log: 'ribotaper_chrN/orf_pred_default/{sample}_{mpsf}/log'
  benchmark: 'time_benchmarks/ribotaper_chrN/orf_pred_default/{sample}_{mpsf}.txt'
  script:
    os.path.join(SCRIPTS_DIR, "ribotaper_analysis.sh")

rule filter_result:
  input: 'ribotaper_chrN/orf_pred_default/{sample}_{mpsf}/ORFs_max'
  output: 'ribotaper_chrN/orf_pred_default/{sample}_{mpsf}/ribotaper_filtered_005.txt'
  conda: RIBOTAPER_ENV
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')

rule raw_table_result:
  input: 'ribotaper_chrN/orf_pred_default/{sample}_{mpsf}/ribotaper_filtered_005.txt'
  output: 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_RiboTaper_raw.txt'
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')
