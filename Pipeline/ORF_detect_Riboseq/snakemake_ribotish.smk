# snakefile_ribotish.smk 

import os
import pandas as pd 

SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_specific.tsv"

try:
    SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t")['sample'].tolist()
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含一个名为 'sample' 的列头。错误: {e}")

MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']

OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribotish"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/merge_chrN"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"

workdir: OUT_DIR

rule all:
  input:
    expand('ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_para.py', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_nom0_50_para.py', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('ribotish_chrN/orf_pred_default/{sample}_{mpsf}/ribotish_filtered_005_longest.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_Ribo-TISH-longest_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule P_site_determination:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    gtf = GTF
  output:
    length_profile = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_qual.txt',
    offsets = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_para.py',
    figure = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_qual.pdf'
  conda:'ribotish_env'
  params:
    # {merge} -> {sample}
    ribotish_prefix = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}'
  benchmark: 'time_benchmarks/ribotish_chrN/P_site_determination/{sample}_{mpsf}.txt'
  threads: 1
  shell:
    r'''ribotish quality -p {threads} -l 25,36 \
    -b {input.ribo_bam} -g {input.gtf} \
    -f {output.figure} \
    -o {output.length_profile} \
    -r {output.offsets}'''

rule offsets_extraction:
  input: 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_para.py'
  output: 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt'
  conda:'ribotish_env'
  script: os.path.join(SCRIPTS_DIR, 'ribotish_offsets_extraction.py')

rule P_site_determination_50_nom0:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    gtf = GTF
  output:
    length_profile = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_nom0_50_qual.txt',
    offsets = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_nom0_50_para.py',
    figure = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_nom0_50_qual.pdf'
  conda:'ribotish_env'
  params:
    # {merge} -> {sample}
    ribotish_prefix = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_nom0_50'
  threads: 1
  shell:
    r'''ribotish quality -p {threads} -l 25,36 --nom0 -d-50,50 \
    -b {input.ribo_bam} -g {input.gtf} \
    -f {output.figure} \
    -o {output.length_profile} \
    -r {output.offsets}'''

rule ribotish_pred_default:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    gtf = GTF,
    fa = FA,
    offsets = 'ribotish_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_para.py'
  output:
    prediction = 'ribotish_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_pred.txt'
  conda:'ribotish_env'
  params:
    min_aa_len = 8
  benchmark: 'time_benchmarks/ribotish_chrN/orf_pred_default/{sample}_{mpsf}.txt'
  threads: 1
  shell:
    r'''ribotish predict \
    -p {threads} \
    --longest \
    --seq \
    --aaseq \
    --blocks \
    --minaalen {params.min_aa_len} \
    --altcodons CTG,GTG,TTG \
    -b {input.ribo_bam} \
    -g {input.gtf} -f {input.fa} \
    --ribopara {input.offsets} \
    -o {output.prediction}'''

rule filter_result:
  input:
    longest_prediction = 'ribotish_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_pred.txt',
  output:
    longest_result = 'ribotish_chrN/orf_pred_default/{sample}_{mpsf}/ribotish_filtered_005_longest.txt'
  conda:'ribotish_env'
  script: os.path.join(SCRIPTS_DIR, 'ribotish_result_filtering.py')

rule raw_table_result:
  input:
    longest_result = 'ribotish_chrN/orf_pred_default/{sample}_{mpsf}/ribotish_filtered_005_longest.txt',
  output:
    longest_result = 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_Ribo-TISH-longest_raw.txt'
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')
