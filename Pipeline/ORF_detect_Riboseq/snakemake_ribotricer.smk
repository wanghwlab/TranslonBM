# snakefile_ribotricer.smk 

import os
import pandas as pd 

SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_specific.tsv"

try:
    SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t")['sample'].tolist()
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含一个名为 'sample' 的列头。错误: {e}")

SPE = 'Human'
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']

RIBOTRICER_ENV = "ribotricer_env"

OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribotricer"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/merge_chrN"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"


workdir: OUT_DIR

rule all:
  input:
    'ribotricer_chrN/index/' + SPE + '_candidate_orfs.tsv',
    'ribotricer_chrN/index/' + SPE + '.fasta',
    expand('ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_translating_ORFs.tsv', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_Ribotricer_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('ribotricer_chrN/count/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_gene_level.tsv', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule ribotricer_index:
  input:
    gtf = GTF,
    fa = FA
  output: 'ribotricer_chrN/index/' + SPE + '_candidate_orfs.tsv'
  conda: RIBOTRICER_ENV
  params:
    ribotricer_index = 'ribotricer_chrN/index/' + SPE,
    start = 'ATG,CTG,GTG,TTG',
    stop = 'TAG,TAA,TGA',
    min_orf_length = 24
  benchmark: 'time_benchmarks/ribotricer_chrN/ribotricer_index/' + SPE + '.txt'
  shell:
    r'''ribotricer prepare-orfs \
    --gtf {input.gtf} \
    --fasta {input.fa} \
    --prefix {params.ribotricer_index} \
    --start_codons {params.start} \
    --stop_codons {params.stop} \
    --min_orf_length {params.min_orf_length}
    '''

rule ribotricer_orfs_seq:
  input:
    fa = FA,
    ribotricer_index = 'ribotricer_chrN/index/' + SPE + '_candidate_orfs.tsv',
  output: 'ribotricer_chrN/index/' + SPE + '.fasta'
  conda: RIBOTRICER_ENV
  benchmark: 'time_benchmarks/ribotricer_chrN/ribotricer_orfs_seq/ribotricer_orfs_seq.txt'
  shell:
    r'''ribotricer orfs-seq \
    --ribotricer_index {input.ribotricer_index} \
    --fasta {input.fa} \
    --saveto {output}''' 

rule ribotricer_pred_chrN_default:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    ribotricer_index = 'ribotricer_chrN/index/' + SPE + '_candidate_orfs.tsv'
  output: 'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_translating_ORFs.tsv'
  conda: RIBOTRICER_ENV
  params:
    # {merge} -> {sample}
    ribotricer_prefix = 'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}'
  log: 'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}.time.log'
  benchmark: 'time_benchmarks/ribotricer_chrN/orf_pred_default/{sample}_{mpsf}.txt'
  shell:
    r'''ribotricer detect-orfs \
    --bam {input.ribo_bam} \
    --ribotricer_index {input.ribotricer_index} \
    --prefix {params.ribotricer_prefix} >{log} ''' 

rule time_calc:
  input: 
    'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}.time.log', 
    'time_benchmarks/ribotricer_chrN/orf_pred_default/{sample}_{mpsf}.txt'
  output: 'time_benchmarks/ribotricer_chrN/P_site_determination/{sample}_{mpsf}.txt'
  conda: RIBOTRICER_ENV
  script: os.path.join(SCRIPTS_DIR, 'ribotricer_psite_time_calc.py')

rule filter_result:
  input: 'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_translating_ORFs.tsv'
  output: 'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_ribotricer_filtered.txt'
  conda: RIBOTRICER_ENV
  script: os.path.join(SCRIPTS_DIR, 'ribotricer_result_filtering.py')

rule raw_table_result:
  input: 'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_ribotricer_filtered.txt'
  output: 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_Ribotricer_raw.txt'
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')

rule ribotricer_quantify_gene_level:
  input:
    detected_orfs = 'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_translating_ORFs.tsv',
    ribotricer_index = 'ribotricer_chrN/index/' + SPE + '_candidate_orfs.tsv'
  output: 'ribotricer_chrN/count/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_gene_level.tsv'
  conda: RIBOTRICER_ENV
  params: 
    orf_types = 'annotated,super_uORF,super_dORF,uORF,dORF,overlap_uORF,overlap_dORF,novel'
  benchmark: 'time_benchmarks/ribotricer_chrN/orf_pred_default/ribotricer_quantify_gene_level/{sample}_{mpsf}.txt'
  shell:
    r'''ribotricer count-orfs \
    --ribotricer_index {input.ribotricer_index} \
    --detected_orfs {input.detected_orfs} \
    --features {params.orf_types} \
    --out {output} \
    --report_all''' 

rule ribotricer_quantify_codon_level:
  input:
    detected_orfs = 'ribotricer_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_translating_ORFs.tsv',
    ribotricer_index = 'ribotricer_chrN/index/' + SPE + '_candidate_orfs.tsv',
    orf_seq_file = 'ribotricer_chrN/index/' + SPE + '.fasta'
  output: 'ribotricer_chrN/count/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_codonwise.tsv'
  conda: RIBOTRICER_ENV
  params: 
    orf_types = 'annotated,super_uORF,super_dORF,uORF,dORF,overlap_uORF,overlap_dORF,novel',
    prefix = 'ribotricer_chrN/count/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}',
  benchmark: 'time_benchmarks/ribotricer_chrN/orf_pred_default/ribotricer_quantify_codon_level/{sample}_{mpsf}.txt'
  shell:
    r'''ribotricer count-orfs-codon \
    --ribotricer_index {input.ribotricer_index} \
    --detected_orfs {input.detected_orfs} \
    --features {params.orf_types} \
    --ribotricer_index_fasta {input.orf_seq_file} \
    --prefix {params.prefix} \
    --report_all'''
