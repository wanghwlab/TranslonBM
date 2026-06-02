# snakefile_orfquant.smk 
import os
import pandas as pd 

#################################
# 1. 设置全局参数和样本列表
#################################
SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_specific.tsv"
try:
    SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t")['sample'].tolist()
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含一个名为 'sample' 的列头。错误: {e}")

SPE = 'Human' 
MAPPING_SOFTWARE = ['tophat2','STAR', 'hisat2'] #'STAR', 'hisat2',
SCIENTIFIC_NAME_LIST = ['Homo.sapiens']
ANNO_NAME_LIST = ['gencodev43']

RIBOSEQ_ENV = 'riboseqc_env'
ORFQUANT_ENV = 'orfquant_env'

#################################
# 2. 路径定义
#################################
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_rerun/orfquant"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/merge_chrN"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"
GENOME_FILE = "/home/tangyuewen/ORF_benchmark/Ref/genome_file" # 修改: 使用您自己的路径

R_ANNO = os.path.basename(GTF)
if SPE == 'Human':
    SCIENTIFIC_NAME = SCIENTIFIC_NAME_LIST[0]
    ANNO_NAME = ANNO_NAME_LIST[0]
# ... (保留其他物种的 elif 判断) ...
else:
    print('please check whether the right scientific name in this snakefile!!!')
    os._exit(1)

#################################
# --- 脚本主体 ---
#################################
workdir: OUT_DIR

rule all:
  input:
    FA + '.2bit',
    #'ORFquant_chrN/index/' + SPE + '/' + R_ANNO + '_Rannot',
    #expand('merge_chrN/{sample}_{mpsf}_mapq_adjust.bam', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    #expand('riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_result.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_ORFquant_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)
    
rule fa2twobit:
  input: fa = FA
  output: twobit = FA + '.2bit'
  shell: r'''/home/tangyuewen/software/faToTwoBit {input} {output}'''

#rule prepare_anno:
#  input:
#    gtf = GTF,
#    twobit = FA + '.2bit',
#    fa = FA
#  output: 
#    annotation_file = 'ORFquant_chrN/index/' + SPE + '/' + R_ANNO + '_Rannot'
#  conda: 'orfquant_env'
#  params:
#    index_dir = 'ORFquant_chrN/index/' + SPE,
#    scientific_name = SCIENTIFIC_NAME,
#    annotation_name = ANNO_NAME
#  benchmark: 'time_benchmarks/ORFquant_chrN/ORFquant_index/'+ SPE + '.txt'
#  threads: 1
#  script: os.path.join(SCRIPTS_DIR, 'ORFquant_prepare_anno.R')

#rule adjust_mapq_tophat2:
#  input:
#    raw_bam = os.path.join(BAM_DIR, '{sample}_tophat2.bam'),
#    genome_file = GENOME_FILE
#  output:
    # 明确定义两个输出文件
#    adjusted_bam = 'merge_chrN/{sample}_tophat2_mapq_adjust.bam',
    # 使用 temp() 标记这个BED文件，它是一个临时的中间文件
#    bed = temp('merge_chrN/{sample}_tophat2.bed') 
#  conda: 'base'
#  threads: 6
#  script:
#    os.path.join(SCRIPTS_DIR, 'adjust_mapq.sh')

rule adjust_mapq_others:
    input:
        raw_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam')
    output:
        adjusted_bam = 'merge_chrN/{sample}_{mpsf}_mapq_adjust.bam'
    wildcard_constraints:
        mpsf="STAR|hisat2"
    conda: 'py3.7'
    threads: 24
    shell: r'''~/miniconda3/bin/samtools view -@ 6 -bh -o {output.adjusted_bam} {input.raw_bam}'''


rule Riboseqc_P_site_determination_chrN:
  input:
    ribo_bam = 'merge_chrN/{sample}_{mpsf}_mapq_adjust.bam',
    annotation_file = 'RiboseQC_chrN/index/' + SPE + '/' + R_ANNO + '_Rannot',
    fa = FA
  output:
    result = 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_results_RiboseQC',
    P_sites_calcs = 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_P_sites_calcs'
  conda: 'riboseqc_env'
  params:
    prefix = 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}'
  threads: 24
  benchmark: 'time_benchmarks/riboseqc_chrN/P_site_determination/{sample}_{mpsf}.txt'
  script: os.path.join(SCRIPTS_DIR, 'Riboseqc_analysis.R')

rule offsets_extraction:
  input:
    result = 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_results_RiboseQC',
    P_sites_calcs = 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_P_sites_calcs'
  output:
    offsets = 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt'
  script: os.path.join(SCRIPTS_DIR, 'Riboseqc_offsets_extraction.py')

rule extract_selected_read_length:
  input:
    riboseqc_P_sites_calcs = 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_P_sites_calcs'
  output:
    read_lengths_file_for_ORFquant = 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_rl_file'
  script: os.path.join(SCRIPTS_DIR, 'ORFquant_extract_rl.py')

rule prepare_for_ORFquant_default:
  input:
    annotation_file = 'ORFquant_chrN/index/' + SPE + '/' + R_ANNO + '_Rannot',
    bam_file = 'merge_chrN/{sample}_{mpsf}_mapq_adjust.bam',
    read_lengths_file = 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_rl_file'
  output:
    for_ORFquant = 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_for_ORFquant'
  params:
    prefix = 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}'
  conda: 'orfquant_env'
  benchmark: 'time_benchmarks/riboseqc_chrN/prepare_for_ORFquant/{sample}_{mpsf}.txt'
  script: os.path.join(SCRIPTS_DIR, 'ORFquant_prepare_for_ORFquant.R')

rule run_ORFquant_default:
  input:
    for_ORFquant = 'riboseqc_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_for_ORFquant',
    annotation_file = 'ORFquant_chrN/index/' + SPE + '/' + R_ANNO + '_Rannot'
  output:
    result = 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_final_ORFquant_results'
  params:
    prefix = 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}' 
  conda: 'orfquant_env'
  benchmark: 'time_benchmarks/ORFquant_chrN/orf_pred_default/{sample}_{mpsf}.txt'
  threads: 8
  script: os.path.join(SCRIPTS_DIR, 'ORFquant_run.R')

rule ORFquant_result_extraction:
  input:
    result = 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_final_ORFquant_results'
  output:
    result_table = 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_result.txt'
  conda: 'orfquant_env'
  script: os.path.join(SCRIPTS_DIR, 'ORFquant_result_extraction.R')
    
rule raw_table_result:
  input: 'ORFquant_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_result.txt' 
  output: 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_ORFquant_raw.txt'
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')
