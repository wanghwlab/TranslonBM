# snakefile_ribocode.smk 

import os
import pandas as pd 

#################################
# 1. 设置全局参数和样本列表
#################################
SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names.tsv"

try:
    #SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t").sample.tolist()
    SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t")['sample'].tolist()
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含一个名为 'sample' 的列头。错误: {e}")


SPE = 'Human'
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2'] 

#################################
# 2. 路径定义
#################################
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribocode"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/merge_chrN_Tx"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
RIBOCODE_BIN_PATH = "~/miniconda3/envs/ribocode_env/bin/"

#################################
# --- 脚本主体 ---
#################################

workdir: OUT_DIR

rule all:
  input:
    # 使用 expand 函数，为 SAMPLES 列表中的每个样本生成最终的 raw_table_result 目标
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_RiboCode_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule ribocode_index:
  input:
    gtf = GTF,
    fa = FA
  output: 'ribocode_chrN/index/' + SPE + '/transcripts.pickle'
  params:
    ribocode_index = 'ribocode_chrN/index/' + SPE
  benchmark: 'time_benchmarks/ribocode_chrN/ribocode_index/' + SPE + '.txt'
  shell:
    r'''{RIBOCODE_BIN_PATH}prepare_transcripts \
    -g {input.gtf} \
    -f {input.fa} \
    -o {params.ribocode_index}'''

rule P_site_determination:
  input:
    ribo_Tx_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    ribocode_transcripts = 'ribocode_chrN/index/' + SPE + '/transcripts.pickle'
  output:
    'ribocode_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_pre_config.txt'
  params:
    ribocode_prefix = 'ribocode_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}',
    ribocode_index = 'ribocode_chrN/index/' + SPE
  benchmark:
    'time_benchmarks/ribocode_chrN/P_site_determination/{sample}_{mpsf}.txt'
  shell:
    r'''{RIBOCODE_BIN_PATH}metaplots -s yes \
    -a {params.ribocode_index} \
    -r {input.ribo_Tx_bam} \
    -o {params.ribocode_prefix}'''

rule offsets_extraction:
  input:
    config = 'ribocode_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_pre_config.txt'
  output:
    offsets = 'ribocode_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt'
  params:
    prefix = '{sample}_{mpsf}'
  script:
    '/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts/ribocode_offsets_extraction.py'

rule ribocode_pred_default:
  input:
    index = 'ribocode_chrN/index/' + SPE + '/transcripts.pickle',
    offset_config = 'ribocode_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_pre_config.txt'
  output:
    'ribocode_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_collapsed.txt',
    temp('{sample}_{mpsf}_psites.hd5')
  params:
    ribocode_index = 'ribocode_chrN/index/' + SPE,
    ribocode_prefix = 'ribocode_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}',
    min_aa_len = 8
  benchmark:
    'time_benchmarks/ribocode_chrN/orf_pred_default/{sample}_{mpsf}.txt'
  shell:
    r'''{RIBOCODE_BIN_PATH}RiboCode -l yes -b \
    -a {params.ribocode_index} \
    -c {input.offset_config} \
    -A CTG,GTG,TTG \
    -m {params.min_aa_len} \
    -o {params.ribocode_prefix}'''

rule raw_table_result:
  input:
    # {merge} -> {sample}
    'ribocode_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_collapsed.txt'
  output:
    # {merge} -> {sample}
    'orf_result/raw_prediction_result_default/{sample}_{mpsf}_RiboCode_raw.txt'
  script:
    '/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts/copyfile.py'
