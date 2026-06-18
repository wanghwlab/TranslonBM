# snakefile_gedi.smk 
import os
import pandas as pd 

SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_specific.tsv"

try:
    SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t")['sample'].tolist()
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含一个名为 'sample' 的列头。错误: {e}")

SPE = 'Human'
MAPPING_SOFTWARE = ['hisat2', 'STAR', 'tophat2']
GEDI = '/home/tangyuewen/software/Gedi/gedi-Gedi_1.0.6a/Gedi/gedi'

OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/gedi"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/merge_chrN"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"

workdir: OUT_DIR

rule all:
  input:
    expand('gedi_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}.orfs.tsv', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('gedi_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_FDR01.orfs.bed', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_PRICE_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule gedi_price_index:
  input:
    gtf = GTF,
    fa = FA
  output:
    oml = 'gedi_chrN/index/' + SPE + '.oml'
  conda:
    "gedi_env"
  params:
    index_dir = 'gedi_chrN/index',
    spe = SPE
  benchmark: 'time_benchmarks/gedi_chrN/index/' + SPE + '.txt'
  shell:
    r'''{GEDI} -e IndexGenome -nobowtie -nostar -nokallisto \
      -s {input.fa} \
      -a {input.gtf} \
      -n {params.spe} \
      -o {output.oml} \
      -f {params.index_dir}'''

rule gedi_price_pred_default:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    oml = 'gedi_chrN/index/' + SPE + '.oml'
  output:
    prediction = 'gedi_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}.orfs.tsv',
    prediction_orf = 'gedi_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}.orfs.cit'
  conda:
    "gedi_env"
  params:
    prefix = 'gedi_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}'
  benchmark: 'time_benchmarks/gedi_chrN/orf_pred_default/{sample}_{mpsf}.txt'
  threads: 1
  shell:
    r'''{GEDI} -e Price \
      -nthreads {threads} \
      -genomic {input.oml} \
      -prefix {params.prefix} \
      -reads {input.ribo_bam}'''

rule gedi_price_pred_FDR01_bed12:
  input: 'gedi_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}.orfs.cit'
  output: 'gedi_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_FDR01.orfs.bed'
  conda:
    "gedi_env"
  benchmark: 'time_benchmarks/gedi_chrN/orf_pred_default/ViewCIT/{sample}_{mpsf}.txt'
  shell:
    r'''{GEDI} -e ViewCIT \
      -m bed -name 'd.transcript+"_"+d.type+"_"+d.orfid' \
      {input} > {output}'''

rule raw_table_result:
  input: 'gedi_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_FDR01.orfs.bed'
  output: 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_PRICE_raw.txt'
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')
