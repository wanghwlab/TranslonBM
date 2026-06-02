# snakefile_ribohmm.smk

import os
import pandas as pd
import glob

#################################
# 1. 设置全局参数和样本列表
#################################
SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_match.tsv"
try:
    SAMPLES_DF = pd.read_csv(SAMPLE_SHEET, sep="\t").set_index("ribo_sample", drop=False)
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含 'ribo_sample' 和 'rna_sample' 列头。错误: {e}")

SAMPLES = SAMPLES_DF.index.tolist()

SPE = 'Human'
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']
READ_LENGTHS = ['28','29','30','31']
RIBOHMM_ENV = "riboHMM"

#################################
# 2. 路径定义
#################################
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribohmm_match"
RIBO_BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/merge_chrN"
RNA_BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.5/Mapping/rnamerge_chrN"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"
RIBOHMM_SRC_DIR = '/home/tangyuewen/ORF_benchmark/Ref/ORFtools/Ribohmm/scripts/riboHMM-py3'
MAPPABILITY_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribohmm/mappability"

#################################
# --- 脚本主体 ---
#################################

workdir: OUT_DIR

rule all:
    input:
        'ribohmm_chrN/index/GenePred',
        expand('ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/infer_CDS.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
        expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_riboHMM_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule gtfTogpd:
    input:
        gtf = GTF
    output: 
        gpd = 'ribohmm_chrN/index/GenePred'
    conda: RIBOHMM_ENV
    benchmark: 'time_benchmarks/ribohmm_chrN/gtfToGenePred/gtfToGenePred.txt'
    shell:
        r'''/home/tangyuewen/software/gtfToGenePred {input.gtf} {output.gpd}''' 

rule ribo_data_transmission:
    input:
        ribo_bam = os.path.join(RIBO_BAM_DIR, '{sample}_{mpsf}.bam'),
        ribo_bam_bai = os.path.join(RIBO_BAM_DIR, '{sample}_{mpsf}.bam.bai')
    output:
        ribo_bam = temp('ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/ribo_bam_tbi/{sample}_{mpsf}.bam'),
        ribo_bam_bai = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/ribo_bam_tbi/{sample}_{mpsf}.bam.bai'
    shell:
        r'''
            cp -f {input.ribo_bam} {output.ribo_bam} && \
            cp -f {input.ribo_bam_bai} {output.ribo_bam_bai}
        '''

rule rna_data_transmission:
    input:
        rna_bam = lambda wildcards: os.path.join(RNA_BAM_DIR, SAMPLES_DF.loc[wildcards.sample, "rna_sample"] + f'_{wildcards.mpsf}.bam'),
        rna_bam_bai = lambda wildcards: os.path.join(RNA_BAM_DIR, SAMPLES_DF.loc[wildcards.sample, "rna_sample"] + f'_{wildcards.mpsf}.bam.bai')
    output:
        rna_bam = temp('ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/rna_bam_tbi/rna_{mpsf}.bam'),
        rna_bam_bai = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/rna_bam_tbi/rna_{mpsf}.bam.bai'
    shell:
        r'''
            cp -f {input.rna_bam} {output.rna_bam} && \
            cp -f {input.rna_bam_bai} {output.rna_bam_bai}
        '''

rule ribo_data_convert:
    input:
        ribo_bam = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/ribo_bam_tbi/{sample}_{mpsf}.bam',
    output:
        fwd_tbi = expand('ribohmm_chrN/ORF_detecting_default/{{sample}}_{{mpsf}}/ribo_bam_tbi/{{sample}}_{{mpsf}}_fwd.{read_length}.gz.tbi', read_length=READ_LENGTHS),
        rev_tbi = expand('ribohmm_chrN/ORF_detecting_default/{{sample}}_{{mpsf}}/ribo_bam_tbi/{{sample}}_{{mpsf}}_rev.{read_length}.gz.tbi', read_length=READ_LENGTHS)
    conda: RIBOHMM_ENV
    benchmark: 'time_benchmarks/ribohmm_chrN/ORF_detecting_default/bam_to_tbi/RIBO_{sample}_{mpsf}.txt'
    shell:
        r'python {RIBOHMM_SRC_DIR}/bam_to_tbi.py --dtype riboseq {input.ribo_bam}'

rule rna_data_convert:
    input:
        rna_bam = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/rna_bam_tbi/rna_{mpsf}.bam',
    output: 
        tbi = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/rna_bam_tbi/rna_{mpsf}.gz.tbi'
    benchmark: 'time_benchmarks/ribohmm_chrN/ORF_detecting_default/bam_to_tbi/RNA_{sample}_{mpsf}.txt'
    conda: RIBOHMM_ENV
    shell:
        r'python {RIBOHMM_SRC_DIR}/bam_to_tbi.py --dtype rnaseq {input.rna_bam}'

rule learn_model_default:
    input:
        fasta = FA,
        gtf = GTF,
        ribo_fwd_tbi = expand('ribohmm_chrN/ORF_detecting_default/{{sample}}_{{mpsf}}/ribo_bam_tbi/{{sample}}_{{mpsf}}_fwd.{read_length}.gz.tbi', read_length=READ_LENGTHS),
        rna_tbi = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/rna_bam_tbi/rna_{mpsf}.gz.tbi'
    output:
        model = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/learn_model.txt'
    params:
        riboseq_prefix = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/ribo_bam_tbi/{sample}_{mpsf}',
        rnaseq_prefix = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/rna_bam_tbi/rna_{mpsf}',
        mappability_prefix = os.path.join(MAPPABILITY_DIR, 'mappability_{mpsf}')
    log: 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/{sample}_{mpsf}_learn_model.log'
    conda: RIBOHMM_ENV
    benchmark: 'time_benchmarks/ribohmm_chrN/ORF_detecting_default/learn_model_default/{sample}_{mpsf}.txt'
    shell:
        r'''
            python {RIBOHMM_SRC_DIR}/learn_model.py \
            --rnaseq_file {params.rnaseq_prefix} \
            --mappability_file {params.mappability_prefix} \
            --log_file {log} \
            --model_file {output.model} \
            {input.fasta} {input.gtf} {params.riboseq_prefix}
        '''

rule infer_CDS_default:
    input:
        fasta = FA,
        gtf = GTF,
        model = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/learn_model.txt',
        rna_tbi = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/rna_bam_tbi/rna_{mpsf}.gz.tbi',
        ribo_fwd_tbi = expand('ribohmm_chrN/ORF_detecting_default/{{sample}}_{{mpsf}}/ribo_bam_tbi/{{sample}}_{{mpsf}}_fwd.{read_length}.gz.tbi', read_length=READ_LENGTHS)
    output:
         infer = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/infer_CDS.txt'
    params:
        riboseq_prefix = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/ribo_bam_tbi/{sample}_{mpsf}',
        rnaseq_prefix = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/rna_bam_tbi/rna_{mpsf}',
        mappability_prefix = os.path.join(MAPPABILITY_DIR, 'mappability_{mpsf}')
    benchmark: 'time_benchmarks/ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}.txt'
    conda: RIBOHMM_ENV
    shell:
        r'''
            python {RIBOHMM_SRC_DIR}/infer_CDS.py \
            --output_file {output.infer} \
            --rnaseq_file {params.rnaseq_prefix} \
            --mappability_file {params.mappability_prefix} \
            {input.model} {input.fasta} {input.gtf} {params.riboseq_prefix}
        '''
        
rule filter_result:
    input:
        infer_result = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/infer_CDS.txt',
        gpd = 'ribohmm_chrN/index/GenePred'
    output:
        filterd_result = 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/ribohmm_filtered.txt'
    conda: RIBOHMM_ENV
    script: os.path.join(SCRIPTS_DIR, 'ribohmm_result_filtering.py')

rule raw_table_result:
    input: 'ribohmm_chrN/ORF_detecting_default/{sample}_{mpsf}/ribohmm_filtered.txt'
    output: 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_riboHMM_raw.txt'
    script: os.path.join(SCRIPTS_DIR, 'copyfile.py')
