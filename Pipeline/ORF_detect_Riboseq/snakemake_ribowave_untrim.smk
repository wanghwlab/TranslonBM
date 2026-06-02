# snakefile_ribowave.smk 

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
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']

RIBOWAVE_ENV = "qc4"


#################################
# 2. 路径定义
#################################
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_untrim/ribowave"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping/merge_chrN"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"
GTFTOGENEPRED = '/home/tangyuewen/software/gtfToGenePred'
RIBOWAVE_SCRIPTS_DIR = '/home/tangyuewen/software/Ribowave-master/scripts'

#################################
# --- 脚本主体 ---
#################################

workdir: OUT_DIR

rule all:
    input:
        GTF + '.gpd',
        'ribowave_chrN/index/' + SPE + '/exons.gtf',
        expand('ribowave_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_result.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
        expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_RiboWave_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule samtools_faidx:
    input: FA
    output: FA + '.fai'
    conda: RIBOWAVE_ENV
    shell: "samtools faidx {input}"

rule create_annotation:
    input:
        gtf = GTF,
        fasta = FA,
    output:
        exon_gtf = 'ribowave_chrN/index/' + SPE + '/exons.gtf',
        start_codon = 'ribowave_chrN/index/' + SPE + '/start_codon.bed',
        orfs = 'ribowave_chrN/index/' + SPE + '/final.ORFs'
    conda: RIBOWAVE_ENV
    params: 'ribowave_chrN/index/' + SPE
    benchmark: 'time_benchmarks/ribowave_chrN/create_annotation/annotation.txt'
    shell:
        r'''bash {RIBOWAVE_SCRIPTS_DIR}/create_annotation.sh -G {input.gtf} -f {input.fasta} -o {params} -s {RIBOWAVE_SCRIPTS_DIR}'''

rule P_site_determination_chrN:
    input:
        ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
        start_codon = 'ribowave_chrN/index/' + SPE + '/start_codon.bed'
    output: 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/P-site/{sample}_{mpsf}.psite1nt.txt'
    conda: RIBOWAVE_ENV
    benchmark: 'time_benchmarks/ribowave_chrN/P_site_determination/{sample}_{mpsf}.txt'
    params:
        out_dir = 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}',
        out_name = '{sample}_{mpsf}'
    shell:
        r'''bash {RIBOWAVE_SCRIPTS_DIR}/P-site_determination.sh -i {input.ribo_bam} -S {input.start_codon} -o {params.out_dir} -n {params.out_name} -s {RIBOWAVE_SCRIPTS_DIR}'''

rule offsets_extraction:
  input:
    psite1nt = 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/P-site/{sample}_{mpsf}.psite1nt.txt'
  output:
    offsets = 'ribowave_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt'
  conda: RIBOWAVE_ENV
  script:
    os.path.join(SCRIPTS_DIR, 'ribowave_offsets_extraction.py')

rule chr_sizes:
    input: FA + '.fai'
    output: 'ribowave_chrN/index/' + SPE + '/sizes.genome'
    benchmark: 'time_benchmarks/ribowave_chrN/sizes_genome/sizes_genome.txt'
    shell:
        r'''cut -f1,2 {input} > {output}'''

rule create_psite_track_chrN_default:
    input:
        ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
        exon_gtf = 'ribowave_chrN/index/' + SPE + '/exons.gtf',
        genome_gtf = GTF,
        chr_sizes = 'ribowave_chrN/index/' + SPE + '/sizes.genome',
        psite = 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/P-site/{sample}_{mpsf}.psite1nt.txt'
    output: 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/bedgraph/{sample}_{mpsf}/final.psite'
    conda: RIBOWAVE_ENV
    benchmark: 'time_benchmarks/ribowave_chrN/create_psite_track/{sample}_{mpsf}.txt'
    params:
        out_dir = 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}',
        out_name = '{sample}_{mpsf}'
    shell:
        r'''bash {RIBOWAVE_SCRIPTS_DIR}/create_track_Ribo.sh -i {input.ribo_bam} -G {input.exon_gtf} -g {input.chr_sizes} -P {input.psite} -o {params.out_dir} -n {params.out_name} -s {RIBOWAVE_SCRIPTS_DIR}'''

rule ribowave_chrN_default:
    input:
        psite = 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/bedgraph/{sample}_{mpsf}/final.psite',
        orfs = 'ribowave_chrN/index/' + SPE + '/final.ORFs'
    output:
        'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/Ribowave/{sample}_{mpsf}.PF_psite',
        'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/Ribowave/{sample}_{mpsf}.feats1',
        'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/Ribowave/result/{sample}_{mpsf}.95%.mx'
    conda: RIBOWAVE_ENV
    benchmark: 'time_benchmarks/ribowave_chrN/orf_pred_default/{sample}_{mpsf}.txt'
    params:
        out_dir = 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/Ribowave',
        out_name = '{sample}_{mpsf}'
    threads: 24
    shell:
        r'''{RIBOWAVE_SCRIPTS_DIR}/Ribowave -PD -a {input.psite} -b {input.orfs} -o {params.out_dir} -n {params.out_name} -s {RIBOWAVE_SCRIPTS_DIR} -p {threads}'''

rule gtfTogpd:
  input:
    gtf = GTF
  output: 
    gpd = GTF + '.gpd'
  conda: RIBOWAVE_ENV
  shell:
    r'''{GTFTOGENEPRED} {input.gtf} {output.gpd}''' 

rule ribowave_type_chrN_default:
    input: 
        result = 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/Ribowave/result/{sample}_{mpsf}.95%.mx',
        gpd = GTF + '.gpd'
    output: 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_result.txt'
    conda: RIBOWAVE_ENV
    script: os.path.join(SCRIPTS_DIR, 'ribowave_mergetype.py')

rule raw_table_result:
  input: 'ribowave_chrN/orf_pred_default/{sample}_{mpsf}/{sample}_{mpsf}_result.txt'
  output: 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_RiboWave_raw.txt'
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')
