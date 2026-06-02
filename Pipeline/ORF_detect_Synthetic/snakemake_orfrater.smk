# snakefile_orfrater.smk 

import os
import pandas as pd # 新增: 导入pandas库

#################################
# 1. 设置全局参数和样本列表
#################################
SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_simu.tsv"

try:
    SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t")['sample'].tolist()
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含一个名为 'sample' 的列头。错误: {e}")

SPE = 'Human'
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']

ORFRATER_ENV = "python27"

#################################
# 2. 路径定义
#################################
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_simu/orfrater"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping/simulation"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect_simu/scripts"
ORFRATER_DIR = '/home/tangyuewen/ORF_benchmark/Ref/ORFtools/ORF-RATER'

#################################
# --- 脚本主体 ---
#################################
workdir: OUT_DIR

rule all:
  input:
    'orfrater_chrN/index/' + SPE + '.bed',
    expand('orfrater_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orfrater_chrN/orf_pred_default/{sample}_{mpsf}/tfams.txt',sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orfrater_chrN/orf_pred_default/{sample}_{mpsf}/ratedorfs.bed', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orfrater_filtered_08.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_ORFRATER_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule gtf2bed:
  input: GTF
  output: 'orfrater_chrN/index/' + SPE + '.bed'
  conda: ORFRATER_ENV
  benchmark: 'time_benchmarks/orfrater_chrN/gtf2bed/gtf2bed.txt'
  shell:
    r'''
    /home/tangyuewen/software/gtfToGenePred -ignoreGroupsWithoutExons {input} stdout | \
    /home/tangyuewen/software/genePredToBed stdin {output}
    '''

rule prune_transcripts:
  input:
    bed = 'orfrater_chrN/index/' + SPE + '.bed',
    fa = FA,
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam')
  output:
    transcripts_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/transcript.bed',
    summarytable = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/tid_removal_summary.txt'
  params:
    workdir = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}',
    min_len = 28,
    max_len = 31,
    inbed = lambda wildcards, input: os.path.abspath(input.bed),
    fa = lambda wildcards, input: os.path.abspath(input.fa),
    ribo_bam = lambda wildcards, input: os.path.abspath(input.ribo_bam),
    summarytable = lambda wildcards, output: os.path.basename(str(output.summarytable)),
    outbed = lambda wildcards, output: os.path.basename(str(output.transcripts_bed))
  log: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/1prune.log'
  benchmark: 'time_benchmarks/orfrater_chrN/prune_transcripts/{sample}_{mpsf}.txt'
  threads: 1
  conda: 'python27'
  shell:
    r"""
    (cd {params.workdir} && python {ORFRATER_DIR}/prune_transcripts.py \
    --inbed {params.inbed} \
    --summarytable {params.summarytable} \
    --outbed {params.outbed} \
    --minlen {params.min_len} --maxlen {params.max_len} \
    -v -p {threads} -f \
    {params.fa} {params.ribo_bam})
    """

rule make_tfams:
  input:
    transcripts_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/transcript.bed',
    genename_looktable = os.path.join(os.path.dirname(GTF), 'looktable')
  output: 
    tfams_txt = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/tfams.txt',
    tfams_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/tfams.bed'
  params:
    workdir = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}',
    genename_looktable = lambda wildcards, input: os.path.abspath(input.genename_looktable),
    transcripts_bed = lambda wildcards, input: os.path.basename(str(input.transcripts_bed))
  log: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/2make_tfams.log'
  benchmark: 'time_benchmarks/orfrater_chrN/make_tfams/{sample}_{mpsf}.txt'
  conda: 'python27'
  shell:
    r"""
    (cd {params.workdir} && python {ORFRATER_DIR}/make_tfams.py \
    -g {params.genename_looktable} \
    --inbed {params.transcripts_bed} \
    -v -f)
    """    

rule find_orfs_and_types:
  input:
    fa = FA,
    transcripts_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/transcript.bed',
    tfams_txt = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/tfams.txt'
  output:
    orf_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orf.h5'
  params:
    workdir = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}',
    start_codons = 'NTG',
    fa = lambda wildcards, input: os.path.abspath(input.fa),
    orfstore = lambda wildcards, output: os.path.basename(str(output.orf_h5)),
    inbed = lambda wildcards, input: os.path.basename(str(input.transcripts_bed))
  threads: 1
  log: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/3find_ORF.log'
  benchmark: 'time_benchmarks/orfrater_chrN/find_orfs_and_types/{sample}_{mpsf}.txt'
  conda: ORFRATER_ENV
  shell:
    r"""
    (cd {params.workdir} && python {ORFRATER_DIR}/find_orfs_and_types.py {params.fa} \
    --orfstore {params.orfstore} \
    --inbed {params.inbed} \
    --codons {params.start_codons} \
    -v -p {threads} -f)
    """

rule psite_trimmed:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    source_bed = 'orfrater_chrN/index/' + SPE + '.bed',
    orf_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orf.h5'
  output:
    offsets = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/Ribo/offsets.txt',
    tallies = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/Ribo/tallies.txt'
  conda: ORFRATER_ENV
  params:
    workdir = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}',
    sub_dir = 'Ribo',
    min_len = 28,
    max_len = 31,
    ribo_bam = lambda wildcards, input: os.path.abspath(input.ribo_bam),
    offsetfile = lambda wildcards, output: os.path.basename(str(output.offsets)),
    cdsbed = lambda wildcards, input: os.path.abspath(input.source_bed),
    tallyfile = lambda wildcards, output: os.path.basename(str(output.tallies))
  log: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/4psite.log'
  benchmark: 'time_benchmarks/orfrater_chrN/P_site_determination/{sample}_{mpsf}.txt'
  threads: 1
  shell:
    r"""
    (cd {params.workdir} && python {ORFRATER_DIR}/psite_trimmed.py {params.ribo_bam} \
    --subdir {params.sub_dir} \
    --offsetfile {params.offsetfile} \
    --cdsbed {params.cdsbed} \
    --minrdlen {params.min_len} \
    --maxrdlen {params.max_len} \
    --tallyfile {params.tallyfile} \
    -v -p {threads} -f)
    """

rule regress_orfs:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    orf_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orf.h5',
    transcripts_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/transcript.bed',
    offsets = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/Ribo/offsets.txt'
  output:
    regression_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/Ribo/regression.h5'
  conda: ORFRATER_ENV
  benchmark: 'time_benchmarks/orfrater_chrN/orf_pred_default/{sample}_{mpsf}_regress_orfs.txt'
  params:
    workdir = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}',
    sub_dir = 'Ribo',
    orfstore = lambda wildcards, input: os.path.basename(str(input.orf_h5)),
    inbed = lambda wildcards, input: os.path.basename(str(input.transcripts_bed)),
    offsetfile = lambda wildcards, input: os.path.basename(str(input.offsets)),
    ribo_bam = lambda wildcards, input: os.path.abspath(input.ribo_bam)
  log: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/5regress_orfs.log'
  threads: 1
  shell:
    r"""
    (cd {params.workdir} && python {ORFRATER_DIR}/regress_orfs.py \
    --subdir {params.sub_dir} \
    --orfstore {params.orfstore} \
    --inbed {params.inbed} \
    --offsetfile {params.offsetfile} \
    -v -p {threads} {params.ribo_bam})
    """ 

rule rate_regression_output:
  input:
    regression_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/Ribo/regression.h5',
    orf_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orf.h5'
  output:
    orfratings_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orfratings.h5',
    csv = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/rate_regression.csv'
  conda: ORFRATER_ENV
  log: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/6rate_regression_output.log'
  benchmark: 'time_benchmarks/orfrater_chrN/orf_pred_default/{sample}_{mpsf}_rate_regression_output.txt'
  threads: 1
  params:
    workdir = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}',
    orfstore = lambda wildcards, input: os.path.basename(str(input.orf_h5)),
    ratingsfile = lambda wildcards, output: os.path.basename(str(output.orfratings_h5)),
    csv = lambda wildcards, output: os.path.basename(str(output.csv)),
    regression_h5 = lambda wildcards, input: os.path.basename(str(input.regression_h5))
  shell:
    r"""
    (cd {params.workdir} && python {ORFRATER_DIR}/rate_regression_output.py \
    --orfstore {params.orfstore} \
    --ratingsfile {params.ratingsfile} \
    --CSV {params.csv} \
    -v -p {threads} -f {params.regression_h5})
    """ 

rule make_orf_bed:
  input:
    transcripts_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/transcript.bed',
    orfratings_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orfratings.h5'
  output:
    ratedorfs_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/ratedorfs.bed'
  params:
    min_len = 8,
    min_rating = 0.8,
    workdir = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}',
    inbed = lambda wildcards, input: os.path.basename(str(input.transcripts_bed)),
    outbed = lambda wildcards, output: os.path.basename(str(output.ratedorfs_bed)),
    ratingsfile = lambda wildcards, input: os.path.basename(str(input.orfratings_h5))
  conda: ORFRATER_ENV
  log: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/7make_orf_bed.log'
  benchmark: 'time_benchmarks/orfrater_chrN/make_orf_bed/{sample}_{mpsf}.txt'
  shell:
    r"""
    (cd {params.workdir} && python {ORFRATER_DIR}/make_orf_bed.py \
    --inbed {params.inbed} \
    --outbed {params.outbed} \
    --ratingsfile {params.ratingsfile} \
    --minlen {params.min_len} \
    --minrating {params.min_rating})
    """ 

rule quantify_orfs:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam'),
    transcripts_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/transcript.bed',
    offsets = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/Ribo/offsets.txt',
    orfratings_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orfratings.h5',
    ratedorfs_bed = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/ratedorfs.bed'
  output:
    quant_h5 = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/quant.h5',
    csv = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/quant.csv'
  conda: ORFRATER_ENV
  params:
    workdir = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}',
    sub_dir = 'Ribo',
    min_len = 8,
    min_rating = 0.8,
    inbed = lambda wildcards, input: os.path.basename(str(input.transcripts_bed)),
    offsetfile = lambda wildcards, input: os.path.basename(str(input.offsets)),
    ratingsfile = lambda wildcards, input: os.path.basename(str(input.orfratings_h5)),
    quantfile = lambda wildcards, output: os.path.basename(str(output.quant_h5)),
    csv = lambda wildcards, output: os.path.basename(str(output.csv)),
    ribo_bam = lambda wildcards, input: os.path.abspath(input.ribo_bam)
  threads: 1
  log: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/8quantify_orfs.log' 
  benchmark: 'time_benchmarks/orfrater_chrN/orf_pred_default/{sample}_{mpsf}_quantify_orfs.txt'
  shell:
    r"""
    (cd {params.workdir} && python {ORFRATER_DIR}/quantify_orfs.py \
    --subdir {params.sub_dir} \
    --inbed {params.inbed} \
    --offsetfile {params.offsetfile} \
    --ratingsfile {params.ratingsfile} \
    --minrating {params.min_rating} \
    --minlen {params.min_len} \
    --quantfile {params.quantfile} \
    --CSV {params.csv} \
    -v -p {threads} -f {params.ribo_bam})
    """ 
    
rule offsets_extraction:
  input:
    orfrater_offsets = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/Ribo/offsets.txt'
  output:
    standard_offsets = 'orfrater_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt'
  script: 
    os.path.join(SCRIPTS_DIR, 'copyfile.py')

rule filter_result:
  input:
    all_result = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/quant.csv'
  output: 
    result = 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orfrater_filtered_08.txt'
  conda: ORFRATER_ENV
  script: os.path.join(SCRIPTS_DIR, 'orfrater_result_filtering.py')

rule raw_table_result:
  input: 'orfrater_chrN/orf_pred_default/{sample}_{mpsf}/orfrater_filtered_08.txt'
  output: 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_ORFRATER_raw.txt'
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')
