# snakefile_riborf.smk 
import os
import pandas as pd 
import glob 

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
RIBORFSCRIPTS_DIR = '/home/tangyuewen/software/RibORF.2.0'


#################################
# 2. 路径定义
#################################
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/riborf"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/merge_chrN"
GTF = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"

#################################
# --- 脚本主体 ---
#################################

workdir: OUT_DIR

rule all:
  input:
    'riborf_chrN/index/' + SPE + '/candidateORF.genepred.txt',
    expand('riborf_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('riborf_chrN/P_site_determination/{sample}_{mpsf}/corrected/readDist.plot.1.R', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('riborf_chrN/orf_pred_default/{sample}_{mpsf}/corrected.{sample}_{mpsf}.sam', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('riborf_chrN/orf_pred_default/{sample}_{mpsf}/repre.valid.pred.pvalue.parameters.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('riborf_chrN/orf_pred_default/{sample}_{mpsf}/riborf_filtered_07.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_RibORF_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)

rule gtfTogpd:
  input:
    gtf = GTF
  output: 
    gpd = 'riborf_chrN/index/' + SPE + '.gpd'
  benchmark: 'time_benchmarks/riborf_chrN/gtfTogpd/' + SPE + '.txt'
  conda:'qc_snakemake'
  shell:
    r'''/home/tangyuewen/software/gtfToGenePred {input.gtf} {output.gpd}''' 

rule riborf_index:
  input:
    gpd = 'riborf_chrN/index/' + SPE + '.gpd',
    fa = FA
  output: 'riborf_chrN/index/' + SPE + '/candidateORF.genepred.txt'
  conda:'qc_snakemake'
  params:
    riborf_index = 'riborf_chrN/index/' + SPE
  benchmark: 'time_benchmarks/riborf_chrN/riborf_index/' + SPE + '.txt'
  shell:
    r'''perl {RIBORFSCRIPTS_DIR}/ORFannotate.pl \
    -g {input.fa} \
    -t {input.gpd} \
    -o {params.riborf_index} \
    -s ATG\/CTG\/GTG\/TTG \
    -l 27
    '''

rule riborf_convert_to_sam_chrN:
  input:
    ribo_bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam')
  output:
    ribo_sam = temp('riborf_chrN/temp_sam/{sample}_{mpsf}.sam') 
  threads: 8
  conda:'base'
  shell:
    r'''samtools view -@ {threads} {input.ribo_bam} -h -o {output.ribo_sam}'''

rule generate_bam_stats:
  input:
    bam = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam')
  output:
    stats = os.path.join(BAM_DIR, '{sample}_{mpsf}.stats')
  threads: 8 
  conda:'base'
  shell:
    r'''samtools stats -@ {threads} {input.bam} > {output.stats}'''

rule get_reads_length_dist:
  input: os.path.join(BAM_DIR, '{sample}_{mpsf}.stats')
  output: temp('riborf_chrN/P_site_determination/{sample}_{mpsf}/reads_length.txt')
  script:
    os.path.join(SCRIPTS_DIR, 'riborf_get_reads_lengths.sh')

def read_content(file_path):
    with open(file_path, 'r') as f:
        content = f.read().strip()
    return content

checkpoint riborf_metagene_profile_chrN:
  input:
    ribo_sam = 'riborf_chrN/temp_sam/{sample}_{mpsf}.sam',
    gpd = 'riborf_chrN/index/' + SPE + '.gpd',
    reads_length_file = 'riborf_chrN/P_site_determination/{sample}_{mpsf}/reads_length.txt'
  output: 
    directory('riborf_chrN/checkpoint/{sample}_{mpsf}')
  conda:'qc_snakemake'
  params:
    out_prefix = 'riborf_chrN/checkpoint/{sample}_{mpsf}',
    reads_length = lambda wildcards,input: read_content(input[2])
  benchmark: 'time_benchmarks/riborf_chrN/P_site_determination/{sample}_{mpsf}.txt'
  shell:
    r''' mkdir -p {output} && 
    perl {RIBORFSCRIPTS_DIR}/readDist.pl \
    -f {input.ribo_sam} \
    -g {input.gpd} \
    -o {params.out_prefix} \
    -d '{params.reads_length}' ''' 

rule rename_sta_read_dist:
  input:
    sta_read_dist_file = lambda wildcards:glob.glob(os.path.join(checkpoints.riborf_metagene_profile_chrN.get(**wildcards).output[0], "sta.*.txt"))[0]
  output:
    'riborf_chrN/P_site_determination/{sample}_{mpsf}/read_dist_for_correct.txt'
  shell:
    r'''cp {input.sta_read_dist_file} {output}'''

rule riborf_select_reads:
  input: 'riborf_chrN/P_site_determination/{sample}_{mpsf}/read_dist_for_correct.txt'
  output: 'riborf_chrN/P_site_determination/{sample}_{mpsf}/offset.corretion.parameters.txt'
  conda:'qc_snakemake'
  benchmark: 'time_benchmarks/riborf_chrN/parameterOffset/{sample}_{mpsf}.txt'
  shell:
    r'''perl {RIBORFSCRIPTS_DIR}/parameterOffsetModified.pl -f {input} -o {output}'''

rule offsets_extraction:
  input: 'riborf_chrN/P_site_determination/{sample}_{mpsf}/offset.corretion.parameters.txt'
  output: 
    offsets_for_correct = 'riborf_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_offsets_for_correct.txt',
    standard_offsets = 'riborf_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_standard_offsets.txt'
  script: 
    os.path.join(SCRIPTS_DIR, 'riborf_offsets_extraction.py')

rule riborf_offsets_correct_default:
  input:
    ribo_sam = 'riborf_chrN/temp_sam/{sample}_{mpsf}.sam',
    offsets = 'riborf_chrN/P_site_determination/{sample}_{mpsf}/{sample}_{mpsf}_offsets_for_correct.txt'
  output:
    corrected_sam = 'riborf_chrN/orf_pred_default/{sample}_{mpsf}/corrected.{sample}_{mpsf}.sam'
  conda:'qc_snakemake'
  benchmark: 'time_benchmarks/riborf_chrN/offsetCorrect/{sample}_{mpsf}.txt'
  shell:
    r'''perl {RIBORFSCRIPTS_DIR}/offsetCorrect.pl \
    -r {input.ribo_sam} \
    -p {input.offsets} \
    -o {output.corrected_sam}'''

rule riborf_metagene_profile_corrected_chrN:
  input:
    corrected_sam = 'riborf_chrN/orf_pred_default/{sample}_{mpsf}/corrected.{sample}_{mpsf}.sam',
    gpd = 'riborf_chrN/index/' + SPE + '.gpd'
  output: 'riborf_chrN/P_site_determination/{sample}_{mpsf}/corrected/readDist.plot.1.R'
  conda:'qc_snakemake'
  benchmark: 'time_benchmarks/riborf_chrN/readDist_Correct/{sample}_{mpsf}.txt'
  params:
    out_prefix = 'riborf_chrN/P_site_determination/{sample}_{mpsf}/corrected'
  shell:
    r'''perl {RIBORFSCRIPTS_DIR}/readDist.pl \
    -f {input.corrected_sam} \
    -g {input.gpd} \
    -o {params.out_prefix} \
    -d 1'''

rule riborf_pred_chrN_default:
  input:
    corrected_sam = 'riborf_chrN/orf_pred_default/{sample}_{mpsf}/corrected.{sample}_{mpsf}.sam',
    can_orf = 'riborf_chrN/index/' + SPE + '/candidateORF.genepred.txt'
  output: 'riborf_chrN/orf_pred_default/{sample}_{mpsf}/repre.valid.pred.pvalue.parameters.txt'
  conda:'qc_snakemake'
  params:
    pred_out_prefix = 'riborf_chrN/orf_pred_default/{sample}_{mpsf}'
  benchmark: 'time_benchmarks/riborf_chrN/orf_pred_default/{sample}_{mpsf}.txt'
  shell:
    r'''perl {RIBORFSCRIPTS_DIR}/ribORF.pl \
    -f {input.corrected_sam} \
    -c {input.can_orf} \
    -o {params.pred_out_prefix} \
    -l 27 \
    -r 11 \
    -p 0.7
    '''

rule filter_result:
  input: 'riborf_chrN/orf_pred_default/{sample}_{mpsf}/repre.valid.pred.pvalue.parameters.txt'
  output: 'riborf_chrN/orf_pred_default/{sample}_{mpsf}/riborf_filtered_07.txt'
  script: os.path.join(SCRIPTS_DIR, 'riborf_result_filtering.py')

rule raw_table_result:
  input: 'riborf_chrN/orf_pred_default/{sample}_{mpsf}/riborf_filtered_07.txt'
  output: 'orf_result/raw_prediction_result_default/{sample}_{mpsf}_RibORF_raw.txt'
  script: os.path.join(SCRIPTS_DIR, 'copyfile.py')
