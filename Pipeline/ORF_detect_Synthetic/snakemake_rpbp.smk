import pandas as pd
import os

#################################
# 1. 设置全局参数和样本列表 (Global Configuration)
#################################

SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_simu.tsv"
TEMPLATE_YAML = "/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp/Human.alignments-only.yaml"
PREPARE_GENOME_YAML = "/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp/prepare_genome.yaml"


try:
    SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t")['sample'].tolist()
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含一个名为 'sample' 的列头。错误: {e}")

SPE = 'Human'
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']
GENOME_NAME = "GRCh38"

#################################
# 2. 路径定义 (Path Definitions)
#################################
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_simu/rpbp"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping/simulation"
GTF_FILE = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.no_underscores.gtf"
FASTA_FILE = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect_simu/scripts"
GENOME_BASE_PATH = "/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_index"
BENCHMARK_DIR = os.path.join(OUT_DIR, "time_benchmarks")


PLACEHOLDER_GTF = "/home/chengennong/Ribo_benchmark/ref/Human/gencode.v43.annotation.gtf"
PLACEHOLDER_FASTA = "/home/chengennong/Ribo_benchmark/ref/Human/GRCh38.primary_assembly.genome.fa"
PLACEHOLDER_RIBOSEQ_DATA = "/path/to/your/c-elegans-example"


#################################
# 3. Snakemake 工作流主体
#################################

workdir: OUT_DIR

rule all:
  input:
    expand('rpbp_chrN/configs_default/{sample}_{mpsf}.alignments-only.yaml', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('rpbp_chrN/orf_pred_default/{sample}_{mpsf}/without-rrna-mapping/{sample}_{mpsf}Aligned.sortedByCoord.out.bam', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('rpbp_chrN/orf_pred_default/{sample}_{mpsf}/metagene-profiles/{sample}_{mpsf}.periodic-offsets.csv.gz', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('rpbp_chrN/orf_pred_default/{sample}_{mpsf}/orf-predictions/{sample}_{mpsf}_rpbp_filtered.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('orf_result/raw_prediction_result_default/{sample}_{mpsf}_Rp-Bp_raw.txt', sample=SAMPLES, mpsf=MAPPING_SOFTWARE)
	
rule prepare_genome_indices:
    input:
        config = PREPARE_GENOME_YAML,
        gtf = GTF_FILE,
        fasta = FASTA_FILE
    output:
        annotated_bed = os.path.join(GENOME_BASE_PATH, f"{GENOME_NAME}.annotated.bed.gz")
    benchmark:
        os.path.join(BENCHMARK_DIR, "prepare_genome_indices/prepare_genome.txt")
    log:
        os.path.join(OUT_DIR, "logs/prepare_genome.log")
    conda:'rpbp_env'
    threads: 1 
    shell:
        # --overwrite 确保可以覆盖不完整的文件
        #"prepare-rpbp-genome {input.config} --num-cpus {threads} --log-file {log} --overwrite"
        "prepare-rpbp-genome {input.config} --num-cpus {threads} --log-file {log} --mem 32G"


rule create_default_yaml:
  input:
    template=TEMPLATE_YAML
  params:
    sample_name = '{sample}_{mpsf}',
    riboseq_data_path = os.path.join(OUT_DIR, 'rpbp_chrN/orf_pred_default/{sample}_{mpsf}'),
    new_gtf = GTF_FILE,
    new_fasta = FASTA_FILE
  output:
    'rpbp_chrN/configs_default/{sample}_{mpsf}.alignments-only.yaml'
  benchmark:
    os.path.join(BENCHMARK_DIR, 'create_default_yaml/{sample}_{mpsf}.txt')
  conda:'rpbp_env'
  shell:
    r'''
    sed -e "s/INPUT_SAMPLE/{params.sample_name}/g" \
        -e "s|{PLACEHOLDER_RIBOSEQ_DATA}|{params.riboseq_data_path}|g" \
        -e "s|{PLACEHOLDER_GTF}|{params.new_gtf}|g" \
        -e "s|{PLACEHOLDER_FASTA}|{params.new_fasta}|g" \
        {input.template} > {output}
    '''

rule create_symlink_default:
  input:
    bam_file = os.path.join(BAM_DIR, '{sample}_{mpsf}.bam')
  output:
    'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/without-rrna-mapping/{sample}_{mpsf}Aligned.sortedByCoord.out.bam'
  run:
    if not os.path.exists(input.bam_file):
        raise FileNotFoundError(f"输入BAM文件未找到: {input.bam_file}")
    os.symlink(os.path.abspath(input.bam_file), output[0])


rule run_rpbp_chrN_default:
  input:
    bam = 'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/without-rrna-mapping/{sample}_{mpsf}Aligned.sortedByCoord.out.bam',
    config = 'rpbp_chrN/configs_default/{sample}_{mpsf}.alignments-only.yaml',
    genome_bed = rules.prepare_genome_indices.output.annotated_bed
  output:
    offset = 'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/metagene-profiles/{sample}_{mpsf}.periodic-offsets.csv.gz'
  conda: 'rpbp_env'
  benchmark:
    os.path.join(BENCHMARK_DIR, 'orf_pred_default/{sample}_{mpsf}.txt')
  log: 'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/orf-predictions/{sample}_{mpsf}_rpbp_log.txt'
  threads: 1
  shell:
    r'''run-all-rpbp-instances {input.config} --num-cpus {threads} --logging-level INFO --log-file {log} --keep-intermediate-files --mem 120G'''


rule rename_result:
  input:
    'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/metagene-profiles/{sample}_{mpsf}.periodic-offsets.csv.gz'
  output:
    'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/orf-predictions/{sample}_{mpsf}_rpbp_filtered_raw.txt'
  conda: 'rpbp_env'
  params:
    prediction_dir = 'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/orf-predictions'
  shell:
    r'''zcat {params.prediction_dir}/*.filtered.predicted-orfs.bed.gz > {output}'''

rule filter_result:
  input:
    'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/orf-predictions/{sample}_{mpsf}_rpbp_filtered_raw.txt'
  output:
    'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/orf-predictions/{sample}_{mpsf}_rpbp_filtered.txt'
  conda: 'rpbp_env'
  script:
    os.path.join(SCRIPTS_DIR, 'rpbp_result_filtering.py')

rule raw_table_result:
  input:
    'rpbp_chrN/orf_pred_default/{sample}_{mpsf}/orf-predictions/{sample}_{mpsf}_rpbp_filtered.txt'
  output:
    'orf_result/raw_prediction_result_default/{sample}_{mpsf}_Rp-Bp_raw.txt'
  script:
    os.path.join(SCRIPTS_DIR, 'copyfile.py')
