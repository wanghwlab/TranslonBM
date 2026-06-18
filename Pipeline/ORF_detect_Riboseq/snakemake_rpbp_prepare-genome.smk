
import pandas as pd
import os

SAMPLE_SHEET = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/configs/sample_names_specific.tsv"
TEMPLATE_YAML = "/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp/Human.alignments-only.yaml"
PREPARE_GENOME_YAML = "/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/prepare_genome.yaml"


try:
    SAMPLES = pd.read_csv(SAMPLE_SHEET, sep="\t")['sample'].tolist()
except Exception as e:
    raise ValueError(f"无法读取或解析样本文件: {SAMPLE_SHEET}. 请确保它是一个tab分割的文件，且包含一个名为 'sample' 的列头。错误: {e}")

SPE = 'Human'
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']
GENOME_NAME = "GRCh38"


OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/rpbp"
BAM_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping_trim5prime/merge_chrN"
GTF_FILE = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.no_underscores.gtf"
FASTA_FILE = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
SCRIPTS_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_detect/scripts"
GENOME_BASE_PATH = "/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome"
BENCHMARK_DIR = os.path.join(OUT_DIR, "time_benchmarks")

PLACEHOLDER_GTF = "/home/chengennong/Ribo_benchmark/ref/Human/gencode.v43.annotation.gtf"
PLACEHOLDER_FASTA = "/home/chengennong/Ribo_benchmark/ref/Human/GRCh38.primary_assembly.genome.fa"
PLACEHOLDER_RIBOSEQ_DATA = "/path/to/your/c-elegans-example"

workdir: OUT_DIR

rule all:
  input:
    expand('rpbp_chrN/configs_default/{sample}_{mpsf}.alignments-only.yaml', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/GRCh38.annotated.bed.gz', sample=SAMPLES, mpsf=MAPPING_SOFTWARE),
    expand('/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/transcript-index/GRCh38.orfs-exons.bed.gz'),
    expand('/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/transcript-index/GRCh38.orfs-genomic.bed.gz'),
    expand('/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/transcript-index/GRCh38.orfs-labels.tab.gz'),
    expand('/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/transcript-index/GRCh38.transcripts.annotated.fa')
	
rule prepare_genome_indices:
    input:
        config = PREPARE_GENOME_YAML,
        gtf = GTF_FILE,
        fasta = FASTA_FILE
    output:
        annotated_bed = os.path.join(GENOME_BASE_PATH, f"{GENOME_NAME}.annotated.bed.gz"),
        trans_bed = "/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/transcript-index/GRCh38.orfs-exons.bed.gz",
        labels_bed = '/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/transcript-index/GRCh38.transcripts.annotated.fa',
        tab_bed = '/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/transcript-index/GRCh38.orfs-labels.tab.gz',
        orf_bed = '/home/tangyuewen/ORF_benchmark/Ref/ORFtools/rpbp_genome/transcript-index/GRCh38.orfs-genomic.bed.gz'
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


