import os

# ================= 变量定义区域 =================
SAMPLES = [
    "simulation_6M_T1",
    "simulation_6M_T3",
    "simulation_60M_T1",
    "simulation_60M_T3"
]

MERGE = SAMPLES
MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']

# 路径配置
FA = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa"
GPPY_PATH = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_reformat/scripts/gtf.py"
REF_GTF_FOR_GPPY = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation.gtf"

INPUT_ROOT = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_simu_untrim"
OUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/untrim"

workdir: OUT_DIR
SCRIPT_CONVERT = "./scripts/ORFquant_convert.R"
SCRIPT_GTF = "./scripts/gtf.py"
SCRIPT_MERGE = "./scripts/ORFquant_ORFmerge.R"

# ================= 规则定义区域 =================

rule all:
  input:
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_gedi_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotricer_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_riborf_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_ribowave_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_ORFquant_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    #expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_ribohmm_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    #expand('final_ORFs/orf_pred_default/reformatted/{merge}_{mpsf}_ribohmm_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_rpbp_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_orfrater_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotish_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_ribocode_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
    #expand('final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotaper_gcoor.tsv.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE)

# ----------------- GEDI -----------------
rule gedi_addchr:
    input:
        raw_bed = os.path.join(INPUT_ROOT, 'gedi/gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_FDR01.orfs.bed')
    output:
        addchr_bed = 'gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_addchr.bed'
    conda:'r_deseq2'
    threads:1
    script:
        'scripts/gedi_addchr.R'

rule gedi_getfasta:
    input:
        addchr_bed = 'gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_addchr.bed'
    output:
        ORF_sequence_fa = 'gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_ORF_sequence.fa',
        ORF_sequence_tsv = 'gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_ORF_sequence.tsv'
    params:
        ref = FA
    conda:'base'
    threads:1
    shell:
        r'''
            bedtools getfasta -fi {params.ref} \
                -bed {input.addchr_bed} \
                -split -nameOnly -s >{output.ORF_sequence_fa} && \
            seqkit fx2tab {output.ORF_sequence_fa} > {output.ORF_sequence_tsv}
        '''

#rule gedi_ORFmerge:
#    input:
#        addchr_bed = 'gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_addchr.bed',
#        ORF_sequence_tsv = 'gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_ORF_sequence.tsv'
#    output:
#        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_gedi_gcoor.tsv.gz',
#        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_gedi_merged_gcoor.tsv.gz'
#    conda:'r_deseq2'
#    threads:24
#    script:
#        'scripts/gedi_ORFmerge.R'

rule gedi_ORFmerge:
    input:
        addchr_bed = 'gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_addchr.bed',
        ORF_sequence_tsv = 'gedi_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_ORF_sequence.tsv'
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_gedi_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_gedi_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    shell:  # 注意这里由 script 改为了 shell
        '''
        Rscript /home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_reformat/scripts/gedi_ORFmerge.R \
        {input.addchr_bed} \
        {input.ORF_sequence_tsv} \
        {output.formatted_ORF} \
        {output.merged_ORF}
        '''
        
# ----------------- RIBOTRICER -----------------
rule ribotricer_convert:
    input:
        raw_file = os.path.join(INPUT_ROOT, 'ribotricer/ribotricer_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_translating_ORFs.tsv')
    output:
        withseq = 'ribotricer_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_withseq.tsv',
    conda:'r_deseq2'
    threads:1
    script:
        'scripts/ribotricer_convert.R'

rule ribotricer_ORFmerge:
    input:
        withseq = 'ribotricer_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_withseq.tsv',
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotricer_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotricer_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribotricer_ORFmerge.R'

# ----------------- RIBORF -----------------
rule riborf_convert:
    input:
        raw_file = os.path.join(INPUT_ROOT, 'riborf/riborf_chrN/orf_pred_default/{merge}_{mpsf}/repre.valid.pred.pvalue.parameters.txt')
    output:
        withseq = 'riborf_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_withseq.tsv'
    conda:'r_deseq2'
    threads:1
    script:
        'scripts/riborf_convert_modify.R'

rule riborf_tid2gid:
    input:
        withseq = 'riborf_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_withseq.tsv',
    output:
        tid_withseq = 'riborf_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_tid_withseq.tsv',
        tid_noseq = 'riborf_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_tid_noseq.tsv'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/riborf_tid2gid.R'

rule riborf_gppy:
    input:
        tid_noseq = 'riborf_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_tid_noseq.tsv'
    output:
        block = 'riborf_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv',
    conda:'py3.7'
    threads:1
    shell:
        r'''
            tail -n +2 {input.tid_noseq} | \
            python {GPPY_PATH} tiv2giv \
            -g {REF_GTF_FOR_GPPY} \
            -i /dev/stdin -a > {output.block}
        '''
        
rule riborf_ORFmerge:
    input:
        block = 'riborf_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv',
        tid_withseq = 'riborf_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_tid_withseq.tsv'
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_riborf_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_riborf_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/riborf_ORFmerge.R'

# ----------------- RIBOWAVE -----------------
rule ribowave_convert:
    input:
        raw_file = os.path.join(INPUT_ROOT, 'ribowave/ribowave_chrN/orf_pred_default/{merge}_{mpsf}/Ribowave/result/{merge}_{mpsf}.95%.mx')
    output:
        gppy = 'ribowave_chrN/orf_pred_default/{merge}_{mpsf}/Ribowave/result/{merge}_{mpsf}_gppy.tsv'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribowave_convert.R'
#    shell:  # 注意：这里改成了 shell
#        '''
#        # 显式创建目录，防止 R 脚本因为目录不存在而失败
#        mkdir -p $(dirname {output.gppy})
#        
#        Rscript /home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_reformat/scripts/ribowave_convert.R \
#        {input.raw_file} \
#        {output.gppy}
#        '''

rule ribowave_gppy:
    input:
        gppy = 'ribowave_chrN/orf_pred_default/{merge}_{mpsf}/Ribowave/result/{merge}_{mpsf}_gppy.tsv'
    output:
        block = 'ribowave_chrN/orf_pred_default/{merge}_{mpsf}/Ribowave/result/{merge}_{mpsf}_block.tsv'
    conda:'py3.7'
    threads:1
    shell:
        r'''
            python {GPPY_PATH} tiv2giv \
            -g {REF_GTF_FOR_GPPY} \
            -i {input.gppy} -a > {output.block}
        '''

rule ribowave_ORFmerge:
    input:
        block = 'ribowave_chrN/orf_pred_default/{merge}_{mpsf}/Ribowave/result/{merge}_{mpsf}_block.tsv'
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribowave_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribowave_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribowave_ORFmerge.R'

# ----------------- ORFQUANT -----------------
rule ORFquant_convert:
    input:
        raw_file = os.path.join(INPUT_ROOT, 'orfquant/ORFquant_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_result.txt')
    output:
        gppy = 'ORFquant_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ORFquant_convert.R'

rule ORFquant_gppy:
    input:
        gppy = 'ORFquant_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    output:
        block = 'ORFquant_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    conda:'py3.7'
    threads:1
    shell:
        r'''
            python {GPPY_PATH} tiv2giv \
            -g {REF_GTF_FOR_GPPY} \
            -i {input.gppy} -a > {output.block}
        '''

rule ORFquant_ORFmerge:
    input:
        block = 'ORFquant_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ORFquant_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ORFquant_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ORFquant_ORFmerge.R'
    
# ----------------- ORFRATER -----------------
rule orfrater_convert:
    input:
        raw_file = os.path.join(INPUT_ROOT, 'orfrater/orfrater_chrN/orf_pred_default/{merge}_{mpsf}/rate_regression.csv')
    output:
        gppy = 'orfrater_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/orfrater_convert.R'

rule orfrater_gppy:
    input:
        gppy = 'orfrater_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    output:
        block = 'orfrater_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    conda:'py3.7'
    threads:1
    shell:
        r'''
            python {GPPY_PATH} tiv2giv \
            -g {REF_GTF_FOR_GPPY} \
            -i {input.gppy} -a > {output.block}
        '''

rule orfrater_ORFmerge:
    input:
        block = 'orfrater_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_orfrater_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_orfrater_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/orfrater_ORFmerge.R'

# ----------------- RIBOHMM -----------------
rule ribohmm_ORFmerge:
    input:
        infer_CDS = os.path.join(INPUT_ROOT, 'ribohmm/ribohmm_chrN/ORF_detecting_default/{merge}_{mpsf}/infer_CDS.txt')
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribohmm_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribohmm_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribohmm_ORFmerge.R'

rule ribohmm_correct_seq:
    input:
        table = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_{type}_gcoor.tsv.gz'
    output:
        table = 'final_ORFs/orf_pred_default/reformatted/{merge}_{mpsf}_{type}_gcoor.tsv.gz'
    params:
        genome = "/home/tangyuewen/ORF_benchmark/Ref/GRCh38.primary_assembly.genome.fa",
        script = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/scripts/ORF_reformat/scripts/ribohmm_modify.sh"
    conda:'base'
    threads: 4 
    shell:
        """
        bash {params.script} {input.table} {output.table} {params.genome}
        """

# ----------------- RPBP -----------------
#rule rpbp_ORFmerge:
#    input:
#        rpbp = os.path.join(INPUT_ROOT, 'rpbp/rpbp_chrN/orf_pred_default/{merge}_{mpsf}/orf-predictions/{merge}_{mpsf}_rpbp_filtered_raw.txt'),
#        fa = os.path.join(INPUT_ROOT, 'rpbp/rpbp_chrN/orf_pred_default/{merge}_{mpsf}/orf-predictions/orfs.fa')
#    output:
#        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_rpbp_gcoor.tsv.gz',
##        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_rpbp_merged_gcoor.tsv.gz'
#    conda:'r_deseq2'
#    threads:24
#    script:
#        'scripts/rpbp_ORFmerge.R'

rule rpbp_ORFmerge:
    input:
        rpbp = os.path.join(INPUT_ROOT, 'rpbp/rpbp_chrN/orf_pred_default/{merge}_{mpsf}/orf-predictions/{merge}_{mpsf}_rpbp_filtered_raw.txt')
        #fa = os.path.join(INPUT_ROOT, 'rpbp/rpbp_chrN/orf_pred_default/{merge}_{mpsf}/orf-predictions/orfs.fa')
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_rpbp_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_rpbp_merged_gcoor.tsv.gz'
    params:
        fa_dir = os.path.join(INPUT_ROOT, 'rpbp/rpbp_chrN/orf_pred_default/{merge}_{mpsf}/orf-predictions/')
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/rpbp_ORFmerge.R'
        
# ----------------- RIBOTISH -----------------
rule ribotish_convert:
    input:
        raw_file = os.path.join(INPUT_ROOT, 'ribotish/ribotish_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_pred.txt')
    output:
        withseq = 'ribotish_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_withseq.tsv'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribotish_convert.R'

rule ribotish_gppy_pre:
    input:
        withseq = 'ribotish_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_withseq.tsv'        
    output:
        gppy = 'ribotish_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    conda:'r_deseq2'
    threads:8
    script:
        'scripts/ribotish_gppy.R' 

rule ribotish_gppy:
    input:
        gppy = 'ribotish_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    output:
        block = 'ribotish_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    conda:'py3.7'
    threads:1
    shell:
        r'''
            python {GPPY_PATH} tiv2giv \
            -g {REF_GTF_FOR_GPPY} \
            -i {input.gppy} -a > {output.block}
        '''

rule ribotish_ORFmerge:
    input:
        block = 'ribotish_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotish_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotish_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribotish_ORFmerge.R' 

# ----------------- RIBOTAPER -----------------
rule ribotaper_convert:
    input:
        raw_file = os.path.join(INPUT_ROOT, 'ribotaper/ribotaper_chrN/orf_pred_default/{merge}_{mpsf}/ORFs_max')
    output:
        gppy = 'ribotaper_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribotaper_convert.R' 

rule ribotaper_gppy:
    input:
        gppy = 'ribotaper_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    output:
        block = 'ribotaper_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    conda:'py3.7'
    threads:1
    shell:
        r'''
            python {GPPY_PATH} tiv2giv \
            -g {REF_GTF_FOR_GPPY} \
            -i {input.gppy} -a > {output.block}
        '''

rule ribotaper_ORFmerge:
    input:
        block = 'ribotaper_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotaper_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribotaper_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribotaper_ORFmerge.R' 

# ----------------- RIBOCODE -----------------
rule ribocode_convert:
    input:
        raw_file = os.path.join(INPUT_ROOT, 'ribocode/ribocode_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}.txt')
    output:
        withseq = 'ribocode_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_withseq.tsv'
    conda:'r_deseq2'
    threads:1
    script:
        'scripts/ribocode_convert.R'

rule ribocode_gppy_pre:
    input:
        withseq = 'ribocode_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_withseq.tsv'        
    output:
        gppy = 'ribocode_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    conda:'r_deseq2'
    threads:8
    script:
        'scripts/ribocode_gppy.R' 

rule ribocode_gppy:
    input:
        gppy = 'ribocode_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_gppy.tsv'
    output:
        block = 'ribocode_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    conda:'py3.7'
    threads:1
    shell:
        r'''
            python {GPPY_PATH} tiv2giv \
            -g {REF_GTF_FOR_GPPY} \
            -i {input.gppy} -a > {output.block}
        '''

rule ribocode_ORFmerge:
    input:
        block = 'ribocode_chrN/orf_pred_default/{merge}_{mpsf}/{merge}_{mpsf}_block.tsv'
    output:
        formatted_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribocode_gcoor.tsv.gz',
        merged_ORF = 'final_ORFs/orf_pred_default/{merge}_{mpsf}_ribocode_merged_gcoor.tsv.gz'
    conda:'r_deseq2'
    threads:24
    script:
        'scripts/ribocode_ORFmerge.R'
