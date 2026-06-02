include:
    config['config_path']

UBU = '/home/chengennong/Ribo_benchmark/deps/ubu-1.2b'
## UNC Bioinformatics Utilities,sam-xlate
workdir: OUT_DIR

MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']

## this workflow will trim 5'end mismatch
## ambiguous nucleotide (N) will not be trimmed
## trim N in sequence in fastq, and then realign them
## futher more,we need compare position with previous alignment to 
## correct positions

rule all:
    input:
        expand('merge_chrN/{merge}_{mpsf}.bam', merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/{merge}_{mpsf}_plus_trimmed.fq', merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/{merge}_{mpsf}_minus_trimmed.fq', merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/{merge}_{mpsf}_trimmed.fq.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/map_tophat2/{merge}/accepted_hits.bam',merge=MERGE),
        expand('trim_five_prime_mismatch/map_STAR/{merge}/{merge}.Aligned.sortedByCoord.out.bam',merge=MERGE),
        expand('trim_five_prime_mismatch/map_hisat2/{merge}/{merge}.sam',merge=MERGE),
        expand('trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}_trimmed.bam',merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}.bam',merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}.bam.bai',merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/merge_chrN_Tx/{merge}_{mpsf}.bam',merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/merge_chrN_Tx/{merge}_{mpsf}.bam.bai',merge=MERGE, mpsf=MAPPING_SOFTWARE)

## step1:distinguish 5'end mismatched reads 
## Note: trimmed all 5'end nucleotide makes no sense
## it can not help with correcting real offsets

## perfect bam means no mismatch in this alignment
rule extract_five_prime_mismatch_plus:
    input:
        raw_bam = 'merge_chrN/{merge}_{mpsf}.bam'
    output:
        perfect_bam_unsort = temp('trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus_unsort.bam'),
        perfect_bam = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus.bam',
        perfect_bam_index = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus.bam.bai',
        mismatch_bam = temp('trim_five_prime_mismatch/{merge}_{mpsf}_plus.bam')
    conda:'Ribo_benchmark'
    threads: 6
    script:
        'scripts/extract_five_prime_mismatch_plus.sh'

rule extract_five_prime_mismatch_minus:
    input:
        raw_bam = 'merge_chrN/{merge}_{mpsf}.bam'
    output:
        perfect_bam_unsort = temp('trim_five_prime_mismatch/perfect/{merge}_{mpsf}_minus_unsort.bam'),
        perfect_bam = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_minus.bam',
        perfect_bam_index = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_minus.bam.bai',
        mismatch_bam = temp('trim_five_prime_mismatch/{merge}_{mpsf}_minus.bam')
    conda:'Ribo_benchmark'
    threads: 6
    script:
        'scripts/extract_five_prime_mismatch_minus.sh'

## trim_five_prime_mismatch will selecet any alignments with 5'end mismatch

rule bam_to_raw_fastq_plus:
    input:
        bam = 'trim_five_prime_mismatch/{merge}_{mpsf}_plus.bam'
    output:
        raw_fq = temp('trim_five_prime_mismatch/{merge}_{mpsf}_plus.fq')
    conda:'Ribo_benchmark'
    threads: 6
    shell:
        r'''
        samtools fastq -@ {threads} -T MD {input.bam} >{output.raw_fq}
        '''

rule bam_to_raw_fastq_minus:
    input:
        bam = 'trim_five_prime_mismatch/{merge}_{mpsf}_minus.bam'
    output:
        raw_fq = temp('trim_five_prime_mismatch/{merge}_{mpsf}_minus.fq')
    conda:'Ribo_benchmark'
    threads: 6
    shell:
        r'''
        samtools fastq -@ {threads} -T MD {input.bam} >{output.raw_fq}
        '''

## R script for trim mismatched nucleotide
## not trim N, like ANCTG will keep NCTG
## we have to keep N,if not PRICE will throw exception

rule trim_mismatch_plus:
    input:
        raw_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_plus.fq'
    output:
        trimmed_fq = temp('trim_five_prime_mismatch/{merge}_{mpsf}_plus_trimmed.fq')
    conda:'R_plot'
    threads: 6
    script:
        'scripts/trim_mismatch_plus.R'

rule trim_mismatch_minus:
    input:
        raw_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_minus.fq'
    output:
        trimmed_fq = temp('trim_five_prime_mismatch/{merge}_{mpsf}_minus_trimmed.fq')
    conda:'R_plot'
    threads: 6
    script:
        'scripts/trim_mismatch_minus.R'

rule combine_fq:
    input:
        plus = 'trim_five_prime_mismatch/{merge}_{mpsf}_plus_trimmed.fq',
        minus = 'trim_five_prime_mismatch/{merge}_{mpsf}_minus_trimmed.fq'
    output:
        final_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_trimmed.fq.gz'
    shell:
        r'''
            gzip -c {input.plus} {input.minus} > {output.final_fq}
        '''

## need to consider setting chunk
## control memory usage

## realign is the same as previous protocol

rule STAR_mapping:
  input: 
    STARindex = INDEXSTAR,
    retained_reads = 'trim_five_prime_mismatch/{merge}_STAR_trimmed.fq.gz'
  output: 
    STAR_bam = 'trim_five_prime_mismatch/map_STAR/{merge}/{merge}.Aligned.sortedByCoord.out.bam',
    STAR_Tx_bam = 'trim_five_prime_mismatch/map_STAR/{merge}/{merge}.Aligned.toTranscriptome.out.bam'
  params:
    STAR_prefix = 'trim_five_prime_mismatch/map_STAR/{merge}/{merge}.',
    outdir = 'trim_five_prime_mismatch/map_STAR/{merge}'
  conda: 'STAR'
  threads: 12
  shell:
    r''' mkdir -p {params.outdir} && STAR \
    --runThreadN {threads} \
    --genomeDir {input.STARindex} \
    --readFilesCommand zcat \
    --readFilesIn {input.retained_reads} \
    --outFileNamePrefix {params.STAR_prefix} \
    --alignEndsType EndToEnd \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --quantMode TranscriptomeSAM GeneCounts \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 10 '''

rule hisat2_mapping:
  input: 
    retained_reads = 'trim_five_prime_mismatch/{merge}_hisat2_trimmed.fq.gz'
  output: 
    hisat2_sam = 'trim_five_prime_mismatch/map_hisat2/{merge}/{merge}.sam'
  params:
    hisat2index = INDEXHISAT2,
    hisat2_prefix = 'trim_five_prime_mismatch/map_hisat2/{merge}/{merge}'
  log: 'trim_five_prime_mismatch/map_hisat2/{merge}/{merge}.log'
  conda: 'Hisat2'
  threads: 12
  shell:
    r'''hisat2 -t -p {threads} \
    --score-min C,-12 \
    --no-softclip \
    --rna-strandness R \
    --un-gz {params.hisat2_prefix}_unmapped.gz \
    -k 10 \
    -x {params.hisat2index} \
    -U {input.retained_reads} \
    -S {output.hisat2_sam} > {log} 2>&1 '''

## --score-min C,-12 for 2 mismatch,-k 10 for mutilmap 10 maxinum
## --score-min C,-5, no mismatch allowed

rule tophat2_mapping:
  input:
    retained_reads = 'trim_five_prime_mismatch/{merge}_tophat2_trimmed.fq.gz'
  output: 
    tophat2_bam = 'trim_five_prime_mismatch/map_tophat2/{merge}/accepted_hits.bam'
  params:
    tophat2index = INDEXTOPHAT2,
    tophat2index_cache = INDEXTOPHAT2_CACHE,
    tophat2_prefix = 'trim_five_prime_mismatch/map_tophat2/{merge}'
  log: 'trim_five_prime_mismatch/map_tophat2/{merge}/{merge}.Log.final.out'
  conda: 'Tophat2'
  threads: 12
  shell:
    r'''tophat2 -p {threads} \
    --no-coverage-search \
    --read-mismatches 2 \
    --max-multihits 10 \
    --library-type fr-firststrand \
    -o {params.tophat2_prefix} \
    --transcriptome-index={params.tophat2index_cache} \
    {params.tophat2index} \
    {input.retained_reads} > {log} 2>&1'''

## --no-coverage-search to save time and memory

## extract uniquely mapped result of re-alignments

rule STAR_extract_uniquely_mapped:
  input:
    raw_bam = 'trim_five_prime_mismatch/map_STAR/{merge}/{merge}.Aligned.sortedByCoord.out.bam'
  output:
    unique_bam_unsort = temp('trim_five_prime_mismatch/merge_chrN/{merge}_STAR_unsort_trimmed.bam'),
    unique_bam = 'trim_five_prime_mismatch/merge_chrN/{merge}_STAR_trimmed.bam',
    unique_bam_index = 'trim_five_prime_mismatch/merge_chrN/{merge}_STAR_trimmed.bam.bai'
  threads: 6
  conda: 'Ribo_benchmark'
  script:
    'scripts/extract_uniquely_mapped.sh'

rule hisat2_extract_uniquely_mapped:
  input:
    raw_bam = 'trim_five_prime_mismatch/map_hisat2/{merge}/{merge}.sam'
  output:
    unique_bam_unsort = temp('trim_five_prime_mismatch/merge_chrN/{merge}_hisat2_unsort_trimmed.bam'),
    unique_bam = 'trim_five_prime_mismatch/merge_chrN/{merge}_hisat2_trimmed.bam',
    unique_bam_index = 'trim_five_prime_mismatch/merge_chrN/{merge}_hisat2_trimmed.bam.bai'
  threads: 6
  conda: 'Ribo_benchmark'
  script:
    'scripts/extract_uniquely_mapped.sh'
  
rule tophat2_extract_uniquely_mapped:
  input:
    raw_bam = 'trim_five_prime_mismatch/map_tophat2/{merge}/accepted_hits.bam'
  output:
    unique_bam_unsort = temp('trim_five_prime_mismatch/merge_chrN/{merge}_tophat2_unsort_trimmed.bam'),
    unique_bam = 'trim_five_prime_mismatch/merge_chrN/{merge}_tophat2_trimmed.bam',
    unique_bam_index = 'trim_five_prime_mismatch/merge_chrN/{merge}_tophat2_trimmed.bam.bai'
  threads: 6
  conda: 'Ribo_benchmark'
  script:
    'scripts/extract_uniquely_mapped.sh'

rule merge_trimmed_bam:
  input:
    trimmed_bam = 'trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}_trimmed.bam',
    perfect_bam_plus = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus.bam',
    perfect_bam_minus = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_minus.bam'
  output:
    unsort_bam = temp('trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}_unsort.bam'),
    merged_bam = 'trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}.bam',
    merged_bam_index = 'trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}.bam.bai'
  threads: 6
  conda: 'Ribo_benchmark'
  shell:
    r'''
      samtools merge -@ {threads} -o {output.unsort_bam} {input.trimmed_bam} {input.perfect_bam_plus} {input.perfect_bam_minus} && \
      samtools sort -@ {threads} -o {output.merged_bam} {output.unsort_bam} && \
      samtools index -@ {threads} {output.merged_bam}
    '''

# before merge,we need sort and index bam

## we do not have to mapping reads to transcriptome
## use sam-xlate to convert coord after extract uniquely mapped reads 

rule convert_transcriptome_based:
  input:
    gcoor_bam = 'trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}.bam'
  output:
    tcoor_bam_unsort = temp('trim_five_prime_mismatch/merge_chrN_Tx/{merge}_{mpsf}_unsort.bam'),
    tcoor_bam = 'trim_five_prime_mismatch/merge_chrN_Tx/{merge}_{mpsf}.bam',
    tcoor_bam_index = 'trim_five_prime_mismatch/merge_chrN_Tx/{merge}_{mpsf}.bam.bai'
  resources:
    mem_mb = 35840
  params:
    ubu = UBU,
    bed = '/home/chengennong/Ribo_benchmark/ref/Human/gencode.v43.annotation.bed12'
  conda: 'Ribo_benchmark'
  threads: 6
  shell:
    r'''
      java -Xmx35G -jar {params.ubu}/ubu-1.2b-SNAPSHOT-jar-with-dependencies.jar sam-xlate \
        --single \
        --reverse \
        --bed {params.bed} \
        --in {input.gcoor_bam} \
        --out {output.tcoor_bam_unsort} && \
      samtools sort -@ {threads} -o {output.tcoor_bam} {output.tcoor_bam_unsort} && \
      samtools index -@ {threads} -o {output.tcoor_bam_index} {output.tcoor_bam}
    '''
  
## absolute path templiy
## convert gtf to bed, needed to polish
## 35Gb for sam-xlate
## maybe we can warp it later
