MAPPING_SOFTWARE = ['STAR', 'hisat2', 'tophat2']
MERGE = ["SRX1254413","SRX740748","SRX5256553_SRX5256554","SRX5256543_SRX5256555","SRX5887328_SRX5887329_SRX5887330","SRX11812007_SRX11812008_SRX11812009","SRX7666669-73","SRX7666674-78","SRX7666679-83","SRX7666684-88","SRX7666689-93","SRX7666694-98"] 
#"SRX7666679-83" "SRX7666684-88" "SRX7666689-93" "SRX7666694-98")  # 可以添加更多SRX号
#"SRX1447296" "SRX1254413" "SRX740748" "SRX876063_SRX876069" "SRX5256553_SRX5256554" "SRX5887328_SRX5887329_SRX5887330" "SRX11812007_SRX11812008_SRX11812009"  "SRX7666669-73" "SRX7666674-78"
merge_chrN = "/home/tangyuewen/ORF_benchmark/Bam/ignoreNbase/UniqueBam_MDdelNnase"

## this workflow will trim 5'end mismatch
## ambiguous nucleotide (N) will not be trimmed
## trim N in sequence in fastq, and then realign them
## futher more,we need compare position with previous alignment to 
## correct positions
rule all:
    input:
        expand(merge_chrN + '/{merge}_{mpsf}.bam', merge=MERGE, mpsf=MAPPING_SOFTWARE),
	expand('trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus.bam', merge=MERGE, mpsf=MAPPING_SOFTWARE),
	expand('trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus.bam.bai', merge=MERGE, mpsf=MAPPING_SOFTWARE),
	expand('trim_five_prime_mismatch/perfect/{merge}_{mpsf}_minus.bam.bai', merge=MERGE, mpsf=MAPPING_SOFTWARE),
	expand('trim_five_prime_mismatch/{merge}_{mpsf}_plus.fq', merge=MERGE, mpsf=MAPPING_SOFTWARE),
	expand('trim_five_prime_mismatch/{merge}_{mpsf}_minus.fq', merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/{merge}_{mpsf}_plus_trimmed.fq', merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/{merge}_{mpsf}_minus_trimmed.fq', merge=MERGE, mpsf=MAPPING_SOFTWARE),
        expand('trim_five_prime_mismatch/{merge}_{mpsf}_trimmed.fq.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE)
        #expand('trim_five_prime_mismatch/{merge}_{mpsf}_trimmed.fq.gz', merge=MERGE, mpsf=MAPPING_SOFTWARE),
        #expand('trim_five_prime_mismatch/map_tophat2/{merge}/accepted_hits.bam',merge=MERGE),
        #expand('trim_five_prime_mismatch/map_STAR/{merge}/{merge}.Aligned.sortedByCoord.out.bam',merge=MERGE),
        #expand('trim_five_prime_mismatch/map_hisat2/{merge}/{merge}.sam',merge=MERGE),
        #expand('trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}_trimmed.bam',merge=MERGE, mpsf=MAPPING_SOFTWARE),
        #expand('trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}.bam',merge=MERGE, mpsf=MAPPING_SOFTWARE),
        #expand('trim_five_prime_mismatch/merge_chrN/{merge}_{mpsf}.bam.bai',merge=MERGE, mpsf=MAPPING_SOFTWARE),
        #expand('trim_five_prime_mismatch/merge_chrN_Tx/{merge}_{mpsf}.bam',merge=MERGE, mpsf=MAPPING_SOFTWARE),
        #expand('trim_five_prime_mismatch/merge_chrN_Tx/{merge}_{mpsf}.bam.bai',merge=MERGE, mpsf=MAPPING_SOFTWARE)

## step1:distinguish 5'end mismatched reads 
## Note: trimmed all 5'end nucleotide makes no sense
## it can not help with correcting real offsets

## perfect bam means no mismatch in this alignment
rule extract_five_prime_mismatch_plus:
    input:
        raw_bam = merge_chrN + '/{merge}_{mpsf}.bam'
    output:
        perfect_bam_unsort = temp('trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus_unsort.bam'),
        perfect_bam = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus.bam',
        perfect_bam_index = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_plus.bam.bai',
        mismatch_bam = 'trim_five_prime_mismatch/{merge}_{mpsf}_plus.bam'
    conda:'STAR2.7.10b'
    threads: 6
    shell:
        """
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ {threads} -bh -F 16 \
        -e '[MD] =~ "^0[ATCG].*"' \
        -o {output.mismatch_bam} {input.raw_bam} && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ {threads} -bh -F 16 \
        -e '[MD] !~ "^0[ATCG].*"' \
        -o {output.perfect_bam_unsort} {input.raw_bam} && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ {threads} \
        -o {output.perfect_bam} {output.perfect_bam_unsort} && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools index {output.perfect_bam}
        """

rule extract_five_prime_mismatch_minus:
    input:
        raw_bam = merge_chrN + '/{merge}_{mpsf}.bam'
    output:
        perfect_bam_unsort = temp('trim_five_prime_mismatch/perfect/{merge}_{mpsf}_minus_unsort.bam'),
        perfect_bam = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_minus.bam',
        perfect_bam_index = 'trim_five_prime_mismatch/perfect/{merge}_{mpsf}_minus.bam.bai',
        mismatch_bam = 'trim_five_prime_mismatch/{merge}_{mpsf}_minus.bam'
    conda:'STAR2.7.10b'
    threads: 6
    shell:
        """
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ {threads} -bh -f 16 \
        -e '[MD] =~ ".*[ATCG]0$"' \
        -o {output.mismatch_bam} {input.raw_bam} && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ {threads} -bh -f 16 \
        -e '[MD] !~ ".*[ATCG]0$"' \
        -o {output.perfect_bam_unsort} {input.raw_bam} && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ {threads} \
        -o {output.perfect_bam} {output.perfect_bam_unsort} && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools index {output.perfect_bam}
        """



## trim_five_prime_mismatch will selecet any alignments with 5'end mismatch

rule bam_to_raw_fastq_plus:
    input:
        bam = 'trim_five_prime_mismatch/{merge}_{mpsf}_plus.bam'
    output:
        raw_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_plus.fq'
    conda:'STAR2.7.10b'
    threads: 6
    shell:
        r'''
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools fastq -@ {threads} -T MD {input.bam} >{output.raw_fq}
        '''

rule bam_to_raw_fastq_minus:
    input:
        bam = 'trim_five_prime_mismatch/{merge}_{mpsf}_minus.bam'
    output:
        raw_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_minus.fq'
    conda:'STAR2.7.10b'
    threads: 6
    shell:
        r'''
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools fastq -@ {threads} -T MD {input.bam} >{output.raw_fq}
        '''

## R script for trim mismatched nucleotide
## not trim N, like ANCTG will keep NCTG
## we have to keep N,if not PRICE will throw exception

rule trim_mismatch_plus:
    input:
        raw_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_plus.fq'
    output:
        trimmed_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_plus_trimmed.fq'
    conda:'R_plot'
    threads: 6
    shell:
        '~/miniconda3/envs/r_deseq2/bin/Rscript scripts/trim_mismatch_plus.R'

rule trim_mismatch_minus:
    input:
        raw_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_minus.fq'
    output:
        trimmed_fq = 'trim_five_prime_mismatch/{merge}_{mpsf}_minus_trimmed.fq'
    conda:'R_plot'
    threads: 6
    shell:
        '~/miniconda3/envs/r_deseq2/bin/Rscript scripts/trim_mismatch_minus.R'

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
