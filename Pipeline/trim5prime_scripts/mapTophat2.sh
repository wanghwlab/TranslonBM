#!/bin/bash
TopHat2index="/home/tangyuewen/ORF_benchmark/Ref/index_tophat2/GRCh38"  
TopHat2index_cache="/home/tangyuewen/ORF_benchmark/Ref/index_tophat2/transcriptome_data/known"  
threads=12
input_dir="trim_five_prime_mismatch"

trimmed_files=($(ls -d ${input_dir}/S*_tophat2_trimmed.fq.gz 2>/dev/null))

if [ ${#trimmed_files[@]} -eq 0 ]; then
    echo "Error: No tophat2_trimmed files found in ${input_dir}/"
    exit 1
fi

echo "Found ${#trimmed_files[@]} tophat2_trimmed files to process"

for trimmed_file in "${trimmed_files[@]}"; do
    merge=$(basename ${trimmed_file} | sed 's/_tophat2_trimmed.fq.gz//')
    outdir="${input_dir}/map_tophat2/${merge}"
    merge_chrN_dir="${input_dir}/merge_chrN"
    
    echo "Processing sample: ${merge}"
    echo "Input file: ${trimmed_file}"
    echo "Output directory: ${outdir}"
    
    mkdir -p ${outdir}
    mkdir -p ${merge_chrN_dir}

    /home/tangyuewen/software/ORFbenchmark/tophat-2.1.1.Linux_x86_64/tophat2 -p ${threads} \
        --no-coverage-search \
        --read-mismatches 2 \
        --max-multihits 10 \
        --library-type fr-firststrand \
        -o ${outdir} \
        --transcriptome-index=${TopHat2index_cache} \
        ${TopHat2index} \
        ${trimmed_file} > ${outdir}/${merge}.Log.final.out 2>&1

    if [ $? -ne 0 ]; then
        echo "Error: TopHat2 mapping failed for ${merge}"
        continue
    fi

    raw_bam="${outdir}/accepted_hits.bam"
    sorted_bam="${outdir}/${merge}.sorted.bam"
    
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ ${threads} -o ${sorted_bam} ${raw_bam}

    unique_bam_unsort="${outdir}/${merge}.unique.unsorted.bam"
    unique_bam="${outdir}/${merge}.unique.bam"
    unique_bam_index="${unique_bam}.bai"

    echo "Extracting unique mapped reads for ${merge}"
    
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ ${threads} -bh \
        -e '[NH]==1 && rname =~ "chr[^M]+"' \
        -o "${unique_bam_unsort}" "${sorted_bam}" && \
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ ${threads} -o "${unique_bam}" "${unique_bam_unsort}" && \
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools index -@ ${threads} -o "${unique_bam_index}" "${unique_bam}"
 
    if [ $? -ne 0 ]; then
        echo "Error: Unique read extraction failed for ${merge}"
        continue
    fi

    echo "Merging BAM files for ${merge}"
    
    mpsf="tophat2"
    
    trimmed_bam="${unique_bam}"
    perfect_bam_plus="${input_dir}/perfect/${merge}_${mpsf}_plus.bam"
    perfect_bam_minus="${input_dir}/perfect/${merge}_${mpsf}_minus.bam"
    
    unsort_bam="${merge_chrN_dir}/${merge}_${mpsf}_unsort.bam"
    merged_bam="${merge_chrN_dir}/${merge}_${mpsf}.bam"
    merged_bam_index="${merged_bam}.bai"
    
    if [[ ! -f ${perfect_bam_plus} || ! -f ${perfect_bam_minus} ]]; then
        echo "Warning: Perfect BAM files not found for ${merge}, skipping merge step"
        echo "Expected files:"
        echo "  ${perfect_bam_plus}"
        echo "  ${perfect_bam_minus}"
        continue
    fi

    ~/miniconda3/envs/STAR2.7.10b/bin/samtools merge -@ ${threads} -o ${unsort_bam} \
        ${trimmed_bam} ${perfect_bam_plus} ${perfect_bam_minus} && \
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ ${threads} -bh \
        -e '[NH]==1 && rname =~ "chr[^M]+"' \
        -o ${merged_bam}.tmp.bam ${unsort_bam} && \
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ ${threads} \
        -o ${merged_bam} ${merged_bam}.tmp.bam && \
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools index -@ ${threads} ${merged_bam} && \
    rm ${merged_bam}.tmp.bam ${unsort_bam}

    if [ $? -eq 0 ]; then
        echo "Successfully merged BAM files for ${merge}"
    else
        echo "Error: BAM merging failed for ${merge}"
    fi
    
    echo "----------------------------------------"
done

echo "All TopHat2 processing completed"
