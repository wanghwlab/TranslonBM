#!/bin/bash
STARindex="/home/tangyuewen/ORF_benchmark/Ref/index_STAR/"   
threads=12   
input_dir="trim_five_prime_mismatch"   
merge_chrN_dir="${input_dir}/merge_chrN" 

trimmed_files=($(ls -d ${input_dir}/S*_STAR_trimmed.fq.gz 2>/dev/null)) # 

if [ ${#trimmed_files[@]} -eq 0 ]; then #
    echo "Error: No STAR_trimmed files found in ${input_dir}/" # 
    exit 1 # 
fi

echo "Found ${#trimmed_files[@]} STAR_trimmed files to process" # 

for trimmed_file in "${trimmed_files[@]}"; do # 
    merge=$(basename ${trimmed_file} | sed 's/_STAR_trimmed.fq.gz//') # 
    outdir="${input_dir}/map_STAR/${merge}" # 
    
    echo "Processing sample: ${merge}" # 
    echo "Input file: ${trimmed_file}" # 
    echo "Output directory: ${outdir}" # 
    
    mkdir -p ${outdir} # 
    mkdir -p ${merge_chrN_dir} 

    raw_bam="${outdir}/${merge}.Aligned.sortedByCoord.out.bam"  
    unique_bam="${outdir}/${merge}.unique.bam" 
    unique_bam_index="${unique_bam}.bai" 

    if [[ -f "${raw_bam}" && -f "${unique_bam}" && -f "${unique_bam_index}" ]]; then
        echo "Core output BAM files already exist for ${merge}. Skipping STAR mapping and unique read extraction."
    else
        echo "Core output BAM files not found for ${merge}. Performing STAR mapping and unique read extraction."
        ~/miniconda3/envs/STAR2.7.10b/bin/STAR \
            --runThreadN ${threads} \
            --genomeDir ${STARindex} \
            --readFilesCommand zcat \
            --readFilesIn ${trimmed_file} \
            --outFileNamePrefix ${outdir}/${merge}. \
            --alignEndsType EndToEnd \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --quantMode TranscriptomeSAM GeneCounts \
            --outReadsUnmapped Fastx \
            --outFilterMismatchNmax 2 \
            --outFilterMultimapNmax 10 # 

        if [ $? -ne 0 ]; then # 
            echo "Error: STAR mapping failed for ${merge}" # 
            continue
        fi

        echo "Extracting unique mapped reads for ${merge}" # 
        
        unique_bam_unsort="${outdir}/${merge}.unique.unsorted.bam" 
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ ${threads} -bh \
            -e '[NH]==1 && rname =~ "chr[^M]+"' \
            -o "${unique_bam_unsort}" "${raw_bam}" && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ ${threads} -o "${unique_bam}" "${unique_bam_unsort}" && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools index -@ ${threads} -o "${unique_bam_index}" "${unique_bam}" # 

        if [ $? -eq 0 ]; then # 
            echo "Successfully processed ${merge}" # 
            rm "${unique_bam_unsort}" 
        else
            echo "Error: Unique read extraction failed for ${merge}" # 
            continue
        fi
    fi

    echo "Merging BAM files for ${merge}" # 
    
    mpsf="STAR" 
    
    trimmed_bam="${unique_bam}"   
    perfect_bam_plus="${input_dir}/perfect/${merge}_${mpsf}_plus.bam" # 
    perfect_bam_minus="${input_dir}/perfect/${merge}_${mpsf}_minus.bam" # 
    
    unsort_bam="${merge_chrN_dir}/${merge}_${mpsf}_unsort.bam" # 
    merged_bam="${merge_chrN_dir}/${merge}_${mpsf}.bam" # 
    merged_bam_index="${merged_bam}.bai" # 
    
    if [[ ! -f ${perfect_bam_plus} || ! -f ${perfect_bam_minus} ]]; then # 
        echo "Warning: Perfect BAM files not found for ${merge}, skipping merge step" # 
        echo "Expected files:" # 
        echo "  ${perfect_bam_plus}" # 
        echo "  ${perfect_bam_minus}" # 
        echo "----------------------------------------"
        continue
    fi

    if [[ -f "${merged_bam}" && -f "${merged_bam_index}" ]]; then
        echo "Merged BAM file and index already exist for ${merge}. Skipping merge step."
    else
        echo "Performing merge operation for ${merge}."
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools merge -@ ${threads} -o ${unsort_bam} \
        ${trimmed_bam} ${perfect_bam_plus} ${perfect_bam_minus} && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ ${threads} -bh \
        -e '[NH]==1 && rname =~ "chr[^M]+"' \
        -o ${merged_bam}.tmp.bam ${unsort_bam} && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ ${threads} \
        -o ${merged_bam} ${merged_bam}.tmp.bam && \
        ~/miniconda3/envs/STAR2.7.10b/bin/samtools index -@ ${threads} ${merged_bam} && \
        rm ${merged_bam}.tmp.bam ${unsort_bam} # 

        if [ $? -eq 0 ]; then # 
            echo "Successfully merged BAM files for ${merge}" # 
        else
            echo "Error: BAM merging failed for ${merge}" # 
        fi
    fi
    
    echo "----------------------------------------"
done

echo "All STAR processing completed" #
