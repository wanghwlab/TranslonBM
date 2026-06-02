#!/bin/bash
HISAT2index="/home/tangyuewen/ORF_benchmark/Ref/index_hisat2/GRCh38"  
threads=12
input_dir="trim_five_prime_mismatch"

trimmed_files=($(ls -d ${input_dir}/SRX125441*_hisat2_trimmed.fq.gz 2>/dev/null))

if [ ${#trimmed_files[@]} -eq 0 ]; then
    echo "Error: No hisat2_trimmed files found in ${input_dir}/"
    exit 1
fi

echo "Found ${#trimmed_files[@]} hisat2_trimmed files to process"

for trimmed_file in "${trimmed_files[@]}"; do
    merge=$(basename ${trimmed_file} | sed 's/_hisat2_trimmed.fq.gz//')
    outdir="${input_dir}/map_hisat2/${merge}"
    merge_chrN_dir="${input_dir}/merge_chrN" 
    
    echo "Processing sample: ${merge}"
    echo "Input file: ${trimmed_file}"
    echo "Output directory: ${outdir}"
    
    mkdir -p ${outdir}
    mkdir -p ${merge_chrN_dir}  

    /home/tangyuewen/software/ORFbenchmark/hisat2-2.2.1/hisat2 -t -p ${threads} \
	--mp 6,6 --score-min C,-12,0 \
        --no-softclip \
        --rna-strandness R \
        --un-gz ${outdir}/${merge}_unmapped.gz \
        -k 10 \
        -x ${HISAT2index} \
        -U ${trimmed_file} \
        -S ${outdir}/${merge}.sam > ${outdir}/${merge}.log 2>&1

    if [ $? -ne 0 ]; then
        echo "Error: HISAT2 mapping failed for ${merge}"
        continue
    fi

    raw_bam="${outdir}/${merge}.sorted.bam"
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ ${threads} -bh -o ${outdir}/${merge}.bam ${outdir}/${merge}.sam
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ ${threads} -o ${raw_bam} ${outdir}/${merge}.bam

    unique_bam_unsort="${outdir}/${merge}.unique.unsorted.bam"
    unique_bam="${outdir}/${merge}.unique.bam"
    unique_bam_index="${unique_bam}.bai"

    echo "Extracting unique mapped reads for ${merge}"
    
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools view -@ ${threads} -bh \
        -e '[NH]==1 && rname =~ "chr[^M]+"' \
        -o "${unique_bam_unsort}" "${raw_bam}" && \
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools sort -@ ${threads} -o "${unique_bam}" "${unique_bam_unsort}" && \
    ~/miniconda3/envs/STAR2.7.10b/bin/samtools index -@ ${threads} -o "${unique_bam_index}" "${unique_bam}"
 
    
    if [ $? -ne 0 ]; then
        echo "Error: Unique read extraction failed for ${merge}"
        continue
    fi

    rm ${outdir}/${merge}.sam ${outdir}/${merge}.bam

    echo "Merging BAM files for ${merge}"
    
    mpsf="hisat2"
    
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
        rm ${unsort_bam}
    else
        echo "Error: BAM merging failed for ${merge}"
    fi
    
    echo "----------------------------------------"
done

echo "All HISAT2 processing completed"
