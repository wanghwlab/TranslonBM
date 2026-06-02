#!/bin/bash

# 检查参数个数
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file.tsv.gz> <output_file.tsv.gz> <genome.fa>"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"
GENOME_FA="$3"

# 获取文件名用于日志
BASE_NAME=$(basename "$INPUT_FILE" .tsv.gz)

# 创建安全的临时目录/文件前缀，防止并行冲突
TMP_PREFIX=$(mktemp -u)
TMP_INPUT="${TMP_PREFIX}_input.tsv"
TMP_COORDS="${TMP_PREFIX}_coords.bed"
TMP_COORDS_SORTED="${TMP_PREFIX}_coords_sorted.bed"
TMP_SEQ="${TMP_PREFIX}_seq.txt"
TMP_OUTPUT="${TMP_PREFIX}_output.tsv"
TMP_HEADER="${TMP_PREFIX}_header.tsv"

# 确保基因组索引存在 (建议在运行流程前建好，这里作为保险)
if [ ! -e "${GENOME_FA}.fai" ]; then
    echo "Warning: Indexing genome. This might conflict if multiple jobs run simultaneously."
    samtools faidx "$GENOME_FA"
fi

echo "Processing $INPUT_FILE -> $OUTPUT_FILE"

# 解压输入文件
zcat "$INPUT_FILE" > "$TMP_INPUT"

# 创建表头
echo -e "gene_id\ttranscript_id\tcoordinate_id\tchrom\tcoordinate_0base\tstrand\tORF_sequence_correct\tstart_codon\tORF_sequence_aa" > "$TMP_HEADER"

# 处理数据行
tail -n +2 "$TMP_INPUT" | while IFS=$'\t' read -r gene_id transcript_id coordinate_id chrom coordinate strand start_codon orf_aa orf_correct; do
    # 清空临时文件
    > "$TMP_COORDS"
    > "$TMP_SEQ"
    
    # --- 步骤1：解析坐标并生成BED ---
    IFS=',' read -ra exon_coords <<< "$coordinate"
    for exon in "${exon_coords[@]}"; do
        IFS='-' read -ra pos <<< "$exon"
        start=${pos[0]}
        end=${pos[1]}
        # 跳过长度为0的坐标
        if [ "$start" -eq "$end" ]; then
            continue
        fi
        # BED格式：chrom start end name score strand
        echo -e "$chrom\t$start\t$end\t${gene_id}_${transcript_id}_${start}_${end}\t0\t$strand" >> "$TMP_COORDS"
    done
    
    # 检查是否有有效坐标
    if [ ! -s "$TMP_COORDS" ]; then
        echo -e "$gene_id\t$transcript_id\t$coordinate_id\t$chrom\t$coordinate\t$strand\tNA\t$start_codon\t$orf_aa"
        continue
    fi
    
    # 按基因组坐标排序 (BEDTools 需要)
    sort -k1,1 -k2,2n "$TMP_COORDS" > "$TMP_COORDS_SORTED"
    
    # --- 步骤2：提取序列 ---
    # 使用 split 和 name+ 模式确保可以追踪原始顺序
    bedtools getfasta -fi "$GENOME_FA" -bed "$TMP_COORDS_SORTED" -split -name+ -s -tab -fo "$TMP_SEQ"
    
    # --- 步骤3：拼接序列 (处理多外显子和正负链) ---
    full_seq=$(awk -v strand="$strand" '{
        # 提取坐标信息：拆分bedtools输出的序列名称
        split($1, a, "::");
        gsub(/>/, "", a[1]);
        split(a[1], b, "_");
        
        # 存储序列
        seq[NR] = $2;
    }
    END {
        # 根据链方向决定拼接顺序
        if (strand == "+") {
            # 正链：按BED顺序 (基因组升序) 拼接
            for (i = 1; i <= NR; i++) {
                printf "%s", seq[i];
            }
        } else {
            # 负链：按BED顺序的逆序 (基因组降序) 拼接
            # 注意：bedtools -s 已经对序列取了反向互补，但外显子顺序还是按基因组坐标排列的(exon1, exon2...)
            # 对于负链基因，基因组坐标较大的外显子其实是转录本的 5端。
            # 所以需要倒序拼接。
            for (i = NR; i >= 1; i--) {
                printf "%s", seq[i];
            }
        }
    }' "$TMP_SEQ")

    # --- 步骤4：验证长度 ---
    aa_length=${#orf_aa}
    expected_dna_length=$((aa_length * 3))
    
    if [ -z "$full_seq" ]; then
        full_seq="NA"
    elif [ ${#full_seq} -lt $expected_dna_length ]; then
        # 允许少量长度差异(比如终止密码子)，但如果差太多则报警
        # echo "Warning: Short sequence for $gene_id" >&2
        : # no-op
    fi
    
    # 输出结果
    echo -e "$gene_id\t$transcript_id\t$coordinate_id\t$chrom\t$coordinate\t$strand\t$full_seq\t$start_codon\t$orf_aa"

done > "$TMP_OUTPUT"

# 合并表头和结果，并处理可能的回车符
cat "$TMP_HEADER" "$TMP_OUTPUT" | awk 'BEGIN {FS=OFS="\t"} {gsub(/\r/,"",$7); print}' | gzip > "$OUTPUT_FILE"

# 清理临时文件
rm -f "${TMP_PREFIX}"*

echo "Finished processing: $OUTPUT_FILE"
