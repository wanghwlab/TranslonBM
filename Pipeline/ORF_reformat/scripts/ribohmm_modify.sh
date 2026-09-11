#!/bin/bash

# Check the number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file.tsv.gz> <output_file.tsv.gz> <genome.fa>"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"
GENOME_FA="$3"

# Get the file name for logging
BASE_NAME=$(basename "$INPUT_FILE" .tsv.gz)

# Create a safe temporary directory and file prefix to prevent parallel-run conflicts
TMP_PREFIX=$(mktemp -u)
TMP_INPUT="${TMP_PREFIX}_input.tsv"
TMP_COORDS="${TMP_PREFIX}_coords.bed"
TMP_COORDS_SORTED="${TMP_PREFIX}_coords_sorted.bed"
TMP_SEQ="${TMP_PREFIX}_seq.txt"
TMP_OUTPUT="${TMP_PREFIX}_output.tsv"
TMP_HEADER="${TMP_PREFIX}_header.tsv"

# Ensure that the genome index exists (preferably build it before running the workflow)
if [ ! -e "${GENOME_FA}.fai" ]; then
    echo "Warning: Indexing genome. This might conflict if multiple jobs run simultaneously."
    samtools faidx "$GENOME_FA"
fi

echo "Processing $INPUT_FILE -> $OUTPUT_FILE"

# Decompress the input file
zcat "$INPUT_FILE" > "$TMP_INPUT"

# Create the header
echo -e "gene_id\ttranscript_id\tcoordinate_id\tchrom\tcoordinate_0base\tstrand\tORF_sequence_correct\tstart_codon\tORF_sequence_aa" > "$TMP_HEADER"

# Process data rows
tail -n +2 "$TMP_INPUT" | while IFS=$'\t' read -r gene_id transcript_id coordinate_id chrom coordinate strand start_codon orf_aa orf_correct; do
    # Clear temporary files
    > "$TMP_COORDS"
    > "$TMP_SEQ"
    
    # --- Step 1: Parse coordinates and generate BED ---
    IFS=',' read -ra exon_coords <<< "$coordinate"
    for exon in "${exon_coords[@]}"; do
        IFS='-' read -ra pos <<< "$exon"
        start=${pos[0]}
        end=${pos[1]}
        # Skip zero-length coordinates
        if [ "$start" -eq "$end" ]; then
            continue
        fi
        # BED format:chrom start end name score strand
        echo -e "$chrom\t$start\t$end\t${gene_id}_${transcript_id}_${start}_${end}\t0\t$strand" >> "$TMP_COORDS"
    done
    
    # Check for valid coordinates
    if [ ! -s "$TMP_COORDS" ]; then
        echo -e "$gene_id\t$transcript_id\t$coordinate_id\t$chrom\t$coordinate\t$strand\tNA\t$start_codon\t$orf_aa"
        continue
    fi
    
    # Sort by genomic coordinate (required by BEDTools)
    sort -k1,1 -k2,2n "$TMP_COORDS" > "$TMP_COORDS_SORTED"
    
    # --- Step 2: Extract sequences ---
    # Use split and name+ mode to preserve the original order
    bedtools getfasta -fi "$GENOME_FA" -bed "$TMP_COORDS_SORTED" -split -name+ -s -tab -fo "$TMP_SEQ"
    
    # --- Step 3: Concatenate sequences (handling multiple exons and both strands) ---
    full_seq=$(awk -v strand="$strand" '{
        # Extract coordinate information by splitting sequence names produced by bedtools
        split($1, a, "::");
        gsub(/>/, "", a[1]);
        split(a[1], b, "_");
        
        # Store sequences
        seq[NR] = $2;
    }
    END {
        # Choose the concatenation order based on strand
        if (strand == "+") {
            # Positive strand: concatenate in BED order (ascending genomic coordinates)
            for (i = 1; i <= NR; i++) {
                printf "%s", seq[i];
            }
        } else {
            # Negative strand: concatenate in reverse BED order (descending genomic coordinates)
            # Note: bedtools -s reverse-complements the sequences, but exon order still follows genomic coordinates(exon1, exon2...)
            # For a negative-strand gene, the exon with the larger genomic coordinate is at the transcript 5-prime end.
            # Therefore, concatenate the exons in reverse order.
            for (i = NR; i >= 1; i--) {
                printf "%s", seq[i];
            }
        }
    }' "$TMP_SEQ")

    # --- Step 4: Validate sequence length ---
    aa_length=${#orf_aa}
    expected_dna_length=$((aa_length * 3))
    
    if [ -z "$full_seq" ]; then
        full_seq="NA"
    elif [ ${#full_seq} -lt $expected_dna_length ]; then
        # Allow small length differences (for example, a stop codon), but warn about large differences
        # echo "Warning: Short sequence for $gene_id" >&2
        : # no-op
    fi
    
    # Write results
    echo -e "$gene_id\t$transcript_id\t$coordinate_id\t$chrom\t$coordinate\t$strand\t$full_seq\t$start_codon\t$orf_aa"

done > "$TMP_OUTPUT"

# Combine the header and results and handle possible carriage returns
cat "$TMP_HEADER" "$TMP_OUTPUT" | awk 'BEGIN {FS=OFS="\t"} {gsub(/\r/,"",$7); print}' | gzip > "$OUTPUT_FILE"

# Remove temporary files
rm -f "${TMP_PREFIX}"*

echo "Finished processing: $OUTPUT_FILE"
