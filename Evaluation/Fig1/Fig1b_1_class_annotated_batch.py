import os
import pandas as pd
import glob

def split_file_by_annotation_status(file1, file2, annotated_dir, non_annotated_dir):
    # Extract base filename without path and extension
    base_filename = os.path.basename(file1)
    
    # Create output file paths (keep original filename)
    annotated_output = os.path.join(annotated_dir, base_filename)
    non_annotated_output = os.path.join(non_annotated_dir, base_filename)

    # Load file1 (gzipped TSV, all columns)
    df1 = pd.read_csv(file1, sep="\t", compression="gzip", header=None)
    df1.columns = [f"col{i+1}" for i in range(df1.shape[1])]

    # Load file2 (TXT, only first 4 columns)
    df2 = pd.read_csv(file2, sep="\t", header=None, usecols=[0, 1, 2, 3])
    df2.columns = ["col1", "col2", "col3", "col4"]

    # Mark rows as annotated or non_annotated
    is_annotated = df1[["col1", "col3", "col4", "col5"]].apply(tuple, axis=1).isin(
        df2.apply(tuple, axis=1)
    )
    # Split and save
    annotated = df1[is_annotated]
    non_annotated = df1[~is_annotated]
    
    # Save with gzip compression
    annotated.to_csv(annotated_output, sep="\t", index=False, header=False, compression='gzip')
    non_annotated.to_csv(non_annotated_output, sep="\t", index=False, header=False, compression='gzip')

    return annotated.shape[0], non_annotated.shape[0]

# File paths
file2 = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.block_Trans_Mod.txt"
#input_dir = "/home/tangyuewen/ORF_benchmark/final_ORFs/trim/orf_pred_default/"
input_dir = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_unique/orf_pred_default_trim/"
annotated_dir = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_canonical/orf_pred_default_trim/"
non_annotated_dir = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_non_canonical/orf_pred_default_trim/"

# Ensure output directories exist
os.makedirs(annotated_dir, exist_ok=True)
os.makedirs(non_annotated_dir, exist_ok=True)

# Get all _gcoor.tsv.gz files (excluding _merged_gcoor.tsv.gz)
file_pattern = os.path.join(input_dir, "*_gcoor.tsv.gz")
file1_list = [f for f in glob.glob(file_pattern) if "_merged_gcoor.tsv.gz" not in f]

print(f"Found {len(file1_list)} files to process")

# Process each file
for file1 in file1_list:
    try:
        num_annotated, num_non_annotated = split_file_by_annotation_status(
            file1, file2, annotated_dir, non_annotated_dir)
        print(f"Processed {os.path.basename(file1)}")
        print(f"  Annotated: {num_annotated}, Non-annotated: {num_non_annotated}")
    except Exception as e:
        print(f"Error processing {file1}: {str(e)}")

print("Processing completed.")
