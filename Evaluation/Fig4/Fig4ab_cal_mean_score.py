import os
import glob
import pandas as pd

# --- 1. Configuration ---

BASE_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/"


INPUT_SUBDIRS = [
    #"TIS_score_trim",
    #"TIS_score_intersect_trim",
    #"TIS_score_intersect_untrim",
    #"TIS_score_union_trim",
    #"TIS_score_union_untrim",
    #"Maxquant_score_union_trim",
    #"Maxquant_score_union_untrim",
    #"Maxquant_score_intersect_trim",
    #"Maxquant_score_intersect_untrim",
    #"TISeq_score_intersect_trim",
    #"TISeq_score_intersect_untrim",
    #"TISeq_score_union_trim",
    #"TISeq_score_union_untrim",
    "merged_gain"
]

OUTPUT_SUBDIR = "Mean_score"


def process_and_average_files():
    """
    Main function: find, read, aggregate, and save CSV files.
    """
    full_output_path = os.path.join(BASE_DIR, OUTPUT_SUBDIR)
    
    os.makedirs(full_output_path, exist_ok=True)
    
    for subdir in INPUT_SUBDIRS:
        full_input_path = os.path.join(BASE_DIR, subdir)
        print("-" * 50)
        print(f"Processing subdirectory: {subdir}")

        search_pattern = os.path.join(full_input_path, '*combination_gains*_peptide.csv')
        csv_files = glob.glob(search_pattern)

        if not csv_files:
            print(f"  [Warning] No CSV files were found in directory '{full_input_path}'.")
            continue
            
        print(f"  Found {len(csv_files)} CSV files; starting processing...")
        
        for file_path in csv_files:
            basename = os.path.basename(file_path)
            print(f"    - Processing: {basename}")
            
            try:

                df = pd.read_csv(file_path)
                
                grouping_columns = ['aligner', 'tool_a', 'tool_b']
                
                if not all(col in df.columns for col in grouping_columns):
                    print(f"      [Warning] File {basename} is missing required grouping columns and was skipped.")
                    continue

                averaged_df = df.groupby(grouping_columns).mean().reset_index()
                
                averaged_df.insert(1, 'sample', 'Average7sample')
                
                new_filename = f"{subdir}_{basename}"
                output_filepath = os.path.join(full_output_path, new_filename)
                
                averaged_df.to_csv(output_filepath, index=False)
                print(f"      -> Saved averaged results to: {new_filename}")

            except Exception as e:
                print(f"      [Error] Processing file {basename} encountered an error: {e}")

    print("\nAll files processed.")


if __name__ == '__main__':
    process_and_average_files()
