import pandas as pd
import glob
import os

def extract_rpbp_offsets():
    # 1. Set the search path
    input_pattern = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/rpbp/rpbp_chrN/orf_pred_default/*/metagene-profiles/*.periodic-offsets.csv.gz"
    
    print(f"Searching for files: {input_pattern} ...")
    files = glob.glob(input_pattern)
    
    if not files:
        print("Error: No .csv.gz files were found; check the path.")
        return

    print(f"Found {len(files)} files; starting processing...")
    
    all_data = []

    for file_path in files:
        try:
            # --- A. Parse the file name to extract the sample ---
            basename = os.path.basename(file_path)
            
            sample_name = basename.replace(".periodic-offsets.csv.gz", "")
            
            # --- B. Read the gzipped CSV ---
            df = pd.read_csv(file_path, compression='gzip')
            
            if 'length' not in df.columns or 'highest_peak_offset' not in df.columns:
                print(f"Skipping file {basename}: column names do not match")
                continue

            # --- C. Clean and transform data ---
            temp_df = df[['length', 'highest_peak_offset']].copy()
            
            temp_df['read_length'] = temp_df['length'].astype(int)
            
            temp_df['offset'] = temp_df['highest_peak_offset'].abs().astype(int)
            temp_df['sample'] = sample_name
            temp_df['ORFtools'] = 'rpbp'
            
            final_cols_df = temp_df[['sample', 'ORFtools', 'read_length', 'offset']]
            
            all_data.append(final_cols_df)
            
        except Exception as e:
            print(f"Processing file {basename} failed: {e}")

    # --- D. Merge and save ---
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        final_df = final_df[['sample', 'ORFtools', 'read_length', 'offset']]
        final_df = final_df.sort_values(by=['sample', 'read_length'])
        
        output_file = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/rpbp_offsets_extracted.csv"
        final_df.to_csv(output_file, index=False)
        
        print("-" * 30)
        print("Processing completed!")
        print(f"Number of samples: {final_df['sample'].nunique()}")
        print(f"File saved as: {os.path.abspath(output_file)}")
        print("\nData preview:")
        print(final_df.head())
    else:
        print("No valid data were extracted.")

if __name__ == "__main__":
    extract_rpbp_offsets()
