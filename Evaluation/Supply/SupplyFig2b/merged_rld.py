import pandas as pd
import glob
import os

def extract_rld_from_stats():
    # 1. Set file paths
    input_pattern = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/Mapping/merge_chrN/*.stats"
    
    print(f"Searching for files: {input_pattern} ...")
    files = glob.glob(input_pattern)
    
    if not files:
        print("Error: No .stats files were found; check the path.")
        return

    print(f"Found {len(files)} files; extracting RL lines...")
    
    all_data = []

    for file_path in files:
        try:
            # --- A. Parse file name ---
            basename = os.path.basename(file_path)
            filename_no_ext = basename.replace('.stats', '')
            
            if '_' in filename_no_ext:
                sample_name, tool_name = filename_no_ext.rsplit('_', 1)
            else:
                sample_name = filename_no_ext
                tool_name = "Unknown"

            # --- B. Extract content (core modification) ---
            file_data = []
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('RL'):
                        parts = line.strip().split()
                        
                        if len(parts) >= 3:
                            try:
                                r_len = int(parts[1])
                                r_count = int(parts[2])
                                
                                file_data.append({
                                    'sample': sample_name,
                                    'ORFtools': tool_name, 
                                    'read_length': r_len,
                                    'read_count': r_count
                                })
                            except ValueError:
                                continue 

            if file_data:
                all_data.extend(file_data)
            else:
                print(f"Warning: File {basename} contains no line beginning with RL.")
            
        except Exception as e:
            print(f"Processing file {basename} failed: {e}")

    # --- C. Save results ---
    if all_data:
        final_df = pd.DataFrame(all_data)
        
        cols = ['sample', 'ORFtools', 'read_length', 'read_count']
        final_df = final_df[cols]
        
        final_df = final_df.sort_values(by=['sample', 'ORFtools', 'read_length'])
        
        output_file = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/untrim/merged_RLD_stats_final.csv"
        final_df.to_csv(output_file, index=False)
        
        print("-" * 30)
        print("Extraction completed!")
        print(f"Total rows: {len(final_df)}")
        print(f"Saved as: {os.path.abspath(output_file)}")
        print("\nData preview:")
        print(final_df.head())
    else:
        print("No data were extracted.")

if __name__ == "__main__":
    extract_rld_from_stats()
