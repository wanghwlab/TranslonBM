import pandas as pd
import glob
import os

def merge_offset_files():
    # 1. Define the search path
    search_pattern = "/home/tangyuewen/ORF_benchmark/rerun_2025.12/ORFdetect_simu/*/*chrN/P_site_determination/*/simulation_*_standard_offsets.txt"
    
    print(f"Searching for files: {search_pattern} ...")
    files = glob.glob(search_pattern)
    
    if not files:
        print("Error: No matching files were found; check the path.")
        return

    print(f"Found {len(files)} files; starting processing...")

    all_data = []

    for file_path in files:
        try:
            # 2. Extract metadata (sample and ORFtools)
            file_name = os.path.basename(file_path)
            sample_name = file_name.replace("_standard_offsets.txt", "")
            
            path_parts = file_path.split(os.sep)
            try:
                base_index = path_parts.index("ORFdetect_simu")
                orf_tool = path_parts[base_index + 1] 
            except ValueError:
                print(f"Warning: ORFdetect_simu was not found in the path; unable to infer the tool name: {file_path}")
                orf_tool = "Unknown"

            # 3. Read file contents
            df = pd.read_csv(file_path, sep=r'\s+', header=None, names=['read_length', 'offset'], engine='python')
            
            df['sample'] = sample_name
            df['ORFtools'] = orf_tool
            
            all_data.append(df)
            
        except Exception as e:
            print(f"Processing file {file_path} failed: {e}")

    # 4. Merge and save
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        
        cols = ['sample', 'ORFtools', 'read_length', 'offset']
        final_df = final_df[cols]
        
        output_file = "/home/tangyuewen/ORF_benchmark/rerun_2025.12/plots/offset_plots/merged_offsets_all.csv"
        final_df.to_csv(output_file, index=False)
        
        print("-" * 30)
        print(f"Processing completed!")
        print(f"Total rows: {len(final_df)}")
        print(f"Number of samples: {final_df['sample'].nunique()}")
        print(f"Number of ORFtools: {final_df['ORFtools'].nunique()}")
        print(f"File saved as: {os.path.abspath(output_file)}")
        
        print("\nData preview:")
        print(final_df.head())
    else:
        print("No data were extracted.")

if __name__ == "__main__":
    merge_offset_files()
