import os
import pandas as pd
import re 

def extract_parameters(path_string, keywords):
    """
    Helper function that finds and joins matching keywords in a path string.
    """
    found_keywords = [k for k in keywords if k in path_string]
    if found_keywords:
        return '_'.join(found_keywords)
    return '-'

def process_benchmark_files(root_directory, output_csv_path):
    """
    Traverse all .txt files under the specified path, extract information, and merge it into a CSV file.
    """
    all_data_rows = []

    print(f"Searching for .txt files under '{root_directory}' ...")

    for root, dirs, files in os.walk(root_directory):
        for filename in files:
            if filename.endswith(".txt"):
                
                if 'time_benchmarks' not in root.split(os.sep):
                    continue 

                full_path = os.path.join(root, filename)
                
                base_name = os.path.splitext(filename)[0]
                try:
                    parts = base_name.rsplit('_', 1)
                    sample = parts[0]
                    aligner = parts[1]
                except IndexError:
                    print(f"Warning: File name has an unexpected format and was skipped: {filename}")
                    continue

                try:
                    data_df = pd.read_csv(full_path, sep='\t')
                    if data_df.empty:
                        print(f"Warning: File is empty and was skipped: {full_path}")
                        continue
                    
                    data_df['sample'] = sample
                    data_df['aligner'] = aligner
                    data_df['doc_path'] = full_path

                    all_data_rows.append(data_df)
                except Exception as e:
                    print(f"Error: Processing file {full_path} failed: {e}")

    if not all_data_rows:
        print("No .txt files were found under a time_benchmarks subdirectory.")
        return

    print("Merging all data...")
    final_df = pd.concat(all_data_rows, ignore_index=True)


    print("Extracting software, chr, and parameter information from paths...")

    regex_pattern = r'(\w+)_(chr[NM])'
    extracted_info = final_df['doc_path'].str.extract(regex_pattern)
    
    final_df['software'] = extracted_info[0].fillna('-')
    final_df['chr'] = extracted_info[1].fillna('-')

    param_keywords = [
        'ribotaper', 'orfrater', 'rpbp', 'orfquant', 'ribowave', 'gedi', 
        'ribotricer', 'ribohmm', 'ribotish', 'ribocode', 'riborf'
    ]
    final_df['parameter'] = final_df['doc_path'].apply(extract_parameters, keywords=param_keywords)

    desired_order = [
        'sample', 'aligner', 'software', 'chr', 'parameter',
        's', 'h:m:s', 'max_rss', 'max_vms', 'max_uss', 
        'max_pss', 'io_in', 'io_out', 'mean_load', 'cpu_time', 'doc_path'
    ]
    
    for col in desired_order:
        if col not in final_df.columns:
            final_df[col] = None
            
    final_df = final_df[desired_order]

    final_df.to_csv(output_csv_path, index=False)
    print(f"\nProcessing completed. Results saved to: {output_csv_path}")
    print(f"Processed {len(final_df)} file records.")


if __name__ == '__main__':

    SEARCH_PATH = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_simu/"
    OUTPUT_FILE = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/time_benchmark/ORF_simu_time_benchmark_summary.csv"

    process_benchmark_files(SEARCH_PATH, OUTPUT_FILE)
