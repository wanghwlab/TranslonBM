import pandas as pd
import glob
import os
import re

def parse_ribotricer_offsets():
    # 1. Define the Ribotricer file search path
    input_pattern = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribotricer/ribotricer_chrN/orf_pred_default/*/*_psite_offsets.txt"
    
    # 2. Define the base offset
    DEFAULT_BASE_OFFSET = 12 

    files = glob.glob(input_pattern)
    print(f"Found {len(files)} Ribotricer offset files; starting conversion...")

    all_data = []

    for file_path in files:
        try:
            # --- Extract metadata (sample and ORFtools) ---
            filename = os.path.basename(file_path)
            
            # Remove the suffix to obtain the sample name
            sample_name = filename.replace("_psite_offsets.txt", "")
            orf_tool = "ribotricer"

            # --- Parse file contents ---
            base_length = None
            offsets_dict = {}

            with open(file_path, 'r') as f:
                lines = f.readlines()

            for line in lines:
                line = line.strip()
                
                # 1. Extract the base length (relative lag to base: 29)
                base_match = re.search(r"relative lag to base:\s*(\d+)", line)
                if base_match:
                    base_length = int(base_match.group(1))
                    offsets_dict[base_length] = DEFAULT_BASE_OFFSET
                    continue

                # 2. Extract the lag for each length (lag of 28: 0)
                lag_match = re.search(r"lag of\s*(\d+):\s*(-?\d+)", line)
                if lag_match:
                    read_len = int(lag_match.group(1))
                    lag = int(lag_match.group(2))
                    
                    # Core calculation formula
                    abs_offset = DEFAULT_BASE_OFFSET + lag
                    offsets_dict[read_len] = abs_offset

            # --- Convert to a DataFrame ---
            if offsets_dict:
                for r_len, off_val in offsets_dict.items():
                    all_data.append({
                        'sample': sample_name,
                        'ORFtools': orf_tool,
                        'read_length': r_len,
                        'offset': off_val
                    })
            else:
                print(f"Warning: File {filename} contained no valid parsed data")

        except Exception as e:
            print(f"Processing file {file_path} failed: {e}")

    # --- Save results ---
    if all_data:
        df = pd.DataFrame(all_data)
        
        # Sort by sample and length
        df = df.sort_values(by=['sample', 'read_length'])
        
        output_file = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/ribotricer_converted_offsets.csv"
        df.to_csv(output_file, index=False)
        
        print("-" * 30)
        print(f"Conversion completed!")
        print(f"Saved file: {os.path.abspath(output_file)}")
        print("\nData preview (first five rows):")
        print(df.head())
        
        print("\n[Validation] Check the 26-nt offset for simulation_6M_T3 (expected: 9):")
        check = df[(df['sample'].str.contains("simulation_6M_T3")) & (df['read_length'] == 26)]
        print(check)
    else:
        print("No data were extracted.")

if __name__ == "__main__":
    parse_ribotricer_offsets()
