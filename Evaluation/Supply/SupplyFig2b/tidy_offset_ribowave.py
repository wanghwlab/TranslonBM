import glob
import os
import csv

# ================= Configuration =================

# 1. Input-file pattern
INPUT_PATTERN = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribowave/ribowave_chrN/orf_pred_default/*/P-site/*.psite1nt.txt"

# 2. Output CSV path
OUTPUT_FILE = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/ribowave_converted_offsets.csv"

# ================= Main program =================
def process_ribowave_offsets():
    files = glob.glob(INPUT_PATTERN)
    
    if not files:
        print("No files were found; check the path configuration.")
        return

    print(f"Found {len(files)} files; starting processing...")

    with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile, delimiter=',') 
        writer.writerow(["sample", "ORFtools", "read_length", "offset"])

        for file_path in files:
            try:
                # --- 1. Get sample name ---
                filename = os.path.basename(file_path)
                sample_name = filename.replace(".psite1nt.txt", "")
                
                # --- 2. Read file contents ---
                with open(file_path, 'r') as infile:
                    for line in infile:
                        line = line.strip()
                        if not line: continue
                        
                        parts = line.split()
                        
                        if len(parts) >= 2:
                            read_len_str = parts[0]
                            psite_pos_str = parts[1]
                            
                            # --- 3. Correct the offset ---
                            try:
                                position = int(psite_pos_str)
                                offset = position - 1
                                
                                writer.writerow([sample_name, "ribowave", read_len_str, offset])
                                
                            except ValueError:
                                print(f"Warning: File {filename} contains a nonnumeric line: {line}")

            except Exception as e:
                print(f"Processing file {file_path} failed: {e}")

    print(f"Processing completed. Results saved to:  {os.path.abspath(OUTPUT_FILE)}")

if __name__ == "__main__":
    process_ribowave_offsets()
