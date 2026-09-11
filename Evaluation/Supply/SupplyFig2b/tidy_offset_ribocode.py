import glob
import os
import csv

# ================= Configuration =================

# 1. Input-file glob pattern
INPUT_PATTERN = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/ribocode/ribocode_chrN/P_site_determination/SR*/SRX*_pre_config.txt"

# 2. Output file name
OUTPUT_FILE = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/ribocode_extracted_offsets.csv"

# ================= Main program =================

def extract_ribocode_offsets():
    # Find all matching files
    files = glob.glob(INPUT_PATTERN)
    
    if not files:
        print("No files were found; check the path configuration.")
        return

    print(f"Found {len(files)} files; starting extraction...")

    # Open the CSV file for writing
    with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile) # comma-delimited by default
        
        # Write the header
        writer.writerow(["sample", "ORFtools", "read_length", "offset"])

        for file_path in files:
            try:
                # Read all lines
                with open(file_path, 'r') as infile:
                    lines = infile.readlines()
                
                # Find a non-comment line (nonempty and not beginning with #)
                config_line = None
                for line in lines:
                    stripped_line = line.strip()
                    if stripped_line and not stripped_line.startswith('#'):
                        config_line = stripped_line
                        break
                
                # Parse the configuration line if one was found
                if config_line:
                    # Split on whitespace (tab or space)
                    parts = config_line.split()
                    
                    if len(parts) >= 5:
                        # parts[0] -> SampleName (for example SRX876063_SRX876069_tophat2)
                        # parts[3] -> ReadLengths (for example 25,28,29)
                        # parts[4] -> Offsets (for example 9,12,12)
                        
                        sample_name = parts[0]
                        lengths_str = parts[3]
                        offsets_str = parts[4]
                        
                        # Convert comma-delimited strings to lists
                        lengths = lengths_str.split(',')
                        offsets = offsets_str.split(',')
                        
                        # Ensure that lengths and offsets have equal counts, then write paired values
                        if len(lengths) == len(offsets):
                            for l, o in zip(lengths, offsets):
                                # Write one row: sample, ribocode, length, offset
                                writer.writerow([sample_name, "ribocode", l, o])
                        else:
                            print(f"Warning: {sample_name} has different numbers of lengths and offsets")
                    else:
                        print(f"Warning: Unexpected file format (insufficient columns): {file_path}")
                else:
                    print(f"Warning: No valid configuration line was found in the file: {file_path}")

            except Exception as e:
                print(f"Error processing file {file_path}: {e}")

    print(f"Processing completed. Results saved to:  {os.path.abspath(OUTPUT_FILE)}")

if __name__ == "__main__":
    extract_ribocode_offsets()
