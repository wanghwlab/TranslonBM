import glob
import os
import csv
import re

# 1. Define the search path
search_path = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_simu_untrim/rpbp/rpbp_chrN/orf_pred_default/*/orf-profiles/*.profiles.mtx.gz"

# 2. Define the output file name
output_csv = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_simu/rpbp_converted_offsets.csv"

def parse_rpbp_filename(file_path):
    """
    Parse the file name to extract sample name, lengths, and offsets
    Example file name: SRX...tophat2.length-25-28-29-33.offset-9-12-12-0.profiles.mtx.gz
    """
    filename = os.path.basename(file_path)
    
    # Use a regular expression to extract key fields
    pattern = re.compile(r"^(.*?)\.length-([\d-]+)\.offset-([\d-]+)\.profiles\.mtx\.gz$")
    match = pattern.match(filename)
    
    if match:
        sample_name = match.group(1)
        lengths_str = match.group(2)
        offsets_str = match.group(3)
        
        lengths = [x for x in lengths_str.split('-') if x]
        offsets = [x for x in offsets_str.split('-') if x]
        
        return sample_name, lengths, offsets
    else:
        print(f"Warning: Unable to parse the file-name format -> {filename}")
        return None, None, None

def main():
    files = glob.glob(search_path)
    print(f"Found {len(files)} files; starting processing...")
    
    results = []
    
    for f in files:
        sample_name, lengths, offsets = parse_rpbp_filename(f)
        
        if sample_name and lengths and offsets:
            if len(lengths) != len(offsets):
                print(f"Error: {sample_name} has different numbers of lengths and offsets.")
                continue
            
            for length, offset in zip(lengths, offsets):
                row = [sample_name, 'rpbp', length, offset]
                results.append(row)

    # Write the CSV file
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=',') 
        writer.writerows(results)
        
    print(f"Processing completed. Results saved to:  {os.path.abspath(output_csv)}")

if __name__ == "__main__":
    main()
