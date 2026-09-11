import gzip
import csv

# --- Configuration (Configuration) ---
# Path to file A
#FILE_A = "/home/tangyuewen/ORF_benchmark/final_ORFs/trim/orf_pred_default/SRX876063_SRX876069_hisat2_ribohmm_gcoor.tsv.gz"
FILE_B = "/home/tangyuewen/ORF_benchmark/rerun_2025.5/ORFblock/final_ORFs/orf_pred_default/SRX876063_SRX876069_hisat2_ORFquant_gcoor.tsv.gz"
# Path to file B
FILE_A = "/home/tangyuewen/ORF_benchmark/rerun_2025.5/ORFblock/final_ORFs/orf_pred_default/SRX876063_SRX876069_hisat2_ribohmm_gcoor.tsv.gz"

# Column index to extract (the third column has index 2)
# The column index to extract (column 3 has an index of 2)
COLUMN_INDEX = 2 
#0 is the gene name, 1 is the transcript name, and 2 is gcoor_id

# --- Function definition (Function Definition) ---

def get_unique_ids_from_file(filepath, col_index):
    """
    Read the specified column from a gzipped TSV file and return the set of unique IDs.
    Reads a specified column from a gzipped TSV file and returns a set of unique IDs.
    """
    unique_ids = set()
    try:
        # Use rt mode to read text
        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            # Use csv.reader for robust TSV parsing
            reader = csv.reader(f, delimiter='\t')
            
            # Skip the header row
            # Skip the header line
            header = next(reader, None)
            if header is None:
                print(f"Warning: File is empty or has no header - {filepath}")
                return unique_ids

            # Read each remaining row
            for row in reader:
                # Ensure that the row contains enough columns
                if len(row) > col_index:
                    unique_ids.add(row[col_index])
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        # Return an empty set when the file is not found
        return set()
    except Exception as e:
        print(f"An error occurred while reading {filepath}: {e}")
        return set()
        
    return unique_ids

# --- Main program (Main Program) ---

if __name__ == "__main__":
    print("Step 1: Extracting unique IDs from File A...")
    ids_A = get_unique_ids_from_file(FILE_A, COLUMN_INDEX)
    
    print("Step 2: Extracting unique IDs from File B...")
    ids_B = get_unique_ids_from_file(FILE_B, COLUMN_INDEX)
    
    # --- Calculation (Calculation) ---
    
    # Get the number of unique IDs in each file
    total_A = len(ids_A)
    total_B = len(ids_B)
    
    # Use set.intersection() to calculate the overlap
    common_ids = ids_A.intersection(ids_B)
    common_count = len(common_ids)
    
    # --- Print results (Print Results) ---
    
    print("----------------------------------------")
    print("Comparison Results:")
    print("----------------------------------------")
    print(f"Total unique coordinate_ids in File A: {total_A}")
    print(f"Total unique coordinate_ids in File B: {total_B}")
    print(f"Number of common coordinate_ids: {common_count}")
    print("----------------------------------------")
    
    # Calculate percentages and handle division by zero
    if total_A > 0:
        percent_A = (common_count / total_A) * 100
        print(f"Overlap: {percent_A:.2f}% of File A's IDs are present in File B.")
    
    if total_B > 0:
        percent_B = (common_count / total_B) * 100
        print(f"Overlap: {percent_B:.2f}% of File B's IDs are present in File A.")
    
    print("----------------------------------------")
