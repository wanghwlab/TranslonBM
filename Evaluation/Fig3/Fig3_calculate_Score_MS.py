import pandas as pd
import os
import itertools
import re
from multiprocessing import Pool
import sys

# --- Configuration (adjust for your environment) ---

WORKDIR = "/home/tangyuewen/ORF_benchmark/Maxquant_2026.1"
BASE_OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_group/PR_recall_Maxquant/"

REF_BLOCK_FILE = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.block_Trans_Mod.txt"

SAMPLES = [ "SRX1254413", "SRX11812007_SRX11812008_SRX11812009", "SRX5887328_SRX5887329_SRX5887330"]
ALIGNERS = ["tophat2", "STAR", "hisat2"]
PRED_FLAGS = ["orf_pred_default"] 
TRIM_FLAGS = ["trim","untrim"]
FDR = [""]

NUM_PROCESSES = 16

ANNOTATED_CDS_SET = set()

# --- Helper functions ---
def normalize_coords_string(coords_str: str) -> str:
    try:
        blocks = coords_str.strip().split(',')
        block_tuples = []
        for b in blocks:
            start, end = map(int, b.split('-'))
            block_tuples.append((start, end))
        block_tuples.sort(key=lambda x: x[0])
        return ",".join([f"{s}-{e}" for s, e in block_tuples])
    except Exception:
        return coords_str

def load_annotated_cds(ref_path: str):
    global ANNOTATED_CDS_SET
    print(f"Loading annotated CDS blocks from: {ref_path} ...")
    try:
        with open(ref_path, 'r') as f:
            for line in f:
                if not line.strip(): continue
                fields = line.strip().split()
                if len(fields) >= 5:
                    chrom, raw_blocks, strand = fields[1].strip(), fields[2].strip(), fields[3].strip()
                    norm_blocks = normalize_coords_string(raw_blocks)
                    ANNOTATED_CDS_SET.add((chrom, strand, norm_blocks))
        print(f"Successfully loaded {len(ANNOTATED_CDS_SET)} annotated CDS records.")
    except Exception as e:
        print(f"[ERROR] Failed to load annotated CDS file: {e}")
        sys.exit(1)

def is_all_seq(protein_id_str):
    if not isinstance(protein_id_str, str): return False
    main_ids = protein_id_str.split('|')[0]
    return all(id_val.startswith('seq_') for id_val in main_ids.split(';'))

def has_any_unique_peptides(count_str):
    if not isinstance(count_str, str): return False
    try:
        return any(int(c) > 0 for c in count_str.split(';'))
    except (ValueError, TypeError):
        return False

def extract_seq_ids(id_string):
    if not isinstance(id_string, str): return []
    return re.findall(r'seq_\d+', id_string)

# --- Core worker function ---
def process_level(params):
    sample, trim_flag, pred_flag, aligner, fdr, level = params
    print(f"Starting [{level.capitalize()} Level]: Sample={sample}, Aligner={aligner}")

    try:
        output_subdir = os.path.join(BASE_OUTPUT_DIR, trim_flag, pred_flag)
        stats_path = os.path.join(output_subdir, f"{sample}_MS_valid_{level}_{aligner}_score.tsv")
        
        ref_path = os.path.join(WORKDIR, "maxquant_ref", f"{pred_flag}_{trim_flag}", f"{sample}_ref.tsv")
        rerun_subdir = f"{pred_flag}_{trim_flag}"
        
        if level == 'peptide':
            ms_data_path = os.path.join(WORKDIR, "maxquant_rerun", rerun_subdir, fdr, aligner, sample, "combined", "txt", "peptides.txt")
            ms_df = pd.read_csv(ms_data_path, sep='\t', low_memory=False)
            valid_items = ms_df[(ms_df['Unique (Proteins)'] == 'yes') & (ms_df['Potential contaminant'].isna()) & (ms_df['Reverse'].isna())].copy()
            id_column = 'Proteins'
        else: # protein level
            ms_data_path = os.path.join(WORKDIR, "maxquant_rerun", rerun_subdir, fdr, aligner, sample, "combined", "txt", "proteinGroups.txt")
            ms_df = pd.read_csv(ms_data_path, sep='\t', low_memory=False)
            has_mask = ms_df['Peptide counts (unique)'].apply(has_any_unique_peptides)
            valid_items = ms_df[(ms_df['Potential contaminant'].isna()) & (ms_df['Reverse'].isna()) & (has_mask)].copy()
            valid_items = valid_items[valid_items['Protein IDs'].apply(is_all_seq)]
            id_column = 'Protein IDs'

        # 1. Load the complete reference file and trim whitespace from column headers
        full_ref_df = pd.read_csv(ref_path, sep='\t')
        full_ref_df.rename(columns=lambda x: x.strip(), inplace=True)
        
        ref_df = full_ref_df[(full_ref_df['samples'] == sample) & (full_ref_df['aligner'] == aligner)].copy()
        
        if ref_df.empty:
            print(f"  -> [Warning] No matching data were found in the reference file for Sample={sample}, Aligner={aligner} ; skipping.")
            return

        # ★★★ Core fix 1: extract the complete set of valid IDs for the current aligner (the universe) ★★★
        aligner_seq_ids = set(ref_df['sequence_id'].astype(str).str.strip())

        # 2. Define the Novel ORF set
        novel_seq_ids = set()
        if 'chrom' in ref_df.columns and 'coordinate_0base' in ref_df.columns:
            for _, row in ref_df.iterrows():
                try:
                    chrom = str(row['chrom']).strip()
                    strand = str(row['strand']).strip()
                    raw_coords = str(row['coordinate_0base']).strip()
                    norm_blocks = normalize_coords_string(raw_coords)
                    
                    if (chrom, strand, norm_blocks) not in ANNOTATED_CDS_SET:
                        novel_seq_ids.add(str(row['sequence_id']).strip())
                except Exception:
                    continue
        else:
            print(f"  -> [Error] The reference file lacks a coordinate column. Current columns: {list(ref_df.columns)}")
            return
            
        print(f"  -> [Statistics] {sample} ({aligner}) Novel ORFs in _ref.tsv: {len(novel_seq_ids)}")

        results = []
        
        # 3. Extract all validated seq_id values from the mass-spectrometry results
        raw_validated_seq_ids = set(valid_items[id_column].apply(extract_seq_ids).explode().dropna())
        
        # ★★★ Core fix 2: intersect the global MS-detected set with the current aligner set ★★★
        validated_seq_ids = raw_validated_seq_ids.intersection(aligner_seq_ids)
        ms_valid_orfs_num = len(validated_seq_ids)
        
        # Calculate MS_novel (MS set intersected with the Novel set)
        ms_novel_seq_ids = validated_seq_ids.intersection(novel_seq_ids)
        ms_novel_orfs_num = len(ms_novel_seq_ids)

        # 4. Calculate by orf_detector group
        for detector, group in ref_df.groupby('orf_detector'):
            predicted_seq_ids = set(group['sequence_id'].astype(str).str.strip().unique())
            pred_orfs_num = len(predicted_seq_ids)
            
            # Overall metrics
            valid_orfs_num = len(predicted_seq_ids.intersection(validated_seq_ids))
            validation_rate = valid_orfs_num / pred_orfs_num if pred_orfs_num > 0 else 0
            recall = valid_orfs_num / ms_valid_orfs_num if ms_valid_orfs_num > 0 else 0
            f1 = (2 * validation_rate * recall) / (validation_rate + recall) if (validation_rate + recall) > 0 else 0
            
            # Novel metrics
            # ★★★ Define pred_novel_seq_ids first (predicted sequences intersected with the full Novel set) ★★★
            pred_novel_seq_ids = predicted_seq_ids.intersection(novel_seq_ids)
            # 1. Total predicted Novel ORFs (denominator of Validation Rate)
            pred_novel_orfs_num = len(pred_novel_seq_ids)

            # 2. Number of MS-validated Novel ORFs (numerator)
            valid_novel_orfs_num = len(pred_novel_seq_ids.intersection(ms_novel_seq_ids))

            # 3. Calculate Novel Recall
            novel_recall = valid_novel_orfs_num / ms_novel_orfs_num if ms_novel_orfs_num > 0 else 0

            # 4. Calculate Novel Validation Rate (added)
            novel_validation_rate = valid_novel_orfs_num / pred_novel_orfs_num if pred_novel_orfs_num > 0 else 0

            results.append({
                'aligner': aligner,
                'orf_detector': detector,
                'MS_valid_ORFs_num': ms_valid_orfs_num,
                'pred_ORFs_num': pred_orfs_num,
                'valid_ORFs_num': valid_orfs_num,
                'Validation_Rate': validation_rate,
                'Recall': recall,
                'F1_score': f1,
                'MS_novel_ORFs_num': ms_novel_orfs_num,
                'pred_novel_ORFs_num': pred_novel_orfs_num,         # Output the number of predicted Novel ORFs
                'valid_novel_ORFs_num': valid_novel_orfs_num,
                'Novel_Validation_Rate': novel_validation_rate,     # Output Novel Validation Rate
                'Novel_Recall': novel_recall
            })

        output_df = pd.DataFrame(results)
        os.makedirs(os.path.dirname(stats_path), exist_ok=True)
        output_df.to_csv(stats_path, sep='\t', index=False)
        print(f"Completed [{level.capitalize()} Level]: Sample={sample}, Aligner={aligner}")

    except FileNotFoundError as e:
        print(f"[Error] File not found for {sample}/{aligner}: {e}. Skip.")
    except Exception as e:
        print(f"[Error] Processing failed for {sample}/{aligner}: {e}. Skip.")

def process_peptide_wrapper(params):
    return process_level(params + ('peptide',))

def process_protein_wrapper(params):
    return process_level(params + ('protein',))

# --- Main entry point ---
if __name__ == '__main__':
    load_annotated_cds(REF_BLOCK_FILE)

    task_params = list(itertools.product(SAMPLES, TRIM_FLAGS, PRED_FLAGS, ALIGNERS, FDR))
    print(f"Generated {len(task_params) * 2} analysis tasks (peptide/protein), using {NUM_PROCESSES} cores in parallel.")

    print("\n--- Starting peptide-level analysis ---")
    with Pool(processes=NUM_PROCESSES) as pool:
        pool.map(process_peptide_wrapper, task_params)
    
    print("\n--- Starting protein-level analysis ---")
    with Pool(processes=NUM_PROCESSES) as pool:
        pool.map(process_protein_wrapper, task_params)

    print("\nAll analysis tasks completed.")
