#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import csv
import multiprocessing
import sys
from typing import List, Dict, Set, Tuple

# --- 配置區 ---
BASE_INPUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/"
BASE_OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/PR_recall_Blocks/merged_ATG/" 
REF_BLOCK_FILE = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.block_Trans_Mod.txt"

NUM_PROCESSES = 32

DIRECTORIES_TO_PROCESS = [
    "orf_pred_default_trim",
    "orf_pred_default_untrim",
    # "orf_pred_undisputed_trim",
    # "orf_pred_undisputed_untrim",
]


def normalize_coords_string(coords_str: str) -> str:
    """
    將坐標字符串標準化，確保格式一致。
    例如輸入 "100-200,50-80" 會被重排為 "50-80,100-200"。
    這能防止因exon順序不同導致的匹配失敗。
    """
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

def load_ground_truth(ref_path: str) -> Set[Tuple]:
    """
    從新的 TXT 文件中加載真實的 CDS Block 坐標。
    文件格式: GeneID <tab> Chrom <tab> BlockCoords <tab> Strand <tab> TranscriptID
    """
    ground_truth_set = set()
    print(f"Loading ground truth Blocks from: {ref_path} ...")
    
    try:
        with open(ref_path, 'r') as f:
            for line in f:
                if not line.strip(): continue
                fields = line.strip().split()
                
                if len(fields) >= 5:
                    chrom = fields[1]
                    raw_blocks = fields[2]
                    strand = fields[3]
                    
                    norm_blocks = normalize_coords_string(raw_blocks)
                    
                    ground_truth_set.add((chrom, strand, norm_blocks))
                    
    except FileNotFoundError:
        print(f"[ERROR] Reference file not found at: {ref_path}. Exiting.", file=sys.stderr)
        sys.exit(1)
        
    print(f"Loaded {len(ground_truth_set)} unique ORF Block records.")
    return ground_truth_set

def parse_filename(filename_str: str) -> tuple:
    """從文件名中解析出 sample, soft, tools。"""
    base = filename_str.replace('_gcoor.tsv.gz', '')
    parts = base.split('_')
    if len(parts) >= 3:
        tools, soft, sample = parts[-1], parts[-2], '_'.join(parts[:-2])
        return sample, soft, tools
    return base, 'unknown', 'unknown'

def extract_orf_info_from_row(row: pd.Series) -> tuple:
    """
    從 DataFrame 的一行中提取 ORF 完整坐標信息。
    返回: (chrom, strand, normalized_block_string)
    """
    try:
        chrom = row['chrom']
        strand = row['strand']
        raw_coords = str(row['coordinate_0base'])
        
        norm_coords = normalize_coords_string(raw_coords)
        
        return (chrom, strand, norm_coords)
    except (ValueError, IndexError, KeyError, AttributeError):
        return None

def process_file_worker(task_info: tuple) -> Dict:
    """處理單個預測文件的工作函數。"""
    filepath, ground_truth_set, filter_type = task_info
    
    filename = os.path.basename(filepath)
    sample, soft, tools = parse_filename(filename)
    total_gt_items = len(ground_truth_set) 

    try:
        df = pd.read_csv(filepath, sep='\t', compression='gzip', low_memory=False)
    except Exception:
        df = pd.DataFrame()

    if df.empty or 'start_codon' not in df.columns or 'coordinate_0base' not in df.columns:
        return {'sample': sample, 'soft': soft, 'tools': tools, 'filter_type': filter_type, 
                'precision': 0.0, 'recall': 0.0, 'F1-score': 0.0,
                'TP': 0, 'FP': 0, 'FN': total_gt_items, 'Hit_Count': 0, 'totalORF': 0, 'GT_Total': total_gt_items}
    

    all_predictions = [res for res in df.apply(extract_orf_info_from_row, axis=1) if res is not None]
    
    unique_predicted_set = set(all_predictions)
    hits_in_gt = unique_predicted_set.intersection(ground_truth_set)
    hit_count = len(hits_in_gt)
    
    fn = total_gt_items - hit_count

    if filter_type == 'ATG':
        df_to_process = df[df['start_codon'] == 'ATG']
    else: 
        df_to_process = df[df['start_codon'] != 'ATG']

    subset_predictions = [res for res in df_to_process.apply(extract_orf_info_from_row, axis=1) if res is not None]
    totalORF = len(subset_predictions)
    
    if totalORF == 0:
        return {'sample': sample, 'soft': soft, 'tools': tools, 'filter_type': filter_type, 
                'precision': 0.0, 'recall': 0.0, 'F1-score': 0.0,
                'TP': 0, 'FP': 0, 'FN': fn, 'Hit_Count': hit_count, 'totalORF': 0, 'GT_Total': total_gt_items}

    tp = sum(1 for pred in subset_predictions if pred in ground_truth_set)
    fp = totalORF - tp
    
    precision = tp / totalORF

    recall_denominator = tp + fn 
    
    recall = tp / total_gt_items if total_gt_items > 0 else 0.0
    f1_score = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0.0
    
    return {'sample': sample, 'soft': soft, 'tools': tools, 'filter_type': filter_type, 
            'precision': precision, 'recall': recall, 'F1-score': f1_score,
            'TP': tp, 'FP': fp, 'FN': fn, 'Hit_Count': hit_count, 'totalORF': totalORF, 'GT_Total': total_gt_items}


def run_benchmark_for_directory(prediction_dir: str, output_basename: str, ground_truth_set: Set[Tuple]):
    """對單個目錄執行完整的評分。"""
    print("-" * 80)
    print(f"Starting Block-level benchmark for directory: {prediction_dir}")
    if not os.path.isdir(prediction_dir):
        print(f"[WARNING] Directory does not exist. Skipping: {prediction_dir}")
        return

    search_path = os.path.join(prediction_dir, "*_gcoor.tsv.gz")
    prediction_files = sorted(glob.glob(search_path))
    if not prediction_files:
        print(f"[WARNING] No '.tsv.gz' files found. Skipping: {prediction_dir}")
        return
        
    print(f"Found {len(prediction_files)} prediction files to process.")

    tasks = []
    for filepath in prediction_files:
        tasks.append((filepath, ground_truth_set, 'ATG'))
        tasks.append((filepath, ground_truth_set, 'NTG'))

    pool_size = min(NUM_PROCESSES, len(tasks))
    results_data = []
    
    if pool_size > 0:
        with multiprocessing.Pool(processes=pool_size) as pool:
            results_data = pool.map(process_file_worker, tasks)
    
    print(f"Processing complete. Writing results...")

    results_atg = [res for res in results_data if res and res['filter_type'] == 'ATG']
    results_ntg = [res for res in results_data if res and res['filter_type'] == 'NTG']

    output_csv_atg = os.path.join(BASE_OUTPUT_DIR, f"{output_basename}_ATG_BlockMatch.csv")
    output_csv_ntg = os.path.join(BASE_OUTPUT_DIR, f"{output_basename}_NTG_BlockMatch.csv")

    header = ['sample', 'soft', 'tools', 'precision', 'recall', 'F1-score', 'TP', 'FP', 'FN', 'Hit_Count', 'totalORF', 'GT_Total']
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    if results_atg:
        with open(output_csv_atg, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            for r in results_atg: del r['filter_type']
            writer.writerows(sorted(results_atg, key=lambda x: (x['sample'], x['soft'], x['tools'])))

    if results_ntg:
        with open(output_csv_ntg, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            for r in results_ntg: del r['filter_type']
            writer.writerows(sorted(results_ntg, key=lambda x: (x['sample'], x['soft'], x['tools'])))
            
    print(f"Saved results to {BASE_OUTPUT_DIR}")


def main():
    print("--- Starting Batch ORF Prediction Benchmark (Exact Block Match) ---")
    
    gt_set = load_ground_truth(REF_BLOCK_FILE)
    if not gt_set: return

    for dir_name in DIRECTORIES_TO_PROCESS:
        full_prediction_dir = os.path.join(BASE_INPUT_DIR, dir_name)
        output_filename_base = f"{dir_name}_score"
        run_benchmark_for_directory(full_prediction_dir, output_filename_base, gt_set)

    print("-" * 80)
    print("--- Script Finished ---")

if __name__ == "__main__":
    main()
