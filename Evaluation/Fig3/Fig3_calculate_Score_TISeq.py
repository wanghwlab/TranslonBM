#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import csv
import multiprocessing
import sys
from typing import Dict, Set, Tuple


BASE_INPUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/"
GROUND_TRUTH_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/TISeq/" 
BASE_OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_group/PR_recall_TISeq/"

REF_BLOCK_FILE = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.block_Trans_Mod.txt"

DIRECTORIES_TO_PROCESS = [
    "orf_pred_default_trim",
    "orf_pred_default_untrim",
    #"orf_pred_undisputed_trim",
    #"orf_pred_undisputed_untrim",
]

GROUND_TRUTH_TOOL = 'ribotish'
RIBO_TO_TI_MAPPING = {
    "SRX876063_SRX876069": "SRX876066_SRX876072",
    "SRX740748": "SRX740745",
    "SRX5256543_SRX5256555": "SRX5256564_SRX5256565",
    "SRX5887328_SRX5887329_SRX5887330": "SRX5887326",
}

NUM_PROCESSES = 32

ANNOTATED_CDS_SET = set()

def normalize_coords_string(coords_str: str) -> str:
    """將坐標字符串標準化（按起始位置排序），確保能夠精準匹配 Block。"""
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
    """加載註釋的 CDS Block 用於區分 Annotated / Non-canonical ORF"""
    global ANNOTATED_CDS_SET
    print(f"Loading annotated CDS blocks from: {ref_path} ...")
    try:
        with open(ref_path, 'r') as f:
            for line in f:
                if not line.strip(): continue
                fields = line.strip().split()
                if len(fields) >= 5:
                    chrom, raw_blocks, strand = fields[1], fields[2], fields[3]
                    norm_blocks = normalize_coords_string(raw_blocks)
                    ANNOTATED_CDS_SET.add((chrom, strand, norm_blocks))
        print(f"Successfully loaded {len(ANNOTATED_CDS_SET)} annotated CDS records.")
    except Exception as e:
        print(f"[ERROR] Failed to load annotated CDS file: {e}", file=sys.stderr)
        sys.exit(1)

def parse_filename_v2(filename_str: str) -> tuple:
    """从复杂文件名中解析出 (样本名, 比对软件, 预测工具)。"""
    base = filename_str.replace('_gcoor.tsv.gz', '')
    parts = base.split('_')
    if len(parts) >= 3:
        tools, aligner, sample = parts[-1], parts[-2], '_'.join(parts[:-2])
        return sample, aligner, tools
    return None, None, None

def preload_ground_truth_library_final(gt_directory_path: str) -> Dict:
    """
    预加载金标准库。
    改为返回 dict: {(chrom, predicted_start, strand): orf_type}
    """
    ground_truth_library = {}
    print(f"--- Pre-loading Ground Truth Library (Source Tool: {GROUND_TRUTH_TOOL}) ---")
    
    for dir_name in DIRECTORIES_TO_PROCESS:
        search_path = os.path.join(gt_directory_path, dir_name, f"*_{GROUND_TRUTH_TOOL}_gcoor.tsv.gz")
        gt_files = glob.glob(search_path)
        
        for filepath in gt_files:
            filename = os.path.basename(filepath)
            sample, aligner, _ = parse_filename_v2(filename)
            if not all([sample, aligner]): continue
            
            gt_key = (sample, aligner)
            try:
                df = pd.read_csv(filepath, sep='\t', compression='gzip', low_memory=False)
                df = df[df['start_codon'] == 'ATG']
                
                gt_collection = {}
                for _, row in df.iterrows():
                    try:
                        chrom, strand, coords_str = row['chrom'], row['strand'], str(row['coordinate_0base'])
                        
                        norm_blocks = normalize_coords_string(coords_str)
                        orf_type = 'annotated_CDS' if (chrom, strand, norm_blocks) in ANNOTATED_CDS_SET else 'non_canonical'
                        
                        exon_coords = [(int(p.split('-')[0]), int(p.split('-')[1])) for p in coords_str.split(',')]
                        valid_exon_coords = [c for c in exon_coords if c[0] != c[1]]
                        if not valid_exon_coords: continue
                        
                        predicted_start = min(c[0] for c in valid_exon_coords) if strand == '+' else max(c[1] for c in valid_exon_coords) - 3
                        
                        existing_type = gt_collection.get((chrom, predicted_start, strand))
                        if existing_type != 'annotated_CDS':
                            gt_collection[(chrom, predicted_start, strand)] = orf_type
                            
                    except (ValueError, IndexError, KeyError):
                        continue
                ground_truth_library[gt_key] = gt_collection
            except Exception as e:
                print(f"[WARNING] Could not process ground truth file {filename} due to an error: {e}")
                continue

    if not ground_truth_library:
        print(f"[ERROR] No ground truth files from tool '{GROUND_TRUTH_TOOL}' found in '{gt_directory_path}'. Exiting.")
        return None

    print(f"--- Ground Truth Library pre-loading complete. Loaded {len(ground_truth_library)} entries. ---")
    return ground_truth_library

def process_file_worker(task_info: tuple) -> list:
    """处理单个预测文件的工作函数，分別計算 annotated_CDS 和 non_canonical 兩組的結果。"""
    filepath, ground_truth_dict = task_info
    
    filename = os.path.basename(filepath)
    sample, aligner, tools = parse_filename_v2(filename)

    try:
        df = pd.read_csv(filepath, sep='\t', compression='gzip', low_memory=False)
        df = df[df['start_codon'] == 'ATG']
    except Exception:
        df = pd.DataFrame()

    orf_types = ['annotated_CDS', 'non_canonical']
    pred_by_type = {t: set() for t in orf_types}
    gt_by_type = {t: set() for t in orf_types}

    # 1. 劃分金標準 (Ground Truth) 分組
    for tis, orf_type in ground_truth_dict.items():
        gt_by_type[orf_type].add(tis)

    # 2. 劃分當前工具的預測結果分組
    if not df.empty:
        for _, row in df.iterrows():
            try:
                chrom, strand, coords_str = row['chrom'], row['strand'], str(row['coordinate_0base'])
                
                norm_blocks = normalize_coords_string(coords_str)
                orf_type = 'annotated_CDS' if (chrom, strand, norm_blocks) in ANNOTATED_CDS_SET else 'non_canonical'

                exon_coords = [(int(p.split('-')[0]), int(p.split('-')[1])) for p in coords_str.split(',')]
                valid_exon_coords = [c for c in exon_coords if c[0] != c[1]]
                if not valid_exon_coords: continue
                
                predicted_start = min(c[0] for c in valid_exon_coords) if strand == '+' else max(c[1] for c in valid_exon_coords) - 3
                pred_by_type[orf_type].add((chrom, predicted_start, strand))
            except (ValueError, IndexError, KeyError):
                continue

    # 3. 分別計算兩組的 recall 和 validationRate
    results = []
    for orf_type in orf_types:
        pred_set = pred_by_type[orf_type]
        gt_set = gt_by_type[orf_type]
        
        totalORF = len(pred_set)     
        total_gt_tis = len(gt_set)   
        
        unique_hits_set = pred_set.intersection(gt_set)
        tp = len(unique_hits_set)
        
        fp = totalORF - tp
        fn = total_gt_tis - tp
        
        validation_rate = tp / totalORF if totalORF > 0 else 0.0
        recall = tp / total_gt_tis if total_gt_tis > 0 else 0.0
        
        f1_score = 2 * (validation_rate * recall) / (validation_rate + recall) if (validation_rate + recall) > 0 else 0.0
        
        results.append({
            'sample': sample, 
            'soft': aligner, 
            'tools': tools, 
            'ORF_type': orf_type,  
            'validationRate': validation_rate, 
            'recall': recall, 
            'F1-score': f1_score,
            'TP': tp, 
            'FP': fp, 
            'FN': fn, 
            'TIS_hit': tp, 
            'totalORF': totalORF, 
            'TIS': total_gt_tis
        })

    return results


def run_benchmark_for_directory(prediction_dir: str, output_csv: str, ground_truth_library: Dict):
    print("-" * 80)
    print(f"Starting benchmark for directory: {prediction_dir}")
    if not os.path.isdir(prediction_dir): return

    prediction_files = sorted(glob.glob(os.path.join(prediction_dir, "*_gcoor.tsv.gz")))
    if not prediction_files: return
        
    print(f"Found {len(prediction_files)} Ribo-seq files to process. Preparing tasks...")
    
    tasks = []
    for filepath in prediction_files:
        filename = os.path.basename(filepath)
        ribo_sample, ribo_aligner, _ = parse_filename_v2(filename)
        if not all([ribo_sample, ribo_aligner]): continue
        
        ti_sample = None
        for key, value in RIBO_TO_TI_MAPPING.items():
            if key in ribo_sample:
                ti_sample = value
                break
        if ti_sample is None: continue
            
        gt_key = (ti_sample, ribo_aligner)
        ground_truth_dict = ground_truth_library.get(gt_key)
        if ground_truth_dict is None: continue
        
        tasks.append((filepath, ground_truth_dict))

    if not tasks: 
        print(f"No matching TI-seq data found for files in this directory. Skipping.")
        return
        
    pool_size = min(NUM_PROCESSES, len(tasks))
    print(f"Initializing multiprocessing pool with {pool_size} workers for {len(tasks)} tasks...")
    
    results_data = []
    with multiprocessing.Pool(processes=pool_size) as pool:
        worker_results = pool.map(process_file_worker, tasks)
        for res_list in worker_results:
            results_data.extend(res_list)
            
    if not results_data: return

    results_data.sort(key=lambda x: (x['sample'], x['tools'], x['ORF_type']))
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    
    header = ['sample', 'soft', 'tools', 'ORF_type', 'validationRate', 'recall', 'F1-score', 'TP', 'FP', 'FN', 'TIS_hit', 'totalORF', 'TIS']
    
    print(f"Writing {len(results_data)} results to {output_csv}...")
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(results_data)
    print(f"Benchmark for this directory finished successfully.")


def main():
    """主執行函數。"""
    print("--- Starting Batch ORF Prediction Benchmark (TI-seq mode, Grouped by ORF Type) ---")
    
    load_annotated_cds(REF_BLOCK_FILE)
    
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    gt_library = preload_ground_truth_library_final(GROUND_TRUTH_DIR)
    if not gt_library: return

    for dir_name in DIRECTORIES_TO_PROCESS:
        full_prediction_dir = os.path.join(BASE_INPUT_DIR, dir_name)
        
        output_filename = f"{dir_name}_score_TISeq.csv"
        full_output_csv_path = os.path.join(BASE_OUTPUT_DIR, output_filename)
        
        run_benchmark_for_directory(full_prediction_dir, full_output_csv_path, gt_library)

    print("\n--- All specified directories have been processed. ---")


if __name__ == "__main__":
    main()
