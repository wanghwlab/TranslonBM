#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import sys
import pandas as pd
import multiprocessing
from collections import defaultdict
import itertools

# --- 1. 配置区---
BASE_INPUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/"
BASE_OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TIS_score/"
GTF_FILE = "/home/tangyuewen/ORF_benchmark/Ref/gencode.v43.annotation_TIS.gtf"

SAMPLES_TO_PROCESS = [
    #"SRX1447296",
    "SRX876063_SRX876069",
    "simulation_6M_T1",
    "SRX1254413",  
    "SRX5256543_SRX5256555",
"SRX740748",
"SRX5887328_SRX5887329_SRX5887330",
"SRX11812007_SRX11812008_SRX11812009",
"simulation_6M_T3",
"simulation_60M_T1",
"simulation_60M_T3",
#"SRX7666669-73",
#"SRX7666674-78",
#"SRX7666679-83",
#"SRX7666684-88",
#"SRX7666689-93",
#"SRX7666694-98",
]

DIRECTORIES_TO_PROCESS = [
    #"orf_pred_default_trim",
    "orf_pred_default_untrim"
    
]

NUM_PROCESSES = 32

# --- 2. 核心功能函数 ---
def load_ground_truth(gtf_path: str) -> set:
    """從GTF文件中加載真實的TIS坐標。"""
    ground_truth_set = set()
    print("Loading ground truth TIS from GTF file...")
    try:
        with open(gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\t')
                if len(fields) > 3 and fields[2] == 'start_codon':
                    chrom, start_1based, strand = fields[0], int(fields[3]), fields[6]
                    ground_truth_set.add((chrom, start_1based - 1, strand))
    except FileNotFoundError:
        print(f"[ERROR] GTF file not found at: {gtf_path}. Exiting.", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(ground_truth_set)} unique TIS records.")
    return ground_truth_set

def parse_filename(filename_str: str) -> tuple:
    """從文件名中解析出 sample, soft (aligner), tools。"""
    base = filename_str.replace('_gcoor.tsv.gz', '')
    parts = base.split('_')
    if len(parts) >= 3:
        tools, soft, sample = parts[-1], parts[-2], '_'.join(parts[:-2])
        return sample, soft, tools
    return base, 'unknown', 'unknown'

def calculate_orf_metrics(predicted_orfs: list, truth_set: set) -> dict:
    """
    为给定的ORF列表计算所有指标。
    - TP/FP 基于ORF的数量。
    - FN 基于去重后TIS的覆盖情况，采用 FN = |Truth Set| - |TIS_hit|。
    """
    if not predicted_orfs:
        return {'tp': 0, 'fp': 0, 'fn': len(truth_set), 'tis_hit': 0, 'precision': 0.0, 'recall': 0.0, 'fscore': 0.0}

    total_predicted_orfs = len(predicted_orfs)

    tp = sum(1 for orf_id, tis in predicted_orfs if tis in truth_set)

    fp = total_predicted_orfs - tp

    unique_predicted_tis = {tis for orf_id, tis in predicted_orfs}
    tis_hit_count = len(unique_predicted_tis.intersection(truth_set))
    fn = len(truth_set) - tis_hit_count

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    return {'tp': tp, 'fp': fp, 'fn': fn, 'tis_hit': tis_hit_count, 'fscore': f_score}
    #return {'tp': tp, 'fp': fp, 'fn': fn, 'fscore': f_score}


def load_orfs_from_file(filepath: str) -> tuple:
    """从单个文件中加载ORF列表，每个ORF包含其标识和TIS坐标。"""
    try:
        sample, aligner, tool = parse_filename(os.path.basename(filepath))
        if not all([sample, aligner, tool]): return None
        
        df = pd.read_csv(filepath, sep='\t', compression='gzip', low_memory=False)
        
        atg_orfs = []
        ntg_orfs = []

        for row_tuple in df.itertuples(index=False):
            try:
                row_dict = row_tuple._asdict()
                strand = row_dict['strand']
                coords_str = str(row_dict['coordinate_0base'])
                start_codon = row_dict['start_codon']
                
                exon_coords = [(int(p.split('-')[0]), int(p.split('-')[1])) for p in coords_str.split(',')]
                valid_exon_coords = [c for c in exon_coords if c[0] != c[1]]
                if not valid_exon_coords: continue
                
                start = min(c[0] for c in valid_exon_coords) if strand == '+' else max(c[1] for c in valid_exon_coords) - 3
                tis = (row_dict['chrom'], start, strand)
                
                orf_identifier = tuple(row_tuple)

                if start_codon == 'ATG':
                    atg_orfs.append((orf_identifier, tis))
                else:
                    ntg_orfs.append((orf_identifier, tis))

            except (ValueError, IndexError, KeyError):
                continue
        return (sample, aligner, tool, atg_orfs, ntg_orfs)
    except Exception as e:
        print(f"Warning: Could not process {filepath}. Error: {e}")
        return None

def main():
    """主执行函数 (内存优化版，按样本循环处理)。"""
    print("--- Starting ORF Combination Performance Analysis (Memory-Optimized Workflow) ---")

    ground_truth_set = load_ground_truth(GTF_FILE)
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    for sample_prefix in SAMPLES_TO_PROCESS:
        print(f"\n{'='*25}\nProcessing Sample: {sample_prefix}\n{'='*25}")

        print(f"  [Phase 1] Finding and loading files for '{sample_prefix}'...")
        
        files_for_this_sample = []
        for dir_name in DIRECTORIES_TO_PROCESS:
            search_path = os.path.join(BASE_INPUT_DIR, dir_name, f"{sample_prefix}*_gcoor.tsv.gz")
            found_files = glob.glob(search_path)
            files_for_this_sample.extend(found_files)

        if not files_for_this_sample:
            print(f"  Warning: No files found for sample prefix '{sample_prefix}'. Skipping.")
            continue  

        print(f"  Found {len(files_for_this_sample)} files. Loading ORFs using {NUM_PROCESSES} cores...")
        
        all_predictions = defaultdict(lambda: defaultdict(dict))
        with multiprocessing.Pool(NUM_PROCESSES) as pool:
            results = pool.map(load_orfs_from_file, files_for_this_sample)
            for res in results:
                if res:
                    sample, aligner, tool, atg_orfs_list, ntg_orfs_list = res
                    if sample.startswith(sample_prefix):
                        all_predictions[aligner][sample][tool] = {
                            'ATG': atg_orfs_list,
                            'NTG': ntg_orfs_list
                        }
        
        print(f"  Data loaded for '{sample_prefix}'.")

        print(f"  [Phase 2] Calculating metrics for '{sample_prefix}'...")

        sample_results = defaultdict(lambda: {'ATG': [], 'NTG': []})
        for filter_type in ['ATG', 'NTG']:
            for aligner, samples in all_predictions.items():
                for sample, tools_data in samples.items():
                    tools = sorted(tools_data.keys())
                    if len(tools) < 2: continue
                    
                    for tool_a, tool_b in itertools.combinations(tools, 2):
                        orfs_a = tools_data[tool_a][filter_type]
                        orfs_b = tools_data[tool_b][filter_type]
                        metrics_a = calculate_orf_metrics(orfs_a, ground_truth_set)
                        metrics_b = calculate_orf_metrics(orfs_b, ground_truth_set)
                        combined_orfs_dict = {orf_id: tis for orf_id, tis in orfs_a + orfs_b}
                        orfs_combined_deduplicated = list(combined_orfs_dict.items())
                        metrics_combined = calculate_orf_metrics(orfs_combined_deduplicated, ground_truth_set)
                        fscore_max = max(metrics_a['fscore'], metrics_b['fscore'])
                        fscore_change = metrics_combined['fscore'] - fscore_max
                        set_a = {orf_id for orf_id, tis in orfs_a}
                        set_b = {orf_id for orf_id, tis in orfs_b}
                        intersection = len(set_a.intersection(set_b))
                        union = len(set_a.union(set_b))
                        jaccard_index = intersection / union if union > 0 else 0
                        result_row = {
                            "aligner": aligner, "sample": sample, "tool_a": tool_a, "tool_b": tool_b,
                            "tp_a": metrics_a['tp'], "fp_a": metrics_a['fp'], "fn_a": metrics_a['fn'], 
                            "tis_hit_a": metrics_a['tis_hit'], "fscore_a": metrics_a['fscore'],
                            "tp_b": metrics_b['tp'], "fp_b": metrics_b['fp'], "fn_b": metrics_b['fn'], 
                            "tis_hit_b": metrics_b['tis_hit'], "fscore_b": metrics_b['fscore'],
                            "tp_combined": metrics_combined['tp'], "fp_combined": metrics_combined['fp'], "fn_combined": metrics_combined['fn'], 
                            "tis_hit_combined": metrics_combined['tis_hit'], "fscore_combined": metrics_combined['fscore'],
                            "fscore_max": fscore_max, "fscore_change": fscore_change, 
                            "extent_of_overlap_jaccard": jaccard_index
                        }
                        sample_results[sample][filter_type].append(result_row)
        
        print(f"  [Phase 3] Saving results for '{sample_prefix}'...")
        if not sample_results:
            print("  No results were generated for this sample.")
        
        for sample_name, results in sample_results.items():
            if results['ATG']:
                atg_df = pd.DataFrame(results['ATG'])
                output_filename = f"{sample_name}_combination_metrics_ATG.csv"
                output_path = os.path.join(BASE_OUTPUT_DIR, output_filename)
                atg_df.to_csv(output_path, index=False)
                print(f"  Saved ATG results for sample '{sample_name}' to {output_path}")
            
            if results['NTG']:
                ntg_df = pd.DataFrame(results['NTG'])
                output_filename = f"{sample_name}_combination_metrics_NTG.csv"
                output_path = os.path.join(BASE_OUTPUT_DIR, output_filename)
                ntg_df.to_csv(output_path, index=False)
                print(f"  Saved NTG results for sample '{sample_name}' to {output_path}")

        print(f"--- Finished processing {sample_prefix} ---")

    print("\n--- All specified samples have been processed. Script Finished Successfully ---")


if __name__ == "__main__":
    main()

