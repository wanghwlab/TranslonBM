#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import csv
import multiprocessing
from collections import defaultdict
import itertools
import sys

# --- 1. 配置区  ---
BASE_INPUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/"
TI_SEQ_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/TISeq/"
BASE_OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TISeq_score_union_untrim/"


SAMPLES_TO_PROCESS = [
    "SRX876063_SRX876069",
    "SRX740748",
    "SRX5256543_SRX5256555",
    "SRX5887328_SRX5887329_SRX5887330",
]

DIRECTORIES_TO_PROCESS = [
    #"orf_pred_default_trim",
    "orf_pred_default_untrim"
]

GROUND_TRUTH_TOOL = 'ribotish'

RIBO_TO_TI_MAPPING = {
    "SRX876063_SRX876069": "SRX876066_SRX876072",
    "SRX740748": "SRX740745",
    "SRX5256543_SRX5256555": "SRX5256564_SRX5256565",
    "SRX5887328_SRX5887329_SRX5887330": "SRX5887326",
}

NUM_PROCESSES = 32


# --- 2. 核心功能函数 ---
def parse_filename(filename_str: str) -> tuple:
    """从复杂文件名中解析出 (样本名, 比对软件, 预测工具)。"""
    base = filename_str.replace('_gcoor.tsv.gz', '')
    parts = base.split('_')
    if len(parts) >= 3:
        tools, aligner, sample = parts[-1], parts[-2], '_'.join(parts[:-2])
        return sample, aligner, tools
    return base, 'unknown', 'unknown'

def preload_ti_seq_gold_standards(gt_directory_path: str) -> dict:
    """
    预加载TI-seq金标准库。
    返回一个字典，键是 (样本名, 比对软件)，值是该样本的TIS坐标集合 (set)。
    """
    ground_truth_library = {}
    print(f"--- Pre-loading TI-seq Ground Truth Library (Source Tool: {GROUND_TRUTH_TOOL}) ---")
    
    for dir_name in DIRECTORIES_TO_PROCESS:
        search_path = os.path.join(gt_directory_path, dir_name, f"*_{GROUND_TRUTH_TOOL}_gcoor.tsv.gz")
        gt_files = glob.glob(search_path)
        
        for filepath in gt_files:
            filename = os.path.basename(filepath)
            sample, aligner, _ = parse_filename(filename)
            if not all([sample, aligner]):
                continue
            
            gt_key = (sample, aligner)
            try:
                df = pd.read_csv(filepath, sep='\t', compression='gzip', low_memory=False)
                gt_collection = set()
                for _, row in df.iterrows():
                    try:
                        chrom, strand, coords_str = row['chrom'], row['strand'], str(row['coordinate_0base'])
                        exon_coords = [(int(p.split('-')[0]), int(p.split('-')[1])) for p in coords_str.split(',')]
                        valid_exon_coords = [c for c in exon_coords if c[0] != c[1]]
                        if not valid_exon_coords: continue
                        
                        predicted_start = min(c[0] for c in valid_exon_coords) if strand == '+' else max(c[1] for c in valid_exon_coords) - 3
                        gt_collection.add((chrom, predicted_start, strand))
                    except (ValueError, IndexError, KeyError):
                        continue
                ground_truth_library[gt_key] = gt_collection
            except Exception as e:
                print(f"[WARNING] Could not process ground truth file {filename} due to an error: {e}")
                continue

    if not ground_truth_library:
        print(f"[ERROR] No ground truth files from tool '{GROUND_TRUTH_TOOL}' found in '{gt_directory_path}'. Exiting.")
        sys.exit(1)

    print(f"--- Ground Truth Library pre-loading complete. Loaded {len(ground_truth_library)} entries. ---")
    return ground_truth_library


def load_orfs_from_file(filepath: str) -> tuple:
    """
    从单个文件中加载ORF列表，每个ORF包含其标识和TIS坐标。
    同时根据起始密码子区分为 ATG 和 NTG 列表。
    """
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
                start_codon = row_dict.get('start_codon', 'ATG')
                
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

def calculate_orf_metrics(predicted_orfs: list, truth_set: set) -> dict:
    """
    为给定的ORF列表计算所有指标，FN基于TIS覆盖情况。
    """
    if not predicted_orfs:
        return {'tp': 0, 'fp': 0, 'fn': len(truth_set), 'tis_hit': 0, 'precision': 0.0, 'recall': 0.0, 'fscore': 0.0, 'total_predicted': 0}

    total_predicted_orfs = len(predicted_orfs)
    
    tp = sum(1 for orf_id, tis in predicted_orfs if tis in truth_set)
    fp = total_predicted_orfs - tp
    
    unique_predicted_tis = {tis for orf_id, tis in predicted_orfs}
    tis_hit_count = len(unique_predicted_tis.intersection(truth_set))
    fn = len(truth_set) - tis_hit_count
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    return {'tp': tp, 'fp': fp, 'fn': fn, 'tis_hit': tis_hit_count, 'fscore': f_score, 'precision': precision, 'recall': recall, 'total_predicted': total_predicted_orfs}


def main():
    """主执行函数。"""
    print("--- Starting ORF Combination Performance Analysis (using TI-seq as Ground Truth) ---")

    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    ti_seq_library = preload_ti_seq_gold_standards(TI_SEQ_DIR)
    if not ti_seq_library:
        return

    for sample_prefix in SAMPLES_TO_PROCESS:
        print(f"\n{'='*25}\nProcessing Ribo-seq Sample: {sample_prefix}\n{'='*25}")

        print(f"  [Phase 1] Finding and loading prediction files for '{sample_prefix}'...")
        
        files_for_this_sample = []
        for dir_name in DIRECTORIES_TO_PROCESS:
            search_path = os.path.join(BASE_INPUT_DIR, dir_name, f"{sample_prefix}*_gcoor.tsv.gz")
            files_for_this_sample.extend(glob.glob(search_path))

        if not files_for_this_sample:
            print(f"  Warning: No prediction files found for sample prefix '{sample_prefix}'. Skipping.")
            continue

        print(f"  Found {len(files_for_this_sample)} files. Loading ORFs using {NUM_PROCESSES} cores...")
        
        all_predictions = defaultdict(lambda: defaultdict(dict))
        with multiprocessing.Pool(NUM_PROCESSES) as pool:
            results = pool.map(load_orfs_from_file, files_for_this_sample)
            for res in results:
                if res:
                    sample, aligner, tool, atg_orfs, ntg_orfs = res
                    if sample.startswith(sample_prefix):
                        all_predictions[aligner][sample][tool] = {'ATG': atg_orfs, 'NTG': ntg_orfs}
        
        if not all_predictions:
            print(f"  Warning: Data could not be loaded for '{sample_prefix}'. Skipping.")
            continue
        print(f"  Data loaded for '{sample_prefix}'.")

        print(f"  [Phase 2] Calculating combination metrics for '{sample_prefix}'...")
        
        aggregated_results_for_sample = {'ATG': [], 'NTG': []}

        for aligner, samples_data in all_predictions.items():
            ti_sample_name = RIBO_TO_TI_MAPPING.get(sample_prefix)
            if not ti_sample_name:
                print(f"  Warning: No TI-seq mapping found for Ribo-seq sample '{sample_prefix}'. Skipping aligner '{aligner}'.")
                continue
            
            ground_truth_key = (ti_sample_name, aligner)
            ground_truth_set = ti_seq_library.get(ground_truth_key)

            if not ground_truth_set:
                print(f"  Warning: No TI-seq gold standard found for key ('{ti_sample_name}', '{aligner}'). Skipping analysis for this aligner.")
                continue
            
            print(f"  Processing aligner '{aligner}' with TI-seq gold standard from sample '{ti_sample_name}' ({len(ground_truth_set)} unique TIS).")
            
            for sample, tools_data in samples_data.items():
                tools = sorted(tools_data.keys())
                if len(tools) < 2:
                    print(f"  Info: Less than 2 tools found for sample '{sample}', aligner '{aligner}'. Skipping combination analysis.")
                    continue
                
                for filter_type in ['ATG', 'NTG']:
                    print(f"    Analyzing combinations for {filter_type} ORFs...")
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
                            "total_a": metrics_a['total_predicted'], "tp_a": metrics_a['tp'], "fp_a": metrics_a['fp'], "fn_a": metrics_a['fn'], "tis_hit_a": metrics_a['tis_hit'], "fscore_a": metrics_a['fscore'],
                            "total_b": metrics_b['total_predicted'], "tp_b": metrics_b['tp'], "fp_b": metrics_b['fp'], "fn_b": metrics_b['fn'], "tis_hit_b": metrics_b['tis_hit'], "fscore_b": metrics_b['fscore'],
                            "total_combined": metrics_combined['total_predicted'], "tp_combined": metrics_combined['tp'], "fp_combined": metrics_combined['fp'], "fn_combined": metrics_combined['fn'], "tis_hit_combined": metrics_combined['tis_hit'], "fscore_combined": metrics_combined['fscore'],
                            "fscore_max": fscore_max, "fscore_change": fscore_change, 
                            "extent_of_overlap_jaccard": jaccard_index,
                            "total_TIS_in_gold_standard": len(ground_truth_set)
                        }
                        aggregated_results_for_sample[filter_type].append(result_row)
        
        print(f"  [Phase 3] Aggregating and saving results for '{sample_prefix}'...")
        
        if aggregated_results_for_sample['ATG']:
            atg_df = pd.DataFrame(aggregated_results_for_sample['ATG'])
            output_filename = f"{sample_prefix}_combination_metrics_ATG_vs_TISeq.csv"
            output_path = os.path.join(BASE_OUTPUT_DIR, output_filename)
            atg_df.to_csv(output_path, index=False)
            print(f"  Saved aggregated ATG results for sample '{sample_prefix}' to {output_path}")
        
        if aggregated_results_for_sample['NTG']:
            ntg_df = pd.DataFrame(aggregated_results_for_sample['NTG'])
            output_filename = f"{sample_prefix}_combination_metrics_NTG_vs_TISeq.csv"
            output_path = os.path.join(BASE_OUTPUT_DIR, output_filename)
            ntg_df.to_csv(output_path, index=False)
            print(f"  Saved aggregated NTG results for sample '{sample_prefix}' to {output_path}")

        print(f"--- Finished processing {sample_prefix} ---")

    print("\n--- All specified samples have been processed. Script Finished Successfully ---")


if __name__ == "__main__":
    main()
