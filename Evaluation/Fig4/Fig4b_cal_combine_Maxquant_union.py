#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
import itertools
import re
from collections import defaultdict

# --- 1. 配置区 ---

WORKDIR = "/home/tangyuewen/ORF_benchmark/Maxquant_2026.1"
BASE_OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/Maxquant_score_union_untrim/"

SAMPLES = ["SRX1254413", "SRX11812007_SRX11812008_SRX11812009", "SRX5887328_SRX5887329_SRX5887330"]#, "SRX1254413", "SRX11812007_SRX11812008_SRX11812009", "SRX5887328_SRX5887329_SRX5887330"


ALIGNERS = ["tophat2", "STAR", "hisat2"]

PRED_FLAG = "orf_pred_default"
TRIM_FLAG = "untrim"

FDR = ""


# --- 2. 辅助功能函数 ---
def is_all_seq(protein_id_str: str) -> bool:
    """检查蛋白质ID是否全部为预测的seq_ ID"""
    if not isinstance(protein_id_str, str): return False
    main_ids = protein_id_str.split('|')[0]
    return all(id_val.startswith('seq_') for id_val in main_ids.split(';'))

def has_any_unique_peptides(count_str: str) -> bool:
    """检查蛋白质组是否含有任何独特肽段"""
    if not isinstance(count_str, str): return False
    try:
        return any(int(c) > 0 for c in count_str.split(';'))
    except (ValueError, TypeError):
        return False

def extract_seq_ids(id_string: str) -> list:
    """从字符串中提取所有seq_ ID"""
    if not isinstance(id_string, str): return []
    return re.findall(r'seq_\d+', id_string)


# --- 3. 核心功能模块 ---
def load_all_predictions_from_ref(sample: str, pred_flag: str, trim_flag: str) -> dict:
    """
    ★ 新的核心函数: 从 ref_index.tsv 文件一次性加载所有工具的预测结果。
    返回一个按 aligner 和 tool 组织好的嵌套字典。
    结构: {aligner: {tool: {set_of_sequence_ids}}}
    """
    ref_path = os.path.join(WORKDIR, "maxquant_ref", f"{pred_flag}_{trim_flag}", f"{sample}_ref_index.tsv")
    print(f"Loading all predictions from reference file: {ref_path}")
    try:
        full_ref_df = pd.read_csv(ref_path, sep='\t')
        
        all_predictions = defaultdict(lambda: defaultdict(set))

        for (aligner, tool), group in full_ref_df.groupby(['aligner', 'orf_detector']):
            all_predictions[aligner][tool].update(group['sequence_id'].dropna())
        
        print(f"-> Successfully loaded predictions for {len(all_predictions)} aligners.")
        return all_predictions

    except FileNotFoundError:
        print(f"[ERROR] Reference file not found: {ref_path}. Skipping sample {sample}.")
        return None
    except Exception as e:
        print(f"[ERROR] Failed to load reference file {ref_path}: {e}")
        return None

def load_ms_gold_standard(sample: str, aligner: str, level: str, pred_flag: str, trim_flag: str) -> tuple:
    """
    加载MaxQuant的质谱验证结果作为金标准 (逻辑保持不变)。
    返回一个包含(金标准ID集合, 金标准总数)的元组。
    """
    print(f"  Loading MS Gold Standard for: {sample} | {aligner} | {level.capitalize()} Level")
    rerun_subdir = f"{pred_flag}_{trim_flag}"
    
    try:
        if level == 'peptide':
            ms_data_path = os.path.join(WORKDIR, "maxquant_rerun", rerun_subdir, FDR, aligner, sample, "combined", "txt", "peptides.txt")
            ms_df = pd.read_csv(ms_data_path, sep='\t', low_memory=False)
            valid_items = ms_df[(ms_df['Unique (Proteins)'] == 'yes') & (ms_df['Potential contaminant'].isna()) & (ms_df['Reverse'].isna())].copy()
            id_column = 'Proteins'
        else:  # protein level
            ms_data_path = os.path.join(WORKDIR, "maxquant_rerun", rerun_subdir, FDR, aligner, sample, "combined", "txt", "proteinGroups.txt")
            ms_df = pd.read_csv(ms_data_path, sep='\t', low_memory=False)
            has_mask = ms_df['Peptide counts (unique)'].apply(has_any_unique_peptides)
            valid_items = ms_df[(ms_df['Potential contaminant'].isna()) & (ms_df['Reverse'].isna()) & (has_mask)].copy()
            valid_items = valid_items[valid_items['Protein IDs'].apply(is_all_seq)]
            id_column = 'Protein IDs'
        
        if valid_items.empty:
            print(f"  [Warning] No valid MS items found for {sample}/{aligner}/{level}. Gold standard is empty.")
            return set(), 0
            
        validated_seq_ids = set(valid_items[id_column].apply(extract_seq_ids).explode().dropna())
        ms_valid_orfs_num = len(validated_seq_ids)
        
        print(f"  -> Found {ms_valid_orfs_num} unique validated ORF IDs in gold standard.")
        return validated_seq_ids, ms_valid_orfs_num

    except FileNotFoundError as e:
        print(f"  [Error] Gold standard file not found: {e}. Skipping.")
        return None, None
    except Exception as e:
        print(f"  [Error] Failed to load gold standard for {sample}/{aligner}/{level}: {e}")
        return None, None

def calculate_metrics(predicted_set: set, gold_standard_set: set, gold_standard_total: int) -> dict:
    """基于预测集和金标准集计算性能指标 (逻辑保持不变)。"""
    tp = len(predicted_set.intersection(gold_standard_set))
    fp = len(predicted_set) - tp
    fn = gold_standard_total - tp
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    return {'tp': tp, 'fp': fp, 'fn': fn, 'fscore': f_score, 'precision': precision, 'recall': recall}


# --- 4. 主程序  ---
def main():
    """主执行函数，采用与原始脚本一致的数据加载和统计逻辑。"""
    print("--- Starting ORF Combination Performance Analysis (Corrected Workflow) ---")
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    for sample in SAMPLES:
        print(f"\n{'='*30}\nProcessing Sample: {sample}\n{'='*30}")
        
        predictions_data = load_all_predictions_from_ref(sample, PRED_FLAG, TRIM_FLAG)
        if not predictions_data:
            continue

        for level in ['peptide', 'protein']:
            print(f"\n--- Analyzing at {level.upper()} Level ---")
            
            all_results_for_sample = [] 

            for aligner, tool_predictions in predictions_data.items():
                print(f"\n-- Aligner: {aligner} --")
                
                gold_standard_set, gold_standard_total = load_ms_gold_standard(sample, aligner, level, PRED_FLAG, TRIM_FLAG)
                if gold_standard_set is None or gold_standard_total == 0:
                    print(f"Skipping combination analysis for {aligner} due to missing or empty gold standard.")
                    continue
                
                if len(tool_predictions) < 2:
                    print(f"  Warning: Need at least 2 ORF tools for comparison. Found {len(tool_predictions)}. Skipping aligner.")
                    continue
                
                tools = sorted(tool_predictions.keys())
                for tool_a, tool_b in itertools.combinations(tools, 2):
                    orfs_a = tool_predictions[tool_a]
                    orfs_b = tool_predictions[tool_b]
                    
                    metrics_a = calculate_metrics(orfs_a, gold_standard_set, gold_standard_total)
                    metrics_b = calculate_metrics(orfs_b, gold_standard_set, gold_standard_total)
                    
                    orfs_combined = orfs_a.union(orfs_b)
                    metrics_combined = calculate_metrics(orfs_combined, gold_standard_set, gold_standard_total)
                    
                    fscore_max = max(metrics_a['fscore'], metrics_b['fscore'])
                    fscore_change = metrics_combined['fscore'] - fscore_max
                    
                    intersection = len(orfs_a.intersection(orfs_b))
                    union = len(orfs_combined)
                    jaccard_index = intersection / union if union > 0 else 0
                    
                    result_row = {
                        "aligner": aligner, "sample": sample, "tool_a": tool_a, "tool_b": tool_b,
                        "tp_a": metrics_a['tp'], "fp_a": metrics_a['fp'], "fn_a": metrics_a['fn'], "precision_a": metrics_a['precision'], "recall_a": metrics_a['recall'], "fscore_a": metrics_a['fscore'],
                        "tp_b": metrics_b['tp'], "fp_b": metrics_b['fp'], "fn_b": metrics_b['fn'], "precision_b": metrics_b['precision'], "recall_b": metrics_b['recall'], "fscore_b": metrics_b['fscore'],
                        "tp_combined": metrics_combined['tp'], "fp_combined": metrics_combined['fp'], "fn_combined": metrics_combined['fn'], "precision_combined": metrics_combined['precision'], "recall_combined": metrics_combined['recall'], "fscore_combined": metrics_combined['fscore'],
                        "fscore_max": fscore_max, "fscore_change": fscore_change, 
                        "extent_of_overlap_jaccard": jaccard_index,
                        "gold_standard_total_orfs": gold_standard_total,
                    }
                    all_results_for_sample.append(result_row)
            
            if not all_results_for_sample:
                print(f"No results generated for sample '{sample}' at {level} level.")
                continue

            output_df = pd.DataFrame(all_results_for_sample)
            output_filename = f"{sample}_combination_metrics_{FDR}_{level}.csv"
            output_path = os.path.join(BASE_OUTPUT_DIR, output_filename)
            output_df.to_csv(output_path, index=False)
            print(f"\n>>> Successfully saved results for sample '{sample}' at {level} level to: {output_path}")

    print("\n--- All specified samples have been processed. Script Finished Successfully ---")

if __name__ == '__main__':
    main()
