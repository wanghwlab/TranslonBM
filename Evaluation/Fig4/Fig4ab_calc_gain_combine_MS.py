#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import numpy as np

# --- 1. 配置区 ---
UNION_CSV_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/Maxquant_score_union_trim/" 
INTERSECT_CSV_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/Maxquant_score_intersect_trim/"
OUTPUT_FILE = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/merged_gain/Maxquant_combination_gains_trim_peptide.csv"

def calculate_gains():
    print("正在读取 CSV 文件...")
    
    union_files = glob.glob(os.path.join(UNION_CSV_DIR, "*ombined_metrics*peptide.csv"))
    intersect_files = glob.glob(os.path.join(INTERSECT_CSV_DIR, "combined_metrics*peptide.csv"))

    if not union_files or not intersect_files:
        print("错误: 未找到足够的 CSV 文件，请检查路径是否正确。")
        return

    df_union = pd.concat([pd.read_csv(f) for f in union_files], ignore_index=True)
    df_intersect = pd.concat([pd.read_csv(f) for f in intersect_files], ignore_index=True)

    print(f"共加载了 {len(df_union)} 条 Union 记录，{len(df_intersect)} 条 Intersect 记录。")

    merge_keys = ['aligner', 'sample', 'tool_a', 'tool_b']

    df_intersect_subset = df_intersect[merge_keys + ['precision_combined', 'recall_combined']].copy()

    df_intersect_subset.rename(columns={
        'precision_combined': 'precision_combined_intersect',
        'recall_combined': 'recall_combined_intersect'
    }, inplace=True)

    df_union.rename(columns={
        'precision_combined': 'precision_combined_union',
        'recall_combined': 'recall_combined_union',
        'fscore_combined': 'fscore_combined_union'
    }, inplace=True)

    print("正在合并两组数据...")
    df_merged = pd.merge(df_union, df_intersect_subset, on=merge_keys, how='inner')

    if df_merged.empty:
        print("警告: 合并后数据为空！请确保两组数据的 aligner, sample, tool_a, tool_b 完全匹配。")
        return

    print("正在计算 Combination Gains...")

    df_merged['max_recall'] = np.maximum(df_merged['recall_a'], df_merged['recall_b'])
    df_merged['max_precision'] = np.maximum(df_merged['precision_a'], df_merged['precision_b'])

    df_merged['Recall_gain_union'] = df_merged['recall_combined_union'] - df_merged['max_recall']

    df_merged['Recall_gain_intersect'] = df_merged['recall_combined_intersect'] - df_merged['max_recall']

    df_merged['Precision_gain_union'] = df_merged['precision_combined_union'] - df_merged['max_precision']

    df_merged['Precision_gain_intersect'] = df_merged['precision_combined_intersect'] - df_merged['max_precision']

    front_columns = [
        'aligner', 'sample', 'tool_a', 'tool_b',
        'precision_a', 'precision_b', 'max_precision',
        'precision_combined_union', 'Precision_gain_union',
        'precision_combined_intersect', 'Precision_gain_intersect',
        'recall_a', 'recall_b', 'max_recall',
        'recall_combined_union', 'Recall_gain_union',
        'recall_combined_intersect', 'Recall_gain_intersect',
        'extent_of_overlap_jaccard'
    ]

    other_columns = [col for col in df_merged.columns if col not in front_columns]
    final_columns = front_columns + other_columns

    df_merged = df_merged[final_columns]

    df_merged.to_csv(OUTPUT_FILE, index=False)
    print(f"\n计算完成！")
    print(f"有效配对数据共 {len(df_merged)} 条。")
    print(f"结果已成功保存至: {OUTPUT_FILE}")

if __name__ == "__main__":
    calculate_gains()
