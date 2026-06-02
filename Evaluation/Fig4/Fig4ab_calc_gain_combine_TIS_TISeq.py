#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import numpy as np

# --- 1. 配置区 ---
TASKS = [
    {
        "name": "TISeq_untrim",
        "union_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TISeq_score_union_untrim/",
        "intersect_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TISeq_score_intersect_untrim/",
        "output_file": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/merged_gain/TISeq_combination_gains_untrim.csv"
    },
    {
        "name": "TIS_trim",
        "union_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TIS_score_union_trim/",
        "intersect_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TIS_score_intersect_trim/",
        "output_file": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/merged_gain/TIS_combination_gains_trim.csv"
    },
    {
        "name": "TISeq_trim",
        "union_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TISeq_score_union_trim/",
        "intersect_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TISeq_score_intersect_trim/",
        "output_file": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/merged_gain/TISeq_combination_gains_trim.csv"
    },
    {
        "name": "TIS_untrim",
        "union_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TIS_score_union_untrim/",
        "intersect_dir": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/TIS_score_intersect_untrim/",
        "output_file": "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/merged_gain/TIS_combination_gains_untrim.csv"
    }
]

# --- 2. 核心计算函数 ---

def calculate_precision_recall(df):
    """
    根据 tp, fp, fn 动态计算 precision, recall 和 f1-score，确保数值一致性。
    """
    df['precision_a'] = np.where((df['tp_a'] + df['fp_a']) > 0, df['tp_a'] / (df['tp_a'] + df['fp_a']), 0.0)
    df['recall_a'] = np.where((df['tp_a'] + df['fn_a']) > 0, df['tp_a'] / (df['tp_a'] + df['fn_a']), 0.0)
    df['fscore_a'] = np.where((df['precision_a'] + df['recall_a']) > 0, 
                              2 * (df['precision_a'] * df['recall_a']) / (df['precision_a'] + df['recall_a']), 0.0)
    
    df['precision_b'] = np.where((df['tp_b'] + df['fp_b']) > 0, df['tp_b'] / (df['tp_b'] + df['fp_b']), 0.0)
    df['recall_b'] = np.where((df['tp_b'] + df['fn_b']) > 0, df['tp_b'] / (df['tp_b'] + df['fn_b']), 0.0)
    df['fscore_b'] = np.where((df['precision_b'] + df['recall_b']) > 0, 
                              2 * (df['precision_b'] * df['recall_b']) / (df['precision_b'] + df['recall_b']), 0.0)
    
    df['precision_combined'] = np.where((df['tp_combined'] + df['fp_combined']) > 0, df['tp_combined'] / (df['tp_combined'] + df['fp_combined']), 0.0)
    df['recall_combined'] = np.where((df['tp_combined'] + df['fn_combined']) > 0, df['tp_combined'] / (df['tp_combined'] + df['fn_combined']), 0.0)
    df['fscore_combined'] = np.where((df['precision_combined'] + df['recall_combined']) > 0, 
                                     2 * (df['precision_combined'] * df['recall_combined']) / (df['precision_combined'] + df['recall_combined']), 0.0)
    return df

def process_single_task(task):
    """处理单个数据集的合并与增益计算"""
    print(f"{'='*40}")
    print(f"正在处理任务: {task['name']}")
    
    union_files = glob.glob(os.path.join(task['union_dir'], "*combined_metrics*.csv"))
    intersect_files = glob.glob(os.path.join(task['intersect_dir'], "*combined_metrics*.csv"))

    if not union_files or not intersect_files:
        print(f"  [跳过] 未找到足够的 CSV 文件。\n  Union目录: {task['union_dir']}\n  Intersect目录: {task['intersect_dir']}")
        return

    df_union = pd.concat([calculate_precision_recall(pd.read_csv(f)) for f in union_files], ignore_index=True)
    df_intersect = pd.concat([calculate_precision_recall(pd.read_csv(f)) for f in intersect_files], ignore_index=True)

    print(f"  共加载了 {len(df_union)} 条 Union 记录，{len(df_intersect)} 条 Intersect 记录。")

    merge_keys = ['aligner', 'sample', 'tool_a', 'tool_b']

    df_intersect_subset = df_intersect[merge_keys + ['precision_combined', 'recall_combined', 'fscore_combined']].copy()
    df_intersect_subset.rename(columns={
        'precision_combined': 'precision_combined_intersect',
        'recall_combined': 'recall_combined_intersect',
        'fscore_combined': 'fscore_combined_intersect'
    }, inplace=True)

    df_union.rename(columns={
        'precision_combined': 'precision_combined_union',
        'recall_combined': 'recall_combined_union',
        'fscore_combined': 'fscore_combined_union'
    }, inplace=True)

    df_merged = pd.merge(df_union, df_intersect_subset, on=merge_keys, how='inner')

    if df_merged.empty:
        print("  [警告] 合并后数据为空！请确保两组数据的 aligner, sample, tool_a, tool_b 完全匹配。")
        return

    df_merged['max_recall'] = np.maximum(df_merged['recall_a'], df_merged['recall_b'])
    df_merged['max_precision'] = np.maximum(df_merged['precision_a'], df_merged['precision_b'])
    df_merged['max_fscore'] = np.maximum(df_merged['fscore_a'], df_merged['fscore_b'])

    df_merged['Recall_gain_union'] = df_merged['recall_combined_union'] - df_merged['max_recall']
    df_merged['Recall_gain_intersect'] = df_merged['recall_combined_intersect'] - df_merged['max_recall']

    df_merged['Precision_gain_union'] = df_merged['precision_combined_union'] - df_merged['max_precision']
    df_merged['Precision_gain_intersect'] = df_merged['precision_combined_intersect'] - df_merged['max_precision']

    df_merged['F1score_gain_union'] = df_merged['fscore_combined_union'] - df_merged['max_fscore']
    df_merged['F1score_gain_intersect'] = df_merged['fscore_combined_intersect'] - df_merged['max_fscore']

    front_columns = [
        'aligner', 'sample', 'tool_a', 'tool_b',
        'precision_a', 'precision_b', 'max_precision',
        'precision_combined_union', 'Precision_gain_union',
        'precision_combined_intersect', 'Precision_gain_intersect',
        'recall_a', 'recall_b', 'max_recall',
        'recall_combined_union', 'Recall_gain_union',
        'recall_combined_intersect', 'Recall_gain_intersect',
        'fscore_a', 'fscore_b', 'max_fscore',
        'fscore_combined_union', 'F1score_gain_union',
        'fscore_combined_intersect', 'F1score_gain_intersect',
        'extent_of_overlap_jaccard'
    ]

    other_columns = [col for col in df_merged.columns if col not in front_columns]
    final_columns = front_columns + other_columns
    df_merged = df_merged[final_columns]

    os.makedirs(os.path.dirname(task['output_file']), exist_ok=True)
    df_merged.to_csv(task['output_file'], index=False)
    print(f"  计算完成！有效配对数据共 {len(df_merged)} 条。")
    print(f"  结果已保存至: {task['output_file']}")

if __name__ == "__main__":
    print("开始执行组合增益批量计算...")
    for task in TASKS:
        process_single_task(task)
    print("\n所有任务处理完毕。")
