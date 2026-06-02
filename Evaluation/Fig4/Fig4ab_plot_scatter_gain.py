#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use('Agg') 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob
import os
import multiprocessing
from matplotlib.backends.backend_pdf import PdfPages

# --- 1. 配置区 ---
INPUT_DIR = '/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/Mean_score'
OUTPUT_DIR = '/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/Mean_plots_Scatter'

NUM_PROCESSES = 8

TOOL_ORDER = [
    'ribotaper', 'ribohmm', 'orfrater', 'ribotish', 'gedi', 
    'ribocode', 'riborf', 'ribowave', 'rpbp', 'ribotricer', 'ORFquant'
]
LETTERS = "ABCDEFGHIJK"
TOOL_MAP = dict(zip(TOOL_ORDER, LETTERS))

COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#393b79'
]
TOOL_COLOR_MAP = dict(zip(TOOL_ORDER, COLORS))
LETTER_COLOR_MAP = {TOOL_MAP[tool]: TOOL_COLOR_MAP[tool] for tool in TOOL_ORDER}

# --- 2. 核心函数 ---
def get_combo_label(t1, t2):
    """根据两个软件的名称，返回按字母顺序排列的组合标签（例如 'AK'）"""
    if t1 not in TOOL_MAP or t2 not in TOOL_MAP:
        return ""
    l1, l2 = TOOL_MAP[t1], TOOL_MAP[t2]
    return "".join(sorted([l1, l2]))

def process_aligner(args):
    """为单个 Aligner 绘制 3x2 的散点图并保存为 PDF"""
    (aligner, trim_status), df_aligner = args
    print(f"正在处理: {aligner} [{trim_status}]")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pdf_path = os.path.join(OUTPUT_DIR, f"{aligner}_Combination_Scatter_Gains_{trim_status}.pdf")

    fig, axes = plt.subplots(3, 2, figsize=(14, 18))
    levels = ['TIS', 'TISeq', 'MS']
    strats = ['union', 'intersect']
    
    for i, level in enumerate(levels):
        df_level = df_aligner[df_aligner['level'] == level]
        
        for j, strat in enumerate(strats):
            ax = axes[i, j]
            strat_Display = strat.capitalize()
            
            if df_level.empty:
                ax.set_title(f"{level} - {strat_Display}\n(No Data Found)", fontsize=14)
                ax.axis('off')
                continue
                
            x_col = f'Precision_gain_{strat}'
            y_col = f'Recall_gain_{strat}'
            
            if x_col not in df_level.columns or y_col not in df_level.columns:
                ax.set_title(f"{level} - {strat_Display}\n(Missing Gain Columns)", fontsize=14)
                ax.axis('off')
                continue
            
            grouped = df_level.groupby(['tool_a', 'tool_b'])[[x_col, y_col]].mean().reset_index()
            
            ax.scatter(grouped[x_col], grouped[y_col], alpha=0)
            
            ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
            ax.axvline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
            
            for _, row in grouped.iterrows():
                x_val = row[x_col]
                y_val = row[y_col]
                label = get_combo_label(row['tool_a'], row['tool_b'])
                
                if not label or pd.isna(x_val) or pd.isna(y_val):
                    continue
                
                l1, l2 = label[0], label[1]
                c1 = LETTER_COLOR_MAP[l1]
                c2 = LETTER_COLOR_MAP[l2]
                
                ax.text(x_val, y_val, l1, color=c1, fontsize=10, ha='right', va='center', weight='bold', alpha=0.9)
                ax.text(x_val, y_val, l2, color=c2, fontsize=10, ha='left', va='center', weight='bold', alpha=0.9)
                
            ax.set_title(f"{level} Level - {strat_Display} Strategy", fontsize=14, pad=10, weight='bold')
            ax.set_xlabel('Precision Gain', fontsize=12)
            ax.set_ylabel('Recall Gain', fontsize=12)
            ax.grid(True, linestyle=':', alpha=0.5)
            
    legend_patches = [
        mpatches.Patch(color=TOOL_COLOR_MAP[tool], label=f"{TOOL_MAP[tool]} : {tool}") 
        for tool in TOOL_ORDER
    ]
    
    fig.legend(handles=legend_patches, title="Tool Map (A-K)", 
               loc='center right', bbox_to_anchor=(0.98, 0.5), 
               fontsize=12, title_fontsize=14,
               facecolor='whitesmoke', edgecolor='lightgray', framealpha=0.9)
             
    plt.tight_layout(rect=[0, 0, 0.88, 1])
    fig.suptitle(f"Combination Gain Scatter Analysis: {aligner} ({trim_status})", fontsize=20, y=0.98, weight='bold')
    plt.subplots_adjust(top=0.92, hspace=0.3, wspace=0.2)

    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
    print(f"  - 完成保存: {os.path.basename(pdf_path)}")


def main():
    csv_files = glob.glob(os.path.join(INPUT_DIR, 'merged_gain*.csv'))
    if not csv_files:
        print(f"错误：在目录 '{INPUT_DIR}' 中找不到 CSV 文件。")
        return
        
    df_list = []
    for f in csv_files:
        try:
            temp_df = pd.read_csv(f)
            base_name = os.path.basename(f).upper()
            if 'TISEQ' in base_name:
                temp_df['level'] = 'TISeq'
            elif 'TIS' in base_name:
                temp_df['level'] = 'TIS'
            elif 'MS' in base_name:
                temp_df['level'] = 'MS'
            else:
                temp_df['level'] = 'Unknown'
                
            if 'UNTRIM' in base_name:
                temp_df['trim_status'] = 'untrim'
            elif 'TRIM' in base_name:
                temp_df['trim_status'] = 'trim'
            else:
                temp_df['trim_status'] = 'unknown'
                
            df_list.append(temp_df)
        except Exception as e:
            print(f"读取 {f} 时出错: {e}")
            
    if not df_list: return
    
    full_df = pd.concat(df_list, ignore_index=True)
    
    tasks = []
    for (al, ts), group in full_df.groupby(['aligner', 'trim_status']):
        tasks.append(((al, ts), group.copy()))
        
    print(f"找到 {len(tasks)} 个 Aligner/Trim 组合进行绘图...")
    
    with multiprocessing.Pool(min(NUM_PROCESSES, len(tasks))) as pool:
        pool.map(process_aligner, tasks)
        
    print("\n所有 Scatter 绘图任务处理完成！图表保存在:", OUTPUT_DIR)
if __name__ == '__main__':
    main()
