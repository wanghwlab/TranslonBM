#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- 1. 配置区域 ---
INPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/tools_overlap/merged_ATG/orf_pred_default_untrim/"
OUTPUT_CSV = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap_study/overlap_summary_stats_real_untrim.csv"
OUTPUT_PDF = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap_study/overlap_donut_charts_real_untrim.pdf"

DATASET_MAPPING = {
    "SRX876063_SRX876069": "Ji et al. (2015)",
    "SRX740748": "Gao et al. (2015)",
    "SRX1254413": "Calviello et al. (2016)",
    "SRX5256543_SRX5256555": "Martinez et al. (2020)",
    "SRX5887328_SRX5887329_SRX5887330": "Chen et al. (2020)",
    "SRX11812007_SRX11812008_SRX11812009": "Chothani et al. (2022)"
}

ALIGNERS = ['hisat2', 'STAR', 'tophat2']

COLORS = plt.cm.tab20(np.linspace(0, 1, 10))


# --- 2. 数据收集与统计 ---
def collect_data(input_dir):
    print("开始读取和统计数据...")
    all_stats = []
    
    file_pattern = os.path.join(input_dir, "*_overlap_count.txt")
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"警告：在 {input_dir} 未找到任何 .txt 文件！")
        return pd.DataFrame()

    for file in files:
        filename = os.path.basename(file)
        
        parts = filename.replace("_overlap_count.txt", "").split("_")
        
        aligner = parts[-1]
        dataset_id = "_".join(parts[:-1])
        
        if dataset_id not in DATASET_MAPPING:
            continue  
            
        try:
            df = pd.read_csv(file, sep='\t')
            counts = df['count'].value_counts().to_dict()
            total_orfs = sum(counts.values())
            
            row = {
                'Dataset_ID': dataset_id,
                'Author_Year': DATASET_MAPPING[dataset_id],
                'Aligner': aligner,
                'Total_ORFs': total_orfs
            }
            
            for i in range(1, 12):  
                num = counts.get(i, 0)
                category_name = "Unique" if i == 1 else f"{i}_tools"
                row[f"{category_name}_count"] = num
                row[f"{category_name}_percent"] = (num / total_orfs) * 100 if total_orfs > 0 else 0
                
            all_stats.append(row)
            
        except Exception as e:
            print(f"读取文件出错 {filename}: {e}")

    df_stats = pd.DataFrame(all_stats)
    return df_stats


# --- 3. 绘图函数 ---
def plot_nested_donuts(df_stats, output_pdf):
    print("开始绘制环形饼图...")
    
    datasets = df_stats['Author_Year'].unique()
    num_datasets = len(datasets)

    fig, axes = plt.subplots(2, 3, figsize=(20, 14))
    axes = axes.flatten()
    
    for idx, dataset_name in enumerate(datasets):
        ax = axes[idx]
        df_sub = df_stats[df_stats['Author_Year'] == dataset_name]
        
        if df_sub.empty:
            ax.axis('off')
            continue
            
        ax.set_title(dataset_name, fontsize=16, fontweight='bold', pad=20)
        
        size = 0.3
        radius_start = 1.0
        
        center_texts = []
        
        for ring_idx, aligner in enumerate(ALIGNERS):
            aligner_data = df_sub[df_sub['Aligner'] == aligner]
            if aligner_data.empty:
                continue
            
            row = aligner_data.iloc[0]
            
            counts_list = []
            labels_list = []
            
            for i in range(2, 12):
                c = row.get(f"{i}_tools_count", 0)
                counts_list.append(c)
                labels_list.append(f"{i}")
            
            total_shared = sum(counts_list)
            unique_pct = row.get('Unique_percent', 0)
            
            center_texts.append(f"{aligner}: Unique {unique_pct:.1f}%")
            
            if total_shared > 0:
                wedges, texts, autotexts = ax.pie(
                    counts_list, 
                    radius=radius_start - (ring_idx * size),
                    colors=COLORS,
                    wedgeprops=dict(width=size, edgecolor='w', linewidth=1),
                    autopct=lambda pct: f"{pct:.1f}%" if pct > 3 else "", 
                    pctdistance=1.0 - (size/2) - (ring_idx * size * 0.3),
                    startangle=90
                )
                
                plt.setp(autotexts, size=8, weight="bold", color="black")

        center_info = "Shared ORFs Breakdown\n(Inner to Outer:\n" + ", ".join(ALIGNERS) + ")\n\n"
        center_info += "\n".join(center_texts)
        ax.text(0, 0, center_info, ha='center', va='center', fontsize=10, fontweight='bold', 
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))
        
        ax.set(aspect="equal")
    
    for idx in range(num_datasets, len(axes)):
        axes[idx].axis('off')

    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS[i], 
                                  markersize=10, label=f"{i+2} tools") for i in range(10)]
    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.05), 
               ncol=10, fontsize=12, title="Number of Concordant Tools", title_fontsize=14)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)
    
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"饼状图已保存至 {output_pdf}")


# --- 4. 主函数 ---
def main():
    # 1. 提取并统计数据
    df_stats = collect_data(INPUT_DIR)
    
    if df_stats.empty:
        return

    # 2. 保存 CSV
    df_stats.to_csv(OUTPUT_CSV, index=False)
    print(f"统计数据已保存至 {OUTPUT_CSV}")
    
    # 3. 绘制环形图并保存为 PDF
    plot_nested_donuts(df_stats, OUTPUT_PDF)
    print("任务完成！")

if __name__ == "__main__":
    main()
