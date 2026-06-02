import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import re
import matplotlib


matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# 1. 定义文件路径
csv1_path = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap_study/1tools_percent/1_tools_proportion_summary_simu_untrim.csv"
csv2_path = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap_study/alltools_percent/overlap_summary_stats_simu_untrim.csv"
output_pdf = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap_study/1tools_percent/ORF_overlap_distribution_simu_untrim.pdf"

# 2. 读取数据
df1 = pd.read_csv(csv1_path)
df2 = pd.read_csv(csv2_path)

df1 = df1.rename(columns={'Dataset': 'Dataset_ID'})
df = pd.merge(df2, df1, on=['Dataset_ID', 'Aligner'], how='inner')

# --- 数据清洗 ---
df['Author_Year_Clean'] = df['Author_Year'].str.replace(r'\(|\)', '', regex=True)

# 3. 计算数值
df['Common_count'] = df['11_tools_count']
df['Overlap_count'] = df['Total_ORFs'] - df['Unique_count'] - df['Common_count']

df['Unique_pct_total'] = (df['Unique_count'] / df['Total_ORFs']) * 100
df['Overlap_pct_total'] = (df['Overlap_count'] / df['Total_ORFs']) * 100
df['Common_pct_total'] = (df['Common_count'] / df['Total_ORFs']) * 100

# 4. 绘图配置
dataset_order = [
    #"Ji et al. 2015", 
    #"Gao et al. 2015", 
    #"Calviello et al. 2016", 
    #"Martinez et al. 2020", 
    #"Chen et al. 2020", 
    #"Chothani et al. 2022",
    "T1,6MB",
    "T3,6MB",
    "T1,60MB",
    "T3,60MB"

]

aligners = ['tophat2', 'hisat2', 'STAR']
colors = ['#89bbe9', '#d2c2d7', '#b2cfa5']

print(f"开始绘制横向柱状图并保存至 {output_pdf} ...")

with PdfPages(output_pdf) as pdf:
    fig, axes = plt.subplots(3, 1, figsize=(8, 12))
    
    for ax, aligner in zip(axes, aligners):
        df_sub = df[df['Aligner'] == aligner].copy()
        
        if df_sub.empty:
            ax.set_title(f'Aligner: {aligner} (No Data)')
            ax.axis('off')
            continue
            
        df_sub['Author_Year_Clean'] = pd.Categorical(df_sub['Author_Year_Clean'], categories=dataset_order, ordered=True)
        df_sub = df_sub.sort_values('Author_Year_Clean', ascending=False) # ascending=False 让 Ji et al. 在最上面
        
        y_labels = df_sub['Author_Year_Clean']
        y = np.arange(len(y_labels))
        height = 0.6 
        
        x_unique = df_sub['Unique_count'].values
        x_overlap = df_sub['Overlap_count'].values
        x_common = df_sub['Common_count'].values
        
        ax.barh(y, x_unique, height, label='Unique (1 tool)', color=colors[0], edgecolor='white')
        ax.barh(y, x_overlap, height, left=x_unique, label='Overlap (2-10 tools)', color=colors[1], edgecolor='white')
        ax.barh(y, x_common, height, left=x_unique + x_overlap, label='Common (11 tools)', color=colors[2], edgecolor='white')
        
        # 5. 添加标注
        for i in range(len(y)):
            ribo_pct_val = df_sub['ribotricer_percent(%)'].iloc[i]
            x_ribo_pos = x_unique[i] * (ribo_pct_val / 100)
            
            ax.vlines(x=x_ribo_pos, ymin=y[i]-height/2, ymax=y[i]+height/2, 
                      colors='black', linestyles='--', linewidth=1.2, alpha=0.8)
            
            # --- 文本标注 ---
            # 1. Unique 标注 (位于第一个柱子中间)
            u_pct = df_sub['Unique_pct_total'].iloc[i]
            ax.text(x_unique[i] / 2, y[i], f"{u_pct:.2f}%\n(RiboTricer:{ribo_pct_val:.2f}%)", 
                    ha='center', va='center', fontsize=8.5, fontweight='bold')
            
            # 2. Overlap 标注 (位于第二个柱子中间)
            o_pct = df_sub['Overlap_pct_total'].iloc[i]
            if o_pct > 0.1:
                ax.text(x_unique[i] + x_overlap[i]/2, y[i], f"{o_pct:.2f}%", 
                        ha='center', va='center', fontsize=9)
            
            # 3. Common 标注 (位于第三个柱子中间)
            c_pct = df_sub['Common_pct_total'].iloc[i]
            x_common_pos = x_unique[i] + x_overlap[i] + x_common[i]/2
            h_align = 'left' if c_pct < 2.0 else 'center'
            ax.text(x_common_pos, y[i], f"{c_pct:.2f}%", 
                    ha=h_align, va='center', fontsize=9, color='black')

        ax.set_title(f'Aligner: {aligner}', fontsize=16, fontweight='bold', pad=15)
        ax.set_xlabel('Total Predicted ORFs', fontsize=12)
        ax.set_yticks(y)
        ax.set_yticklabels(y_labels, fontsize=11)
        ax.legend(loc='lower right', framealpha=0.5)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

print("横向柱状图绘制完毕！请检查 PDF 文件。")
