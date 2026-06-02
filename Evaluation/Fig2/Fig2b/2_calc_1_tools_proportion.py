import os
import pandas as pd

# 1. 定义数据集、比对软件和ORF预测软件
datasets = [
    "SRX876063_SRX876069", "SRX740748", "SRX1254413", 
    "SRX5256543_SRX5256555", "SRX5887328_SRX5887329_SRX5887330", 
    "SRX11812007_SRX11812008_SRX11812009"
]
aligners = ["tophat2", "hisat2", "STAR"]
tools = [
    "riborf", "ribocode", "ribotricer", "ribohmm", "ribotish", 
    "ribowave", "orfrater", "orfquant", "gedi", "ribotaper", "rpbp"
]

# 2. 定义文件路径
overlap_dir = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/tools_overlap/merged_ATG/orf_pred_default_untrim/"
pred_dir = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_ATG/orf_pred_default_untrim/"
output_csv = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots/tools_overlap_study/1tools_percent/1_tools_proportion_summary_real_untrim.csv"

results_list = []

print("开始统计 1_tools 的软件来源比例...")

# 3. 遍历所有数据集和比对软件的组合
for ds in datasets:
    for al in aligners:
        overlap_file = os.path.join(overlap_dir, f"{ds}_{al}_overlap_count.txt")
        
        if not os.path.exists(overlap_file):
            print(f"[跳过] 找不到文件: {overlap_file}")
            continue
            
        print(f"正在处理: {ds} - {al}")
        
        df_overlap = pd.read_csv(overlap_file, sep=r'\s+')
        
        if 'tools' not in df_overlap.columns or 'coordinate_id' not in df_overlap.columns:
            print(f"[警告] 文件 {overlap_file} 缺少必要的列，请检查格式。")
            continue
            
        set_1_tools = set(df_overlap[df_overlap['tools'] == '1_tools']['coordinate_id'])
        total_1_tools = len(set_1_tools)
        
        row_data = {
            'Dataset': ds,
            'Aligner': al,
            'Total_1_tools_ORFs': total_1_tools
        }
        
        if total_1_tools == 0:
            print(f"  -> {ds} {al} 中没有 1_tools 的 ORF。")
            for tool in tools:
                row_data[f'{tool}_count'] = 0
                row_data[f'{tool}_percent(%)'] = 0.0
            results_list.append(row_data)
            continue
            
        # 4. 遍历 11 种 ORF 预测软件，寻找交集
        for tool in tools:
            pred_file = os.path.join(pred_dir, f"{ds}_{al}_{tool}_gcoor.tsv.gz")
            
            if os.path.exists(pred_file):
                try:
                    df_pred = pd.read_csv(pred_file, sep='\t', usecols=['coordinate_id'], compression='gzip')
                    tool_coords = set(df_pred['coordinate_id'])
                    
                    overlap_count = len(tool_coords.intersection(set_1_tools))
                    
                except Exception as e:
                    print(f"  [错误] 读取文件 {pred_file} 失败: {e}")
                    overlap_count = 0
            else:
                overlap_count = 0
                
            percentage = (overlap_count / total_1_tools) * 100
            
            row_data[f'{tool}_count'] = overlap_count
            row_data[f'{tool}_percent(%)'] = round(percentage, 2)
            
        results_list.append(row_data)

# 5. 将结果保存为 CSV
df_results = pd.DataFrame(results_list)
df_results.to_csv(output_csv, index=False)

print(f"\n统计完成！结果已成功保存至当前目录下的: {output_csv}")
