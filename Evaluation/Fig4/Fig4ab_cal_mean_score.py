import os
import glob
import pandas as pd

# --- 1. 配置区 ---

BASE_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_combine/"


INPUT_SUBDIRS = [
    #"TIS_score_trim",
    #"TIS_score_intersect_trim",
    #"TIS_score_intersect_untrim",
    #"TIS_score_union_trim",
    #"TIS_score_union_untrim",
    #"Maxquant_score_union_trim",
    #"Maxquant_score_union_untrim",
    #"Maxquant_score_intersect_trim",
    #"Maxquant_score_intersect_untrim",
    #"TISeq_score_intersect_trim",
    #"TISeq_score_intersect_untrim",
    #"TISeq_score_union_trim",
    #"TISeq_score_union_untrim",
    "merged_gain"
]

OUTPUT_SUBDIR = "Mean_score"


def process_and_average_files():
    """
    主函数：查找、读取、聚合和保存CSV文件。
    """
    full_output_path = os.path.join(BASE_DIR, OUTPUT_SUBDIR)
    
    os.makedirs(full_output_path, exist_ok=True)
    
    for subdir in INPUT_SUBDIRS:
        full_input_path = os.path.join(BASE_DIR, subdir)
        print("-" * 50)
        print(f"开始处理子文件夹: {subdir}")

        search_pattern = os.path.join(full_input_path, '*combination_gains*_peptide.csv')
        csv_files = glob.glob(search_pattern)

        if not csv_files:
            print(f"  [警告] 在目录 '{full_input_path}' 中找不到任何CSV文件。")
            continue
            
        print(f"  找到 {len(csv_files)} 个CSV文件，开始处理...")
        
        for file_path in csv_files:
            basename = os.path.basename(file_path)
            print(f"    - 正在处理: {basename}")
            
            try:

                df = pd.read_csv(file_path)
                
                grouping_columns = ['aligner', 'tool_a', 'tool_b']
                
                if not all(col in df.columns for col in grouping_columns):
                    print(f"      [警告] 文件 {basename} 缺少必要的分组列，已跳过。")
                    continue

                averaged_df = df.groupby(grouping_columns).mean().reset_index()
                
                averaged_df.insert(1, 'sample', 'Average7sample')
                
                new_filename = f"{subdir}_{basename}"
                output_filepath = os.path.join(full_output_path, new_filename)
                
                averaged_df.to_csv(output_filepath, index=False)
                print(f"      -> 已保存平均化结果至: {new_filename}")

            except Exception as e:
                print(f"      [错误] 处理文件 {basename} 时发生错误: {e}")

    print("\n所有文件处理完毕。")


if __name__ == '__main__':
    process_and_average_files()
