import pandas as pd
import glob
import os

def extract_rpbp_offsets():
    # 1. 设置搜索路径
    input_pattern = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/rpbp/rpbp_chrN/orf_pred_default/*/metagene-profiles/*.periodic-offsets.csv.gz"
    
    print(f"正在搜索文件: {input_pattern} ...")
    files = glob.glob(input_pattern)
    
    if not files:
        print("错误：未找到任何 .csv.gz 文件，请检查路径。")
        return

    print(f"共找到 {len(files)} 个文件，开始处理...")
    
    all_data = []

    for file_path in files:
        try:
            # --- A. 解析文件名提取 Sample ---
            basename = os.path.basename(file_path)
            
            sample_name = basename.replace(".periodic-offsets.csv.gz", "")
            
            # --- B. 读取 Gzip CSV ---
            df = pd.read_csv(file_path, compression='gzip')
            
            if 'length' not in df.columns or 'highest_peak_offset' not in df.columns:
                print(f"跳过文件 {basename}: 列名不匹配")
                continue

            # --- C. 数据清洗与转换 ---
            temp_df = df[['length', 'highest_peak_offset']].copy()
            
            temp_df['read_length'] = temp_df['length'].astype(int)
            
            temp_df['offset'] = temp_df['highest_peak_offset'].abs().astype(int)
            temp_df['sample'] = sample_name
            temp_df['ORFtools'] = 'rpbp'
            
            final_cols_df = temp_df[['sample', 'ORFtools', 'read_length', 'offset']]
            
            all_data.append(final_cols_df)
            
        except Exception as e:
            print(f"处理文件 {basename} 时出错: {e}")

    # --- D. 合并并保存 ---
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        final_df = final_df[['sample', 'ORFtools', 'read_length', 'offset']]
        final_df = final_df.sort_values(by=['sample', 'read_length'])
        
        output_file = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/rpbp_offsets_extracted.csv"
        final_df.to_csv(output_file, index=False)
        
        print("-" * 30)
        print("处理完成！")
        print(f"包含样本数: {final_df['sample'].nunique()}")
        print(f"文件已保存为: {os.path.abspath(output_file)}")
        print("\n数据预览:")
        print(final_df.head())
    else:
        print("未提取到任何有效数据。")

if __name__ == "__main__":
    extract_rpbp_offsets()
