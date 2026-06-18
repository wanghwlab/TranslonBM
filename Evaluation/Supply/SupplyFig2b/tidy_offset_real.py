import pandas as pd
import glob
import os

def merge_offset_files():
    # 1. 定义搜索路径模式
    search_pattern = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect/*/*chrN/P_site_determination/*/*_standard_offsets.txt"
    
    print(f"正在搜索文件: {search_pattern} ...")
    files = glob.glob(search_pattern)
    
    if not files:
        print("错误：未找到任何匹配的文件，请检查路径是否正确。")
        return

    print(f"共找到 {len(files)} 个文件，开始处理...")

    all_data = []

    for file_path in files:
        try:
            # 2. 提取元数据 (Sample 和 ORFtools)
            file_name = os.path.basename(file_path)
            sample_name = file_name.replace("_standard_offsets.txt", "")
            
            path_parts = file_path.split(os.sep)
            try:
                base_index = path_parts.index("ORFdetect")
                orf_tool = path_parts[base_index + 1] 
            except ValueError:
                print(f"警告: 路径中未找到 'ORFdetect_simu'，无法自动提取工具名: {file_path}")
                orf_tool = "Unknown"

            # 3. 读取文件内容
            df = pd.read_csv(file_path, sep=r'\s+', header=None, names=['read_length', 'offset'], engine='python')
            
            df['sample'] = sample_name
            df['ORFtools'] = orf_tool
            
            all_data.append(df)
            
        except Exception as e:
            print(f"处理文件 {file_path} 时出错: {e}")

    # 4. 合并并保存
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        
        cols = ['sample', 'ORFtools', 'read_length', 'offset']
        final_df = final_df[cols]
        
        output_file = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_real/merged_offsets_all.csv"
        final_df.to_csv(output_file, index=False)
        
        print("-" * 30)
        print(f"处理完成！")
        print(f"总行数: {len(final_df)}")
        print(f"包含 Sample 数量: {final_df['sample'].nunique()}")
        print(f"包含 ORFtools 数量: {final_df['ORFtools'].nunique()}")
        print(f"文件已保存为: {os.path.abspath(output_file)}")
        
        print("\n数据预览:")
        print(final_df.head())
    else:
        print("未提取到任何数据。")

if __name__ == "__main__":
    merge_offset_files()
