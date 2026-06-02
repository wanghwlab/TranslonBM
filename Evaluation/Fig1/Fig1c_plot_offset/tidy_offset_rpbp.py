import glob
import os
import csv
import re

# 1. 定义搜索路径
search_path = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_simu_untrim/rpbp/rpbp_chrN/orf_pred_default/*/orf-profiles/*.profiles.mtx.gz"

# 2. 定义输出文件名
output_csv = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/plots/offset_plots_simu/rpbp_converted_offsets.csv"

def parse_rpbp_filename(file_path):
    """
    解析文件名，提取 Sample Name, Lengths, Offsets
    文件名示例: SRX...tophat2.length-25-28-29-33.offset-9-12-12-0.profiles.mtx.gz
    """
    filename = os.path.basename(file_path)
    
    # 使用正则表达式提取关键部分
    pattern = re.compile(r"^(.*?)\.length-([\d-]+)\.offset-([\d-]+)\.profiles\.mtx\.gz$")
    match = pattern.match(filename)
    
    if match:
        sample_name = match.group(1)
        lengths_str = match.group(2)
        offsets_str = match.group(3)
        
        lengths = [x for x in lengths_str.split('-') if x]
        offsets = [x for x in offsets_str.split('-') if x]
        
        return sample_name, lengths, offsets
    else:
        print(f"警告: 无法解析文件名格式 -> {filename}")
        return None, None, None

def main():
    files = glob.glob(search_path)
    print(f"找到 {len(files)} 个文件，开始处理...")
    
    results = []
    
    for f in files:
        sample_name, lengths, offsets = parse_rpbp_filename(f)
        
        if sample_name and lengths and offsets:
            if len(lengths) != len(offsets):
                print(f"错误: {sample_name} 的长度数量与偏移量数量不匹配。")
                continue
            
            for length, offset in zip(lengths, offsets):
                row = [sample_name, 'rpbp', length, offset]
                results.append(row)

    # 写入 CSV 文件
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=',') 
        writer.writerows(results)
        
    print(f"处理完成！结果已保存至: {os.path.abspath(output_csv)}")

if __name__ == "__main__":
    main()
