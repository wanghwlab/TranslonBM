import pandas as pd
import glob
import os
import re
import csv

# ================= 配置区域 =================
INPUT_DIR = "/home/tangyuewen/ORF_benchmark/rerun_2025.9/ORFdetect_simu_untrim"
BASE_DIR = '/home/tangyuewen/ORF_benchmark/rerun_2025.9/' 
OUTPUT_FILE = os.path.join(BASE_DIR, "plots/offset_plots_simu/merged_offsets_all.csv")

# ================= 工具处理函数 =================
def get_tool_name_from_path(file_path, root_dir):
    """
    从路径中动态提取工具名称
    """
    try:
        rel_path = os.path.relpath(file_path, root_dir)
        tool_name = rel_path.split(os.sep)[0]
        return tool_name
    except ValueError:
        return "Unknown"

def get_generic_offsets():
    """
    [1] 通用格式软件 (Standard Offsets)
    包括: orfquant, price, ribotish, ribotaper 等
    匹配: *_standard_offsets.txt
    """
    print(f"\n[1/6] 正在搜索通用格式软件 (Standard Offsets)...")
    search_pattern = os.path.join(INPUT_DIR, "*", "*chrN", "P_site_determination", "*", "*_standard_offsets.txt")
    files = glob.glob(search_pattern)
    print(f"   - 找到 {len(files)} 个文件")
    
    data_list = []
    for file_path in files:
        try:
            file_name = os.path.basename(file_path)
            sample_name = file_name.replace("_standard_offsets.txt", "")
            
            orf_tool = get_tool_name_from_path(file_path, INPUT_DIR)

            df = pd.read_csv(file_path, sep=r'\s+', header=None, names=['read_length', 'offset'], engine='python')
            df['sample'] = sample_name
            df['ORFtools'] = orf_tool
            data_list.append(df)
        except Exception as e:
            print(f"   ! 错误: {file_path} -> {e}")
            
    return pd.concat(data_list, ignore_index=True) if data_list else pd.DataFrame()


def get_ribotricer_offsets():
    """
    [2] Ribotricer
    逻辑: offset = 12 + lag
    """
    print("\n[2/6] 正在处理 Ribotricer...")
    input_pattern = os.path.join(INPUT_DIR, "ribotricer", "ribotricer_chrN", "orf_pred_default", "*", "*_psite_offsets.txt")
    files = glob.glob(input_pattern)
    print(f"   - 找到 {len(files)} 个文件")
    
    DEFAULT_BASE_OFFSET = 12 
    rows = []
    
    for file_path in files:
        try:
            filename = os.path.basename(file_path)
            sample_name = filename.replace("_psite_offsets.txt", "")
            
            with open(file_path, 'r') as f:
                lines = f.readlines()

            for line in lines:
                line = line.strip()
                # Case A: Base offset
                base_match = re.search(r"relative lag to base:\s*(\d+)", line)
                if base_match:
                    rows.append({
                        'sample': sample_name, 'ORFtools': 'ribotricer',
                        'read_length': int(base_match.group(1)), 'offset': DEFAULT_BASE_OFFSET
                    })
                    continue
                # Case B: Lag
                lag_match = re.search(r"lag of\s*(\d+):\s*(-?\d+)", line)
                if lag_match:
                    read_len = int(lag_match.group(1))
                    lag = int(lag_match.group(2))
                    rows.append({
                        'sample': sample_name, 'ORFtools': 'ribotricer',
                        'read_length': read_len, 'offset': DEFAULT_BASE_OFFSET + lag
                    })
        except Exception as e:
            print(f"   ! 错误: {filename} -> {e}")
            
    return pd.DataFrame(rows)


def get_ribocode_offsets():
    """
    [3] Ribocode
    解析: _pre_config.txt
    """
    print("\n[3/6] 正在处理 Ribocode...")
    input_pattern = os.path.join(INPUT_DIR, "ribocode", "ribocode_chrN", "P_site_determination", "SR*", "SRX*_pre_config.txt")
    files = glob.glob(input_pattern)
    print(f"   - 找到 {len(files)} 个文件")
    
    rows = []
    for file_path in files:
        try:
            with open(file_path, 'r') as infile:
                lines = infile.readlines()
            
            config_line = None
            for line in lines:
                if line.strip() and not line.strip().startswith('#'):
                    config_line = line.strip()
                    break
            
            if config_line:
                parts = config_line.split()
                if len(parts) >= 5:
                    sample_name = parts[0]
                    lengths = parts[3].split(',')
                    offsets = parts[4].split(',')
                    
                    if len(lengths) == len(offsets):
                        for l, o in zip(lengths, offsets):
                            rows.append({
                                'sample': sample_name, 'ORFtools': 'ribocode',
                                'read_length': int(l), 'offset': int(o)
                            })
        except Exception as e:
            print(f"   ! 错误: {file_path} -> {e}")
            
    return pd.DataFrame(rows)


def get_ribowave_offsets():
    """
    [4] Ribowave
    逻辑: offset = Position - 1
    """
    print("\n[4/6] 正在处理 Ribowave...")
    input_pattern = os.path.join(INPUT_DIR, "ribowave", "ribowave_chrN", "orf_pred_default", "*", "P-site", "*.psite1nt.txt")
    files = glob.glob(input_pattern)
    print(f"   - 找到 {len(files)} 个文件")
    
    rows = []
    for file_path in files:
        try:
            filename = os.path.basename(file_path)
            sample_name = filename.replace(".psite1nt.txt", "")
            
            with open(file_path, 'r') as infile:
                for line in infile:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            offset = int(parts[1]) - 1
                            rows.append({
                                'sample': sample_name, 'ORFtools': 'ribowave',
                                'read_length': int(parts[0]), 'offset': offset
                            })
                        except ValueError:
                            pass
        except Exception as e:
            print(f"   ! 错误: {file_path} -> {e}")
            
    return pd.DataFrame(rows)


def get_rpbp_offsets():
    """
    [5] Rpbp
    解析: 文件名中的 length-X.offset-Y
    """
    print("\n[5/6] 正在处理 Rpbp...")
    search_path = os.path.join(INPUT_DIR, "rpbp", "rpbp_chrN", "orf_pred_default", "S*", "orf-profiles", "*.profiles.mtx.gz")
    files = glob.glob(search_path)
    print(f"   - 找到 {len(files)} 个文件")
    
    rows = []
    pattern = re.compile(r"^(.*?)\.length-([\d-]+)\.offset-([\d-]+)\.profiles\.mtx\.gz$")
    
    for f in files:
        filename = os.path.basename(f)
        match = pattern.match(filename)
        if match:
            sample_name = match.group(1)
            lengths = [x for x in match.group(2).split('-') if x]
            offsets = [x for x in match.group(3).split('-') if x]
            
            if len(lengths) == len(offsets):
                for l, o in zip(lengths, offsets):
                    rows.append({
                        'sample': sample_name, 'ORFtools': 'rpbp',
                        'read_length': int(l), 'offset': int(o)
                    })
                    
    return pd.DataFrame(rows)


def get_ribohmm_offsets():
    """
    [6] Ribohmm (固定值生成)
    逻辑: 扫描目录获取 Sample -> 生成固定 (28,29,30,31 -> 12)
    """
    print("\n[6/6] 正在处理 Ribohmm (固定值生成)...")
    
    sample_search_pattern = os.path.join(INPUT_DIR, "ribohmm", "ribohmm_chrN", "ORF_detecting_default", "*")
    sample_dirs = glob.glob(sample_search_pattern)
    sample_names = [os.path.basename(d) for d in sample_dirs if os.path.isdir(d)]
    
    print(f"   - 找到 {len(sample_names)} 个 Ribohmm 样本: {sample_names[:3]} ...")
    
    rows = []
    fixed_lengths = [28, 29, 30, 31]
    fixed_offset = 12
    
    for sample in sample_names:
        for length in fixed_lengths:
            rows.append({
                'sample': sample,
                'ORFtools': 'ribohmm',
                'read_length': length,
                'offset': fixed_offset
            })
            
    return pd.DataFrame(rows)


# ================= 主程序执行 =================
def main():
    if not os.path.exists(INPUT_DIR):
        print(f"严重错误: 输入目录不存在 -> {INPUT_DIR}")
        return

    all_dfs = []
    
    all_dfs.append(get_generic_offsets()) 
    all_dfs.append(get_ribotricer_offsets())
    all_dfs.append(get_ribocode_offsets())
    all_dfs.append(get_ribowave_offsets())
    all_dfs.append(get_rpbp_offsets())
    all_dfs.append(get_ribohmm_offsets())
    
    print("\n正在合并所有数据...")
    if any(not df.empty for df in all_dfs):
        final_df = pd.concat(all_dfs, ignore_index=True)
        
        final_df = final_df.sort_values(by=['ORFtools', 'sample', 'read_length'])
        
        final_df = final_df.drop_duplicates()

        cols = ['sample', 'ORFtools', 'read_length', 'offset']
        final_df = final_df[cols]
        
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        final_df.to_csv(OUTPUT_FILE, index=False)
        
        print("=" * 40)
        print(f"全部完成！")
        print(f"数据来源: {INPUT_DIR}")
        print(f"总行数: {len(final_df)}")
        print(f"包含工具: {final_df['ORFtools'].unique()}")
        print(f"结果已保存至: {os.path.abspath(OUTPUT_FILE)}")
        print("=" * 40)
        
        print("\n数据预览 (Top 5):")
        print(final_df.head().to_string(index=False))
    else:
        print("错误: 未提取到任何数据，请检查路径是否正确。")

if __name__ == "__main__":
    main()
