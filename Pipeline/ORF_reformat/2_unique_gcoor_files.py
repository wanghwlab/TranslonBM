import gzip
from pathlib import Path
from multiprocessing import Pool, cpu_count
import sys
from functools import partial
from typing import List, Tuple

# --- 1. 配置---

TARGET_DIRECTORIES = [
    #"/home/tangyuewen/ORF_benchmark/final_ORFs/trim/orf_pred_default/",
    #"/home/tangyuewen/ORF_benchmark/final_ORFs/trim/orf_pred_undisputed/",
    #"/home/tangyuewen/ORF_benchmark/final_ORFs/untrim/orf_pred_default/",
    #"/home/tangyuewen/ORF_benchmark/final_ORFs/untrim/orf_pred_undisputed/",
    "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/trim/final_ORFs/orf_pred_default/",
    "/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/untrim/final_ORFs/orf_pred_default/",
]

BASE_OUTPUT_DIR = Path("/home/tangyuewen/ORF_benchmark/rerun_2025.9/final_ORFs/merged_unique/")

NUM_PROCESSES = min(24, 32)

# --- 2. 辅助函数 ---
def clean_coordinate_string(coord_str: str) -> str:
    """
    过滤掉坐标字符串中长度为0的block (例如 '123-123')。
    """
    if not coord_str or coord_str.lower() == 'na' or '-' not in coord_str:
        return coord_str

    valid_blocks = []
    blocks = coord_str.split(',')
    for block in blocks:
        try:
            start, end = block.split('-')
            if start != end:
                valid_blocks.append(block)
        except ValueError:
            valid_blocks.append(block)
    
    return ','.join(valid_blocks)

# --- 3. 文件处理工作函数 ---

def process_file(input_file: Path, output_dir: Path):
    """
    处理单个 _gcoor.tsv.gz 文件:
    1. 清理 coordinate_id 和 coordinate_0base 列中的零长度block。
    2. 删除 transcript_id 列。
    3. 对处理后的行进行去重。
    4. 将结果写入新的目标目录。
    """
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / input_file.name

        unique_rows = set()
        new_header = ""

        with gzip.open(input_file, 'rt', encoding='utf-8') as infile:
            header_line = infile.readline()
            if not header_line:
                return f"Skipped empty file: {input_file.name}"
            
            header_parts = header_line.strip().split('\t')
            header_parts.pop(1) 
            new_header = '\t'.join(header_parts) + '\n'

            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) < 5: continue 

                parts[2] = clean_coordinate_string(parts[2]) # coordinate_id
                parts[4] = clean_coordinate_string(parts[4]) # coordinate_0base

                parts.pop(1)
                
                unique_rows.add(tuple(parts))

        with gzip.open(output_path, 'wt', encoding='utf-8') as outfile:
            outfile.write(new_header)
            for row_tuple in sorted(list(unique_rows)):
                outfile.write('\t'.join(row_tuple) + '\n')

        return f"Successfully processed: {input_file.name}"

    except Exception as e:
        return f"ERROR processing {input_file.name}: {e}"

# --- 4. 主程序 ---

def main():
    """
    主函数：从指定目录列表查找文件，并并行处理。
    """
    print("Starting file formatting script...")
    print(f"Output directory: {BASE_OUTPUT_DIR}")
    print(f"Using {NUM_PROCESSES} CPU cores.")

    tasks = []
    print("\nLocating files and preparing tasks...")
    
    for dir_path_str in TARGET_DIRECTORIES:
        input_dir = Path(dir_path_str.rstrip('/'))
        if not input_dir.is_dir():
            print(f"Warning: Directory not found, skipping: {dir_path_str}")
            continue
        
        path_parts = input_dir.parts
        if len(path_parts) < 2:
            print(f"Warning: Cannot determine output subdir from path, using parent name: {dir_path_str}")
            sub_dir_name = input_dir.name
        else:
            sub_dir_name = f"{path_parts[-1]}_{path_parts[-2]}_{path_parts[-3]}"
        
        output_subdir = BASE_OUTPUT_DIR / sub_dir_name
        
        found_files = list(input_dir.glob('*_gcoor.tsv.gz'))
        print(f"Found {len(found_files)} files in '{input_dir.name}'. Output will be in '{output_subdir.name}'.")
        for file_path in found_files:
            tasks.append((file_path, output_subdir))

    if not tasks:
        print("\nNo files found to process. Exiting.")
        return

    print(f"\nTotal tasks to process: {len(tasks)}")
    print(f"Initializing parallel processing...\n")
    
    with Pool(processes=NUM_PROCESSES) as pool:
        results = pool.starmap(process_file, tasks)

    print("\n--- Processing Summary ---")
    for res in results:
        print(res)
    print("\nAll tasks completed.")


if __name__ == "__main__":
    main()
