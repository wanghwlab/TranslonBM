import pandas as pd
import numpy as np
import os
import glob
import re
from pathlib import Path

# --- 配置区 ---
BASE_PROCESS_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_group/PR_recall_Maxquant/"
OUTPUT_DIR = "/home/tangyuewen/ORF_benchmark/final_ORFs_2026.1/plots_group/PR_recall_Maxquant/"


def main():
    print("--- 开始聚合Maxquant分析结果 (拆分 Annotated_CDS 和 Novel) ---")
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    search_pattern = os.path.join(BASE_PROCESS_DIR, '**', '*_MS_valid_*_*_score.tsv')
    file_list = glob.glob(search_pattern, recursive=True)
    
    if not file_list:
        print(f"[警告] 在目录 '{BASE_PROCESS_DIR}' 中没有找到任何文件。")
        return
        
    print(f"找到 {len(file_list)} 个结果文件需要处理。")
    
    filename_pattern = re.compile(
        r'^(?P<sample>.+)_MS_valid_(?P<level>protein|peptide)_(?P<aligner>[a-zA-Z0-9]+)(?P<fdr_part>.*)_score\.tsv$'
    )
    
    all_data_list = []
    for filepath in file_list:
        try:
            p = Path(filepath)
            base_p = Path(BASE_PROCESS_DIR)
            relative_p = p.relative_to(base_p)
            
            if len(relative_p.parts) < 3:
                continue
                
            trim_flag = relative_p.parts[0]
            pred_flag = relative_p.parts[1]
            
            filename = os.path.basename(filepath)
            match = filename_pattern.match(filename)
            if not match:
                continue
            
            parsed_info = match.groupdict()
            
            fdr_raw = parsed_info['fdr_part'].strip('_')
            if fdr_raw.startswith('FDR_'):
                fdr_raw = fdr_raw[4:]
            fdr_val = fdr_raw if fdr_raw else "NA"
            
            df = pd.read_csv(filepath, sep='\t')
            
            df['sample'] = parsed_info['sample']
            df['trim_flag'] = trim_flag
            df['pred_flag'] = pred_flag
            df['fdr'] = fdr_val
            df['level'] = parsed_info['level']
            
            all_data_list.append(df)
            
        except Exception as e:
            print(f"[错误] 处理文件失败: {filepath}. 错误信息: {e}")
            
    if not all_data_list:
        print("没有成功读取任何数据，程序退出。")
        return
        
    master_df = pd.concat(all_data_list, ignore_index=True)
    
    df_novel = master_df.copy()
    df_anno = master_df.copy()

    # --- 1. 处理 Novel 组 ---
    df_novel['Expression_Group'] = 'novel'
    
    df_novel['valid_ORFs_num'] = df_novel['valid_novel_ORFs_num']
    df_novel['MS_valid_ORFs_num'] = df_novel['MS_novel_ORFs_num']        
    df_novel['pred_ORFs_num'] = df_novel['pred_novel_ORFs_num'] 
    
    df_novel['Validation_Rate'] = df_novel['Novel_Validation_Rate']
    df_novel['Recall'] = df_novel['Novel_Recall']
    
    df_novel['F1_score'] = np.where(
        (df_novel['Validation_Rate'] + df_novel['Recall']) > 0,
        2 * df_novel['Validation_Rate'] * df_novel['Recall'] / (df_novel['Validation_Rate'] + df_novel['Recall']), 
        0
    )

    # --- 2. 处理 Annotated CDS 组 ---
    df_anno['Expression_Group'] = 'annotated_CDS'
    
    df_anno['valid_ORFs_num'] = df_anno['valid_ORFs_num'] - df_anno['valid_novel_ORFs_num']
    df_anno['MS_valid_ORFs_num'] = df_anno['MS_valid_ORFs_num'] - df_anno['MS_novel_ORFs_num']
    df_anno['pred_ORFs_num'] = df_anno['pred_ORFs_num'] - df_anno['pred_novel_ORFs_num']
    
    df_anno['Validation_Rate'] = np.where(df_anno['pred_ORFs_num'] > 0, df_anno['valid_ORFs_num'] / df_anno['pred_ORFs_num'], 0)
    df_anno['Recall'] = np.where(df_anno['MS_valid_ORFs_num'] > 0, df_anno['valid_ORFs_num'] / df_anno['MS_valid_ORFs_num'], 0)
    
    df_anno['F1_score'] = np.where(
        (df_anno['Validation_Rate'] + df_anno['Recall']) > 0,
        2 * df_anno['Validation_Rate'] * df_anno['Recall'] / (df_anno['Validation_Rate'] + df_anno['Recall']), 
        0
    )


    final_master_df = pd.concat([df_anno, df_novel], ignore_index=True)
    
    rename_map = {
        'aligner': 'soft', 
        'orf_detector': 'tools'
    }
    final_master_df.rename(columns=rename_map, inplace=True)
    
    final_columns_order = [
        'sample', 'soft', 'tools', 'Expression_Group', 
        'Validation_Rate', 'Recall', 'F1_score', 
        'MS_valid_ORFs_num', 'pred_ORFs_num', 'valid_ORFs_num'
    ]
    
    print("\n--- 开始按分组保存最终结果 ---")
    grouped = final_master_df.groupby(['pred_flag', 'trim_flag', 'fdr', 'level'])
    
    for (pred, trim, fdr, level), group_df in grouped:
        output_filename = f"{pred}_{trim}_fdr{fdr}_score_{level}.csv"
        output_filepath = os.path.join(OUTPUT_DIR, output_filename)
        
        final_df = group_df[final_columns_order].copy()

        final_df = final_df.sort_values(by=['sample', 'tools', 'soft', 'Expression_Group'])
        
        final_df.to_csv(output_filepath, index=False)
        print(f"成功保存文件: {output_filepath} (包含 {len(final_df)} 行)")
        
    print("\n--- 所有任务处理完毕 ---")

if __name__ == '__main__':
    main()
