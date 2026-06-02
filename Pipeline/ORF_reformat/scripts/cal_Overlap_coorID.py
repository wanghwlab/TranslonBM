import gzip
import csv

# --- 設定 (Configuration) ---
# 檔案 A 的路徑
#FILE_A = "/home/tangyuewen/ORF_benchmark/final_ORFs/trim/orf_pred_default/SRX876063_SRX876069_hisat2_ribohmm_gcoor.tsv.gz"
FILE_B = "/home/tangyuewen/ORF_benchmark/rerun_2025.5/ORFblock/final_ORFs/orf_pred_default/SRX876063_SRX876069_hisat2_ORFquant_gcoor.tsv.gz"
# 檔案 B 的路徑
FILE_A = "/home/tangyuewen/ORF_benchmark/rerun_2025.5/ORFblock/final_ORFs/orf_pred_default/SRX876063_SRX876069_hisat2_ribohmm_gcoor.tsv.gz"

# 要提取的欄位索引 (第三欄的索引是 2)
# The column index to extract (column 3 has an index of 2)
COLUMN_INDEX = 2 
#0是基因名，1是转录本名，2是gcoor_id

# --- 函數定義 (Function Definition) ---

def get_unique_ids_from_file(filepath, col_index):
    """
    從一個 gzipped TSV 檔案中讀取指定的欄位，並返回一個包含所有不重複 ID 的集合 (set)。
    Reads a specified column from a gzipped TSV file and returns a set of unique IDs.
    """
    unique_ids = set()
    try:
        # 使用 'rt' 模式以文字模式讀取
        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            # 使用 csv.reader 處理 TSV 格式，更穩健
            reader = csv.reader(f, delimiter='\t')
            
            # 跳過標頭行
            # Skip the header line
            header = next(reader, None)
            if header is None:
                print(f"Warning: File is empty or has no header - {filepath}")
                return unique_ids

            # 讀取剩餘的每一行
            for row in reader:
                # 確保行中有足夠的欄位
                if len(row) > col_index:
                    unique_ids.add(row[col_index])
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        # 在檔案找不到時返回一個空集合
        return set()
    except Exception as e:
        print(f"An error occurred while reading {filepath}: {e}")
        return set()
        
    return unique_ids

# --- 主程式 (Main Program) ---

if __name__ == "__main__":
    print("Step 1: Extracting unique IDs from File A...")
    ids_A = get_unique_ids_from_file(FILE_A, COLUMN_INDEX)
    
    print("Step 2: Extracting unique IDs from File B...")
    ids_B = get_unique_ids_from_file(FILE_B, COLUMN_INDEX)
    
    # --- 計算 (Calculation) ---
    
    # 獲取各個檔案的獨立 ID 總數
    total_A = len(ids_A)
    total_B = len(ids_B)
    
    # 使用集合的 intersection() 方法計算交集
    common_ids = ids_A.intersection(ids_B)
    common_count = len(common_ids)
    
    # --- 輸出結果 (Print Results) ---
    
    print("----------------------------------------")
    print("Comparison Results:")
    print("----------------------------------------")
    print(f"Total unique coordinate_ids in File A: {total_A}")
    print(f"Total unique coordinate_ids in File B: {total_B}")
    print(f"Number of common coordinate_ids: {common_count}")
    print("----------------------------------------")
    
    # 計算百分比並處理除以零的錯誤
    if total_A > 0:
        percent_A = (common_count / total_A) * 100
        print(f"Overlap: {percent_A:.2f}% of File A's IDs are present in File B.")
    
    if total_B > 0:
        percent_B = (common_count / total_B) * 100
        print(f"Overlap: {percent_B:.2f}% of File B's IDs are present in File A.")
    
    print("----------------------------------------")
