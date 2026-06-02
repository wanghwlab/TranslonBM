import pandas as pd
import re
def fdr(prediction_table_file,tsv_01,tsv_005):
    prediction_table = pd.read_csv(prediction_table_file, header=0, sep='\t')
    from scipy.stats import rankdata
    p_vals = prediction_table['p value']
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    prediction_table['fdr'] = fdr
    prediction_table = prediction_table.loc[(prediction_table['Codon'] == 'ATG') | (prediction_table['Codon'] == 'CTG') | (prediction_table['Codon'] == 'GTG') | (prediction_table['Codon'] == 'TTG')]
    prediction_table['transcript_id'] = prediction_table['Id'].str.split('_',expand=True)[0]
    prediction_table['chrom'], prediction_table['strand'], prediction_table['ORF_gstart'], prediction_table['ORF_gstop'] = zip(
        *prediction_table['Location'].apply(lambda x: dissect_location(x)))
    fdr01 = prediction_table.loc[prediction_table['fdr'] <= 0.1]
    fdr005 = prediction_table.loc[prediction_table['fdr'] <= 0.05]                   
    fdr01.to_csv(tsv_01, sep='\t', index=False)
    fdr005.to_csv(tsv_005, sep='\t', index=False)

def dissect_location(Location):
    strand=(Location.split(':')[0])[-1]
    if SPE == "Human" or SPE == "Mouse":
        chrom='chr' + (Location.split(':')[0]).strip("+-")
    if SPE == "Zebrafish":
        chrom=(Location.split(':')[0]).strip("+-")
    block_start=(re.findall(r'(\d+)', (Location.split(':')[1])))[0::2]
    block_start = list(map(int, block_start))
    block_end=(re.findall(r'(\d+)', (Location.split(':')[1])))[1::2]
    block_end = list(map(int, block_end))
    ORF_gstart=min(block_start)
    ORF_gstop=max(block_end)
    return chrom, strand, ORF_gstart, ORF_gstop