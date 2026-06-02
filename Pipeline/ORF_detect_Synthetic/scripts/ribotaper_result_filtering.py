import pandas as pd

def result_filtering(prediction_table_file,ribotaper_filtered_005):
    prediction_table = pd.read_csv(prediction_table_file, header=0, sep='\t')
    prediction_table.drop_duplicates(subset=None, keep='first', inplace=True)
    ribotaper_result = prediction_table.loc[prediction_table['ORF_pval_multi_ribo'] < 0.05]
    ribotaper_result_genome = ribotaper_result.ORF_id_gen.str.split("_", expand=True)
    ribotaper_result_genome.columns = ['chrom','ORF_gstart','ORF_gstop']
    ribotaper_result_total = pd.concat( [ribotaper_result_genome, ribotaper_result], axis=1)
    ribotaper_result_total.sort_values(by=['ORF_pval_multi_ribo'],ascending=True)
    ribotaper_result_total.rename(columns={'start_pos':'ORF_tstart', 'stop_pos':'ORF_tstop'}, inplace = True)
    ribotaper_result_total['ORF_tstop'] = ribotaper_result_total['ORF_tstop'] + 2
    ribotaper_result_total.to_csv(ribotaper_filtered_005, sep='\t', index=False)

result_filtering(snakemake.input[0], snakemake.output[0])