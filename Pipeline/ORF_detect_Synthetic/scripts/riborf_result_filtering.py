import pandas as pd
#def filter_result_longest_transcript(prediction_table,riborf_07):
def result_filtering(prediction_table_file,riborf_07):
    prediction_table = pd.read_csv(prediction_table_file, header=0, sep='\t')
    prediction_table.drop_duplicates(subset=None, keep='first', inplace=True)
    riborf_result = prediction_table.loc[prediction_table['pred.pvalue'] > 0.7]
    riborf_result_tx1 = riborf_result.orfID.str.split("|", expand=True)[0].str.split(":", expand=True)
    riborf_result_tx1.columns = ['transcript_id','chrom','strand']
    riborf_result_tx1 = riborf_result_tx1.drop(['chrom','strand'], axis=1)
    riborf_result_tx2 = riborf_result.orfID.str.split("|", expand=True)[2].str.split(":", expand=True)
    riborf_result_tx2.columns = ['transcript_length','ORF_tstart','ORF_tstop']
    riborf_result_tx2 = riborf_result_tx2.astype({'ORF_tstop': 'int32'})
    riborf_result_tx2['ORF_tstop'] = riborf_result_tx2['ORF_tstop'] - 1
    riborf_result_total = pd.concat( [riborf_result_tx1, riborf_result_tx2, riborf_result], axis=1)
    ###keep same ORF with logest_transcript
    riborf_result_total.sort_values(by=['transcript_length'],ascending=False)
    #riborf_result_total.drop_duplicates(['chrom','strand','codon5','codon3','length','readNum','f1','f2','f3','entropy','MAXentropy', 'PME','codonNum','f1max','pred.pvalue'],  keep='first', inplace=True)
    riborf_result_total.to_csv(riborf_07, sep='\t', index=False)

result_filtering(snakemake.input[0], snakemake.output[0])
