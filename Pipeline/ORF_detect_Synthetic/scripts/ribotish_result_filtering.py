import pandas as pd

#def filter_result_longest_transcript(prediction_table,ribotish_07):
def result_filtering(prediction_table_file,ribotish_filtered_005):
    prediction_table = pd.read_csv(prediction_table_file, header=0, sep='\t')
    prediction_table.drop_duplicates(subset=None, keep='first', inplace=True)
    ribotish_result = prediction_table.loc[prediction_table['FrameQvalue'] < 0.05]
    ribotish_result.sort_values(by=['FrameQvalue'],ascending=True, inplace=True)
    ribotish_result_genome = ribotish_result.GenomePos.str.split(":", expand=True)
    ribotish_result_genome.columns = ['chrom','g_coords','strand']
    ribotish_result_total = pd.concat( [ribotish_result_genome, ribotish_result], axis=1)
    ###keep same ORF with logest_transcript
    #ribotish_result_total.drop_duplicates(['chrom','strand','codon5','codon3','length','readNum','f1','f2','f3','entropy','MAXentropy', 'PME','codonNum','f1max','pred.pvalue'],  keep='first', inplace=True)
    ribotish_result_total.rename(columns={'Tid':'transcript_id', 'Start':'ORF_tstart', 'Stop':'ORF_tstop'}, inplace = True)
    ribotish_result_total.to_csv(ribotish_filtered_005, sep='\t', index=False)

result_filtering(snakemake.input[0], snakemake.output[0])