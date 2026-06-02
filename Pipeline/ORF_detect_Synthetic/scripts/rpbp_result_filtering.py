import pandas as pd

def result_filtering(prediction_table_file,rpbp_filtered_bayes_factor_mean):
    prediction_table = pd.read_csv(prediction_table_file, header=0, sep='\t')
    prediction_table.drop_duplicates(subset=None, keep='first', inplace=True)
    rpbp_result = prediction_table.loc[prediction_table['#bayes_factor_mean'] > 5]
    rpbp_result.sort_values(by=['#bayes_factor_mean'],ascending=False, inplace=True)
    rpbp_result.rename(columns={'#bayes_factor_mean':'bayes_factor_mean', '#id':'id', '#seqname':'chrom', '#start':'ORF_gstart', '#end':'ORF_gstop', '#strand':'strand'}, inplace = True)
    rpbp_result_genome = rpbp_result.id.str.split("_", expand=True)
    rpbp_result_genome.columns = ['transcript_id', 'ORF_gpos']
    rpbp_result_total = pd.concat([rpbp_result_genome, rpbp_result], axis=1)
    rpbp_result_total.to_csv(rpbp_filtered_bayes_factor_mean, sep='\t', index=False)

result_filtering(snakemake.input[0], snakemake.output[0])
