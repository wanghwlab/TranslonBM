import pandas as pd

def result_filtering(prediction_table_file,ribotricer_filtered_translating):
    prediction_table = pd.read_csv(prediction_table_file, header=0, sep='\t')
    prediction_table.drop_duplicates(subset=None, keep='first', inplace=True)
    ribotricer_result = prediction_table.loc[prediction_table['status'] == "translating"]
    ribotricer_result.sort_values(by=['phase_score'],ascending=False, inplace=True)
    ribotricer_result_genome = ribotricer_result.ORF_ID.str.split("_", expand=True)
    ribotricer_result_genome.columns = ['transcript_id', 'ORF_gstart', 'ORF_gstop', 'ORF_length']
    ribotricer_result_genome.drop(['transcript_id'], axis=1, inplace=True)
    ribotricer_result_total = pd.concat( [ribotricer_result_genome, ribotricer_result], axis=1)
    ribotricer_result_total.to_csv(ribotricer_filtered_translating, sep='\t', index=False)

result_filtering(snakemake.input[0], snakemake.output[0])