import pandas as pd

def merge_type(gpd, result, output):
    gpd = pd.read_csv(gpd, sep='\t', header=None, names=[
            'transcript_id', 'chrom', 'strand', 'transcript_gstart',
            'transcript_gstop', 'annotated_gstart',
            'annotated_gstop', 'exon_number',
            'exon_gstart_blocks', 'exon_gstop_blocks', 'score',
            'gene_id', 'cdsStartStat', 'cdsEndStat', 'exonFrames'
        ])
    candidator = pd.read_csv(result, sep='\t', header=0)
    mergeresult = pd.merge(left=candidator, right=gpd, how='inner', on=['transcript_id'])
    mergeresult.to_csv(output, sep='\t', index=False)

merge_type(snakemake.input[1], snakemake.input[0], snakemake.output[0])