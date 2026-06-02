from datetime import datetime
import pandas as pd
import re

def time_calc(time_log,bench_time,qc_resource):
    pattern = re.compile(r'([0-9]+):([0-9]+):([0-9]+)')
    with open(time_log, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line.endswith("started ribotricer detect-orfs"):
                start = pattern.search(line).group(0)
            if line.endswith("started calculating phase scores for each ORF"):
                end = pattern.search(line).group(0) 
    start_time=datetime.strptime(str(start),"%H:%M:%S")
    end_time=datetime.strptime(str(end),"%H:%M:%S")
    time=(end_time-start_time).seconds
    resource=pd.read_csv(bench_time, sep="\t", header=0)
    resource['s']=time
    resource.to_csv(qc_resource, index = None, sep = '\t')
    
time_calc(snakemake.input[0],snakemake.input[1],snakemake.output[0])