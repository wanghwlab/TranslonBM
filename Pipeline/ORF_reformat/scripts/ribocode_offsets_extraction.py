## #!/usr/bin/env python
import re 

def offsets_extraction(ribocode_offsets,standard_offsets,prefix):
    matched_lines = []
    with open(ribocode_offsets, "r") as config:
        for line in config:
            line = line.strip()
            if not line.startswith("#") and re.search(prefix, line):
                matched_lines.append(line)
    with open(standard_offsets, "w") as f:
        for line in matched_lines:
            result = line.split('\t')[-2:]
            result_list = [index.split(',') for index in result]
            table = [list(x) for x in zip(*result_list)]
            for row in table:
                f.write('\t'.join(row) + '\n') ## 

# ============================================================================ #
# need consider that if no psite avaliable

offsets_extraction(snakemake.input[0], snakemake.output[0],snakemake.params[0])
