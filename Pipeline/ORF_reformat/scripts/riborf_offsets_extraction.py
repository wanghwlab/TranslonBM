def offsets_extraction_Asite(riborf_detected_Asite,offset_correct_parameters):
    with open(riborf_detected_Asite, 'r') as file_in, open(offset_correct_parameters, 'w') as file_out:
        for line in file_in:
            columns = line.split()
            col2 = columns[0] 
            last_col = columns[-1]
            file_out.write(f"{col2}\t{last_col}\n")

def offsets_extraction_standard(riborf_detected_Asite,riborf_stardard_offsets):
    with open(riborf_detected_Asite, 'r') as file_in, open(riborf_stardard_offsets,'w') as file_out:
        for line in file_in:
            columns = line.split()
            col2 = columns[0] 
            last_col = str(int(columns[-1]) - 3)
            file_out.write(f"{col2}\t{last_col}\n")   

offsets_extraction_Asite(snakemake.input[0],snakemake.output[0])
offsets_extraction_standard(snakemake.input[0],snakemake.output[1])