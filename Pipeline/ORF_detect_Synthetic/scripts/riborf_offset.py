import numpy as np
import os
def offset_riborf(riborf):
    offsets_list = [READ_LENGTHS, PSITE_OFFSETS]
    offsets_matrix = np.mat(offsets_list).T 
    np.savetxt(riborf, offsets_matrix, delimiter="\t", fmt='%s')