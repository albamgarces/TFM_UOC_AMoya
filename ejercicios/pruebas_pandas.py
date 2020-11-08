'''
Created on 8 nov. 2020

@author: alba
'''

import pandas as pd
import numpy as np
import protein_gbf

# data = np.array([["","col1", "col2"], ["Fila1", 11,22], ["Fila2", 33,44]])
# print(pd.DataFrame(data=data[1:,1:], index=data[1:,0], columns=data[0,1:]))
# 
# series = pd.Series()

gbf_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gbf"
for rec in protein_gbf.gbf_parser(gbf_file):
    series = pd.Series({rec.id:str(rec.seq)})
    #index=np.array[str(rec.id)]
    data= np.array([["","sequence"], [rec.id, str(rec.seq)]])
    df = pd.DataFrame(data=data[1:,1:], index=data[1:,0], columns=data[0,1:])
    print(df)