'''
Created on 8 nov. 2020

@author: alba
'''

import pandas as pd
import numpy as np
import protein_gbf
from Bio import SeqIO, Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser


# data = np.array([["","col1", "col2"], ["Fila1", 11,22], ["Fila2", 33,44]])
# print(pd.DataFrame(data=data[1:,1:], index=data[1:,0], columns=data[0,1:]))
# 
# series = pd.Series()

#gbf_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gbf"
#out_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot-proteins.fa"
# for rec in protein_gbf.gbf_parser(gbf_file):
#     series = pd.Series({rec.id:str(rec.seq)})
#     #index=np.array[str(rec.id)]
#     data= np.array([["","id","sequence"], ["",rec.id, str(rec.seq)]])
#     df = pd.DataFrame(data=data[1:,1:], index=data[1:,0], columns=data[0,1:])
#     print(df)


# for rec in protein_gbf.gbf_parser(gbf_file):
#     datos = {"id":[rec.id], "sequence":[str(rec.seq)]}
#     df=pd.DataFrame(datos, index = [nro)
#     print(df)
# for rec in protein_gbf.gbf_parser(gbf_file):
#     series = pd.Series({"id":rec.id, "sequence":str(rec.seq)})
#     df=pd.DataFrame(series)
#     print(df)

# with open(out_file) as fasta_file:
#         identifiers = []
#         seq = []
#         for title, sequence in SimpleFastaParser(fasta_file):
#                         identifiers.append(title)
#                         seq.append(sequence)
#                         s1 = pd.Series(identifiers, name="ID")
#                         s2 = pd.Series(seq, name="sequence")
#                         gbs_df = pd.DataFrame(dict(ID=s1, sequence=s2)).set_index(["ID"])
#                         print(gbs_df)
                        
with open(out_file) as fasta_file:
        identifiers = []
        seq = []
        for title, sequence in SimpleFastaParser(fasta_file):
                        identifiers.append(title)
                        seq.append(sequence)
                        s1 = pd.Series(identifiers, name="ID")
                        s2 = pd.Series(seq, name="sequence")
                        gbf_df = pd.DataFrame(dict(ID=s1, sequence=s2)).set_index(["ID"])
                        gbf_df.reset_index().to_csv("out_gbf.csv", header=True, index=False)
                                  
                        
                        
                        