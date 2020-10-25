'''
Created on 25 oct. 2020

@author: alba
'''


import Bio
from Bio import SeqIO
import BCBio
from BCBio import GFF
# ref_file= "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_proteins.faa"
# ref_handle = open(ref_file)
# #parse the fasta file with sequences info
# seq_dict= SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))
# ref_handle.close()

# gff_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gff"
# file_handle = open(gff_file)
# for rec in GFF.parse(file_handle, base_dict=seq_dict):
#     print(rec)
# gff_file.close()

ref_file= "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_protein.faa"
ref_handle = open(ref_file)
#parse the fasta file with sequences info
seq_dict= SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))
ref_handle.close()

gff_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.gff"
file_handle = open(gff_file)
for rec in GFF.parse(file_handle, base_dict=seq_dict):
    print(rec)
file_handle.close()
