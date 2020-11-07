'''
Created on 25 oct. 2020

@author: alba
'''
import os.path
import Bio
from Bio import SeqIO
import sys
import protein_gbf

protein_gbf.test()

#extract extension

def extension(file_given):
    compt = {".gbf", ".gb", ".gbk", ".gff", ".faa"}
    #get file absolute path
    file_name_abs_path = os.path.abspath(file_given)
    name_file, extension, path_file = os.path.splitext(file_name_abs_path)
    if extension in compt:
        if extension == ".gbf" or extension == ".gbk" or extension == ".gb":
            protein_gbf(file_given) 
        elif extension == ".gff":
            protein_gff(file_given, ref_file)
    else:
        print("File format not available")
        
def print ():
        
# def protein_gbf (file_given):
#     seq_rec = [rec for rec in SeqIO.parse(file_given, "genbank")]
#     for rec in seq_rec:
#         feats = [feat for feat in rec.features if feat.type == "CDS"]
#         for feat in feats:
#             print ("locus_tag: ", feat.qualifiers["locus_tag"][0])
#             print("product: ", feat.qualifiers["product"])
#             print(feat.qualifiers["translation"])
#             
# def protein_gff(file_given, ref_file):
#     fasta_abs_path = os.path.abspath(ref_file)
#     fasta_name, fasta_extension
    