'''
Created on 25 oct. 2020

@author: alba
'''

import os.path
import sys
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

# ref_file= "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_protein.faa"
# ref_handle = open(ref_file)
# #parse the fasta file with sequences info
# seq_dict= SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))
# ref_handle.close()
# 
# gff_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.gff"
# file_handle = open(gff_file)
# for rec in GFF.parse(file_handle, base_dict=seq_dict):
#     print(rec)
# file_handle.close()

#will needed gff and fasta files

def gff_parser (gff_file, ref_file):
    gff_abs_path = os.path.abspath(gff_file) ## obtiene el path absolute al fichero. te ahorra problemas de path incompletos futuros
    gff_name, gff_extension = os.path.splitext(gff_abs_path)
    fasta_abs_path = os.path.abspath(ref_file) 
    fasta_name, fasta_extension = os.path.splitext(fasta_abs_path)
    if gff_extension == ".gff" and fasta_extension == ".fna":
        with open(ref_file) as in_handle:
            seq_dict= SeqIO.to_dict(SeqIO.parse(ref_file, "fasta"))

        for rec in GFF.parse(gff_file):
            print(rec.features)
    else: 
        print("The file extension is {}". format(gff_extension), "and {}". format(fasta_extension))
             
def main():
    if len(sys.argv)>2:
        print("")
    else:
        print ("python %s annotation_file ref_fasta\n" %os.path.realpath(__file__))
        exit()
    gff_parser(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
    main()       
        
    
    
    
    
    
    
    
