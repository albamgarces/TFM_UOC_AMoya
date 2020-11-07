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

from argparse import ArgumentParser

#extract extension

def extension(file_given):
    compt = {".gbf", ".gb", ".gbk", ".gff", ".faa"}
    #get file absolute path
    file_name_abs_path = os.path.abspath(file_given)
    name_file, extension, path_file = os.path.splitext(file_name_abs_path)
    if extension in compt:
        if extension == ".gbf" or extension == ".gbk" or extension == ".gb":
            gbf_parser= protein_gbf.gbf_parser(file_given, debug)
            protein_gbf.main(gbf_parser, debug)
        
             
        elif extension == ".gff":
            protein_gff(file_given, ref_file)
            
       
    else:
        print("File format not available")
        

    

def parseArgument():
    parser = ArgumentParser(description="Get proteins sequences from annotation file")
    parser.add_argument("-a", "--annot_file", metavar="", help="gb, gff or gtf annotation file")
    parser.add_argument("-r", "--ref_file", metavar="", help="Genome references FASTA file")
    parser.add_argument("-p", "--prot_file", metavar="", help="Protein sequence file")
    parser.add_argument("--debug", metavar="")
    parser.add_argument("-o", "--out_folder", metavar="", help="Results folder")
    return parser.parse_args()

arg= parseArgument()
if arg.annot_file == None:
    print("Please provide an annotation file")

        
