'''
Created on 25 oct. 2020

@author: alba
'''
import os.path
import Bio
from Bio import SeqIO
import sys
import protein_gbf
import protein_gff

#protein_gbf.test()
import argparse
from argparse import ArgumentParser

############
    # 1: identificar annot_file: GFF / Genbank / None
    
    # 2: Si es genbank: gbf_parser
    
    # 3: Si es GFF:
    # 3.1: Comprobar ref_file: Yes/No
    # 3.2: gff_parser(annot_file, ref_file)
    
    # ToDO: create gtf parser
#############

def extension(args_dict):
    compt = {}
    compt["fasta"] = [".fa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]
    compt["genbank"] = [".genbank", ".gb", ".gbf", ".gbff", ".gbk"]
    compt["GFF"] = {".gff"}
    
    if arg.debug:
        print("## DEBUG: Format input file ")
        print(compt)
     
    #get  annot file absolute path
     
    file_name_abs_path = os.path.abspath(arg_dict["annot_file"])
    name_file, extension = os.path.splitext(file_name_abs_path)
    
    if arg.debug:
        print("## Debug: name_file and extension ")
        print(os.path.splitext(file_name_abs_path))
         
    #get protein fasta file and annotation csv file
 
    if extension in compt["genbank"]:
        if not arg_dict["ref_file"] is None:
            print("######")
            print("Sorry, we don't need a ref file to parse a genbank file'")
            print("######")
        else:
            gbf_parser = protein_gbf.gbf_parser(arg_dict["annot_file"])
            for fasta in gbf_parser:
                protein_gbf.main(arg_dict["annot_file"])
    
    elif extension in compt["GFF"]:
        if arg_dict["ref_file"]==None:
            print("######")
            print("Please provide a ref file FASTA format")
            print("######")
        else:
            ref_name_abs_path = os.path.abspath(arg_dict["ref_file"])
            ref_name, ref_extension = os.path.splitext(ref_name_abs_path)
            
            if arg.debug:
                print("## DEBUG: ref_name and extension ")
                print(os.path.splitext(ref_name_abs_path))
            
            if ref_extension in compt["fasta"]:
                    gff_parser = protein_gff.protein_recs(arg_dict["annot_file"], arg_dict["ref_file"])
                    protein_gff.main(arg_dict["annot_file"], arg_dict["ref_file"])
            else:
                print("######")
                print("Compatible ref_file formats: ")
                print(compt["fasta"])
                print("######")
    else:
        print("######")
        print("Compatible file formats: ")
        print(compt)
        print("######")
       

############
####################
## starting
####################
####################
    

parser = ArgumentParser(prog='inputParser',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Get proteins sequences from annotation file")
parser.add_argument("-a", "--annot_file", metavar="", help="Annotation file")
parser.add_argument("-r", "--ref_file", metavar="", help="Genome references FASTA file")
parser.add_argument("-p", "--prot_file", metavar="", help="Protein sequence file")
parser.add_argument("-o", "--out_folder", metavar= "", help="Results folder")
parser.add_argument("--debug", action="store_true", default=False)   

arg = parser.parse_args()
arg_dict = vars(arg)
  


# if not arg.annot_file is None:
#     file_given = arg.annot_file
#     print("file given")
#     print(file_given)
#     print(" ")
# elif not arg.ref_file is None:
#     ref_file = arg.ref_file
#     print("ref file given")
#     print(ref_file)
#     print("")
    
if arg.annot_file is None and arg.out_folder is None:
    print("######")
    print("Please provide either an annotation file or an output path folder")
    print("######")
    print(parser.print_help())
    exit()
    
if arg.debug:
    print("##DEBUG: ##")
    print("arguments dictionary: ")
    print(arg)

    

## Cuando no se importe como un modulo este script, se ejecutara esto
if __name__ == '__main__':
    extension(arg_dict)
         
