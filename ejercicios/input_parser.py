'''
Created on 25 oct. 2020

@author: alba
'''
import os.path
import Bio
from Bio import SeqIO
import sys
import protein_gbf

#protein_gbf.test()

from argparse import ArgumentParser

#extract extension

def extension(file_given):
    compt = {}
    compt["fasta"] = [".fa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]
    compt["genbank"] = [".genbank", ".gb", ".gbf", ".gbff", ".gbk"]
    compt["GFF"] = {".gff"}
    
    #get file absolute path
    file_name_abs_path = os.path.abspath(file_given)
    #name_file, extension, path_file = os.path.splitext(file_name_abs_path)
    name_file, extension = os.path.splitext(file_name_abs_path)
    print(os.path.splitext(file_name_abs_path))
    
    
    if len(sys.argv)==3 and extension in compt["genbank"]:
        gbf_parser = protein_gbf.gbf_parser(file_given)
        
        for fasta in gbf_parser:
            protein_gbf.main(file_given)
#     elif extension in compt[""]        
    
      
    else:
        print("")
        print("Sorry, doesn't work, try again!")
        print("")
        
        
        
#         
#         else:
#             print("File format not available")
    
#         if extension == ".gbf" or extension == ".gbk" or extension == ".gb":
#             gbf_parser= protein_gbf.gbf_parser(file_given, debug)
#             protein_gbf.main(gbf_parser, debug)
        
             
#         elif extension == ".gff":
#             protein_gff(file_given, ref_file)
   
def main():
    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        print ("## Usage input_parser: ")
        print ("python %s annotation_file\n" %os.path.realpath(__file__))
        exit()        

    
    extension(sys.argv[2])   

    ## ATENCION: sys.argv[0] es el propio script de python!
    

## Cuando no se importe como un modulo este script, se ejecutara esto
if __name__ == '__main__':
    main()
        

    

def parseArgument():
    parser = ArgumentParser(description="Get proteins sequences from annotation file")
    parser.add_argument("-a", "--annot_file", metavar="", help="Annotation file")
    parser.add_argument("-r", "--ref_file", metavar="", help="Genome references FASTA file")
    parser.add_argument("-p", "--prot_file", metavar="", help="Protein sequence file")
    parser.add_argument("--debug", metavar="")
    parser.add_argument("-o", "--out_folder", metavar="", help="Results folder")
    return parser.parse_args()
  
arg= parseArgument()

####################
## starting
####################


if not arg.annot_file is None:
    file_given = arg.annot_file
    print("file given")
    print(file_given)
    print(" ")
elif not arg.ref_file is None:
    ref_file = arg.ref_file
    print("ref file given")
    print(ref_file)
    print("")
    
# if arg.annot_file == None:
#     print("Please provide an annotation file")
else:
    print("Please provide an annotation file")
    

         
