'''
Created on 6 nov. 2020

@author: alba
'''

#https://www.youtube.com/watch?v=SbQmQ4T4E8E

import sys

# print(sys.argv)

from argparse import ArgumentParser

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

    
    
    