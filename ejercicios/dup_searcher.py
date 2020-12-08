'''
Created on 3 dic. 2020

@author: alba
'''
import os

import pandas as pd
import numpy as np

import argparse
from argparse import ArgumentParser

from Bio import SearchIO
from Bio.Blast import NCBIXML
from Bio.Blast.Record import Blast
#def filter_data(arg_dict):

#https://biopython-tutorial.readthedocs.io/en/latest/notebooks/07%20-%20Blast.html#Parsing-BLAST-output

## -> evalue col 11
## bitscore -> col 12

  
#test_file = "/home/alba/git/TFM_UOC_AMoya/test_BLAST_raw_results.txt"  
######
#test_file:
#Total:16prots
#1=2=3; 4=5=6; 7=8=9
#10=11=12; 13=14=15; 16

######
def filter_data(arg_dict):
    columns = ["qseqid", "sseqid", "pident", "length",
           "mismatch", "gapopen", "qstart", "qend",
           "sstart", "send", "evalue", "bitscore",
           "qlen", "slen"]
#################################################################
# qseqid -> query (e.g., unknown gene) sequence id
# sseqid -> subject (e.g., reference genome) sequence id
# pident -> percentage of identical matches
# length -> alignment length (sequence overlap)
# mismatch -> number of mismatches
# gapopen -> number of gap openings
# qstart -> start of alignment in query
# qend -> end of alignment in query
# sstart -> start of alignment in subject
# send -> end of alignment in subject
# evalue -> expect value
#            number of expected hits of similar quality (score) that could be found just by chance.
#            Blast results are sorted by E-value by default (best hit in first line).
# bitscore -> bit score
#            requires size of a sequence database in which the current match could be found just by chance.
#            The higher the bit-score, the better the sequence similarity
# qlen -> Query sequence length
# slen -> Subject sequence length
################################################################
    file_name_abs_path = os.path.abspath(arg_dict["text_file"])
    name_file, extension = os.path.splitext(file_name_abs_path)
    if extension == ".txt":
        #dataframe with raw blast results
        raw_blast = pd.read_csv(file_name_abs_path, sep="\t",
                                header = None, names=columns)
        print(raw_blast)
        raw_csv = "%s_csv.csv" % name_file
        raw_blast.to_csv(raw_csv, header=True, index=False)
    else:
        print("#####")
        print("Please provide a .txt file")
        print("#####")
        
    #when Blast_file_results is done, evalue>10 is gone
    filter_evalue = raw_blast["evalue"] <= arg_dict["evalue"]
    filter_bitscore = raw_blast["bitscore"] >= arg_dict["evalue"]
    filter_pident = raw_blast["pident"] >= arg_dict["percentage"]
    filter_id = raw_blast["qseqid"] != raw_blast["sseqid"]
    filtered_results = raw_blast[filter_evalue & filter_bitscore & filter_pident & filter_id]
    #sort by pident (desc), evalue(asc), bitscore(desc)
    by_pident = filtered_results.sort_values(["pident", "evalue", "bitscore"],
                                       ascending=[False, True, False])
    print(by_pident)
    sort_csv = "filtered_results.csv"
    by_pident.to_csv(sort_csv, header=True, index=False)


######
##
######
parser = ArgumentParser(prog='dupSearcher',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Search for genomic duplicated proteins")
parser.add_argument("-t", "--text_file", metavar="", help="Blast raw results text file")
parser.add_argument("-e", "--evalue", type=float, metavar="", default= 1e-05, help="BLAST e-value: number of expected hits of similar quality (score) that could be found just by chance.")
parser.add_argument("-b", "--bitscore", type=float, metavar="", default=50, help="BLAST bit-score: requires size of a sequence database in which the current match could be found just by chance.")
parser.add_argument("-p", "--percentage", type=float, metavar="", default=80, help="Percentage of identical matches.")
parser.add_argument("-o", "--out_folder", metavar= "", help="Results folder")
parser.add_argument("--debug", action="store_true", default=False)   
 
arg = parser.parse_args()
arg_dict = vars(arg)
 
if arg.text_file is None:
    print("#####")
    print("Please provide a BLAST results plain text file")
    print("#####")
    print(parser.print_help())
     
if arg.evalue is None:
    print("#####")
    print("Note e-value = 1e-05 is set by default")
    print("#####")
    print(parser.print_help())
     
if arg.bitscore is None:
    print("#####")
    print("Note bit_score = 50 is set by default")
    print("#####")
    print(parser.print_help())
    
if arg.percentage is None:
    print("#####")
    print("Note match percentage = 80% is set by default")
    print("#####")
    print(parser.print_help())

if arg.debug:
    print("##DEBUG: ##")
    print("arguments dictionary: ")
    print(arg)

filter_data(arg_dict)