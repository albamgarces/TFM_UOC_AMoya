'''
Created on 28 nov. 2020

@author: alba
'''

import os
import sys

import subprocess
import argparse
from argparse import ArgumentParser

from Bio.Blast.Applications import NcbiblastpCommandline

#######
#
compt = {}
compt["fasta"] = [".fa", ".faa", ".mpfa", ".fna", ".fsa", ".fas", ".fasta"]

#blastp -query fasta_file -db name -outfmt '6 std qlen slen' -num_threads X -out name_out
# def makeblast_results(arg_dict):
#     #create a commandline fotr the NCBI BLAST+ program blastp
#     #cline = NcbiblastpCommandline(query=arg_dict["fasta_file"], db=(arg_dict["db_name"] + '.phr'), )
#     file_name_abs_path = os.path.abspath(arg_dict["fasta_file"])
#     base, extension = os.path.splitext(file_name_abs_path)
#     if extension in compt["fasta"]:
#         output_path = "%s_blastp_results.txt" % base
#         db_path = "%s.phr" % arg_dict["db_name"]
#     
#         cmd_makeblast_results = "blastp -query %s -db  %s -outfmt \'6 std qlen slen|' -num_threads 1 -out %s" %(arg_dict["fasta_file"], arg_dict["db_name"], output_path)
#         code = system_call(cmd_makeblast_results)
#         if (code == 'FAIL'):
#                 print ('****ERROR: Some error happened during the blastp command')
#                 print (cmd_makeblast_results)
#                 exit()  
#     else:
#             print("#####")
#             print("Please provide a FASTA file")
#             print(compt)
#             print("#####")

#https://github.com/HCGB-IGTP/HCGB_python_functions/blob/ce1762b709e3b9094c0630f63bfb9d33544ed4c3/HCGB/functions/blast_functions.py
def makeblastdb(arg_dict):
    ## generate blastdb for genome
    #phr is the header file, pin is the index file, psq is the sequence file
    file_name_abs_path = os.path.abspath(arg_dict["fasta_file"])
    name_file, extension = os.path.splitext(file_name_abs_path)
    output_path = "%s_blastp_results.txt" % arg_dict["db_name"]
    if arg.debug:
        print("## Debug: name_file and extension ")
        print(os.path.splitext(file_name_abs_path))
            
    if extension in compt["fasta"]:
            
        if arg.debug:
            print("## DEBUG: Format input file ")
            print(compt)
    
        if (os.path.isfile(arg_dict["db_name"] + '.phr')):
            print ("+ BLAST database is already generated...")
                
        cmd_makeblastdb = "%s -in %s -input_type fasta -dbtype %s -out %s" %(arg_dict["executable"], arg_dict["fasta_file"], 'prot', arg_dict["db_name"])
        code = system_call(cmd_makeblastdb)
            ##
            ## -> ¿distinguir que sea un fasta con aminoácidos y no nucleótidos?
            ##
        if (code == 'FAIL'):
            print ('****ERROR: Some error happened during the makeblastDB command')
            print (cmd_makeblastdb)
            exit()
    
        cmd_makeblast_results = "blastp -query %s -db  %s -outfmt \'6 std qlen slen|' -num_threads 1 -out %s" %(arg_dict["fasta_file"], arg_dict["db_name"], output_path)
        code_results = system_call(cmd_makeblast_results)
        if (code_results == 'FAIL'):
            print ('****ERROR: Some error happened during the blastp command')
            print (cmd_makeblast_results)
            exit()  
    
    else:
        print("#####")
        print("Please provide a FASTA file")
        print(compt)
        print("#####")
            
        

#https://github.com/HCGB-IGTP/HCGB_python_functions/blob/2b3b1132fb885c8cb22f4d10fd4c00c25fa10fb8/HCGB/functions/system_call_functions.py
def system_call(cmd, returned=False, message=True):
    """Generates system call using subprocess.check_output"""
    ## call system
    ## send command
    if (message):
        print ("[** System: %s **]" % cmd)

    try:
        out = subprocess.check_output(cmd, shell = True)
        if (returned):
            return (out)
        return ('OK')
    except subprocess.CalledProcessError as err:
        if (returned):
            return (err.output)
        if (message):
            print ("** ERROR **")
            print (err.output)
            print ("** ERROR **")
        
        return ('FAIL')
    

    
############
####################
## starting
####################
####################
    

parser = ArgumentParser(prog='makeblastDB',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description="Create a BLAST database")
parser.add_argument("-d", "--db_name", metavar="", help="New database name")
parser.add_argument("-f", "--fasta_file", metavar="", help="Proteins sequences FASTA file")
parser.add_argument("-e", "--executable", metavar="", help="BLAST executable: makeblastdb,...")
# parser.add_argument("-o", "--out_folder", metavar= "", help="Results folder")
parser.add_argument("--debug", action="store_true", default=False)   

arg = parser.parse_args()
arg_dict = vars(arg)
if arg.debug:
    print(arg)  
    
if arg.db_name is None:
    print("#####")
    print("Please provide a name for your new database")
    print("#####")
 
if arg.fasta_file is None:
    print("#####")
    print("Please provide a proteins sequences FASTA file")
    print("#####")
     
if arg.executable is None:
    print("#####")
    print("Please provide an executable")
    print("#####")
         
if arg.debug:
    print("##DEBUG: ##")
    print("arguments dictionary: ")
    print(arg)
  
if __name__ == "__main__":
    makeblastdb(arg_dict)
    #makeblast_results(arg_dict)  
  
    