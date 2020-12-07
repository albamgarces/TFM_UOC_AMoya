'''
Created on 29 nov. 2020

@author: alba
'''
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

##############################
##### BLAST over internet#####
##############################

#open the sequence file using IO module
# sequence_Data = open("/home/alba/git/TFM_UOC_AMoya/data/subset_output/gff_proteins.fa").read()
# #print(sequence_Data)
# 
# #call the qblast function passinq seq data as main parameter
# result_handle = NCBIWWW.qblast("blastp", "nr", sequence_Data)
# 
# #create the results file
# with open("results_internet.xml", "w") as save_file:
#     blast_results = result_handle.read()
#     save_file.write(blast_results)
# 
# #get BLAST record objects
# result_handle = open("results_internet.xml")    
# blast_records = NCBIXML.parse(result_handle)
# 
# for blast_record in blast_record:
#     blast_records = list(blast_records)
#     print(blast_records)    
    
##############################
##### BLAST locally#####
##############################

#create a database with our FASTA file
