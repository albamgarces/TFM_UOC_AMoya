'''
Created on 27 oct. 2020

@author: alba
'''
from BCBio import GFF

in_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.gff"
#in_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot2.gff"

##SeqRecord to the IDs referenced
# in_handle = open(in_file)
# for rec in GFF.parse(in_handle):
#     print(rec.features)
#     print(rec)
# in_handle.close()

import pprint
from BCBio.GFF import GFFExaminer

##feature attributes summary with count for the number of times appeared
examiner = GFFExaminer()
in_handle = open(in_file)
#pprint.pprint(examiner.parent_child_map(in_handle))
pprint.pprint(examiner.available_limits(in_handle))
in_handle.close()

## subset features to parse
# limit_info = dict(gff_type=["CDS"])
#  
# in_handle = open(in_file)
# for rec in GFF.parse(in_handle, limit_info=limit_info):
#     print(rec.features) ##all records##
#     print(rec.features[0]) ##first record##
# #for rec in GFF.parse(in_handle, target_lines=1000): ##limit the lines read at once time
# in_handle.close()
# 
# from Bio import SeqIO
##annotation + genomic sequence
# in_seq_file ="/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_genomic.fna"
# in_seq_handle = open(in_seq_file)
# seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
# in_seq_handle.close()
# 
# in_handle = open(in_file)
# for rec in GFF.parse(in_handle, base_dict=seq_dict):
#     print(rec)
# in_handle.close()

# in_seq_file ="/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.fna"
# in_seq_handle = open(in_seq_file)
# limit_info = dict(gff_type=["CDS"])
# in_handle = open(in_file)
# ref_recs = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))    
# for rec in GFF.parse(in_handle, limit_info = limit_info, base_dict=ref_recs):
#     print(rec.features)
# in_handle.close()
