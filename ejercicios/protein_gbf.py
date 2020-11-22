'''
Created on 30 oct. 2020

@author: alba
'''
from Bio import SeqIO, Seq
import os
import sys
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np

from Bio.SeqIO.FastaIO import SimpleFastaParser


def main (gbf_file, debug=False):
    with open(gbf_file, "r") as input_handle:
        sequences = [rec for rec in SeqIO.parse(input_handle, "genbank")]
#     if (debug):
#         print ("## DEBUG:")
#         print (sequences)
#         print ()
        
    base, ext = os.path.splitext(gbf_file)
    out_file = "%s-gbf_proteins.fa" % base

    if (debug):
        print ("## DEBUG:")
        print ("base:" , base)
        print ("ext:" , ext)
        print ("out_file:" , out_file)
        print ()
        
    with open(out_file, "w") as output_handle:
        SeqIO.write(gbf_parser(gbf_file, debug=debug), output_handle, "fasta")
    
    


def gbf_parser(gbf_file, debug=False):
    
    #create an empty dataframe. The most important info as first columns
    #all the other file info will be filled next
    columns = ['locus_tag', 'protein_id', 'gene', 'start',
               'end', 'strand', 'pseudo', 'product', 'db_xref',
               'EC_number', 'old_locus_tag', 'inference']
    annot_df = pd.DataFrame(data=None, columns=columns)
    
    for rec in SeqIO.parse(gbf_file, "genbank"):
        if (debug):
            print("## DEBUG: rec")
            print(rec)
            print()
            
        for feature in rec.features:
            
            #sort by CDS type. Duplicate genes analysis needs coding regions to proteins.
            if feature.type=="CDS":
                if int(feature.strand) > 0:
                    strand = "pos"
                else:
                    strand = "neg"
            
                #we create an ID for each entry     
                protID = feature.type + "_" + rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
#                 if (debug):
#                     print("#DEBUG: protID")
#                     print(protID)
#                     print()
                        
                annot_df.loc[protID, ["start", "end", "strand"]] = [feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
                qualif = feature.qualifiers
#                 if (debug):
#                     print("#DEBUG: Qualifiers")
#                     print(qualif)
#                     print()
                for keys, values in qualif.items():
                 
                        #fill the dataframe info
                    if keys not in columns:
                        continue
                    annot_df.loc[protID,[keys]] = [values[0]]
                    
                    if keys=="pseudo":
                        pseudo = "True"
                        annot_df.loc[protID,["pseudo"]] = [pseudo]
            
                        #we create a FASTA file with protein sequences
            
#                     if (debug):
#                         print ("## DEBUG: feature")
#                         print (feature)
#                         print ()
                    
#                         print (feature.location)
#                         print (feature.location.nofuzzy_start)
#                         print (feature.location.nofuzzy_end)
#                         print (feature.strand)
#                         print (feature.id)
#                         print (feature.qualifiers)   
                   
                if keys=="translation":
                    #pseudogenes have no translation item
                    gene_seq = Seq.Seq(feature.qualifiers["translation"][0])
                else:
                    pass
                yield(SeqRecord(gene_seq, feature.qualifiers["locus_tag"][0],"",""))
        
    
    base, ext = os.path.splitext(gbf_file)
    csv_file = "%s-gbf_df.csv" % base            
    annot_df.to_csv(csv_file, header=True)
    
    if (debug):
        print("## DEBUG: dataframe")
        print(annot_df)
               
               
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__)
        
       
        print ("## Usage protein_gbf")
        print ("python %s gbf_file\n" %sys.argv[0])
       

        sys.exit()

    main(*sys.argv[1:], debug=True)
  
    
    
    

            
            