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
    columns = ['locus_tag', 'type', 'gene', 'product', 'EC_number',
           'db_xref', 'start', 'end', 'strand']
    annot_df = pd.DataFrame(data=None, columns=columns)
    
    for rec in SeqIO.parse(gbf_file, "genbank"):
        if (debug):
            print("## DEBUG: rec")
            print(rec)
            print()
            
        for feature in rec.features:
            #we create a dataframe with all items
            if int(feature.strand) > 0:
                strand = "pos"
            else:
                strand = "neg"
                
            protID = feature.type + "_" + rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
#             if (debug):
#                 print("#DEBUG: protID")
#                 print(protID)
#                 print()
                        
            annot_df.loc[protID, ["type", "start", "end", "strand"]] = [feature.type, feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
            qualif = feature.qualifiers
#             if (debug):
#                 print("#DEBUG: Qualifiers")
#                 print(qualif)
#                 print()
            for keys, values in qualif.items():
                annot_df.loc[protID,[keys]] = [values[0]]
            
            #we create a FASTA file with protein sequences
            if feature.type=="CDS":
#                 if (debug):
#                     print ("## DEBUG: feature")
#                     print (feature)
#                     print ()
                    
#                     print (feature.location)
#                     print (feature.location.nofuzzy_start)
#                     print (feature.location.nofuzzy_end)
#                     print (feature.strand)
#                     print (feature.id)
#                     print (feature.qualifiers)   
                   
                if keys=="translation":
                    #pseudogenes have no translation item
                    gene_seq = Seq.Seq(feature.qualifiers["translation"][0])
                else:
                    pass
                yield(SeqRecord(gene_seq, feature.qualifiers["locus_tag"][0],"",""))
    
    base, ext = os.path.splitext(gbf_file)
    csv_file = "%s-df_gbf.csv" % base            
    annot_df.to_csv(csv_file, header=True)
    
    if (debug):
        print("## DEBUG: dataframe")
        print(annot_df)
    
 
def test():
    print("ejemplo")
#     
# def gbf_frame (gbf_fasta):
#     with open(gbf_fasta) as fasta_file:
#         identifiers = []
#         seq = []
#         for title, sequence in SimpleFastaParser(fasta_file):
#                         identifiers.append(title)
#                         seq.append(sequence)
#                         s1 = pd.Series(identifiers, name="ID")
#                         s2 = pd.Series(seq, name="sequence")
#                         gbf_df = pd.DataFrame(dict(ID=s1, sequence=s2)).set_index(["ID"])
#                         gbf_df.reset_index().to_csv("out_gbf.csv", header=True, index=False)
    
                                                 
    

               
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__)
        
       
        print ("## Usage protein_gff")
        print ("python %s gff_file ref_fasta_file\n" %sys.argv[0])
        test()

        sys.exit()

    main(*sys.argv[1:], debug=True)
  
    
    
    

            
            