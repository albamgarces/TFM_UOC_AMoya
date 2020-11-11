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
    if (debug):
        print ("## DEBUG:")
        print (sequences)
        print ()
        
    base, ext = os.path.splitext(gbf_file)
    out_file = "%s-proteins.fa" % base

    if (debug):
        print ("## DEBUG:")
        print ("base:" , base)
        print ("ext:" , ext)
        print ("out_file:" , out_file)
        print ()
        
    with open(out_file, "w") as output_handle:
        SeqIO.write(gbf_parser(gbf_file, debug=debug), output_handle, "fasta")

def gbf_parser(gbf_file, debug=False):
    for rec in SeqIO.parse(gbf_file, "genbank"):
        if (debug):
            print("## DEBUG: rec")
            print(rec)
            print()
            
        for feature in rec.features:
            if feature.type=="CDS":
                if (debug):
                    print ("## DEBUG: feature")
                    print (feature)
                    print ()
                    
                    print (feature.location)
                    print (feature.location.nofuzzy_start)
                    print (feature.location.nofuzzy_end)
                    print (feature.strand)
                    print (feature.id)
                    print (feature.qualifiers)
                   
         
                gene_seq = Seq.Seq(feature.qualifiers["translation"][0])
                yield(SeqRecord(gene_seq, feature.qualifiers["locus_tag"][0],"",""))
                
                protein_annot = SeqRecord(seq = gene_seq, id = feature.qualifiers["locus_tag"][0],
                                          #name=feature.qualifiers["gene"][0],
                                        description=feature.qualifiers["product"][0])
                identifiers = []
                protein = []
                product = []
#                 for sequence, name, prod in protein_annot:
#                     identifiers.append(name)
#                     protein.append(sequence)
#                     product.append(prod)
#                     s1 = pd.Series(identifiers, name="ID")
#                     s2 = pd.Series(protein, name="sequence")
#                     s3 = pd.Series(product, name="product")
#                     gbf_df = pd.DataFrame(dict(ID=s1, sequence=s2, product=s3)).set_index(["ID"])
#                     gbf_df.reset_index().to_csv("out_gbf.csv", header=True, index=False)
#                     print(gbf_df)
                
                identifiers.append(protein_annot.id)
                protein.append(protein_annot.seq)
                product.append(protein_annot.description)
                s1 = pd.Series(identifiers, name="ID")
                s2 = pd.Series(protein, name="sequence")
                s3 = pd.Series(product, name="product")
                gbf_df = pd.DataFrame(dict(ID=s1, sequence=s2, product=s3)).set_index(["ID"])
                gbf_df.reset_index().to_csv("out_gbf.csv", header=True, index=False)
                print(gbf_df)
            
            #for rec in protein_gbf.gbf_parser(gbf_file):
#     series = pd.Series({"id":rec.id, "sequence":str(rec.seq)})
#     df=pd.DataFrame(series)
#     print(df)   
 
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
  
    
    
    

            
            