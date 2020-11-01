'''
Created on 30 oct. 2020

@author: alba
'''
from Bio import SeqIO, Seq
import os
import sys
from Bio.SeqRecord import SeqRecord


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
         
                gene_seq = Seq.Seq(feature.qualifiers["translation"][0])
                yield(SeqRecord(gene_seq, feature.qualifiers["locus_tag"][0],"",""))
           
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print (__doc__)
        
       
        print ("## Usage gff_parser")
        print ("python %s gff_file ref_fasta_file\n" %sys.argv[0])

        sys.exit()

    main(*sys.argv[1:], debug=True)
  
    
    
    

            
            