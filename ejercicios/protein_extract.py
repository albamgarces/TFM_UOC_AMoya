'''
Created on 21 oct. 2020

@author: alba
'''

if __name__ == '__main__':
    pass
import os.path
import Bio
from Bio import SeqIO

#https://gist.github.com/peterk87/5422267
#https://www.biostars.org/p/151891/

#seq records for the genbank file

#genbank_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gbf"
def extension (file_name):
    compt = {".gbf", ".gff", ".fna", ".faa"}
    name, extension = os.path.splitext(file_name)
    if extension in compt:
        if extension == ".gbf":
            seq_rec = [rec for rec in SeqIO.parse("/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gbf", "genbank")]
            for rec in seq_rec:
                feats = [feat for feat in rec.features if feat.type == "CDS"]
                for feat in feats:
                    print ("locus_tag: ", feat.qualifiers["locus_tag"][0])
                    print("product: ", feat.qualifiers["product"])
                    print (feat.qualifiers["translation"])
        else: 
            print ("The file extension is {}". format(extension))        
    else:
        print("This file is not compatible")
extension("example_annot.gbf")   

