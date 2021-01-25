'''
Created on 21 oct. 2020

@author: alba
'''

if __name__ == '__main__':
    pass
import os.path
import Bio
from Bio import SeqIO
import sys

#https://gist.github.com/peterk87/5422267
#https://www.biostars.org/p/151891/

#seq records for the genbank file

#genbank_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gbf"
# Es mejor pasar los ficheros como variable para automatizarlo y no tener que modificarlo cada vez
# especialmente si lo queremos distribuir y hacer automatico

def extension (file_given):
    compt = {".gbf", ".gbk", ".gff", ".fna", ".faa"}
    file_name_abs_path = os.path.abspath(file_given) ## obtiene el path absolute al fichero. te ahorra problemas de path incompletos futuros
    name, extension = os.path.splitext(file_name_abs_path)
    if extension in compt:
        if extension == ".gbf" or extension == ".gbk":
            ##seq_rec = [rec for rec in SeqIO.parse("/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gbf", "genbank")]
            seq_rec = [rec for rec in SeqIO.parse(file_given, "genbank")]
            for rec in seq_rec:
                feats = [feat for feat in rec.features if feat.type == "CDS"]
                for feat in feats:
                    print ("locus_tag: ", feat.qualifiers["locus_tag"][0])
                    print("product: ", feat.qualifiers["product"])
                    print (feat.qualifiers["translation"])
                    
#         elif extension == ".gff":
#             import BCBio
#             from BCBio import GFF
#             seq_rec= [for rec in GFF.parse(file_given)]
#             for rec in seq_rec:
#                 feats = [feat for feat in rec_]
            
                    
        else: 
            print ("The file extension is {}". format(extension))        
    else:
        print("This file is not compatible")

def main():
    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        print ("python %s annotation_file\n" %os.path.realpath(__file__))
        exit()        

    ## llama a la funcion con el argumento 1
    extension(sys.argv[1])   

    ## ATENCION: sys.argv[0] es el propio script de python!
	

## Cuando no se importe como un modulo este script, se ejecutara esto
if __name__ == '__main__':
    main()



