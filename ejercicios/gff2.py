'''
Created on 28 oct. 2020

@author: alba
'''

gff_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.gff"
ref_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.fna"

import sys
import os
import operator
from functools import reduce

import Bio
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import BCBio
from BCBio import GFF
## from basic_GFF_parsing import in_handle ### JF:  Esto para que? Estas importando un script pero no hace nada

def main (gff_file, ref_file, debug=False):
    with open (ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))
    
    ## JF
    if (debug):
        print ("## DEBUG:")
        print (ref_recs)
        print ()

    base, ext = os.path.splitext(gff_file)
    out_file = "%s-proteins.fa" % base

    ## JF
    if (debug):
        print ("## DEBUG:")
        print ("base:" , base)
        print ("ext:" , ext)
        print ("out_file:" , out_file)
        print ()

    with open(out_file, "w") as out_handle:
        SeqIO.write(protein_recs(gff_file, ref_recs, debug=debug), out_handle, "fasta")

#generate protein records from gff_predct        
def protein_recs(gff_file, ref_recs, debug=False):

    with open(gff_file) as in_handle:
        ##parse the output. Generate SeqRecord and SeqFeatures for predictions
        limit_info = dict(gff_type=["CDS"])    
        for rec in GFF.parse(in_handle, limit_info = limit_info, base_dict=ref_recs):
            ## JF
            if (debug):
                print ("## DEBUG: rec")
                print (rec)
                print ()
        
            for feature in rec.features:
                ## JF
                if (debug):
                    print ("## DEBUG: feature")
                    print (feature)
                    print ()

                    ## check feature qualifiers: some might not be always present: pseudo, gene, Name
                    print (feature.location)
                    print (feature.location.nofuzzy_start)
                    print (feature.location.nofuzzy_end)
                    print (feature.strand)
                    print (feature.id)
                    print (feature.qualifiers)

                    #print (feature.qualifiers["gene"][0])
                    #print (feature.qualifiers["Name"][0])
                    #print (feature.qualifiers["Parent"][0])
                    #print (feature.qualifiers["locus_tag"][0])
        
        ## get gene sequence
                gene_seq = Seq.Seq(str(rec.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end][:-3]))

                ## JF
                if (debug):
                    print ("## DEBUG: gene_seq")
                    print (gene_seq)

                if feature.strand == -1:
                    gene_seq = gene_seq.reverse_complement()
                protein_seq = gene_seq.translate()
                yield(SeqRecord(protein_seq, feature.qualifiers["ID"][0], "", ""))

        
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print (__doc__)
        
        ## JF: Tambien es util mostrar que se necesta
        print ("## Usage gff_parser")
        print ("python %s gff_file ref_fasta_file\n" %sys.argv[0])

        sys.exit()

    ## prueba las dos lineas de codigo y entiende la diferencia.
    main(*sys.argv[1:], debug=True)
    #main(*sys.argv[1:])

    # la variable debug no es obligatoria. tiene un "por defecto definido"
    # fijate donde se define main (linea 22)
    # Se utiliza el "=" para indicar el default.
        
        
        
