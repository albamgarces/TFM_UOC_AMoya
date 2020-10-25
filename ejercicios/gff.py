'''
Created on 24 oct. 2020

@author: alba
'''
from __future__ import with_statement #esto no se para qué es
from _ast import With

if __name__ == '__main__':
    pass

#https://biopython.org/wiki/GFF_Parsing
#https://biopython.org/wiki/Gene_predictions_to_protein_sequences

#GFF records contain annotation data; seq info is in a separate FASTA file
import sys
import os
import operator
from functools import reduce

import Bio
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import BCBio
from BCBio import GFF

def gff_parse (gff_file, ref_file):
    with open (ref_file) as in_handle:
        #parse input FASTA into a python dict
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))
    base, ext = os.path.splitext(gff_file)
    out_file = "%s-proteins.faa" % base
    with open(out_file, "w") as out_handle:
        #write output FASTA file
        SeqIO.write(protein_recs(gff_file, ref_recs), out_handle, "fasta")

#generate protein records from gene predictions
def protein_recs(gff_file, ref_recs):
    with open (gff_file) as in_handle:
        for rec in gff_predct(in_handle, ref_recs):
            for feature in rec.features:
                seq_exons = []
                for cds in feature.sub_features:
                    seq_exons.append(rec.seq[
                        cds.location.nofuzzy_start:
                    cds.location.nofuzzy_end])
                    gene_seq = Seq.Seq(str(reduce(operator.add, seq_exons, "")))
                    if feature.strand == -1:
                        gene_seq = gene_seq.reverse_complement()
                    protein_seq = gene_seq.translate()
                    yield SeqRecord(protein_seq, feature.qualifiers["ID"][0], "", "")

#parse gff output, generating SeqRecord and SeqFeatures for predictions
def gff_predct(in_handle, ref_recs):
    for rec in GFF.parse(in_handle, target_lines=1000, base_dict=ref_recs):
        yield rec
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print (__doc__)
        sys.exit()
    gff_parse(*sys.argv[1:])       

    