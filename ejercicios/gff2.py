'''
Created on 28 oct. 2020

@author: alba
'''

# gff_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.gff"
# ref_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.fna"

import sys
import os
import operator
from functools import reduce

import Bio
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import BCBio
from BCBio import GFF
from basic_GFF_parsing import in_handle

def main (gff_file, ref_file):
    with open (ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))
    base, ext = os.path.splitext(gff_file)
    out_file = "%s-proteins.fa" % base
    with open(out_file, "w") as out_handle:
        SeqIO.write(protein_recs(gff_file, ref_recs), out_handle, "fasta")

##parse the output. Generate SeqRecord and SeqFeatures for predictions
def gff_predct (in_handle, ref_recs):
    limit_info = dict(gff_type=["CDS"])    
    for rec in GFF.parse(in_handle, limit_info = limit_info, base_dict=ref_recs):
        yield (rec)

#generate protein records from gff_predct        
def protein_recs(gff_file, ref_recs):
    with open(gff_file) as in_handle:
        for rec in gff_predct(in_handle, ref_recs):
            for feature in rec.features:
                seq_exons = []
                for CDS in feature.sub_features:
                    seq_exons.append(
                        rec.seq[CDS.location.nofuzzy_start : CDS.location.nofuzzy_end]
                        )
                gene_seq =reduce(operator.add, seq_exons,"")
                if feature.strand == -1:
                    gene_seq = gene_seq.reverse_complement()
                protein_seq = gene_seq.translate()
                yield SeqRecord(protein_seq, feature.qualifiers["ID"],"","")
                    
        
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print (__doc__)
        sys.exit()
    main(*sys.argv[1:])

        

        
        