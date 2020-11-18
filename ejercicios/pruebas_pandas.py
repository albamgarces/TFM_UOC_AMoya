'''
Created on 8 nov. 2020

@author: alba
'''
import os
import pandas as pd
import numpy as np
import protein_gbf
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
import BCBio
from BCBio import GFF

# data = np.array([["","col1", "col2"], ["Fila1", 11,22], ["Fila2", 33,44]])
# print(pd.DataFrame(data=data[1:,1:], index=data[1:,0], columns=data[0,1:]))
# 
# series = pd.Series()

gbf_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gbf"
out_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot-proteins.fa"
# for rec in protein_gbf.gbf_parser(gbf_file):
#     series = pd.Series({rec.id:str(rec.seq)})
#     index=np.array[str(rec.id)]
#     data= np.array([["","id","sequence"], ["",rec.id, str(rec.seq)]])
#     df = pd.DataFrame(data=data[1:,1:], index=data[1:,0], columns=data[0,1:])
#     print(df)


# for rec in protein_gbf.gbf_parser(gbf_file):
#     datos = {"id":[rec.id], "sequence":[str(rec.seq)]}
#     df=pd.DataFrame(datos)
#     print(df)
# for rec in protein_gbf.gbf_parser(gbf_file):
#     series = pd.Series({"id":rec.id, "sequence":str(rec.seq)})
#     df=pd.DataFrame(series)
#     print(df)

# with open(out_file) as fasta_file:
#         identifiers = []
#         seq = []
#         for title, sequence in SimpleFastaParser(fasta_file):
#                         identifiers.append(title)
#                         seq.append(sequence)
#                         s1 = pd.Series(identifiers, name="ID")
#                         s2 = pd.Series(seq, name="sequence")
#                         gbs_df = pd.DataFrame(dict(ID=s1, sequence=s2)).set_index(["ID"])
#                         print(gbs_df)
                        
# with open(out_file) as fasta_file:
#         identifiers = []
#         seq = []
#         for title, sequence in SimpleFastaParser(fasta_file):
#                         identifiers.append(title)
#                         seq.append(sequence)
#                         s1 = pd.Series(identifiers, name="ID")
#                         s2 = pd.Series(seq, name="sequence")
#                         gbf_df = pd.DataFrame(dict(ID=s1, sequence=s2)).set_index(["ID"])
#                         gbf_df.reset_index().to_csv("out_gbf.csv", header=True, index=False)
#        
                           
gbf_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot.gbf"                  
subset_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/subset.gbf"

# annot_df = pd.DataFrame()
# 
# for rec in SeqIO.parse(gbf_file, "genbank"):
#     print("*************")
#     print(rec)
#     for feature in rec.features:
#         print()
#         print("#FEATURE:")
#         print(feature)
#         print()
#         if feature.type=="CDS":
#             if int(feature.strand) > 0:
#                 strand = "pos"
#             else:
#                 strand = "neg"
#                 
#             protID = rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
#             print(protID)
#             annot_df.loc[protID, "strand"] = strand
#             
#             print(feature.qualifiers)
#             
#             if not "gene" in feature.qualifiers:
#                 pass
#             else:
#                 annot_df.loc[protID,["gene"]]=feature.qualifiers["gene"]


# columns = ['locus_tag', 'gene', 'product', 'EC_number',
#            'db_xref', 'start', 'end', 'strand', 'translation']
# annot_df = pd.DataFrame(data=None, columns=columns)
# print(annot_df)
# 
# 
# for rec in SeqIO.parse(gbf_file, "genbank"):
#     print("*************")
#     print(rec)
#     for feature in rec.features:
#         print()
#         print("#FEATURE:")
#         print(feature)
#         print()
#         if feature.type=="CDS":
#             if int(feature.strand) > 0:
#                 strand = "pos"
#             else:
#                 strand = "neg"
#                 
#             protID = rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
#             print(protID)
#             annot_df.loc[protID, ["start", "end", "strand"]] = [feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
#             qualif = feature.qualifiers
#             print("#Qualifiers:")
#             print(qualif)
#             print()
#             for keys, values in qualif.items():
#                 print(keys)
#                 print(values)
#                 if keys in columns:
#                     annot_df.loc[protID,[keys]] = [values[0]]
# #             
#                     
# print(annot_df)
# annot_df.to_csv("out_subset.csv", header=True)

#             columns = ['gene']
#             index=feature.qualifiers["locus_tag"]
#             df =pd.DataFrame(columns=columns, index=index)
#             df.loc[:,["gene"]]=feature.qualifiers["gene"]
#             print(df)
#             df.to_csv("out_subset.csv", header=True)
# 
#         rec.ID_feature.location.nofuzzy_start_feature.location.nofuzzy_end_feature.strand
#         
        
#                 
#                     
#           
#                 gene_seq = Seq.Seq(feature.qualifiers["translation"][0])
#                 yield(SeqRecord(gene_seq, feature.qualifiers["locus_tag"][0],"",""))
#                  
# 
# subset_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/subset.gbf"
# for rec in protein_gbf.gbf_parser(subset_file):
#     print(rec)
# columns = ['locus_tag', 'gene', 'product', 'EC_number',
#            'db_xref', 'start', 'end', 'strand', 'translation']
# df =pd.DataFrame(columns=columns)
# df = df.fillna(0)
# data=np.array([np.arange(10)]*9).T
# df = pd.DataFrame(data, columns=columns).set_index(["locus_tag"])
# print(df)

## dataframe from gbf file ##
# columns = ['locus_tag', 'gene', 'product', 'EC_number',
#            'db_xref', 'start', 'end', 'strand', 'translation']
# annot_df = pd.DataFrame(data=None, columns=columns)
# 
# 
# 
# for rec in SeqIO.parse(gbf_file, "genbank"):
#     print("*************")
#     print(rec)
#     for feature in rec.features:
#         print()
#         print("#FEATURE:")
#         print(feature)
#         print()
#         if feature.type=="CDS":
#             if int(feature.strand) > 0:
#                 strand = "pos"
#             else:
#                 strand = "neg"
#                 
#             protID = rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
#             print(protID)
#             annot_df.loc[protID, ["start", "end", "strand"]] = [feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
#             qualif = feature.qualifiers
#             print("#Qualifiers:")
#             print(qualif)
#             print()
#             for keys, values in qualif.items():
#                 annot_df.loc[protID,[keys]] = [values[0]]
# #             
#                     
# print(annot_df)
# annot_df.to_csv("out_subset.csv", header=True)

## dataframe from gff file ##

#gff_file = "/home/alba/git/TFM_UOC_AMoya/data/example_denovo/example_annot2.gff"
ref_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.fna"
subset_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/subset.gff"
gff_file = "/home/alba/git/TFM_UOC_AMoya/data/example_NCBI/example1_genomic.gff"

columns = ['type', 'locus_tag', 'gene', 'product', 'EC_number',
           'Dbxref', 'start', 'end', 'strand', 'translation']
annot_df = pd.DataFrame(data=None, columns=columns)

with open (ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))
with open(gff_file) as in_handle:
    ##parse the output. Generate SeqRecord and SeqFeatures for predictions
    limit_info = dict(gff_type=["CDS", "pseudogene"])    
    for rec in GFF.parse(in_handle, limit_info=limit_info, base_dict=ref_recs):
     
        for feature in rec.features:
            if feature.strand == -1:
                strand = "neg"
            else:
                strand = "pos"
                
            protID = feature.type + "_" + rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
            annot_df.loc[protID, ["type", "start", "end", "strand"]] = [feature.type, feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
            qualif = feature.qualifiers
            for keys, values in qualif.items():
                annot_df.loc[protID,[keys]] = [values[0]]
        
        ## get gene sequence
            gene_seq = Seq.Seq(str(rec.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end][:-3]))

            if feature.strand == -1:
                gene_seq = gene_seq.reverse_complement()
            protein_seq = gene_seq.translate()
                #yield(SeqRecord(protein_seq, feature.qualifiers["ID"][0], "", ""))
            annot_df.loc[protID,["translation"]] = [protein_seq]
                
base, ext = os.path.splitext(gff_file)
csv_file = "%s-gff_df.csv" % base            
annot_df.to_csv(csv_file, header=True)

print(annot_df)


# for rec in SeqIO.parse(gff_file, "genbank"):
#     print("*************")
#     print(rec)
#     for feature in rec.features:
#         print()
#         print("#FEATURE:")
#         print(feature)
#         print()
#         if feature.type=="CDS":
#             if int(feature.strand) > 0:
#                 strand = "pos"
#             else:
#                 strand = "neg"
#                 
#             protID = rec.id + "_" + str(feature.location.nofuzzy_start) + "_" + str(feature.location.nofuzzy_end) + "_" + strand
#             print(protID)
#             annot_df.loc[protID, ["start", "end", "strand"]] = [feature.location.nofuzzy_start, feature.location.nofuzzy_end, strand]
#             qualif = feature.qualifiers
#             print("#Qualifiers:")
#             print(qualif)
#             print()
#             for keys, values in qualif.items():
#                 annot_df.loc[protID,[keys]] = [values[0]]
# #             
#                     
# print(annot_df)
# annot_df.to_csv("out_subset.csv", header=True)



                     
                        