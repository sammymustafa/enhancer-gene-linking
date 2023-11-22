from setup import *

# filter rna-seq data to have relevant data
rna_data = rna["BSS00385"].var
chr_values = h3k27ac_chmm["Chromosome"].unique()
rna_data = rna_data[rna_data['seqnames'].isin(chr_values)]
rna_data

# make function that produces boolean to test if overlap is present
def check_overlap(row):
    return any((row['Start'] <= rna_data['start']) & (row['End'] >= rna_data['end']))