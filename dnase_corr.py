from setup import *
from correlations import pair_corr

# Predicition 2: DNase-seq * Cor(DNase-seq, RNA)^2

# set correlation threshold at 0.07
corr_2 = pair_corr[pair_corr['corr_2'] > 0.07]

# add column saying if overlap from eqtl data is present within each gene-enhancer pair
corr_2['overlap'] = corr_2.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = corr_2[corr_2['overlap']]

# merging data by gene name
merge = pd.merge(overlap, eqt1_filt, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['gene_ids', 'name', 'enhancer', 'dist', 'tss_distance', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from eqtl
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
filtered_merge

# calculate proportion accuracy
prop_acc = len(filtered_merge) / len(corr_2)
print("The predicitive accuracy with a 0.07 threshold on DNase-seq * Cor(DNase-seq, RNA)^2 when linking on the eQTLs is", prop_acc)