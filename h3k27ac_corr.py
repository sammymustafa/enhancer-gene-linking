from setup import *
from correlations import pair_corr

# Predicition 1: H3K27ac * Cor(H3K27ac, RNA)^2

# set correlation threshold at 0.07
corr_1 = pair_corr[pair_corr['corr_1'] > 0.07]

# add column saying if overlap from eqtl data is present within each gene-enhancer pair
corr_1['overlap'] = corr_1.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = corr_1[corr_1['overlap']]

# merging data by gene name
merge = pd.merge(overlap, eqt1_filt, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['gene_ids', 'name', 'enhancer', 'dist', 'tss_distance', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from eqtl
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
filtered_merge

# calculate proportion accuracy
prop_acc = len(filtered_merge) / len(corr_1)
print("The predicitive accuracy with a 0.07 threshold on H3K27ac * Cor(H3K27ac, RNA)^2 when linking on the eQTLs is", prop_acc)