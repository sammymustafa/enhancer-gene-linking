from setup import *
from correlations import pair_corr

# PREDICTION 5: H3K27ac * Cor(H3K27ac, RNA)^2

# set correlation threshold at 0.07
corr_1 = pair_corr[pair_corr['corr_1'] > 0.07]

# add column saying if overlap from rna data is present within each gene-enhancer pair
corr_1['overlap'] = corr_1.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = corr_1[corr_1['overlap']]

# merging data by gene name
merge = pd.merge(overlap, rna_data, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['symbol', 'gene_ids', 'name', 'enhancer', 'dist', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from rna
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
corr_1_trans = filtered_merge['gene_ids'].unique().tolist()
filtered_merge['gene_ids'].unique().tolist()