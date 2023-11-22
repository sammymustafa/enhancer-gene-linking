from setup import *
from correlations import pair_corr

# PREDICTION 6: DNase-seq * Cor(DNase-seq, RNA)^2

# set correlation threshold at 0.07
corr_2 = pair_corr[pair_corr['corr_2'] > 0.07]

# add column saying if overlap from rna data is present within each gene-enhancer pair
corr_2['overlap'] = corr_2.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = corr_2[corr_2['overlap']]

# merging data by gene name
merge = pd.merge(overlap, rna_data, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['symbol', 'gene_ids', 'name', 'enhancer', 'dist', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from rna
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
corr_2_trans = filtered_merge['gene_ids'].unique().tolist()
filtered_merge['gene_ids'].unique().tolist()