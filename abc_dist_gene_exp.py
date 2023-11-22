from setup import *
from correlations import pair_corr

# PREDICTION 1: abc_dist

# set abc_dist threshold at 0.05
abc = pair_corr[pair_corr['abc_dist'] > 0.05]

# add column saying if overlap from rna data is present within each gene-enhancer pair
abc['overlap'] = abc.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = abc[abc['overlap']]

# merging data by gene name
merge = pd.merge(overlap, rna_data, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['symbol', 'gene_ids', 'name', 'enhancer', 'dist', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from rna
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
abc_trans = filtered_merge['gene_ids'].unique().tolist()
filtered_merge['gene_ids'].unique().tolist()