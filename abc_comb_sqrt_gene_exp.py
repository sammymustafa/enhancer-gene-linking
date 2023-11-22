from setup import *
from correlations import pair_corr

# PREDICTION 4: abc_comb_sqrt

# set abc_comb_sqrt threshold at 0.05
abc_comb_sqrt = pair_corr[pair_corr['abc_comb_sqrt'] > 0.05]

# add column saying if overlap from rna data is present within each gene-enhancer pair
abc_comb_sqrt['overlap'] = abc_comb_sqrt.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = abc_comb_sqrt[abc_comb_sqrt['overlap']]

# merging data by gene name
merge = pd.merge(overlap, rna_data, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['symbol', 'gene_ids', 'name', 'enhancer', 'dist', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from eqtl
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
abc_comb_sqrt_trans = filtered_merge['gene_ids'].unique().tolist()
filtered_merge['gene_ids'].unique().tolist()