from setup import *
from correlations import pair_corr

# PREDICTION 2: abc_h3k27ac

# set abc_h3k27ac threshold at 0.05
abc_h3k27ac = pair_corr[pair_corr['abc_h3k27ac'] > 0.05]

# add column saying if overlap from rna data is present within each gene-enhancer pair
abc_h3k27ac['overlap'] = abc_h3k27ac.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = abc_h3k27ac[abc_h3k27ac['overlap']]

# merging data by gene name
merge = pd.merge(overlap, rna_data, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['symbol', 'gene_ids', 'name', 'enhancer', 'dist', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from rna
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
abc_h3k27ac_trans = filtered_merge['gene_ids'].unique().tolist()
filtered_merge['gene_ids'].unique().tolist()