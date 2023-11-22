from setup import *

# PREDICTION 2: abc_h3k27ac

# set abc_h3k27ac threshold at 0.05
abc_h3k27ac = gen_enh_dist[gen_enh_dist['abc_h3k27ac'] > 0.05]

# add column saying if overlap from eqtl data is present within each gene-enhancer pair
abc_h3k27ac['overlap'] = abc_h3k27ac.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = abc_h3k27ac[abc_h3k27ac['overlap']]

# merging data by gene name
merge = pd.merge(overlap, eqt1_filt, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['gene_ids', 'name', 'enhancer', 'dist', 'tss_distance', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from eqtl
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
filtered_merge

# calculate proportion accuracy
prop_acc = len(filtered_merge) / len(abc_h3k27ac)
print("The predicitive accuracy with a 0.05 threshold on abc_h3k27ac when linking on the eQTLs is", prop_acc)