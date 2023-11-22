from setup import *

# PREDICTION 4: abc_comb_sqrt

# set abc_comb_sqrt threshold at 0.05
abc_comb_sqrt = gen_enh_dist[gen_enh_dist['abc_comb_sqrt'] > 0.05]

# add column saying if overlap from eqtl data is present within each gene-enhancer pair
abc_comb_sqrt['overlap'] = abc_comb_sqrt.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = abc_comb_sqrt[abc_comb_sqrt['overlap']]

# merging data by gene name
merge = pd.merge(overlap, eqt1_filt, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['gene_ids', 'name', 'enhancer', 'dist', 'tss_distance', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from eqtl
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
filtered_merge

# calculate proportion accuracy
prop_acc = len(filtered_merge) / len(abc_comb_sqrt)
print("The predicitive accuracy with a 0.05 threshold on abc_distance when linking on the eQTLs is", prop_acc)