from setup import *

# PREDICTION 1: abc_dist

# set abc_dist threshold at 0.05
abc = gen_enh_dist[gen_enh_dist['abc_dist'] > 0.05]

# filter eqtl values of interest within chromosomes 20-22
chr_values = h3k27ac_chmm["Chromosome"].unique()
eqt1_filt = eqtl[eqtl['chrom'].isin(chr_values)].rename(columns = {'gene_id': 'gene_ids'})
eqt1_filt['tss_distance'] = abs(eqt1_filt['tss_distance'])
eqt1_filt = eqt1_filt[eqt1_filt['tss_distance'] < 100000]

# make function that produces boolean to test if overlap is present (takes ~1.5 mins to run)
def check_overlap(row):
    return any((row['Start'] <= eqt1_filt['start']) & (row['End'] >= eqt1_filt['end']))
# add column saying if overlap from eqtl data is present within each gene-enhancer pair
abc['overlap'] = abc.apply(check_overlap, axis=1)
# filtering to see samples with overlapping regions
overlap = abc[abc['overlap']]

# merging data by gene name
merge = pd.merge(overlap, eqt1_filt, on = 'gene_ids', how = 'inner')
# fixing data for analysis
merge = merge[['gene_ids', 'name', 'enhancer', 'dist', 'tss_distance', 'Start', 'End', 'start', 'end']]
# only keep rows with overlap with the matched gene from eqtl
merge['overlap'] = (merge['start'] <= merge['End']) & (merge['end'] >= merge['Start'])
filtered_merge = merge[merge['overlap']]
filtered_merge

# calculate proportion accuracy
prop_acc = len(filtered_merge) / len(abc)
print("The predicitive accuracy with a 0.05 threshold on abc_distance when linking on the eQTLs is", prop_acc)