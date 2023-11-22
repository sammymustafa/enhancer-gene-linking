from setup import *

def abc_distance(distance, hic_power=0.87):
    scale = -4.80 + 11.63 * hic_power
    offset = np.clip(np.abs(distance), 5000, np.inf)
    return np.exp(scale + -1 * hic_power * np.log1p(offset))


# extracting all the tss and enhancer start/end data needed, using to calculate distance

# get tss values
eg_pairs_tss = pd.merge(enhancer_gene_pairs_df, rna_relevant[['gene_ids', 'tss']], on='gene_ids', how='inner')
# get enhancer start and end values (be aware of repeats)
h3k27ac_chmm_deduped = h3k27ac_chmm.drop_duplicates(subset='name')
eg_pairs_all = pd.merge(eg_pairs_tss, h3k27ac_chmm_deduped[['name', 'Start', 'End']], on='name', how='inner')
# calculate distance from tss to start
eg_pairs_all['tss_start'] = abs(eg_pairs_all['tss'] - eg_pairs_all['Start'])
# calculate distance from tss to end
eg_pairs_all['tss_end'] = abs(eg_pairs_all['tss'] - eg_pairs_all['End'])
# confirm they are below 100,000
eg_pairs_all = eg_pairs_all[(eg_pairs_all['tss_start'] <= 100000) & (eg_pairs_all['tss_end'] <= 100000)]
# pick which one is smaller
eg_pairs_all['dist'] = eg_pairs_all.apply(lambda row: min(row['tss_start'], row['tss_end']), axis=1)
# select columns of interest
gen_enh_dist = eg_pairs_all[['gene_ids', 'enhancer', 'name', 'dist', 'Start', 'End']]
gen_enh_dist


# use abc_distance function with calculated distance, add data to df
gen_enh_dist['abc_dist'] = gen_enh_dist['dist'].apply(abc_distance)
gen_enh_dist



# Adjusting hic_power based on different activities

# H3K27AC ACTIVITY
h3k27ac_act = h3k27ac["BSS00385"].X.mean()
gen_enh_dist['abc_h3k27ac'] = gen_enh_dist['abc_dist'] * h3k27ac_act

# DNASE ACTIVITY
dnase_act = dnase["BSS00385"].X.mean()
gen_enh_dist['abc_dnase'] = gen_enh_dist['abc_dist'] * dnase_act

# COMBINED H3K27AC & DNASE ACTIVITY
import math
comb_act = math.sqrt(h3k27ac_act * dnase_act)
gen_enh_dist['abc_comb_sqrt'] = gen_enh_dist['abc_dist'] * comb_act