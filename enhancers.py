from setup import *

eqtl = pd.read_csv("linking/gtex_frontal_cortex_eQTL_hg19.bed.gz", sep="\t")

### Use the 18-state ChromHMM model
chmm = pd.read_csv("https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg19/CALLS/BSS00369_18_CALLS_segments.bed.gz", sep="\t", header=None).iloc[:, range(4)]
chmm.columns = ["Chromosome", "Start", "End", "state"]
chmm = chmm.loc[chmm["state"].str.startswith("Enh"), :]
print(chmm["state"].value_counts())
## TO IMPLEMENT: convert to PyRanges and intersect with h3k27ac.var
chmm_pyr = pr.PyRanges(chmm)
h3k27ac.var = h3k27ac["BSS00385"].var.rename(columns={'seqnames': 'Chromosome', 'start': 'Start', 'end': 'End'})
h3k27ac_pyr = pr.PyRanges(h3k27ac.var)
h3k27ac_chmm = chmm_pyr.join(h3k27ac_pyr, suffix = "_chmm").df
h3k27ac_chmm