from setup import *

def correlation(A, B, wt=None):
    "Fast pearson correlation for one matrix against another"
    if wt is None:
        wt = np.ones(A.shape[1])
    wt = np.ravel(wt) / np.sum(wt)
    Am = (A.T - np.dot(A, wt)).T
    Bm = (B.T - np.dot(B, wt)).T
    Am_sum_sq = np.dot(Am * Am, wt)
    Bm_sum_sq = np.dot(Bm * Bm, wt)
    numer = np.dot(Am * Bm, wt)
    denom = np.sqrt(Am_sum_sq * Bm_sum_sq)
    cor = np.divide(numer, denom, out=np.zeros_like(numer), where=denom != 0)
    return cor

def compute_correlation(links:pd.DataFrame, rna:anndata.AnnData, dna:anndata.AnnData,
                        batch_size=1000, enh_col="name", gene_col="gene_ids"):
    comm_samples = list(set(rna.obs_names.values) & set(h3k27ac.obs_names.values))
    cor = np.zeros(links.shape[0])
    for start in np.arange(0, links.shape[0], batch_size):
        end = min(start + batch_size, links.shape[0])
        enhancers = links[enh_col].values[start:end]
        genes = links[gene_col].values[start:end]
        genes = list(genes)
        rna_data = rna[comm_samples, genes].X.T
        dna_data = dna[comm_samples, enhancers].X.T
        cor[start:end] = correlation(rna_data, dna_data)
    return cor

def abc_distance(distance, hic_power=0.87):
    scale = -4.80 + 11.63 * hic_power
    offset = np.clip(np.abs(distance), 5000, np.inf)
    return np.exp(scale + -1 * hic_power * np.log1p(offset))




# fix rna complications
rna_df = rna.var
rna.var_names = rna_df['gene_ids'].values

# H3K27ac * Cor(H3K27ac, RNA)^2
corr_1 = compute_correlation(gen_enh_dist, rna, dnase)
corr_1_h3k27ac = corr_1 * h3k27ac_act

# DNase-seq * Cor(DNase-seq, RNA)^2
corr_2 = compute_correlation(gen_enh_dist, rna, dnase)
corr_2_dnase = corr_1 * dnase_act

# sqrt(H3K27ac * DNase-seq) * Cor(sqrt(H3K27ac*DNase-seq), RNA)^2
corr_3 = compute_correlation(gen_enh_dist, rna, dnase)
corr_3_sqrt = corr_3 * comb_act


corr_data = pd.DataFrame({'corr_1': corr_1_h3k27ac, 'corr_2': corr_2_dnase, 'corr_3': corr_3_sqrt})
pair_corr = pd.concat([gen_enh_dist, corr_data], axis=1)
pair_corr = pair_corr.dropna(subset=['gene_ids'])
pair_corr