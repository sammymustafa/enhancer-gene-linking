from setup import *

# IDENTIFYING GENE-ENHANCER PAIRS (takes ~6 min to run, 3100 iterations)

# making RNA data into dataframe and filtering by histone chromosome 20-22 to make process quicker
rna_var_df = rna.var.copy()
chr_values = h3k27ac_chmm["Chromosome"].unique()
rna_relevant = rna_var_df[rna_var_df['seqnames'].isin(chr_values)]
# making empty dictionary for gene-enhancer pairs to be stored
enhancer_gene_pairs = {}

# iterate through gene data
for _, gene in tqdm(rna_relevant.iterrows(), desc='Processing'):
    # extract start and end locations
    tss = gene['tss']
    gene_id = gene['gene_ids']
    # identify enhancers within 100k bps from gene tss from both start and end, merge data
    start = h3k27ac_chmm.loc[abs(tss - h3k27ac_chmm['Start']) < 100000]
    end = h3k27ac_chmm.loc[abs(tss - h3k27ac_chmm['End']) < 100000]
    merged = pd.merge(start, end, on='name', how='outer')['name']
    time.sleep(0.1)
    for i in merged:
      enhancer_gene_pairs[i] = gene_id


# convert dictionary to df, data cleaning and identifying enhancer name
enhancer_chunk_pairs_df = pd.DataFrame(list(enhancer_gene_pairs.items()), columns=['Key', 'Value']).rename(columns={'Key': 'name', 'Value': 'gene_ids'})
enhancer_gene_pairs_df = enhancer_chunk_pairs_df.merge(h3k27ac_chmm[['name', 'state']], on='name', how='left').rename(columns={'value': 'gene_ids', 'state': 'enhancer'})
enhancer_gene_pairs_df = enhancer_gene_pairs_df[['gene_ids', 'enhancer', 'name']]
enhancer_gene_pairs_df