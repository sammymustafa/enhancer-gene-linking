# Enhancer-Gene Linking

This project focuses on enhancer-gene linking, employing Hidden Markov Models and the Activity-by-Contact (ABC) model to unravel the dynamics of gene regulation. Key techniques such as enhancer-gene pair identification within 100Kbp of genes, ABC scoring metrics, and correlating enhancers with gene expressions are utilized. These methods are instrumental in advancing our understanding of the regulatory roles enhancers play in gene expression.

### 1. Finding Enhancer-Gene Pairs
* Objective: Identify enhancer-gene pairs within 100,000 base pairs (100Kbp) from the Transcription Start Site (TSS) of a gene.
* Methodology: Utilized AnnData objects to represent samples and features (peaks or genes), computing enhancer-gene pairs based on "interval" fields in RNA and chromatin/DNase data.
* Specific Task: Focused on finding peaks overlapping Enhancer states in BSS00369 (Brain frontal cortex) data​​.

### 2. Activity-by-Contact Model with ABC Distance and Three Scoring Metrics
* ABC Model: Developed an Activity-by-Contact (ABC) model, which employs a power law formula based on the distance between elements and multiplies it by specific activities, such as the geometric mean of H3K27ac and DNase-seq activities.
* Distance Calculation: Implemented an abc_distance function that uses a power law approach for determining distances between enhancer-promoter pairs​​.
* Scoring Metrics: Calculated three ABC-like scores using the activities of H3K27ac, DNase-seq, and the square root of the product of H3K27ac and DNase-seq activities​​.
* Comparative Analysis: Compared the four predictions (including the ABC distance) with GTEx eQTLs to evaluate the predictive accuracy of the enhancer-gene linking​​.

### 3. Correlating Enhancers to Gene Expression with Three Methods
* Method 1: Utilizing H3K27ac activity data, the model multiplies the ABC distance by the mean H3K27ac activity to generate one set of scores​​.
* Method 2: For DNase-seq activity, a similar approach is applied, multiplying the ABC distance by the DNase-seq activity mean​​.
* Method 3: Combines H3K27ac and DNase-seq activities, using the square root of their product as a factor in generating scores​​.

### 4. Estimating Gene Expression from Chromatin Data
* Integration of Data: The methods involve integrating chromatin data (H3K27ac and DNase-seq) with the computed ABC distances to predict enhancer-gene interactions.
* Predictive Accuracy Assessment: The accuracy of these methods is assessed by their ability to correlate enhancer regions with gene expressions, particularly in comparison to established eQTL data.
