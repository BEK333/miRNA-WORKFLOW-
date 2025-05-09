import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import sys

# Usage: python deseq2.py counts.csv meta.tsv output.csv
counts = pd.read_csv(sys.argv[1], index_col=0)
meta = pd.read_csv(sys.argv[2], sep='\t')

groups = meta['condition'].values
unique_groups = np.unique(groups)
assert len(unique_groups) == 2, "Only two-group comparison supported."
mask1 = groups == unique_groups[0]
mask2 = groups == unique_groups[1]

results = []
for gene in counts.index:
    vals1 = counts.loc[gene, mask1]
    vals2 = counts.loc[gene, mask2]
    stat, pval = ttest_ind(vals1, vals2, equal_var=False)
    log2fc = np.log2(np.mean(vals2) + 1) - np.log2(np.mean(vals1) + 1)
    results.append({'gene': gene, 'log2FoldChange': log2fc, 'pvalue': pval})
res_df = pd.DataFrame(results)
res_df['padj'] = res_df['pvalue'] * len(res_df) / (res_df['pvalue'].rank())  # Bonferroni
res_df.to_csv(sys.argv[3], index=False)
