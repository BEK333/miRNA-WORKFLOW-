import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

counts = pd.read_csv(sys.argv[1], index_col=0)
meta = pd.read_csv(sys.argv[2], sep='\t')
variances = counts.var(axis=1)
top20 = variances.sort_values(ascending=False).head(20).index
subset = counts.loc[top20]
subset_z = (subset - subset.mean(axis=1).values[:,None]) / subset.std(axis=1).values[:,None]
plt.figure(figsize=(10,8))
sns.heatmap(subset_z, cmap='vlag', yticklabels=True, xticklabels=meta['sample'])
plt.title('Top 20 Most Variable miRNAs')
plt.tight_layout()
plt.savefig(sys.argv[3])
