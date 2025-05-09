import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

df = pd.read_csv(sys.argv[1])
plt.figure(figsize=(8,6))
df['-log10_padj'] = -np.log10(df['padj'].replace(0, np.nan))
sig = (df['padj'] < 0.05) & (np.abs(df['log2FoldChange']) > 1)
sns.scatterplot(x='log2FoldChange', y='-log10_padj', data=df, hue=sig, palette={True:'red', False:'grey'}, legend=False)
plt.axhline(-np.log10(0.05), color='blue', linestyle='--')
plt.axvline(-1, color='blue', linestyle='--')
plt.axvline(1, color='blue', linestyle='--')
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10 Adjusted p-value')
plt.title('Volcano Plot')
plt.tight_layout()
plt.savefig(sys.argv[2])
