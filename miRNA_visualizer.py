import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Load your data
data = pd.read_csv('./visualize_grouped_miRNA_folded_pvalues.csv')

# Create a volcano plot assuming your data has 'Fold_Change' and 'p_Value' columns
plt.scatter(data['AF_to_AM_FoldChange_HK'], -np.log10(data['AF_AM_p_value_HK']), c='blue', alpha=0.5)
plt.scatter(data['AF_to_FC_FoldChange_HK'], -np.log10(data['AF_FC_p_value_HK']), c='red', alpha=0.5)
plt.scatter(data['AM_to_MC_FoldChange_HK'], -np.log10(data['AM_MC_p_value_HK']), c='orange', alpha=0.5)
plt.scatter(data['FC_to_MC_FoldChange_HK'], -np.log10(data['FC_MC_p_value_HK']), c='green', alpha=0.5)
plt.scatter(data['AF_to_AM_FoldChange_NEG'], -np.log10(data['AF_AM_p_value_NEG']), c='blue', alpha=0.5)
plt.scatter(data['AF_to_FC_FoldChange_NEG'], -np.log10(data['AF_FC_p_value_NEG']), c='red', alpha=0.5)
plt.scatter(data['AM_to_MC_FoldChange_NEG'], -np.log10(data['AM_MC_p_value_NEG']), c='orange', alpha=0.5)
plt.scatter(data['FC_to_MC_FoldChange_NEG'], -np.log10(data['FC_MC_p_value_NEG']), c='green', alpha=0.5)
plt.title('Volcano Plot')
plt.xlabel('Fold Change')
plt.ylabel('-log10(p-value)')
plt.axhline(y=-np.log10(0.05), color='red', linestyle='--')  # Threshold line for p-value
plt.show()


