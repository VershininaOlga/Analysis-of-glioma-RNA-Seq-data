import numpy as np
import pandas as pd
import scipy.stats as scist


def correlation(subdata):
    n_genes, n_obj = subdata.shape
    corr_matrix = np.zeros((n_genes, n_genes))
    for i in range(0, n_genes):
        for j in range(0, n_genes):
            corr_coef, _ = scist.pearsonr(subdata.loc[i][1:], subdata.loc[j][1:])
            corr_matrix[i][j] = corr_coef
    corr_matrix = pd.DataFrame(corr_matrix, index=subdata['gene_name'].values, columns=subdata['gene_name'].values)
    return corr_matrix
