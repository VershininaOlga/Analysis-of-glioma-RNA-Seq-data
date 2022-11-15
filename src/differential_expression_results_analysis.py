import pandas as pd


def deg_identification(DSeq2_results, pval_thr=None, baseMean_thr=None, log2FC_thr=None):
    deg = DSeq2_results
    if pval_thr is not None:
        deg = deg[deg['padj'] <= pval_thr]
    if baseMean_thr is not None:
        deg = deg[deg['baseMean'] > baseMean_thr]
    if log2FC_thr is not None:
        deg = deg[(deg['log2FoldChange'] >= log2FC_thr) | (deg['log2FoldChange'] <= -log2FC_thr)]
    deg.reset_index(drop=True, inplace=True)
    return deg


def up_down_deg(deg):
    deg_up = deg[deg['log2FoldChange'] > 0]
    deg_up.reset_index(drop=True, inplace=True)
    deg_down = deg[deg['log2FoldChange'] < 0]
    deg_down.reset_index(drop=True, inplace=True)
    return deg_up, deg_down


def common_deg(group1, deg_1, group2, deg_2):
    common = deg_1.merge(deg_2, left_on='gene_name', right_on='gene_name', suffixes=('_' + group1, '_' + group2))
    return common


def get_fold_changes(log2fc):
    fc = [2**x if x >=0 else -2**abs(x) for x in log2fc]
    return fc
