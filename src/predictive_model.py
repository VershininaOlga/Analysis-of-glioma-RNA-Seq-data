import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
import differential_expression_results_analysis as diffexpr


def identification_genes_for_predictive_model(DSeq2_results_DvsA, DSeq2_results_AvsC, DSeq2_results_DvsC, 
                                              DSeq2_results_psvsdc, DSeq2_results_mtxvsdc):
    DSeq2_results_psvsdc['gene_name'] = DSeq2_results_psvsdc['gene_name'].apply(lambda x: x.upper())
    DSeq2_results_mtxvsdc['gene_name'] = DSeq2_results_mtxvsdc['gene_name'].apply(lambda x: x.upper())
    
    # 1 - Dead vs Alive, significant DEG with |FC|>2
    deg_DvsA_fc2 = diffexpr.deg_identification(DSeq2_results_DvsA, pval_thr=0.05, baseMean_thr=50, log2FC_thr=1)
    deg_DvsA_fc2_up, deg_DvsA_fc2_down = diffexpr.up_down_deg(deg_DvsA_fc2)
    print(f'DEG (|FC|>2) for DEAD vs ALIVE:\n\tall = {len(deg_DvsA_fc2)}, up = {len(deg_DvsA_fc2_up)}, down = {len(deg_DvsA_fc2_down)}')
    
    # 2 - Alive vs Control, significant DEG and common with #1
    deg_AvsC = diffexpr.deg_identification(DSeq2_results_AvsC, pval_thr=0.05, baseMean_thr=50)
    deg_AvsC = deg_AvsC[deg_AvsC['gene_name'].isin(deg_DvsA_fc2['gene_name'].values)]
    deg_AvsC_up, deg_AvsC_down = diffexpr.up_down_deg(deg_AvsC)
    print(f'DEG for ALIVE vs CONTROL:\n\tall = {len(deg_AvsC)}, up = {len(deg_AvsC_up)}, down = {len(deg_AvsC_down)}')
    
    # 3 - DEG for DCPS and DCMTX
    deg_ps = diffexpr.deg_identification(DSeq2_results_psvsdc, pval_thr=0.05, baseMean_thr=100)
    deg_ps_up, deg_ps_down = diffexpr.up_down_deg(deg_ps)
    deg_mtx = diffexpr.deg_identification(DSeq2_results_mtxvsdc, pval_thr=0.05, baseMean_thr=100)
    deg_mtx_up, deg_mtx_down = diffexpr.up_down_deg(deg_mtx)
    common_deg_ps_mtx = diffexpr.common_deg('DCPS', deg_ps, 'DCMTX', deg_mtx)
    common_deg_ps_mtx_up = diffexpr.common_deg('DCPS', deg_ps_up, 'DCMTX', deg_mtx_up)
    common_deg_ps_mtx_down = diffexpr.common_deg('DCPS', deg_ps_down, 'DCMTX', deg_mtx_down)
    print(f'common DEG for DCPS and DCMTX:\n\tall = {len(common_deg_ps_mtx)}, \n\tsimultaneously up = {len(common_deg_ps_mtx_up)}, \n\tsimultaneously down = {len(common_deg_ps_mtx_down)}, \n\tbehave differently = {len(common_deg_ps_mtx)-len(common_deg_ps_mtx_up)-len(common_deg_ps_mtx_down)}')

    # 4
    unidirectional_up = deg_AvsC_up[deg_AvsC_up['gene_name'].isin(common_deg_ps_mtx_up['gene_name'])]
    unidirectional_up = unidirectional_up.merge(DSeq2_results_DvsC, left_on='gene_name', right_on='gene_name', suffixes=('_A', '_D'))
    unidirectional_up.reset_index(drop=True, inplace=True)

    unidirectional_up_genes = []
    for k in range(0, len(unidirectional_up)):
        if (unidirectional_up['padj_D'][k] > 0.05) or (abs(unidirectional_up['log2FoldChange_D'][k]) < abs(unidirectional_up['log2FoldChange_A'][k])):
            unidirectional_up_genes.append(unidirectional_up['gene_name'][k])

     # 5
    unidirectional_down = deg_AvsC_down[deg_AvsC_down['gene_name'].isin(common_deg_ps_mtx_down['gene_name'])]
    unidirectional_down = unidirectional_down.merge(DSeq2_results_DvsC, left_on='gene_name', right_on='gene_name', suffixes=('_A', '_D'))
    unidirectional_down.reset_index(drop=True, inplace=True)

    unidirectional_down_genes = []
    for k in range(0, len(unidirectional_down)):
        if (unidirectional_down['padj_D'][k] > 0.05) or (abs(unidirectional_down['log2FoldChange_D'][k]) < abs(unidirectional_down['log2FoldChange_A'][k])):
            unidirectional_down_genes.append(unidirectional_down['gene_name'][k])
    
    return unidirectional_up_genes + unidirectional_down_genes


def take_genes_and_transposition(data, genes):
    data = data[data['gene_name'].isin(genes)]
    data_trans = pd.DataFrame({'sample': data.columns[1:]})
    for gene in genes:
        data_trans[gene] = data[data['gene_name'] == gene].values[0][1:]
    return data_trans


def univatiate(genes, data, attributes, pval_thr=0.05):
    sign_genes = []
    for gene in genes:
        cph = CoxPHFitter()
        cph.fit(data[[gene] + ['duration', 'observed']], duration_col='duration', event_col='observed')
        results = cph.summary
        if results['p'][gene] < pval_thr:
            sign_genes.append(gene)
    return sign_genes


def multivariate(genes, data, attributes):
    cph = CoxPHFitter()
    cph.fit(data[genes + ['duration', 'observed']], duration_col='duration', event_col='observed')
    return cph.concordance_index_, cph.params_


def risk_score(data, coefs):
    rs = []
    genes = coefs.index.values
    for k in range(0, len(data)):
        rs.append(np.sum(coefs.values * data[genes].iloc[k].values))
    return rs


def logrank_test_low_high(data_low, data_high):
    results = logrank_test(data_low['duration'], data_high['duration'], 
                           event_observed_A=data_low['observed'], event_observed_B=data_high['observed'])
    return results.summary['p'][0]
    