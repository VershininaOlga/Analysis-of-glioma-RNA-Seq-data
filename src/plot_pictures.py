import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
from lifelines import KaplanMeierFitter
import scipy.stats as scist
from sklearn.linear_model import LinearRegression
import predictive_model as predm


def plot_venn_diagram_deg(n_deg_ps, n_deg_mtx, n_deg_common, labels, picture_name):
    vd = venn2(subsets=[n_deg_ps, n_deg_mtx, n_deg_common], set_labels=labels, set_colors=("tomato","deepskyblue"))
    
    for text in vd.set_labels:
        text.set_fontsize(18)
    for text in vd.subset_labels:
        text.set_fontsize(14)
    vd.get_patch_by_id('11').set_color((0.77, 0.49, 0.49))
    vd.get_patch_by_id('11').set_alpha(1)
    vd.get_patch_by_id('10').set_alpha(0.7)
    vd.get_patch_by_id('01').set_alpha(0.7)
    
    plt.savefig(f'{picture_name}.png', dpi=600, bbox_inches='tight', facecolor='w')
    
    
def plot_hist_deg(ps_fc, mtx_fc, picture_name):
    fig, ax = plt.subplots(figsize=(12, 5))
    
    fontsize = 20
    ticksize=15
    bins = [-22, -10, -5, -2, 2, 5, 10, 22, 100, 1000, 5000, 10000, 22325]

    hist_ps, bin_edges_ps = np.histogram(ps_fc, bins=bins)
    hist_mtx, bin_edges_mxt = np.histogram(mtx_fc, bins=bins)
    ax.bar(range(0, len(hist_ps)), hist_ps, width=1, color='deepskyblue', alpha=0.7, label='DC+PS-PDT')
    ax.bar(range(0, len(hist_mtx)), hist_mtx, width=1, color='tomato', alpha=0.7, label='DC+MTX')
    
    xlabels = list(map(str, bin_edges_mxt))
    ax.set_xticks(ticks=np.arange(-0.5, len(hist_mtx)-0.5, 1))
    ax.set_xticklabels(xlabels[:-1], fontsize=ticksize)
    ax.yaxis.set_tick_params(labelsize=ticksize)
    ax.set_xlabel('Fold Change', size=fontsize)
    ax.set_ylabel('Number of genes', size=fontsize)
    ax.legend(loc = 'upper right', fontsize=fontsize)

    fig.savefig(f'{picture_name}.png', dpi=600, bbox_inches='tight', facecolor='w')
    

def star_for_pval(pvaladj):
    if pvaladj <= 0.0001:
        return '****'
    elif pvaladj <= 0.001:
        return '***'
    elif pvaladj <= 0.01:
        return '**'
    elif pvaladj <= 0.05:
        return '*'
    else:
        return 'ns'
    

def rename_and_replace_groups(data):
    for k in range(0, len(data)):
        if 'PS' in data['Group'][k]:
            data.loc[k, 'Group'] = 'DC+PS-PDT'
        elif 'MTX' in data['Group'][k]:
            data.loc[k, 'Group'] = 'DC+MTX'
        else:
            data.loc[k, 'Group'] = 'Control'
    data = pd.concat([data[data['Group'] == 'DC+PS-PDT'], data[data['Group'] == 'Control'], data[data['Group'] == 'DC+MTX']],
                     axis=0, ignore_index=True)
    return data


def plot_boxplots_th17cells(genes, expr_levels, DSeq2_results_psvsdc, DSeq2_results_mtxvsdc, picture_name):
    expr_levels = rename_and_replace_groups(expr_levels)
    
    sign_ps = [star_for_pval(DSeq2_results_psvsdc[DSeq2_results_psvsdc['gene_name'] == gene]['padj'].values[0]) for gene in genes]
    sign_mtx = [star_for_pval(DSeq2_results_mtxvsdc[DSeq2_results_mtxvsdc['gene_name'] == gene]['padj'].values[0]) for gene in genes]

    fig, ax = plt.subplots(1, len(genes), figsize=(12, 4.5))
    
    colors = ['deepskyblue', 'lightgreen', 'tomato']
    sns.set_palette(sns.color_palette(colors))
    
    for k in range(0, len(genes)):
        ax[k].set(title = genes[k])    
        ax[k] = sns.boxplot(x='Group', y=genes[k], data=expr_levels,  ax=ax[k], showfliers=False)
        ax[k].axes.set_title(genes[k], fontsize=14, fontweight=540)
        ax[k].set_xlabel("Group", fontsize=14)
        ax[k].set(ylim=(0, 9))
        ax[k].set_ylabel(r"$log_{2} TPM$ gene expression", fontsize=14)
        ax[k].tick_params(labelsize=11)
        x1, x2 = 0, 1
        y, h, col = max(expr_levels[genes[k]]) - 0.5, 0.2, "k"
        ax[k].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        ax[k].text((x1+x2)*.5, y+h, sign_ps[k], ha="center", va="bottom", color=col, fontsize=15)
        x1, x2 = 1, 2
        y, h, col = max(expr_levels[genes[k]]) + 0.3, 0.2, "k"
        ax[k].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        ax[k].text((x1+x2)*.5, y+h, sign_mtx[k], ha="center", va="bottom", color=col, fontsize=15)

    fig.subplots_adjust(wspace=0.4)
    
    fig.savefig(f'{picture_name}.png', dpi=600, bbox_inches='tight', facecolor='w')
    

def overall_surv_for_metagenes(metagene_expr_low, metagene_expr_high, cell_type, picture_name):
    
    logrank_pval = predm.logrank_test_low_high(metagene_expr_low, metagene_expr_high)
    
    MS_low = int(np.median(metagene_expr_low['duration']))
    MS_high = int(np.median(metagene_expr_high['duration']))
    delta_MS = ((MS_high - MS_low) / MS_high) * 100
    delta_MS_str = f'+ {int(round(delta_MS))}' if delta_MS > 0 else f'{int(round(delta_MS))}'
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    kmf = KaplanMeierFitter()
    kmf.fit(metagene_expr_low['duration'], event_observed=metagene_expr_low['observed'], label='Low', alpha=1)
    kmf.plot(ax=ax, c='b', lw=2.5)
    kmf.fit(metagene_expr_high['duration'], event_observed=metagene_expr_high['observed'], label='High', alpha=1)
    kmf.plot(ax=ax, c='r', lw=2.5)

    ax.set_title(fr'$\bf{cell_type}$', size=20)
    ax.set_xlabel('Time (days)', size=18)
    ax.set_ylabel('Overall survival probability', size=18)
    ax.set_xticks([0, 1000, 2000, 3000, 4000, 5000, 6000])
    ax.tick_params(labelsize=15)
    ax.set_xlim([0, 6500])
    ax.set_ylim([0, 1])
    ax.axvline(x=730, c='k', ls='--', lw=2) 
    ax.axvline(x=1825, c='k', ls='--', lw=2) 

    ax.legend(bbox_to_anchor=(0.75, 0.9, 0, 0), edgecolor='w', fontsize=15)
    if logrank_pval <= 10**(-2):
        ax.text(x=4900, y=0.82, s='P = %.2e' % (logrank_pval), fontsize=15)
    else:
        ax.text(x=4900, y=0.82, s='P = %.2f' % (logrank_pval), fontsize=15)

    ax.text(x=3500, y=0.7, s=r'$MS^{High}\;=\;$' + fr'${MS_high}$', fontsize=15, c='r')
    ax.text(x=3550, y=0.65, s=r'$MS^{Low}\;=\;$' + fr'${MS_low}$', fontsize=15, c='b')
    ax.text(x=3540, y=0.60, s=fr'$\% \Delta MS\;=\;{delta_MS_str}\%$', fontsize=15, c='k')
    
    ax.text(x=800, y=0.95, s='2 y.s.', fontsize=15)
    ax.text(x=1900, y=0.95, s='5 y.s.', fontsize=15)
    
    fig.savefig(f'{picture_name}.png', dpi=600, transparent=True, facecolor='w')
    

def create_pairs(cells):
        pairs = []
        for c in cells:
            pairs.append(((c, 'Alive'), (c, 'Dead')))
        return pairs
    
    
def plot_boxplot_deconv_res(cells, mannwhitneyu_res, prepare_deconv_res, picture_name):
    fig, ax = plt.subplots(figsize=(20, 7))
    sns.boxplot(y='Fraction', x='Cell type', data=prepare_deconv_res, hue='patient_group', ax=ax, palette={'Alive': 'g', 'Dead': 'r'},
               flierprops={'marker': 'o', 'markerfacecolor': 'k'})
    
    mannwhitneyu_res['sign'] = [star_for_pval(pval) for pval in mannwhitneyu_res['padj'].values]
    for k in range(0, len(cells)):
        x1, x2 = k-0.2, k+0.2
        y, h, col = prepare_deconv_res[prepare_deconv_res['Cell type'] == cells[k]]['Fraction'].max() + 0.05, 0.02, 'k'
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        ax.text((x1+x2)*.5, y+h, mannwhitneyu_res[mannwhitneyu_res['cell'] == cells[k]]['sign'].values[0], 
                ha='center', va='bottom', color=col, fontsize=16)

    ax.set_xlabel('')
    ax.set_xticklabels(cells, fontsize=17)
    ax.set_ylabel('Fraction', fontsize=17)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8], fontsize=17)
    ax.set_ylim([-0.04,  0.91])

    ax.legend(fontsize=17, loc='upper left')

    fig.savefig(f'{picture_name}.png', dpi=600, bbox_inches='tight', facecolor='w') 

    
def plot_surv_roc_curves_for_predictive_model(data_low, data_high, roc_1y, roc_3y, roc_5y, aucs, picture_name):
    logrank_pval = predm.logrank_test_low_high(data_low, data_high)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw = {'wspace': 0.25})

    kmf = KaplanMeierFitter()
    kmf.fit(data_low['duration'], event_observed=data_low['observed'], label='low-risk', alpha=1)
    kmf.plot(ax=ax1, c='steelblue', lw=2.5)
    kmf.fit(data_high['duration'], event_observed=data_high['observed'], label='high-risk', alpha=1)
    kmf.plot(ax=ax1, c='crimson', lw=2.5)

    ax1.set_xlabel('Time (Days)', fontsize=15)
    ax1.set_ylabel('Overall survival probability', fontsize=15)
    ax1.set_xlim([0, 365*5])
    ax1.set_ylim([0, 1])
    ax1.set_xticks([0, 300, 600, 900, 1200, 1500, 1800])
    ax1.set_xticklabels([0, 300, 600, 900, 1200, 1500, 1800], fontsize=14)
    ax1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax1.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=14)
    ax1.legend(fontsize=14, loc='lower left')
    if logrank_pval <= 10**(-4):
        ax1.text(x=650, y=0.07, s='P = %.2e' % (logrank_pval), fontsize=15)
    else:
        ax1.text(x=650, y=0.07, s='P = %.2f' % (logrank_pval), fontsize=15)


    ax2.plot([0, 1], [0, 1], 'k-')
    ax2.plot(roc_1y['FPR'], roc_1y['TPR'], '-', c='deeppink', 
             label=f"1 year, AUC = {round(aucs[aucs['year'] == 1]['AUC'].values[0], 2)}", lw=2.5)
    ax2.plot(roc_3y['FPR'], roc_3y['TPR'], '-', c='darkslateblue', 
             label=f"3 year, AUC = {round(aucs[aucs['year'] == 3]['AUC'].values[0], 2)}", lw=2.5)
    ax2.plot(roc_5y['FPR'], roc_5y['TPR'], '-', c='forestgreen', 
             label=f"5 year, AUC = {round(aucs[aucs['year'] == 5]['AUC'].values[0], 2)}", lw=2.5)
    ax2.set_xlabel('1 - Specificity', fontsize=15)
    ax2.set_ylabel('Sensitivity', fontsize=15)
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1])
    ax2.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_xticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=14)
    ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=14)
    ax2.legend(fontsize=14)
    
    fig.savefig(f'{picture_name}.png', dpi=600, bbox_inches='tight', facecolor='w')

    
def plot_log2fc_genes_from_predictive_model(sign_genes, DSeq2_results_AvsC, DSeq2_results_DvsC, 
                                            DSeq2_results_psvsdc, DSeq2_results_mtxvsdc, picture_name):
    
    def point_for_gene(x, gene, Dseq2res, label, color, ax, fontsize):
        log2fc = Dseq2res[Dseq2res['gene_name'] == gene]['log2FoldChange']
        ax.plot(x, log2fc, 'o', c=color)
        if gene == 'GALNT3' and label == 'DC+PS-PDT':
            ax.text(x + 0.05, log2fc + 0.02, label, fontsize=fontsize, c=color)
        else:
            ax.text(x + 0.05, log2fc + 0.05, label, fontsize=fontsize, c=color)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    fontsize = 14
    ax.axhline(0.0, c='limegreen', ls='-.')
    ax.text(1.5, 0.05, 'Control', c='limegreen', fontsize=fontsize)
    ax.set_ylim([-2.0, 0.65])
    ax.set_xlim([-0.2, len(sign_genes)])
    
    xlabels = []
    for k in range(0, len(sign_genes)):
        xlabels.append(f'{sign_genes[k]}\n{sign_genes[k].capitalize()}')
        point_for_gene(k, sign_genes[k], DSeq2_results_AvsC, 'Alive', 'orange', ax, fontsize)
        point_for_gene(k, sign_genes[k], DSeq2_results_DvsC, 'Dead', 'black', ax, fontsize)
        point_for_gene(k, sign_genes[k], DSeq2_results_psvsdc, 'DC+PS-PDT', 'deepskyblue', ax, fontsize)
        point_for_gene(k, sign_genes[k], DSeq2_results_mtxvsdc, 'DC+MTX', 'tomato', ax, fontsize)
    ax.set_xticks(range(0, len(sign_genes)))
    ax.set_xticklabels(xlabels, fontsize=fontsize)
    ax.yaxis.set_tick_params(labelsize=fontsize)
    ax.set_xlabel(r'$Gene$', fontsize=fontsize+2)
    ax.set_ylabel(r'$log_{2}(Fold~Change)$', fontsize=fontsize+2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.savefig(f'{picture_name}.png', dpi=600, transparent=True, facecolor='w')

    
def plot_risk_score_expr_level_for_genes_in_predictive_model(sign_genes, data_low, data_high, rs_median, picture_name):
    
    def point_colors(status):
        colors = np.array(status)
        colors[colors == 'Alive'] = 'orange'
        colors[colors == 'Dead'] = 'black'
        return colors  
        
    n_low = len(data_low)
    n_high = len(data_high)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 9), sharex=True, 
                                        gridspec_kw = {'hspace': 0.5, 'height_ratios': [1, 1, 1]})

    ax1.scatter(range(1, n_low + 1), data_low['risk_score'], c='steelblue', label='low-risk', s=12)
    ax1.scatter(range(n_low + 1, n_low + n_high + 1), data_high['risk_score'], c='crimson', label='high-risk', s=12)
    ax1.set_ylabel('Risk Score', fontsize=14, labelpad=46)
    ax1.legend(fontsize=12)
    ax1.set_ylim([0, 6.1])
    ax1.set_yticks([0, 2, 4, 6])
    ax1.set_yticklabels([0, 2, 4, 6], fontsize=14)
    ax1.plot([1, n_low+1], [rs_median, rs_median], c='k', ls='--', lw=2)
    ax1.plot([n_low+1, n_low+1], [np.min([np.min(data_low['risk_score']), np.min(data_high['risk_score'])]), rs_median], 
             c='k', ls='--', lw=2)

    ss = ax2.scatter(range(1, n_low + 1), data_low['duration'], c=point_colors(data_low['condition'].values), s=12)
    ax2.scatter(range(n_low + 1, n_low + n_high + 1), data_high['duration'], c=point_colors(data_high['condition'].values), s=12)
    leg = ax2.legend(['Alive', 'Dead'], fontsize=12, loc='upper right')
    m = 0
    for marker in leg.legendHandles:
        if m == 0:
            marker.set_color('orange')
        else:
            marker.set_color('black')
        m += 1

    ax2.plot([n_low+1, n_low+1], [0, 6500], c='k', ls='--', lw=2, zorder=10)
    ax2.set_ylabel('Time (Days)', fontsize=14, labelpad=19)
    ax2.set_yticks([0, 2000, 4000, 6000])
    ax2.set_yticklabels([0, 2000, 4000, 6000], fontsize=14)
                
    data = data_low[sign_genes].append(data_high[sign_genes])
    data = data.astype('float64')
    zscore = (data - np.mean(data, axis=0)) / np.std(data, axis=0)
    
    im = ax3.pcolormesh(zscore.T, cmap='bwr', vmin=-4.5, vmax=4.5)
    ax3.set_xlim([0, n_low+n_high])
    ax3.set_ylabel('Genes', fontsize=14, labelpad=2)
    ax3.set_yticks([0.5, 1.5, 2.5, 3.5])
    ax3.set_yticklabels(sign_genes, fontsize=13)
    ax3.set_xticks(range(0, n_low+n_high+5, 40))
    ax3.set_xticklabels(range(1, n_low+n_high+6, 40), fontsize=14)
    ax3.set_xlabel('Patients', fontsize=14)
    cbaxes = fig.add_axes([0.2, 0.03, 0.2, 0.02])
    cb = fig.colorbar(im, ax=ax3, cax=cbaxes, orientation='horizontal')
    cb.ax.set_title('Row Z-score', fontsize=13, x=-0.3, y=-0.3)
    cb.ax.tick_params(labelsize=12)

    plt.savefig(f'{picture_name}.png', dpi=600, bbox_inches='tight', facecolor='w')


def overall_surv_for_cells(cells_fraction_low, cells_fraction_high, cell_type, picture_name):
    
    logrank_pval = predm.logrank_test_low_high(cells_fraction_low, cells_fraction_high)
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    kmf = KaplanMeierFitter()
    kmf.fit(cells_fraction_low['duration'], event_observed=cells_fraction_low['observed'], label='Low', alpha=1)
    kmf.plot(ax=ax, c='g', lw=2.5)
    kmf.fit(cells_fraction_high['duration'], event_observed=cells_fraction_high['observed'], label='High', alpha=1)
    kmf.plot(ax=ax, c='r', lw=2.5)

    ax.set_title(f'{cell_type}', size=20)
    ax.set_xlabel('Time (days)', size=18)
    ax.set_ylabel('Overall survival probability', size=18)
    ax.set_xticks([0, 1000, 2000, 3000, 4000, 5000, 6000])
    ax.tick_params(labelsize=15)
    ax.set_xlim([0, 6500])
    ax.set_ylim([0, 1])
    ax.axvline(x=730, c='k', ls='--', lw=2) 
    ax.axvline(x=1825, c='k', ls='--', lw=2) 

    ax.legend(bbox_to_anchor=(0.75, 0.9, 0, 0), edgecolor='w', fontsize=15)
    if logrank_pval <= 10**(-2):
        ax.text(x=4900, y=0.82, s='P = %.2e' % (logrank_pval), fontsize=15)
    else:
        ax.text(x=4900, y=0.82, s='P = %.2f' % (logrank_pval), fontsize=15)

    ax.text(x=800, y=0.95, s='2 y.s.', fontsize=15)
    ax.text(x=1900, y=0.95, s='5 y.s.', fontsize=15)
    
    fig.savefig(f'{picture_name}.png', dpi=600, transparent=True, facecolor='w')
    
    
def correlation_for_cells(data_ct, cell_type, picture_name):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    
    r, pval = scist.pearsonr(data_ct['risk_score'], data_ct[cell_type])
        
    reg = LinearRegression().fit(data_ct[['risk_score']], data_ct[cell_type])
    ox = np.sort(data_ct['risk_score'])
    fit = reg.coef_ * ox + reg.intercept_

    plt.plot(data_ct['risk_score'], data_ct[cell_type], 'ok', ms=5)
    plt.plot(ox, fit, 'r-', lw=3) 

    ax.set_xlabel('Risk Score', size=18)
    ax.set_ylabel(f'{cell_type}', size=18)
    ax.tick_params(labelsize=15)
    ax.set_xlim([0, 6])
    
    if cell_type == 'B cells' or cell_type == 'NK cells':
        ax.set_ylim([-0.001, 0.02])
    elif cell_type == 'CAFs' or cell_type == 'T cells CD4+':
        ax.set_ylim([-0.01, 0.2])
    elif cell_type == 'T cells CD8+':
        ax.set_ylim([-0.01, 0.1])
    elif cell_type == 'Endothelial cells' or cell_type == 'Macrophages':
        ax.set_ylim([-0.01, 0.35])

    if pval < 0.05:
        ax.annotate(text='r=%.4f\np=%.1e' % (r, pval), xy=(0.8, 0.9), xycoords='axes fraction', c='r', fontsize=15)
    else:
        ax.annotate(text='r=%.4f\np=%.1e' % (r, pval), xy=(0.8, 0.9), xycoords='axes fraction', fontsize=15)
        
    
    fig.savefig(f'{picture_name}.png', dpi=600, transparent=True, facecolor='w')