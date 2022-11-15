import pandas as pd
import numpy as np
import scipy.stats as scist
from statsmodels.stats.multitest import fdrcorrection


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


def mannwhitneyu_deconv_res(cells, deconvGroup1, deconvGroup2, alternative):
    pval = []
    for c in cells:
        _, p = scist.mannwhitneyu(deconvGroup1[c],  deconvGroup2[c], alternative=alternative)
        pval.append(p)
    mwt = pd.DataFrame({'cell': cells, 'pval': pval})
    _, padj = fdrcorrection(pval)
    mwt['padj'] = padj
    
    sign = []
    for padj in mwt['padj'].values:
        sign.append(star_for_pval(padj))
    mwt['sign'] = sign
    
    return mwt


def check_infiltrating(deconv_res):
    n = len(deconv_res)
    no_infiltrating = []
    for c in deconv_res.columns[1:-1]:
        if np.sum(deconv_res[c] < 10**(-8)) > (n * 75) / 100:
            no_infiltrating.append(c)
    return no_infiltrating


def prepare_deconv_res(deconv_res):
    prepare_data = pd.DataFrame()
    for c in deconv_res.columns[1:-1]:
        df = pd.DataFrame({'Fraction': deconv_res[c]})
        df['Cell type'] = c  
        df['patient_group'] = deconv_res['patient_group']
        df['sample'] = deconv_res['sample']
        prepare_data = prepare_data.append(df, ignore_index=True)
    return prepare_data
