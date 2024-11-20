# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 14:26:22 2024

@author: piercetf
"""

NORMPROT = 'Normalized_FG_Quantity'
N_COND = 4

import itertools
import time
import multiprocessing as multiproc
import logging
import platform



from data_prepper import DataPreparer
from nparc_model import NPARCModel
import load_monocyte_cetsa_data as load
import cetsa_paths
import t_infl

warnlogname = str(cetsa_paths.get_logging_path('cetsa_debug.log'))
logging.basicConfig(filename=warnlogname, level=logging.DEBUG)
logging.captureWarnings(True)

#from sklearn.preprocessing import PolynomialFeatures
import seaborn
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import pandas
import numpy as np
from scipy import stats

RNG_SEED = time.time_ns()

if platform.system() == 'Windows':
    N_PROCS = 20
elif platform.system() == 'Linux':
    import os
    proc_env_var = 'SLURM_CPUS_PER_TASK'
    if proc_env_var in os.environ:
        procs_str = os.environ[proc_env_var]
        N_PROCS = int(procs_str)
    else:
        N_PROCS = 2
    

def mse(predicteds, actuals):
    diffs = actuals - predicteds
    sqdif = diffs**2
    sse = sqdif.sum()
    return sse / len(predicteds)


def ymax_logit(y_act_max, y_values):
    return np.log(y_values) - np.log((y_act_max + 1) - y_values)


def ymax_logit_mse(y_predicted, y_actuals):
    y_act_max = max(y_actuals)
    y_pred_logit = ymax_logit(y_act_max, y_predicted)
    y_act_logit = ymax_logit(y_act_max, y_actuals)
    return mse(y_pred_logit, y_act_logit)


#def _permutation_test(model, datapart, ident, rng, n=50_000):
def _permutation_test(model, indata, outdata, treatdata, cat1, cat2, ident, rng, n=50_000):
    
    # need to denote tables where we don't have the raw data
    # required to even do permutation testing
    # and skip them to further cut down on execution time
    if indata[:,2].sum() == 0 or indata[:,3].sum() == 0:
        pvalue = 1.0
        interval = (0.0, 1.0)
        return pvalue, *interval, 'data-insufficient', -1, ident, cat1, cat2

    base_pred = model.predict(indata)
    neg_base_mse = -ymax_logit_mse(base_pred, outdata)
    permuted_indata = indata.copy()

    fast_better_count = 0
    fast_mmts = np.empty((500,))
    interact_permuted = np.empty_like(permuted_indata)
    
    for j in range(500):
        cat1_permuted = rng.permutation(indata[:,2])
        cat2_permuted = 1 - cat1_permuted
        permuted_indata[:,2] = cat1_permuted
        permuted_indata[:,3] = cat2_permuted

        interact_permuted[:,0:4] = permuted_indata[:,0:4]
        interact_permuted[:,4] = interact_permuted[:,1] * interact_permuted[:,2]
        interact_permuted[:,5] = interact_permuted[:,1] * interact_permuted[:,3]
        

        permuted_pred = model.predict(interact_permuted)
        neg_perm_mse = -ymax_logit_mse(permuted_pred, outdata)
        
        if neg_perm_mse >= neg_base_mse:
            fast_better_count += 1
        fast_mmts[j] = neg_perm_mse
    
    fast_pvalue = (fast_better_count + 1) / (501)
    if fast_pvalue > 0.1:
        p_variance = fast_pvalue * (1 - fast_pvalue) / 500
        interval = stats.norm.interval(0.95, fast_pvalue, p_variance)
        return fast_pvalue, *interval, 'ecdf-faststop', -1, ident, cat1, cat2
    
    better_count = 0
    mmts = np.empty((n,))
    for i in range(n):
        cat1_permuted = rng.permutation(indata[:,2])
        cat2_permuted = 1 - cat1_permuted
        permuted_indata[:,2] = cat1_permuted
        permuted_indata[:,3] = cat2_permuted
        
        interact_permuted[:,0:4] = permuted_indata[:,0:4]
        interact_permuted[:,4] = interact_permuted[:,1] * interact_permuted[:,2]
        interact_permuted[:,5] = interact_permuted[:,1] * interact_permuted[:,3]
        
        permuted_pred = model.predict(interact_permuted)
        neg_perm_mse = -ymax_logit_mse(permuted_pred, outdata)
        if neg_perm_mse >= neg_base_mse:
            better_count += 1
        mmts[i] = neg_perm_mse

    
    if better_count >= 10:
        p_value = (better_count + 1) / (n + 1)
        p_variance = p_value * (1 - p_value) / n
        interval = stats.norm.interval(0.95, p_value, p_variance)
        return p_value, *interval, 'ecdf', -1, ident, cat1, cat2
    
    else:
        sorted_measures = np.sort(mmts)
        thresh_idx = -250
        largest = sorted_measures[thresh_idx:]
        thresh = (sorted_measures[thresh_idx - 1] + sorted_measures[thresh_idx]) / 2
        exceeds = largest - thresh
        gpd_params = stats.genpareto.fit(exceeds)
        # use the 1-sample Kolmogorov-Smirnov test 
        # to check whether the exceedances are from a GPD
        test = stats.ks_1samp(exceeds, stats.genpareto.cdf, args=gpd_params)
        iterations = 0
        
        # we want to have a GPD, which is KS null hypothesis
        while test.pvalue < 0.05:
            iterations += 1
            thresh_idx += 10
            if thresh_idx >= 0: # if we cannot find a GPD, fallback to basic method
                p_value = (better_count + 1) / (n + 1)
                p_variance = p_value * (1 - p_value) / n
                interval = stats.norm.interval(0.95, p_value, p_variance)
                return (better_count + 1) / (n + 1), *interval, 'ecdf-fallback', iterations, ident, cat1, cat2
            largest = sorted_measures[thresh_idx:]
            thresh = (sorted_measures[thresh_idx - 1] + sorted_measures[thresh_idx]) / 2
            exceeds = largest - thresh
            gpd_params = stats.genpareto.fit(exceeds)
            test = stats.ks_1samp(exceeds, stats.genpareto.cdf, args=gpd_params)
        
        sf = stats.genpareto.sf(neg_base_mse - thresh, *gpd_params)
        p_value = (len(exceeds) / n) * sf
        
        # calculate confidence interval using binomial distribution
        count_low, count_high = stats.binom.interval(0.95, n, p_value)
        low_bound = count_low / n
        high_bound = count_high / n
        
        return p_value, low_bound, high_bound, 'gcd', iterations, ident, cat1, cat2
    


def _permutation_test_task(args, n=50_000):
    rng = np.random.default_rng(RNG_SEED)
    model, indata, outdata, treatdata, cat1, cat2, ident = args
    return _permutation_test(model, indata, outdata, treatdata, cat1, cat2, ident, rng, n)


def permutation_tests(models, intables, outcols, treatinfo, catlefts, catrights, idents):
    arguments = list(zip(models,
                         intables,
                         outcols,
                         treatinfo,
                         catlefts,
                         catrights,
                         idents))
    pvalues = []
    with multiproc.Pool(N_PROCS) as pool:
        pvalues_iter = pool.imap_unordered(_permutation_test_task, arguments, 256)
        for row in pvalues_iter:
            pvalues.append(row)
            print(row)
    return pvalues

def create_subtables(focused_data, dataprep):
    
    logger = logging.Logger('subtable')
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(cetsa_paths.get_logging_path("subtables.log"))
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    
    inputs = []
    outputs = []
    treats = []
    cond1s = []
    cond2s = []
    prot_identities = []
    subtables = []
    
    treatments = set(focused_data['Treatment'])
    treatments.remove('DMSO')
    controls = {'DMSO', 'Myricetin'}
    
    pairs = itertools.product(treatments, controls)
    depairs = filter(lambda x : x[0] != x[1], pairs)
    
    
    for cond1, cond2 in depairs:
        try:
            indatas, outdatas, treatdatas, prot_ids = dataprep.transform(focused_data, cond1, cond2)
            prots = prot_ids.reset_index()
            reset_treat = treatdatas.reset_index()
            reset_focused = focused_data.reset_index()
            for prot in prots.groupby(by=['PG.ProteinAccessions', 'PG.Genes']):
                idx = prot[1].index
                
                indata = indatas[idx, :]
                outdata = outdatas[idx]
                treatdata = reset_treat.loc[idx]
                prot_id = prot[0]
                prot_table = reset_focused.loc[idx]
                
                inputs.append(indata)
                outputs.append(outdata)
                treats.append(treatdata)
                cond1s.append(cond1)
                cond2s.append(cond2)
                prot_identities.append(prot_id)
                subtables.append(prot_table)
                
        except Exception as e:
            raise e
    
    return (inputs,
            outputs,
            treats,
            cond1s,
            cond2s,
            prot_identities,
            subtables)
    
    
    
    # proteins = focused_data.groupby(by=['PG.ProteinAccessions', 'PG.Genes'])
    
    # for prot_id, prot_table in proteins:
        
    #     clean_table = prot_table.dropna(subset=['Temperature', NORMPROT],
    #                                     ignore_index=True)
        
    #     pairs = itertools.product(treatments, controls)
    #     depairs = filter(lambda x : x[0] != x[1], pairs)
    #     for cond1, cond2 in depairs:
    #         try:
    #             indata, outdata, treatdata, _ = dataprep.transform(clean_table, cond1, cond2)
    #         except Exception as e:
    #             logger.warning(f"Could not subtable {prot_id}, {cond1}, {cond2}",
    #                           exc_info=e)
    #             continue
    #         inputs.append(indata)
    #         outputs.append(outdata)
    #         treats.append(treatdata)
    #         cond1s.append(cond1)
    #         cond2s.append(cond2)
    #         prot_identities.append(prot_id)
    #         subtables.append(prot_table)
    
    # return inputs, outputs, treats, cond1s, cond2s, prot_identities, subtables

    

def _model_fit_task(model_inputs):
    interact_in, out = model_inputs
    model = NPARCModel(alpha=0.0)
    model.fit(interact_in, out)
    return model

def pfit_models(datapieces):
    with multiproc.Pool(N_PROCS) as p:
        models = p.map(_model_fit_task, datapieces)
    return models



def display_graphs(filename, sig_table, data_table, dataprep, palette=None):
    with PdfPages(cetsa_paths.get_outdir() / filename, keep_empty=False) as pdf:
        for idx in sig_table.index:
            ax = pyplot.subplot()
            acc = sig_table.loc[idx, 'PG.ProteinAccessions']
            gene = sig_table.loc[idx, 'PG.Genes']
            #model = sig_table.loc[idx, 'model']
            treatment1 = sig_table.loc[idx, 'Treatment 1']
            treatment2 = sig_table.loc[idx, 'Treatment 2']
            
            subdata = data_table.loc[(data_table['PG.ProteinAccessions'] == acc) & (data_table['PG.Genes'] == gene),:]
            subdata = subdata.loc[subdata['Treatment'].isin([treatment1,treatment2])]
            
            ins, outs, treats = dataprep.transform(subdata, treatment1, treatment2)
            
            hypo = pandas.DataFrame({'Temperature' : list(range(37,71)) * 2,
                                     'Treatment': [treatment1] * (71 - 37) + [treatment2] * (71-37),
                                     NORMPROT: [0.0] * (2 * (71 - 37))})
            
            seaborn.scatterplot(subdata,
                                x='Temperature',
                                y='Normalized_FG_Quantity',
                                hue='Treatment',
                                style='Treatment',
                                ax=ax,
                                palette=palette)
            
            hypo_in, _, hypo_treats = dataprep.transform(hypo, treatment1, treatment2)
            
            model = sig_table.loc[idx, 'model']
            
            hypo_preds = model.predict(hypo_in)
            
            seaborn.lineplot(x=hypo['Temperature'],
                             y=hypo_preds,
                             hue=hypo_treats,
                             style=hypo_treats,
                             ax=ax,
                             palette=palette)
            
            ax.set_title(gene)
            ax.set_ylabel("Relative Soluble Protein")
            fig = ax.get_figure()
            pdf.savefig(fig)
            ax.cla()     

def main(datapath=None, candidatepath=None, outdir=None):
    data, candidates = load.prepare_data(False, data_path=datapath, candidate_path=candidatepath)
    # we will reload the candidates when we need them later
    del candidates
    narrow_data = data.loc[:, ['PG.ProteinAccessions',
                                  'PG.Genes',
                                  'R.Replicate',
                                  'Temperature', 
                                  'Treatment', 
                                  NORMPROT]]

    #half the time gets spend on this shit in particular
    # and its critical because otherwise we get _nonsense_
    narrow_data = narrow_data.drop_duplicates()
    dataprep = DataPreparer(narrow_data)
    
    prepped_info = create_subtables(narrow_data, dataprep)
    
    in_datas, out_datas, treat_labels, cond_lefts, cond_rights, prot_idents, subtables = prepped_info
    
    model_inputs = list(zip(in_datas, out_datas))
    
    print("begin model fitting")
    
    models = pfit_models(model_inputs)
    
    print("models ready")
    
    prot_accessions = [p[0] for p in prot_idents]
    gene_ids = [p[1] for p in prot_idents]
    datatable = pandas.DataFrame({'PG.ProteinAccessions': prot_accessions,
                                  'PG.Genes': gene_ids,
                                  'Treatment 1' : cond_lefts,
                                  'Treatment 2': cond_rights,
                                  'model': models})
    
    print('started calculating inflection')
    
    with multiproc.Pool(N_PROCS) as pool:
        t1_inflects = pool.map(t_infl.get_T1_inflection, models)
        t2_inflects = pool.map(t_infl.get_T2_inflection, models)
    
    print('finished calculating inflection')
    
    datatable.loc[:, 'T_infl_Treatment_1'] = t1_inflects
    datatable.loc[:, 'T_infl_Treatment_2'] = t2_inflects
    
    datatable.loc[:,'converged'] = datatable['model'].map(lambda m : m.fit_success_)
    
    print('begin perm tests')
    
    perm_tests = permutation_tests(models,
                                   in_datas,
                                   out_datas,
                                   treat_labels,
                                   cond_lefts,
                                   cond_rights,
                                   prot_idents)
    
    perm_table = pandas.DataFrame(data=perm_tests,
                                  columns=['pvalue',
                                           'pvalue_lowbound',
                                           'pvalue_highbound',
                                           'method',
                                           'gcd_iterations',
                                           'ident',
                                           'Treatment 1',
                                           'Treatment 2'])
    
    perm_table.loc[:, 'PG.ProteinAccessions'] = perm_table['ident'].map(lambda x : x[0])
    perm_table.loc[:, 'PG.Genes'] = perm_table['ident'].map(lambda x : x[1])
    
    table = datatable.merge(perm_table,
                             how='left',
                             on=['PG.ProteinAccessions',
                                 'PG.Genes',
                                 'Treatment 1',
                                 'Treatment 2'])
    
    if outdir is None:
        outdir = cetsa_paths.get_outdir()
    
    outpath = outdir / 'nparc_outputs_Oct2024.csv'
    
    
    fisetin = table.loc[table['Treatment 1'] == 'Fisetin', :].copy()
    
    quercetin = table.loc[table['Treatment 1'] == 'Quercetin', :].copy()
    
    myricetin = table.loc[table['Treatment 1'] == 'Myricetin',:].copy()
    
    table.loc[:, 'bh_pval'] = stats.false_discovery_control(table['pvalue'])
    
    table.to_csv(outpath)
    
    fisetin.loc[:, 'bh_pval'] = stats.false_discovery_control(fisetin['pvalue'])
    
    quercetin.loc[:, 'bh_pval'] = stats.false_discovery_control(quercetin['pvalue'])
    
    myricetin.loc[:, 'bh_pval'] = stats.false_discovery_control(myricetin['pvalue'])
    
    fisetinpath = outdir / 'nparc_fisetin_Oct2024.csv'
    
    quercetinpath = outdir / 'nparc_quercetin_Oct2024.csv'
    
    myricetinpath = outdir / 'nparc_myricetin_Oct2024.csv'
    
    fisetin.to_csv(fisetinpath)
    
    quercetin.to_csv(quercetinpath)
    
    myricetin.to_csv(myricetinpath)
    
    sig_fisetin = fisetin.loc[fisetin['bh_pval'] < 0.05, :]
    
    sig_quercetin = quercetin.loc[quercetin['bh_pval'] < 0.05, :]
    
    sig_myricetin = myricetin.loc[myricetin['bh_pval'] < 0.05, :]
    
    sig_fisetin.to_csv(outdir / 'nparc_sig_fisetin_Oct2024.csv')
    
    sig_quercetin.to_csv(outdir / 'nparc_sig_quercetin_Oct2024.csv')
    
    sig_myricetin.to_csv(outdir / 'nparc_sig_myricetin_Oct2024.csv')
    
    myr_genes = set(sig_myricetin['PG.Genes'])
    
    unshared_fisetin = sig_fisetin.loc[~sig_fisetin['PG.Genes'].isin(myr_genes), :]
    
    unshared_quercetin = sig_quercetin.loc[~sig_quercetin['PG.Genes'].isin(myr_genes), :]
    
    unshared_fisetin.to_csv(outdir / 'nparc_unshared_fisetin_Oct2024.csv')
    
    unshared_quercetin.to_csv(outdir / 'nparc_unshared_quercetin_Oct2024.csv')
    
    display_graphs("nparc_Oct2024_fisetin.pdf", unshared_fisetin, narrow_data, dataprep, dataprep.palette())
    
    display_graphs("nparc_Oct2024_quercetin.pdf", unshared_quercetin, narrow_data, dataprep, dataprep.palette())
    
    
    return table

                  

if __name__ == '__main__':
    import time
    start = time.time()
    alldata = main()
    end = time.time()
    elapsed = end - start
    print(elapsed / 60)