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
import pathlib


from data_prepper import DataPreparer
from nparc_model import ScaledNPARCModel, calc_bayes_factor
import load_monocyte_cetsa_data as load
import cetsa_paths
import t_infl
import permutation_test as ptest

warnlogname = str(cetsa_paths.get_logging_path('cetsa_debug.log'))
logging.basicConfig(filename=warnlogname, level=logging.DEBUG)
logging.captureWarnings(True)

#from sklearn.preprocessing import PolynomialFeatures
import seaborn
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import pandas

from scipy import stats

import toml


RNG_SEED = time.time_ns()

if platform.system() == 'Windows':
    N_PROCS = multiproc.cpu_count() - 1
elif platform.system() == 'Linux':
    import os
    proc_env_var = 'SLURM_CPUS_PER_TASK'
    if proc_env_var in os.environ:
        procs_str = os.environ[proc_env_var]
        N_PROCS = int(procs_str)
    else:
        N_PROCS = 2


def _perm_test_w(arg):
    return _perm_test(*arg)

def _perm_test(model, intab, outcol, treat, cat1, cat2, ident, n=50_000):
    test = ptest.PermutationTest(model, intab, outcol, treat, cat1, cat2, ident)
    return test.permutation_test(n=n)

def permutation_tests(pool, models, intables, outcols, treatinfo, cats1, cats2, idents):
    input_data = list(zip(models,
                          intables,
                          outcols,
                          treatinfo,
                          cats1,
                          cats2,
                          idents))
    
    pvals = []
    pvalues_iter = pool.imap_unordered(
        _perm_test_w, 
        input_data,
        1024)
    for row in pvalues_iter:
        pvals.append(row)
        print(row)
    
    return pvals

def permutation_tests_iter(models, intables, outcols, treatinfo, cats1, cats2, idents):
    input_data = list(zip(models,
                          intables,
                          outcols,
                          treatinfo,
                          cats1,
                          cats2,
                          idents))
    
    pvals = []
    for box in input_data:
        pval = _perm_test_w(box)
        print(pval)
        pvals.append(pval)
    
    return pvals

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
    null_inputs = []
    
    treatments = set(focused_data['Treatment'])
    
    configpath = pathlib.Path(__file__).with_name("cetsa_config.toml")
    with configpath.open("r") as config:
        params = toml.load(config)
    
    v_control = params['controls']['vehicle']
    nonsenolytics = set(params['controls']['nonsenolytic'])
    treatments.remove(v_control)
    controls = {v_control, *nonsenolytics}
    pairs = itertools.product(treatments, controls)
    depairs = filter(lambda x : x[0] != x[1], pairs)
    
    
    # we were able to reduce runtime from 2.5 minutes
    # to 12 seconds by reducing the number of times
    # we have to call the dataprep.transform function
    # by running it across the entire dataset for each
    # comparison and separating the proteins out second
    # rather than runnning it on each protein subpage
    # in separate calls
    for cond1, cond2 in depairs:
        try:
            indatas, outdatas, treatdatas, prot_ids = dataprep.transform(focused_data, cond1, cond2)
            
            null_indatas, _, _, _ = dataprep.null_model_cols_transform(focused_data, cond1, cond2)
            
            #breakpoint()
            
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
                
                null_indata = null_indatas[idx, :]
                
                inputs.append(indata)
                outputs.append(outdata)
                treats.append(treatdata)
                cond1s.append(cond1)
                cond2s.append(cond2)
                prot_identities.append(prot_id)
                subtables.append(prot_table)
                null_inputs.append(null_indata)
                
        except Exception as e:
            raise e
    
    return (inputs,
            outputs,
            treats,
            cond1s,
            cond2s,
            prot_identities,
            subtables,
            null_inputs)
    

def _model_fit_task(model_inputs):
    interact_in, out = model_inputs
    model = ScaledNPARCModel(alpha=0.0)
    model.fit(interact_in, out)
    return model

def pfit_models(pool, datapieces):
    models = pool.map(_model_fit_task, datapieces)
    return models



def display_graphs(filename, sig_table, data_table, dataprep, palette=None, outdir=None):
    if outdir is None:
        outdir = cetsa_paths.get_outdir()
    with PdfPages(outdir / filename, keep_empty=False) as pdf:
        for idx in sig_table.index:
            ax = pyplot.subplot()
            acc = sig_table.loc[idx, 'PG.ProteinAccessions']
            gene = sig_table.loc[idx, 'PG.Genes']
            #model = sig_table.loc[idx, 'model']
            treatment1 = sig_table.loc[idx, 'Treatment 1']
            treatment2 = sig_table.loc[idx, 'Treatment 2']
            
            subdata = data_table.loc[(data_table['PG.ProteinAccessions'] == acc) & (data_table['PG.Genes'] == gene),:]
            subdata = subdata.loc[subdata['Treatment'].isin([treatment1,treatment2])]
            
            ins, outs, treats, _protids = dataprep.transform(subdata, treatment1, treatment2)
            
            hypo = pandas.DataFrame({'Temperature' : list(range(37,71)) * 2,
                                     'Treatment': [treatment1] * (71 - 37) + [treatment2] * (71-37),
                                     NORMPROT: [0.0] * (2 * (71 - 37)),
                                     'PG.ProteinAccessions' : [''] * (2 * (71 - 37)),
                                     'PG.Genes' : [''] * (2 * (71 - 37))})
            
            seaborn.scatterplot(subdata,
                                x='Temperature',
                                y='Normalized_FG_Quantity',
                                hue='Treatment',
                                style='Treatment',
                                ax=ax,
                                palette=palette)
            
            hypo_in, hypo_out, hypo_treats, hypo_protids = dataprep.transform(hypo, treatment1, treatment2)
            
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
    
    in_datas, out_datas, treat_labels, cond_lefts, cond_rights, prot_idents, subtables, null_inputs = prepped_info
    
    model_inputs = list(zip(in_datas, out_datas))
    
    nullmodel_inputs = list(zip(null_inputs, out_datas))
    
    with multiproc.Pool(N_PROCS) as pool:
    
        print("begin model fitting")
        models = pfit_models(pool, model_inputs)
        print("models ready")
        
        print("begin null model fitting")
        null_models = pfit_models(pool, nullmodel_inputs)
        print("null models ready")
        
        print("begin calculating bayes factors")
        bfactors = [calc_bayes_factor(alt, null) for alt, null in zip(models, null_models)]
        print("bayes factors ready")
        
        
        
        print('started calculating inflection')
        t1_inflects = pool.map(t_infl.get_T1_inflection, models)
        t2_inflects = pool.map(t_infl.get_T2_inflection, models)
        #datatable.loc[:, 'T_infl_Treatment_1'] = t1_inflects
        #datatable.loc[:, 'T_infl_Treatment_2'] = t2_inflects
        print('finished calculating inflection')
        
        prot_accessions = [p[0] for p in prot_idents]
        gene_ids = [p[1] for p in prot_idents]
        datatable = pandas.DataFrame({'PG.ProteinAccessions': prot_accessions,
                                      'PG.Genes': gene_ids,
                                      'Treatment 1' : cond_lefts,
                                      'Treatment 2': cond_rights,
                                      'model': models,
                                      'T_infl_Treatment_1' : t1_inflects,
                                      'T_infl_Treatment_2' : t2_inflects,
                                      'Bayes_factors': bfactors})
        
        datatable.loc[:,'converged'] = datatable['model'].map(lambda m : m.fit_success_)
        
        datatable.loc[:, 'T1_Decreasing'] = datatable['model'].map(lambda m : m.treatment1_decreasing)
        datatable.loc[:, 'T2_Decreasing'] = datatable['model'].map(lambda m : m.treatment2_decreasing)
        
        print('begin perm tests')
        
        perm_tests = permutation_tests(pool,
                                       models,
                                       in_datas,
                                       out_datas,
                                       treat_labels,
                                       cond_lefts,
                                       cond_rights,
                                       prot_idents)
        print("perm tests finished")
        
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
    
    #data2, candidates = load.prepare_data
    data2, candidates2 = load.prepare_data(False, datapath, candidatepath)
    del data2
    
    
    if outdir is None:
        outdir = cetsa_paths.get_outdir()
        
    
    
    table.loc[:, 'bh_pval'] = stats.false_discovery_control(table['pvalue'])
    table.to_csv(outdir / 'nparc_outputs_Oct2024.csv')
    
    conditions = table['Treatment 1'].unique()
    
    configparams = cetsa_paths.loadparams()
    
    nonsenolytics = configparams['controls']['nonsenolytic']
    
    # set of genes which interact with nonsenolytic controls
    ns_genes = set()
    
    for non_seno in nonsenolytics:
        nonseno_nospace = non_seno.replace(' ', '_')
        nstable = table.loc[table['Treatment 1'] == non_seno, :].copy()
        nstable.loc[:, 'bh_pval'] = stats.false_discovery_control(nstable['pvalue'])
        filename = 'nparc_nscontrol_{}_Jan2025.csv'.format(nonseno_nospace)
        filepath = outdir / filename
        nstable.to_csv(filepath)
        ns_sig = nstable.loc[nstable['bh_pval'] < 0.05, :]
        sig_genes = ns_sig['PG.Genes']
        ns_genes.update(sig_genes)
        
        
    
    for condition in conditions:
        cond_nospace = condition.replace(" ", "_")
        cond_table = table.loc[table['Treatment 1'] == condition, :].copy()
        cond_table.loc[:, 'bh_pval'] = stats.false_discovery_control(cond_table['pvalue'])
        filename = 'nparc_{}_Jan2025.csv'.format(cond_nospace)
        filepath = outdir / filename
        cond_table.to_csv(filepath)
        sig_table = cond_table.loc[cond_table['bh_pval'] < 0.05, :]
        sigfilename = 'nparc_sig_{}_Jan2025.csv'.format(cond_nospace)
        sigpath = outdir / sigfilename
        sig_table.to_csv(sigpath)
        unshared = sig_table.loc[~sig_table['PG.Genes'].isin(ns_genes), :]
        unshare_filename = 'nparc_unshare_{}_Jan2025.csv'.format(cond_nospace)
        unshare_path = outdir / unshare_filename
        unshared.to_csv(unshare_path)
        
        graphname = "nparc_Jan2025_{}.pdf".format(cond_nospace)
        display_graphs(graphname,
                       unshared,
                       narrow_data,
                       dataprep,
                       dataprep.palette(),
                       outdir)
    

    
    return table

                  

if __name__ == '__main__':
    alldata = main()
    