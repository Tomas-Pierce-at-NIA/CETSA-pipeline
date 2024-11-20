# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:30:58 2024

@author: piercetf
"""

NORMPROT = 'Normalized_FG_Quantity'
N_COND = 4

import cetsa2
from data_prepper import DataPreparer
#from data_prepper2 import DataPreparer
import t_infl
import load_monocyte_cetsa_data as load
import cProfile
import pandas
import numpy as np

def main():
    print("loading")
    data, cand = load.prepare_data(False)
    print("narrowing")
    narrow_data = data.loc[:, ['PG.ProteinAccessions',
                                  'PG.Genes',
                                  'R.Replicate',
                                  'Temperature', 
                                  'Treatment', 
                                  NORMPROT]]
    print("de-duping")
    narrow_data = narrow_data.drop_duplicates()
    print("data prepping")
    dataprep = DataPreparer(narrow_data)
    prepped_subtables = cetsa2.create_subtables(narrow_data, dataprep)
    in_datas, out_datas, treat_labels, cond_lefts, cond_rights, prot_idents, subtables = prepped_subtables
    model_inputs = list(zip(in_datas, out_datas))
    print("model fitting")
    models = []
    i = 0
    for model_train in model_inputs:
        if i == 100:
            break
        model = cetsa2._model_fit_task(model_train)
        models.append(model)
        i += 1
        print(i)
    
    prot_accessions = [p[0] for p in prot_idents]
    gene_ids = [p[1] for p in prot_idents]
    datatable = pandas.DataFrame({'PG.ProteinAccessions': prot_accessions[:100],
                                  'PG.Genes': gene_ids[:100],
                                  'Treatment 1' : cond_lefts[:100],
                                  'Treatment 2': cond_rights[:100],
                                  'model': models[:100]})
    t1_inflects = []
    t2_inflects = []
    j = 0
    for model in models:
        if j == 100:
            break
        t1 = t_infl.get_T1_inflection(model)
        t2 = t_infl.get_T2_inflection(model)
        t1_inflects.append(t1)
        t2_inflects.append(t2)
        j += 1
        print(j)
    
    datatable.loc[:, 'T_infl_Treatment_1'] = t1_inflects
    datatable.loc[:, 'T_infl_Treatment_2'] = t2_inflects
    
    datatable.loc[:,'converged'] = datatable['model'].map(lambda m : m.fit_success_)
    
    
    k = 0
    for model in models:
        if k == 100:
            break
        cetsa2._permutation_test(model, 
                                 in_datas[k],
                                 out_datas[k],
                                 treat_labels[k],
                                 cond_lefts[k],
                                 cond_rights[k],
                                 prot_idents[k],
                                 np.random.default_rng(20242011112),
                                 n=50_000)
        k += 1
        print(k)


if __name__ == '__main__':
    profile = cProfile.Profile()
    profile.enable()
    try:
        main()
    finally:
        profile.disable()
        profile.dump_stats(r"C:\Users\piercetf\Projects\nparc_profile.profile")