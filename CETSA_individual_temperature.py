# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:47:46 2024

@author: piercetf

Implement the individual temperature analysis method in Python
"""

NORMPROT = 'Normalized_FG_Quantity'

from scipy import stats, integrate
import pandas
import numpy as np
import seaborn
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot

import load_monocyte_cetsa_data as load
import cetsa_paths


def get_all_student_tests(data):
    conditions = set(data['Treatment'])
    controls = {'DMSO', 'Myricetin'}
    treatments = list(conditions.difference(controls))
    treatments.append('Myricetin')
    
    groups = data.groupby(by=['PG.ProteinAccessions',
                              'PG.Genes'])
    
    student_table = []
    
    for ident, table in groups:
        temp_gs = table.groupby(by=['Temperature'])
        for temp, subtable in temp_gs:
            for treatment in treatments:
                for control in controls:
                    treat_sample = subtable.loc[subtable['Treatment'] == treatment, NORMPROT]
                    control_sample = subtable.loc[subtable['Treatment'] == control, NORMPROT]
                    t_test = stats.ttest_ind(treat_sample, control_sample)
                    row = (*ident,
                           *temp,
                           treatment,
                           control,
                           t_test.pvalue,
                           t_test.statistic,
                           t_test.df)
                    student_table.append(row)
    
    student_frame = pandas.DataFrame(student_table,
                                     columns=['PG.ProteinAccessions',
                                              'PG.Genes',
                                              'Temperature',
                                              'Treatment',
                                              'Control',
                                              'T-test_pvalue',
                                              'T-test_statistic',
                                              'T-test_df'])
    
    return student_frame


def get_all_dunnett_tests(data):
    conditions = set(data['Treatment'])
    control = 'DMSO'
    treatments = list(conditions.difference({control}))
    
    groups = data.groupby(by=['PG.ProteinAccessions',
                              'PG.Genes'])
    
    dunnet_table = []
    
    for ident, table in groups:
        temp_gs = table.groupby(by=['Temperature'])
        for temp, subtable in temp_gs:
            samples = []
            for treatment in treatments:
                sample = subtable.loc[subtable['Treatment'] == treatment, NORMPROT]
                samples.append(sample)
            control_measures = subtable.loc[subtable['Treatment'] == control, NORMPROT]
            if len(samples[0]) < 2 or len(samples[1]) < 2 or len(samples[2]) < 2 or len(control_measures) < 2:
                continue
            dunnet_res = stats.dunnett(*samples, control=control_measures)
            dunnet_row0 = (*ident, 
                           *temp, 
                           treatments[0], 
                           control, 
                           dunnet_res.pvalue[0], 
                           dunnet_res.statistic[0])
            dunnet_row1 = (*ident, 
                           *temp, 
                           treatments[1], 
                           control, 
                           dunnet_res.pvalue[1], 
                           dunnet_res.statistic[1])
            dunnet_row2 = (*ident, 
                           *temp, 
                           treatments[2], 
                           control, 
                           dunnet_res.pvalue[2], 
                           dunnet_res.statistic[2])
            dunnet_table.append(dunnet_row0)
            dunnet_table.append(dunnet_row1)
            dunnet_table.append(dunnet_row2)
    
    dunnetts_frame = pandas.DataFrame(data=dunnet_table,
                                      columns=['PG.ProteinAccessions',
                                               'PG.Genes',
                                               'Temperature',
                                               'Treatment',
                                               'Control',
                                               'dunnetts_pvalue',
                                               'dunnetts_statistic'])
    
    return dunnetts_frame


def calc_fisher_statistic(pvals):
    return -2 * np.sum(np.log(pvals))

def get_chi_mean(cov_mtx):
    return 2 * len(cov_mtx)

def get_chi_variance(cov_mtx):
    k_conds = len(cov_mtx)
    cov_total = 0
    for i in range(cov_mtx.shape[0]):
        for j in range(i + 1, cov_mtx.shape[1]):
            cov_total += cov_mtx[i,j]
    return (4 * k_conds) + (2 * cov_total)

def get_p_value(deg_free, scale_param, fisher_stat):
    return stats.chi2.sf(fisher_stat, deg_free, scale=scale_param)

def collect_students(students_table):
    
    students_table.loc[:, 'log_pval'] = np.log(students_table['T-test_pvalue']) * (-2)
    
    #breakpoint()
    pivot = students_table.pivot_table(index=['PG.Genes'],
                                       columns="Temperature",
                                       values="log_pval")
    stud_cov = pivot.cov()
    
    stud_cov = np.array(stud_cov)
    
    chi_mean = get_chi_mean(stud_cov)
    
    chi_variance = get_chi_variance(stud_cov)
    
    # degrees of freedom
    f = 2 * (chi_mean**2) / chi_variance
    
    # scale parameter
    c = chi_variance / (2 * chi_mean)
    
    combo_t_pvals = []
    
    for specific, subtable in students_table.groupby(by=['PG.ProteinAccessions',
                                                         'PG.Genes',
                                                         'Treatment',
                                                         'Control']):
        accession, gene, treatment, control = specific
        
        pvals = subtable['T-test_pvalue']
        
        fisher_stat = calc_fisher_statistic(pvals)
        
        brown_pval = get_p_value(f, c, fisher_stat)
        
        combo_t_pvals.append((accession, gene, treatment, control, brown_pval, fisher_stat))
    
    combo_frame = pandas.DataFrame(combo_t_pvals,
                                   columns=['PG.ProteinAccessions',
                                            'PG.Genes',
                                            'Treatment',
                                            'Control',
                                            'combo_student_pvalue',
                                            'Fisher statistic'])
    
    return combo_frame


def collect_dunnetts(dunnetts_table):
    combo_d_pvals = []
    
    dunnetts_table.loc[:, 'log_pval'] = np.log(dunnetts_table.loc[:,'dunnetts_pvalue']) * (-2)
    
    pivot = dunnetts_table.pivot_table(index=["PG.Genes"], 
                                 columns="Temperature",
                                 values="log_pval")
    
    dunn_cov = pivot.cov()
    
    dunn_cov = np.array(dunn_cov)
    
    chi_mean = get_chi_mean(dunn_cov)
    
    chi_variance = get_chi_variance(dunn_cov)
    
    # degrees freedom
    f = 2 * (chi_mean**2) / chi_variance
    
    #scale parameter
    c = chi_variance / (2 * chi_mean)
    
    for specific, subtable in dunnetts_table.groupby(by=['PG.ProteinAccessions',
                                                         'PG.Genes',
                                                         'Treatment']):
        accession, gene, treatment = specific
        
        pvals = subtable['dunnetts_pvalue']
        
        #combo_res = stats.combine_pvalues(pvals)
        #combo_pval = stats.hmean(pvals)
        #combo_pval = asym_hmean_pval(pvals)
        fisher_stat = calc_fisher_statistic(pvals)
        brown_pval = get_p_value(f, c, fisher_stat)
        
        combo_d_pvals.append((accession, gene, treatment, brown_pval, fisher_stat))
    
    combo_frame = pandas.DataFrame(combo_d_pvals,
                                   columns=['PG.ProteinAccessions',
                                             'PG.Genes',
                                             'Treatment',
                                             'combo_dunnett_pvalue',
                                             'Fisher statistic'])
    
    return combo_frame


def graph_all_proteins(data, all_comps, focus, filename, treatment, control):
    
    treatments = ['DMSO', 'Fisetin', 'Quercetin', 'Myricetin']
    
    colors = seaborn.color_palette('hls', len(treatments))
    palette = dict(zip(treatments, colors))
    
    p_tab = focus.loc[:, ['PG.ProteinAccessions',
                          'PG.Genes',
                          'bh_pval']]
    
    subdata = data[data['PG.ProteinAccessions'].isin(focus['PG.ProteinAccessions'])]
    
    subdata = subdata.merge(p_tab,
                            how='inner',
                            on=['PG.ProteinAccessions',
                                'PG.Genes'])
    
    subdata = subdata.sort_values(by=['bh_pval'])
    
    sub_it_comps = all_comps.loc[(all_comps['Treatment'] == treatment) & (all_comps['Control'] == control),
                                 :]
    
    with PdfPages(filename, keep_empty=False) as pdf:
        for prot_name, prot_table in subdata.groupby(by=['PG.ProteinAccessions',
                                                      'PG.Genes'],
                                                     sort=False):
            
            prot_acc, gene_id = prot_name
            
            pvalue = focus.loc[focus['PG.ProteinAccessions'] == prot_acc,
                               'bh_pval'].item()
            
            prot_comps = sub_it_comps.loc[sub_it_comps['PG.ProteinAccessions'] == prot_acc,
                                          :]
            
            if 'T-test_pvalue' in prot_comps.columns:
                comp_table = prot_comps.loc[:, ['Temperature',
                                                'T-test_pvalue']]
                comp_table.loc[:, 'p-value'] = comp_table['T-test_pvalue']
                comp_table.loc[:, '-log(p-value)'] = -np.log(comp_table['p-value'])
            else:
                comp_table = prot_comps.loc[:, ['Temperature',
                                                'dunnetts_pvalue']]
                comp_table.loc[:, 'p-value'] = comp_table['dunnetts_pvalue']
                comp_table.loc[:, '-log(p-value)'] = -np.log(comp_table['p-value'])
                
            
            ptable = prot_table.loc[(prot_table['Treatment'] == treatment) | (prot_table['Treatment'] == control), :]
            
            mean_table = ptable.groupby(
                by=['Treatment',
                    'Temperature']).mean(numeric_only=True).reset_index()
            
            fig, axes = pyplot.subplots(nrows=2, height_ratios=[3,1])
            
            ax = seaborn.scatterplot(ptable, 
                                     x='Temperature',
                                     y=NORMPROT,
                                     hue='Treatment',
                                     style='Treatment',
                                     palette=palette,
                                     ax=axes[0])
            
            seaborn.lineplot(mean_table,
                             x='Temperature',
                             y=NORMPROT,
                             hue='Treatment',
                             ax=ax,
                             palette=palette,
                             errorbar=None)
            
            title = f"{prot_name[1]} {treatment}"
            ax.set_title(title)
            ax.set_ylabel("Relative Soluble Protein")
            ax.set_xlabel("Temperature")
            
            #fig = ax.get_figure()
            
            fig.text(0.9, 
                     0.01, 
                     f"adjusted p-value = {round(pvalue, 5)}", 
                     wrap=True, horizontalalignment='right')
            
            seaborn.barplot(comp_table,
                            x='Temperature',
                            y='p-value',
                            ax=axes[1])
            axes[1].axhline(0.05, color='red', ls='-.')
            
            #axes[1].sharex(axes[0])
            
            pdf.savefig(fig)
            ax.cla()
            
            pyplot.close(fig)
            
            print(prot_name)
    
    return


def auc_all_proteins(data):
    
    table = []
    
    for prot_name, prot_table in data.groupby(by=['PG.ProteinAccessions',
                                                  'PG.Genes']):
        mean_table = prot_table.groupby(by=['Temperature',
                                            'Treatment']
                                        ).mean(numeric_only=True).reset_index()
        
        dmso = mean_table.loc[mean_table['Treatment'] == 'DMSO', :]
        fisetin = mean_table.loc[mean_table['Treatment'] == 'Fisetin', :]
        quercetin = mean_table.loc[mean_table['Treatment'] == 'Quercetin', :]
        myricetin = mean_table.loc[mean_table['Treatment'] == 'Myricetin', :]
        
        dmso_auc = integrate.trapezoid(dmso[NORMPROT], 
                                       x=dmso['Temperature'])
        fisetin_auc = integrate.trapezoid(fisetin[NORMPROT],
                                          x=fisetin['Temperature'])
        quercetin_auc = integrate.trapezoid(quercetin[NORMPROT], 
                                            x=quercetin['Temperature'])
        myricetin_auc = integrate.trapezoid(myricetin[NORMPROT], 
                                            x=myricetin['Temperature'])
        
        fisetin_delta_auc = fisetin_auc - dmso_auc
        quercetin_delta_auc = quercetin_auc - dmso_auc
        myricetin_delta_auc = myricetin_auc - dmso_auc
        
        fisetin_r_auc = fisetin_auc / dmso_auc
        quercetin_r_auc = quercetin_auc / dmso_auc
        myricetin_r_auc = myricetin_auc / dmso_auc
        
        fisetin_delta2_auc = fisetin_auc - myricetin_auc
        quercetin_delta2_auc = quercetin_auc - myricetin_auc
        myricetin_delta2_auc = myricetin_auc - myricetin_auc
        
        fisetin_r2_auc = fisetin_auc / myricetin_auc
        quercetin_r2_auc = quercetin_auc / myricetin_auc
        myricetin_r2_auc = myricetin_auc / myricetin_auc
        
        
        
        
        rowf = (*prot_name, 
                'Fisetin', 
                fisetin_auc, 
                fisetin_delta_auc,
                fisetin_r_auc,
                fisetin_delta2_auc,
                fisetin_r2_auc)
        rowq = (*prot_name, 
                'Quercetin', 
                quercetin_auc, 
                quercetin_delta_auc,
                quercetin_r_auc,
                quercetin_delta2_auc,
                quercetin_r2_auc)
        rowm = (*prot_name, 
                'Myricetin', 
                myricetin_auc, 
                myricetin_delta_auc,
                myricetin_r_auc,
                myricetin_delta2_auc,
                myricetin_r2_auc)
        
        table.append(rowf)
        table.append(rowq)
        table.append(rowm)
    
    frame = pandas.DataFrame(table,
                             columns=['PG.ProteinAccessions',
                                      'PG.Genes',
                                      'Treatment',
                                      'AUC',
                                      'diff AUC vs control',
                                      'fold change AUC vs control',
                                      'diff AUC vs Myricetin',
                                      'fold change AUC vs Myricetin'])
    
    frame.loc[:, 'log_fold_change_auc'] = np.log(frame['fold change AUC vs control'])
    
    frame.loc[:, 'logfold_myricetin'] = np.log(frame['fold change AUC vs Myricetin'])
    
    return frame


def run_analysis(data, candidates, method='dunnett'):
    
    datadir = cetsa_paths.get_outdir()
    
    if method == 'dunnett':
        dunnetts = get_all_dunnett_tests(data)
        #ignore 37 degrees
        dunnetts = dunnetts.loc[dunnetts['Temperature'] > 37, :]
        combo_dunnetts = collect_dunnetts(dunnetts)
        auc_table = auc_all_proteins(data)
        combo_stats = combo_dunnetts.merge(auc_table,
                                           how='left',
                                           on=['PG.ProteinAccessions',
                                               'PG.Genes',
                                               'Treatment'])
        cd_quercetin = combo_stats.loc[combo_dunnetts['Treatment'] == 'Quercetin',:].copy()
        cd_fisetin = combo_stats.loc[combo_dunnetts['Treatment'] == 'Fisetin',:].copy()
        cd_myricetin = combo_stats.loc[combo_dunnetts['Treatment'] == 'Myricetin',:].copy()
        cd_quercetin.loc[:, 'bh_pval'] = stats.false_discovery_control(cd_quercetin['combo_dunnett_pvalue'],
                                                                       method='bh')
        cd_fisetin.loc[:,'bh_pval'] = stats.false_discovery_control(cd_fisetin['combo_dunnett_pvalue'],
                                                                    method='bh')
        cd_myricetin.loc[:,'bh_pval'] = stats.false_discovery_control(cd_myricetin['combo_dunnett_pvalue'],
                                                                      method='bh')
        sig_quercetin = cd_quercetin[cd_quercetin['bh_pval'] < 0.05]
        sig_fisetin = cd_fisetin[cd_fisetin['bh_pval'] < 0.05]
        sig_myricetin = cd_myricetin[cd_myricetin['bh_pval'] < 0.05]
        myr_genes = set(sig_myricetin['PG.Genes'])
        notshared_quercetin = sig_quercetin[~sig_quercetin['PG.Genes'].isin(myr_genes)]
        notshared_fisetin = sig_fisetin[~sig_fisetin['PG.Genes'].isin(myr_genes)]
        candidate_genes = candidates[['Genes', 
                                      'UniProtIds', 
                                      'ProteinDescriptions',
                                      'GO Biological Process',
                                      'GO Molecular Function',
                                      'GO Cellular Component']].drop_duplicates()
        focal_quercetin = notshared_quercetin.merge(candidate_genes,
                                                    how='left',
                                                    left_on=['PG.ProteinAccessions', 'PG.Genes'],
                                                    right_on=['UniProtIds', 'Genes'])
        focal_fisetin = notshared_fisetin.merge(candidate_genes,
                                                how='left',
                                                left_on=['PG.ProteinAccessions', 'PG.Genes'],
                                                right_on=['UniProtIds', 'Genes'])
        all_quercetin = cd_quercetin.merge(candidate_genes,
                                           how='left',
                                           left_on=['PG.ProteinAccessions', 'PG.Genes'],
                                           right_on=['UniProtIds', 'Genes'])
        all_fisetin = cd_fisetin.merge(candidate_genes,
                                       how='left',
                                       left_on=['PG.ProteinAccessions', 'PG.Genes'],
                                       right_on=['UniProtIds', 'Genes'])
        
        focal_quercetin.sort_values(by=['bh_pval'], inplace=True, kind='mergesort')
        focal_fisetin.sort_values(by=['bh_pval'], inplace=True, kind='mergesort')
        all_quercetin.sort_values(by=['bh_pval'], inplace=True, kind='mergesort')
        all_fisetin.sort_values(by=['bh_pval'], inplace=True, kind='mergesort')
        
        combo_stats.sort_values(by=['combo_dunnett_pvalue'], inplace=True, kind='mergesort')
        combo_stats.to_csv(datadir / "ITA_allcomps_Oct2024_Dunnett.csv")
        
        focal_quercetin.to_csv(datadir / "ITA_focal_quercetin_Oct2024_Dunnet.csv")
        focal_fisetin.to_csv(datadir / "ITA_focal_fisetin_Oct2024_Dunnet.csv")
        all_quercetin.to_csv(datadir / "ITA_all_quercetin_Oct2024_Dunnet.csv")
        all_fisetin.to_csv(datadir / "ITA_all_fisetin_Oct2024_Dunnet.csv")
        graph_all_proteins(data, 
                           dunnetts, 
                           focal_fisetin, 
                           datadir / "ITA_Fisetin_signif_graphs_Dunnet.pdf",
                           'Fisetin',
                           'DMSO')
        graph_all_proteins(data,
                           dunnetts,
                           focal_quercetin,
                           datadir / "ITA_Quercetin_signif_graphs_Dunnet.pdf",
                           'Quercetin',
                           'DMSO')
        
        return focal_fisetin, focal_quercetin
        
    elif method == 'student':
        students = get_all_student_tests(data)
        students = students.loc[students['Temperature'] > 37, :]
        combo_students = collect_students(students)
        auc_table = auc_all_proteins(data)
        combo_stats = combo_students.merge(auc_table,
                                           how='left',
                                           on=['PG.ProteinAccessions',
                                               'PG.Genes',
                                               'Treatment'])
        
        good_combo_stats = combo_stats.dropna(subset=['combo_student_pvalue'])
        
        good_combo_stats.loc[:, 'bh_pval'] =stats.false_discovery_control(good_combo_stats['combo_student_pvalue'])
        
        candidate_genes = candidates[['Genes', 
                                      'UniProtIds', 
                                      'ProteinDescriptions',
                                      'GO Biological Process',
                                      'GO Molecular Function',
                                      'GO Cellular Component']].drop_duplicates()
        
        all_info = good_combo_stats.merge(candidate_genes,
                                          how='left',
                                          left_on=['PG.ProteinAccessions', 'PG.Genes'],
                                          right_on=['UniProtIds', 'Genes'])
        
        sigtable = all_info.loc[all_info['bh_pval'] < 0.05, :]
        
        sigtable.sort_values(by=['bh_pval'], inplace=True, kind='mergesort')
        all_info.sort_values(by=['bh_pval'], inplace=True, kind='mergesort')
        
        sigtable.to_csv(datadir / "ITA_signif_comps_Student_Oct2024.csv")
        all_info.to_csv(datadir / "ITA_all_comps_Student_Oct2024.csv")
        
        fisetin1 = sigtable.loc[(sigtable['Treatment'] == 'Fisetin') & (sigtable['Control'] == 'DMSO'),:]
        fisetin2 = sigtable.loc[(sigtable['Treatment'] == 'Fisetin') & (sigtable['Control'] == 'Myricetin'),:]
        
        quercetin1 = sigtable.loc[(sigtable['Treatment'] == 'Quercetin') & (sigtable['Control'] == 'DMSO'),:]
        quercetin2 = sigtable.loc[(sigtable['Treatment'] == 'Quercetin')&(sigtable['Control']=='Myricetin'),:]
        
        myricetin = sigtable.loc[(sigtable['Treatment'] == 'Myricetin'),:]
        
        myr_genes = set(myricetin['PG.Genes'])
        
        fisetin_all = pandas.concat([fisetin1, fisetin2])
        quercetin_all = pandas.concat([quercetin1, quercetin2])
        
        fisetin_all = fisetin_all.loc[~fisetin_all['PG.Genes'].isin(myr_genes),:]
        quercetin_all = quercetin_all.loc[~quercetin_all['PG.Genes'].isin(myr_genes),:]
        
        fisetin_all.to_csv(datadir / "ITA_fisetin_Student_Oct2024.csv")
        quercetin_all.to_csv(datadir / "ITA_quercetin_Student_Oct2024.csv")
        
        graph_all_proteins(data,
                           students,
                           fisetin1,
                           datadir / "ITA_fisetinVDMSO_Student_Oct2024.pdf",
                           'Fisetin',
                           'DMSO')
        
        graph_all_proteins(data,
                           students,
                           fisetin2,
                           datadir / "ITA_fisetinVmyricetin_Student_Oct2024.pdf",
                           'Fisetin',
                           'Myricetin')
        
        graph_all_proteins(data,
                           students,
                           quercetin1,
                           datadir / "ITA_quercetinVDMSO_Student_Oct2024.pdf",
                           'Quercetin',
                           'DMSO')
        
        graph_all_proteins(data,
                           students,
                           quercetin2,
                           datadir / "ITA_quercetinVmyricetin_Student_Oct2024.pdf",
                           'Quercetin',
                           'Myricetin')
        
        return fisetin_all, quercetin_all
        


if __name__ == '__main__':
    
    data, candidates = load.prepare_data(False)
    stud_fis, stud_quer = run_analysis(data, candidates, method='student')
    #dun_fis, dun_quer = run_analysis(data, candidates, method='dunnett')

