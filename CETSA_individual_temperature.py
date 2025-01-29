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


def get_all_student_tests(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Takes as input the normalized data of a CETSA experiment
    such that the treatment, temperature, and relative soluble fraction
    of each sample is known.
    Runs T-tests at each temperature between each pair of conditions
    for every protein identified in the data,
    excepting comparisons between treatment conditions,
    where the treatments are determined from data and the controls
    are hard-coded to be DMSO and Myricetin.
    
    Returns a datatable describing the results of the T-test for each
    protein-temperature-treatment-control combination.
    """
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
                    if treatment == control:
                        continue
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



def calc_fisher_statistic(pvals):
    """
    Calculates test statistic for combinations of pvalues that is
    -2 * sum(log(p-values)) 
    as used by Brown's method for combination p-values
    for a collection of p-values
    """
    return -2 * np.sum(np.log(pvals))

def get_chi_mean(cov_mtx):
    """
    Calculate mean of a chi-squared variable 
    that would be consistent with the input
    covariance matrix
    as used in Brown's method for combination of non-independent p-values
    """
    return 2 * len(cov_mtx)

def get_chi_variance(cov_mtx):
    """
    Calculate variance of a chi-squared variable
    that would be consistent with the input 
    covariance matrix
    as used in Brown's method for combination of non-independent p-values
    """
    k_conds = len(cov_mtx)
    cov_total = 0
    for i in range(cov_mtx.shape[0]):
        for j in range(i + 1, cov_mtx.shape[1]):
            cov_total += cov_mtx[i,j]
    return (4 * k_conds) + (2 * cov_total)

def get_p_value(deg_free, scale_param, fisher_stat):
    """
    Calculate combination p-value using a Chi-squared variable
    with the specified parameters and input Fisher statistic
    ie -2*sum(log(p-values))
    
    Used to implement Brown's method for combining non-independent p-values
    """
    return stats.chi2.sf(fisher_stat, deg_free, scale=scale_param)

def collect_students(students_table: pandas.DataFrame) -> pandas.DataFrame:
    """
    Take as input a datatable describing the results of T-tests
    for
    protein-temperature-treatment-control
    and combine the p-values thereof by Brown's method to produce
    a table of p-values for
    protein-treatment-control
    cases.
    """
    
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




def graph_all_proteins(data, all_comps, focus, filename, treatment, control):
    """
    Graph all proteins where the comparison between treatment and control
    has been detected as significant ( adjusted p < 0.05).
    Displays line connecting averages at each temperature for treatment
    and control. Displays scatter dots for individual measurements.
    Also displays bar graph of individual-temperature p-values beneath
    x-axis and outputs adjusted p-value in bottom-right corner.
    """
    
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


def auc_all_proteins(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Calculates change in area under curve between each pair of
    treatment, control conditions.
    Currently considered suspect for not being consistent
    with same approach implemented via Excel.
    Computes area-under-curve using trapezoidal approximation.
    """
    
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


def run_analysis(data, candidates, datadir=None):
    """
    Responsible for overall running of individual-temperature based analysis.
    """
    
    if datadir is None:
        datadir = cetsa_paths.get_outdir()
        
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
    
    params = cetsa_paths.loadparams()
    
    conditions = all_info['Treatment'].unique()
    
    ns_controls = params['controls']['nonsenolytic']
    v_control = params['controls']['vehicle']
    controls = [v_control, *ns_controls]
    
    ns_genes = set()
    for ns_control in ns_controls:
        ns_table = sigtable.loc[sigtable['Treatment'] == ns_control, :]
        ns_genes.update(ns_table['PG.Genes'])
    
    for condition in conditions:
        cond_nospace = condition.replace(" ", "_")
        cond_table = sigtable.loc[sigtable['Treatment'] == condition, :]
        cond_unshare = cond_table.loc[~cond_table['PG.Genes'].isin(ns_genes),:]
        cond_fname = datadir / "ITA_un_{}_Jan2025.csv".format(cond_nospace)
        cond_unshare.to_csv(cond_fname)
        
        for control in controls:
            pair_tab = cond_table.loc[cond_table['Control'] == control, :]
            fname = datadir / "ITA_{}_v_{}_Jan2025.pdf".format(condition, 
                                                               control)
            graph_all_proteins(data,
                               students,
                               pair_tab,
                               fname,
                               condition,
                               control)
        
    
    return all_info
        


if __name__ == '__main__':
    
    data, candidates = load.prepare_data(False)
    results = run_analysis(data, candidates)

