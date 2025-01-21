# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:02:01 2024

@author: piercetf
"""

import pandas
import seaborn
from matplotlib import pyplot

import cetsa_paths



def load_candidates(candidate_path = None):
    
    CANDIDATE_COLUMNS = ['Condition Numerator',
                         'Condition Denominator',
                         'ProteinGroups',
                         'ProteinNames',
                         'ProteinDescriptions',
                         'Genes',
                         'UniProtIds',
                         '# Unique Total Peptides']
    
    
    if candidate_path is None:
        cached_path = cetsa_paths.get_candidates_filepath(True)
        canon_path = cetsa_paths.get_candidates_filepath(False)
    
    else:
        cached_path = candidate_path
        canon_path = candidate_path
    
    try:
        candidates = pandas.read_csv(cached_path, 
                                     sep='\t',
                                     usecols=CANDIDATE_COLUMNS,
                                     engine='pyarrow')
    except Exception as e:
        print("Trying to load cached version failed, falling back to canonical path")
        print("Relevant error was: {}".format(e))
        candidates = pandas.read_csv(canon_path, 
                                     sep='\t',
                                     usecols=CANDIDATE_COLUMNS)
    
    return candidates


def load_basedata(data_path = None):
    
    DATA_COLUMNS = ['PG.ProteinAccessions',
                    'PG.Genes',
                    'R.Condition',
                    'R.Replicate',
                    'FG.Quantity']
    
    if data_path is None:
        cached_path = cetsa_paths.get_data_filepath(True)
        canon_path = cetsa_paths.get_data_filepath(False)
    
    else:
        cached_path = data_path
        canon_path = data_path
    
    try:
        basedata = pandas.read_csv(cached_path, 
                                   sep='\t',
                                   usecols=DATA_COLUMNS,
                                   engine='pyarrow')
    except Exception as e:
        print("Trying to load cached version failed, falling back to canonical path")
        print("Relevant error was: {}".format(e))
        basedata = pandas.read_csv(canon_path, 
                                   sep='\t',
                                   usecols=DATA_COLUMNS)
    
    # custom string is more accurate and avoids NaN-based pathology
    # later on
    basedata['PG.Genes'] = basedata['PG.Genes'].fillna('GENEMISSING')
    return basedata

def remove_unipeptides(data_table, candidate_table) -> (pandas.DataFrame, pandas.DataFrame):
    """ Remove proteins identified by only 1 unique peptide from both the
    data table and the candidates table.
    Returns the results in the same order as the parameters.
    """
    
    multipeptide_idx = candidate_table['# Unique Total Peptides'] > 1
    multipeptide_candidates = candidate_table.loc[multipeptide_idx, :]
    multipeptide_ids = set(multipeptide_candidates['UniProtIds'])
    multipeptide_data_idx = data_table['PG.ProteinAccessions'].isin(multipeptide_ids)
    multipeptide_data_table = data_table.loc[multipeptide_data_idx, :]
    return multipeptide_data_table, multipeptide_candidates

def remove_deprecated_columns(table):
    colnames = list(table.columns)
    for cname in colnames:
        if cname.startswith('[DEPRECATED]'):
            del table[cname]



def rename_special_columns(candidate_table):
    """take a prior table as input and produce a new table,
    where variables with special characters have been renamed to
    prevent the potential for problems.
    """
    return candidate_table.rename(columns={
        "# of Ratios" : "Number_of_Ratios",
        "% Change" : "Percent_Change",
        "# Unique Total Peptides" : "Number_Unique_Total_Peptides",
        "# Unique Total EG.Id" : "Number_Unique_Total_EGid"
        })



def load_data(data_path = None, candidate_path = None):
    """Load in data and perform filtering and data processing steps
    """
    basedata = load_basedata(data_path)
    candidates = load_candidates(candidate_path)
    
    multipep_data, multipep_candidates = remove_unipeptides(basedata, candidates)
    
    # unneeded because we only load in columns we actually use
    # remove_deprecated_columns(multipep_candidates)
    
    # remove this column by not loading it to begin with
    #del multipep_candidates['Valid'] # don't know, based on template
    
    # add aliases for variables to avoid problems with special characters
    multipep_candidates = rename_special_columns(multipep_candidates)
    

    
    treat_temp_num = multipep_candidates['Condition Numerator'].str.split(
        ' ',
        expand=True
        )
    
    multipep_candidates.loc[:, 'Treatment_Numerator'] = treat_temp_num[0]
    multipep_candidates.loc[:, 'Temperature_Numerator'] = treat_temp_num[1].astype(float)
    
    treat_temp_denom = multipep_candidates['Condition Denominator'].str.split(
        ' ',
        expand=True)
    multipep_candidates.loc[:, 'Treatment_Denominator'] = treat_temp_denom[0]
    multipep_candidates.loc[:, 'Temperature_Denominator'] = treat_temp_denom[1].astype(float)
    
    # only consider comparisons at the same temperature
    sametemp_multipep_candidates = multipep_candidates.loc[multipep_candidates.Temperature_Numerator == multipep_candidates.Temperature_Denominator,:].copy()
    
    
    temp_ints = pandas.to_numeric(sametemp_multipep_candidates["Temperature_Numerator"])
    sametemp_multipep_candidates['Temperature'] = temp_ints
    
    treat_temp_dat = multipep_data['R.Condition'].str.split(
        ' ',
        expand=True)
    
    # split out between substance and temperature 
    # use assign to intentionally make copy
    multipep_data = multipep_data.assign(
        Treatment=treat_temp_dat[0],
        Temperature=treat_temp_dat[1].astype(float)
        )
    
    return multipep_data, sametemp_multipep_candidates

def calc_total_protein_quantity(peptide_data):
    """Total the amount of material observed per-protein identification"""
    grouped = peptide_data.groupby(by=['PG.ProteinAccessions', 
                                     'PG.Genes', 
                                     'R.Replicate',
                                     'Treatment',
                                     'Temperature'])
    
    # I believe this is the total protein quantity
    summed = grouped[['FG.Quantity']].transform("sum")
    peptide_data.loc[:, 'Total_FG_Quantity'] = summed['FG.Quantity']
    


def display_counts(data):
    """display the number of tested temperature and number of detected proteins
    as a function of the replicate and treatment"""
    itemcounts = data.groupby(
        by=["R.Replicate", "Treatment"]
        ).nunique().loc[:,["PG.ProteinAccessions", "Temperature"]].reset_index()
    print(itemcounts)
    itemcounts['R.Replicate'] = itemcounts['R.Replicate'].astype(str)
    seaborn.barplot(data=itemcounts, 
                    x="Treatment", 
                    y="Temperature", 
                    hue="R.Replicate"
                    )
    pyplot.ylabel("Number of Temperatures tested")
    fig = pyplot.gcf()
    pyplot.show(block=False)
    pyplot.pause(2.0)
    pyplot.close(fig)
    
    seaborn.barplot(data=itemcounts,
                    x="Treatment",
                    y="PG.ProteinAccessions",
                    hue="R.Replicate"
                    )
    pyplot.ylabel("Number of detected proteins")
    
    fig = pyplot.gcf()
    pyplot.show(block=False)
    pyplot.pause(2.0)
    pyplot.close(fig)


def norm_protein_mintemp(data):
    """"Normalize the protein quantity against the protein quantity observed
    at the lowest tested temperature
    """
    
    # minimum temperature data
    mt_data = data[data['Temperature'] == min(data['Temperature'])]
    
    mt_groups = mt_data.groupby(by=['PG.Genes',
                                     'PG.ProteinAccessions',
                                     'Treatment'])
    
    mt_avgs = mt_groups.mean(numeric_only=True)
    
    mt_avgs.loc[:, 'Referent_Protein'] = mt_avgs['Total_FG_Quantity']
    del mt_avgs['R.Replicate']
    del mt_avgs['Temperature']
    del mt_avgs['FG.Quantity']
    del mt_avgs['Total_FG_Quantity']
    
    lowtemp = mt_avgs.reset_index()
    
    data2 = data.copy()
    del data2['FG.Quantity']
    data2 = data2.drop_duplicates()
    
    merged = data2.merge(lowtemp, 
                        how='left',
                        on=['PG.Genes',
                            'PG.ProteinAccessions',
                            'Treatment'],
                        validate="m:1")
    
    merged.loc[:, 'Normalized_FG_Quantity'] = merged['Total_FG_Quantity'] / merged['Referent_Protein']
    
    return merged


def prepare_data(display=False, data_path = None, candidate_path = None):
    """General data loading routine, including filtering and normalization,
    loads both data (normalized) and candidates (filtered).
    Primary point of interaction with the loading module, all other components
    are implementation.
    """
    filtered_data, filtered_candidates = load_data(data_path, candidate_path)
    calc_total_protein_quantity(filtered_data)
    #breakpoint()
    if display:
        display_counts(filtered_data)
    
    normalized_data = norm_protein_mintemp(filtered_data)
    
    normalized_data = normalized_data.dropna(subset=['Temperature', 'Normalized_FG_Quantity'])
    
    return normalized_data, filtered_candidates



if __name__ == '__main__':
    d,c = prepare_data()
    