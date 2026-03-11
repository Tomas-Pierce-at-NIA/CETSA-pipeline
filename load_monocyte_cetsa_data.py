# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:02:01 2024

@author: piercetf
"""

import pathlib
import polars as pl
import pandas


def load_data2(data_path: pathlib.Path) -> pl.LazyFrame:
    """loads data file from .tsv or .csv format file with
    certain columns expected to be present.
    Handles some aspects of intial data preparation achievable without
    reference to other files
    """
    # column names expected in data file. others are allowed,
    # but these are required
    DATA_COLUMNS = ['PG.ProteinAccessions',
                    'PG.Genes',
                    'R.Condition',
                    'R.Replicate',
                    'FG.Quantity']
    if data_path.suffix == '.tsv':
        table = pl.scan_csv(data_path,
                            separator='\t'
                            )
    elif data_path.suffix == '.csv':
        table = pl.scan_csv(data_path,
                            separator=','
                            )
    #extentions to support more formats would go here
    else:
        raise RuntimeError("format: {} not currently supported".format(data_path.suffix))
    
    table1 = table.select(pl.col(DATA_COLUMNS))
    # convert condition column to explicitly list temperature and treatment
    # as different columns
    table2 = table1.with_columns(
        pl.col('R.Condition').str.split(' ').list.to_struct(upper_bound=2)
        ).unnest('R.Condition')
    table3 = table2.rename({'field_0': 'Treatment', 'field_1': 'Temperature'})
    # numeric temperature helps later on
    table3 = table3.cast({'Temperature': pl.Int64})
    # mark values with missing gene name with special value
    # missing-value issues could occur later otherwise
    table4 = table3.select(pl.exclude('PG.Genes'), pl.col('PG.Genes').replace(None, pl.lit('GENEMISSING')))
    return table4


def load_candidates2(candidate_path :pathlib.Path) -> pl.LazyFrame:
    """Load candidates file from .csv or .tsv file
       Needed mostly to remove proteins identified from only 1 peptide.
       Handles some aspects of the preparation which don't reference other data
    """
    # these column names are required
    CANDIDATE_COLUMNS = ['Condition Numerator',
                         'Condition Denominator',
                         'ProteinGroups',
                         'ProteinNames',
                         'ProteinDescriptions',
                         'Genes',
                         'UniProtIds',
                         '# Unique Total Peptides']
    if candidate_path.suffix == '.tsv':
        table = pl.scan_csv(candidate_path,
                            separator='\t')
    elif candidate_path.suffix == '.csv':
        table = pl.scan_csv(candidate_path,
                            separator=',')
    else:
        msg = "format: {} not currenlty supported".format(candidate_path.suffix)
        raise RuntimeError(msg)

    table1 = table.select(pl.col(CANDIDATE_COLUMNS))
    # remove proteins identified by only 1 peptide
    table2 = table1.filter(pl.col('# Unique Total Peptides') > 1)
    # special characters could cause problems
    table3 = table2.rename({
        "# Unique Total Peptides" : "Number_Unique_Total_Peptides"
        })

    # getting temperature and treatment information from numerator and
    # denominator. may be a hold-over from the previous project that
    # used a different analysis strategy, which is also why there is
    # an explicity restriction of the candidates to same-temperature
    # comparisons.
    table4 = table3.with_columns(
        pl.col('Condition Numerator').str.split(' ').list.to_struct(upper_bound=2)
        ).unnest(pl.col('Condition Numerator'), separator="::")
    table5 = table4.rename({'Condition Numerator::field_0': 'Treatment_Numerator',
                            'Condition Numerator::field_1': 'Temperature_Numerator'
                            })
    table6 = table5.with_columns(
        pl.col('Condition Denominator').str.split(' ').list.to_struct(upper_bound=2)
        ).unnest(pl.col('Condition Denominator'), separator="::")
    table7 = table6.rename({'Condition Denominator::field_0': 'Treatment_Denominator',
                            'Condition Denominator::field_1': 'Temperature_Denominator'
                            })
    table8 = table7.filter(pl.col('Temperature_Numerator').eq(pl.col('Temperature_Denominator')))
    
    return table8

def filt_data_unipep(lz_data :pl.LazyFrame, lz_candidates :pl.LazyFrame) -> pl.LazyFrame:
    """Uses a candidates frame to remove proteins not in those candidates from the data.
       Utility is to remove proteins which are only identified by 1 peptide from the data,
       hence the name.
    """
    unique_uniprot = lz_candidates.select(pl.col('UniProtIds').unique())
    joined = lz_data.join(
        unique_uniprot,
        how='semi', # return x if present in y is a semi-join
        left_on=pl.col('PG.ProteinAccessions'),
        right_on=pl.col('UniProtIds')
        )
    return joined


def compute_total_prot2(lz_data :pl.LazyFrame) -> pl.LazyFrame:
    """The data file initially has quantifications for individual peptides.
       To conver that into a quantification for each protein, we simply
       sum the quantifications for peptides corresponding to a protein
       within each experimental condition.
    """
    grouped = lz_data.group_by(['PG.ProteinAccessions',
                                'PG.Genes',
                                'R.Replicate',
                                'Treatment',
                                'Temperature']
                               )
    return grouped.sum()


def norm_prot_mintemp2(lz_data :pl.LazyFrame, use_avg=False) -> pl.LazyFrame:
    """Normalize each protein against the level of that protein against the
       lowest temperature quantified. Will either normalize within replicates
       if use_avg == False or normalize against mean of lowest-temperature
       replicates if use_avg == True. 
    """
    mintemp = lz_data.filter(
        pl.col('Temperature') == pl.col('Temperature').min()
        )
    mintemp = mintemp.select(
        pl.exclude(['FG.Quantity', 'Temperature']),
        pl.col('FG.Quantity').alias('MinTempQuant')
        )
    if use_avg:
        gmintemp = mintemp.group_by(
            ['PG.ProteinAccessions',
             'PG.Genes',
             'Treatment']
            )
        mintemp2 = gmintemp.mean()
        mintemp2 = mintemp2.select(
            pl.col('PG.ProteinAccessions'),
            pl.col('PG.Genes'),
            pl.col('Treatment'),
            pl.col('MinTempQuant')
            )
        table = lz_data.join(
            mintemp2,
            on=['PG.ProteinAccessions',
                'PG.Genes',
                'Treatment'
                ],
            how='inner'
            )
        table2 = table.with_columns(
            Normalized_FG_Quant = pl.col('FG.Quantity') / pl.col('MinTempQuant')
            )
        return table2
    else:
        table = lz_data.join(
            mintemp,
            on=['PG.ProteinAccessions',
                'PG.Genes',
                'R.Replicate',
                'Treatment'],
            how='inner'
            )
        table2 = table.with_columns(
            Normalized_FG_Quant = pl.col('FG.Quantity') / pl.col('MinTempQuant')
            )
        return table2


def prep_data2(data_path :pathlib.Path, candidate_path :pathlib.Path) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """Loads data from data file and candidates file into memory as polars lazyframes,
       before doing initial data preparation steps.
    """
    lzdat = load_data2(data_path)
    lzcan = load_candidates2(candidate_path)
    lzdat = filt_data_unipep(lzdat, lzcan)
    lzdat = compute_total_prot2(lzdat)
    lzdat = lzdat.cast({'Temperature': pl.Int64})
    lzdat = norm_prot_mintemp2(lzdat, True)
    return lzdat, lzcan


if __name__ == '__main__':
    import cetsa_paths
    dpath = cetsa_paths.get_data_filepath()
    cpath = cetsa_paths.get_candidates_filepath()
    lzdat, lzcan = prep_data2(dpath, cpath)
    
    
    
    
