# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:47:46 2024

@author: piercetf

Implement the individual temperature analysis method in Python
"""

NORMPROT = 'Normalized_FG_Quantity'

from scipy import stats
import pandas
import numpy as np
import seaborn
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot

import load_monocyte_cetsa_data as load
import cetsa_paths


def get_all_U_tests(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Take as input the normalized data of a CETSA experiment such that
    treatment, temp, relative soluble fraction known.
    Runs Mann-Whitney U-tests at each temperature betweeen each 
    treatment and each control for each protein in the data.
    
    Returns a datatable describing the result of the U-tests for each
    protein-temperature-treatment-control combination
    """
    
    pass

