# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:49:28 2024

@author: piercetf
"""

from nparc_model import NPARCModel

import pandas as pd
cimport pandas as cpd

import numpy as np
cimport numpy as cnp
cnp.import_array()

cpdef enum TestProcedure:
    ECDF = 1
    GPD = 2

cdef class PermutationTestResult:
    cdef readonly TestProcedure proc_used
    cdef readonly double p_value, conf90_low, conf90_high
    
    def __init__(self, str proc, double p_value, double conf90_low, double conf90_high):
        if proc == "ECDF":
            self.proc_used = TestProcedure.ECDF
        else:
            self.proc_used = TestProcedure.GPD
        
        self.p_value = p_value
        self.conf90_low = conf90_low
        self.conf90_high = conf90_high
    
    def __str__(self):
        if self.proc == TestProcedure.ECDF:
            return f"""PermutationTestResult("ECDF", {self.p_value}, {self.conf90_low}, {self.conf90_high})"""
        else:
            return f"""PermutationTestResult("GPD", {self.p_value}, {self.conf90_low}, {self.conf90_high})"""



