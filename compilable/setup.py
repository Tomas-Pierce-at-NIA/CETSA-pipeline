# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 16:13:51 2024

@author: piercetf
"""

# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name='nparc_model2',
    ext_modules=cythonize("nparc_model2.pyx"),
    include_dirs=[np.get_include()]
)