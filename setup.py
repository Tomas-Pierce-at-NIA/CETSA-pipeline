# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 13:57:07 2024

@author: piercetf
"""

from setuptools import setup
from Cython.Build import cythonize

setup(name="data_prepper2",
      ext_modules=cythonize("data_prepper2.pyx",
                            language_level = "3"))

