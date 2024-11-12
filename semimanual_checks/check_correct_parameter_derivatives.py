# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:15:26 2024

@author: piercetf
"""

from sympy import Symbol
import sympy

sympy.init_printing()

x = Symbol('X', real=True)
w = Symbol('w', real=True)

z = w * x

sigmoid = 1 / (1 + sympy.exp(-z))

p = Symbol('p', real=True)

model = sigmoid - (p * sigmoid) + p

y = Symbol('y', real=True)

n = Symbol('n', real=True, integer=True)

#k = Symbol('k', real=True)

logloss = -(1/n) * (y * sympy.log(model) + (1 - y) * sympy.log(1 - model))

logloss_pderiv = logloss.simplify().diff(p)

logloss_wderiv = logloss.simplify().diff(w)

log_model = sympy.log(model)

log_1minus_model = sympy.log(1 - model)