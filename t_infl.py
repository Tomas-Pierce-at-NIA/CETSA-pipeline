# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 15:51:26 2024

@author: piercetf
"""


from scipy import optimize
import sympy
from nparc_model import NPARCModel
from typing import Callable

# analytic derivative is computationally expensive,
# so getting it done once and reusing is more efficient

# I have verified this empirically,
# doing this changed the runtime from ~ 1 minute
# to ~ 1.3 __seconds__ for the test sequence at the bottom
# in if __name__ == '__main__'

t = sympy.Symbol('t', real=True)
c1 = sympy.Symbol('c1', real=True)
c2 = sympy.Symbol('c2', real=True)
p = sympy.Symbol('p', real=True)
b0 = sympy.Symbol('b0', real=True)
w_t = sympy.Symbol('w_t', real=True)
w_c1 = sympy.Symbol('w_c1', real=True)
w_c2 = sympy.Symbol('w_c2', real=True)
w_tc1 = sympy.Symbol('w_tc1', real=True)
w_tc2 = sympy.Symbol('w_tc2', real=True)
_z = b0 + (w_t * t) + (w_c1 * c1) + (w_c2 * c2) + (w_tc1 * c1 * t) + (w_tc2 * c2 * t)
expit = 1 / (1 + sympy.exp(-_z))
analytic_form = expit - (p * expit) + p
analytic_form = analytic_form.simplify()
analytic_derivative = analytic_form.diff(t)

def analytic_deriv(model :NPARCModel, treatment :int) -> Callable[[float], float]:
    
    params = model.params_
    
    substitutes = [(p, params[0]),
                         (b0, params[1]),
                         (w_t, params[2]),
                         (w_c1, params[3]),
                         (w_c2, params[4]),
                         (w_tc1, params[5]),
                         (w_tc2, params[6])]
    
    if treatment == 1:
        substitutes.extend([(c1, 1),
                    (c2, 0)])
    elif treatment == 2:
        substitutes.extend([(c1, 0),
                    (c2, 1)])
    
    withnum_deriv = analytic_derivative.subs(substitutes)
    
    as_func = sympy.lambdify(t, withnum_deriv)
    
    return as_func
    

def find_inflection(analytic_deriv :Callable[[float], float]) -> float:
    res = optimize.minimize_scalar(analytic_deriv, bounds=[0,1])
    return res.x

def get_T_inflection(model :NPARCModel, treatment: int) -> float:
    deriv = analytic_deriv(model, treatment)
    return find_inflection(deriv)

def get_T1_inflection(model :NPARCModel) -> float:
    return get_T_inflection(model, 1)

def get_T2_inflection(model :NPARCModel) -> float:
    return get_T_inflection(model, 2)


if __name__ == '__main__':
    import cProfile
    import load_monocyte_cetsa_data as load
    from data_prepper import DataPreparer
    from nparc_model import NPARCModel
    
    data, can = load.prepare_data()
    profile = cProfile.Profile()
    #profile.enable()
    i = 0
    try:
        for gene_stuff, table in data.groupby(by=['PG.Genes']):
            if i > 200:
                break
            try:
                dprep = DataPreparer(table)
                inputs, outputs, treats = dprep.transform(table, 'Fisetin', 'DMSO')
                model = NPARCModel()
                model.fit(inputs, outputs)
            except ValueError:
                continue
            except KeyError:
                continue
            profile.runcall(get_T1_inflection, model)
            profile.runcall(get_T2_inflection, model)
            i += 1
    finally:
        print("preparing to dump")
        profile.dump_stats(r"C:\Users\piercetf\Documents\t_infl.profile")
        print("dumped")

