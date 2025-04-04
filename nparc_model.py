# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 11:58:30 2024

@author: piercetf
"""

import warnings
import numpy as np
from scipy import optimize
from scipy import special
from sklearn.base import BaseEstimator
from sklearn.utils import validation
from sklearn.exceptions import DataConversionWarning, NotFittedError



def _predict_nparc(p: float, w: np.ndarray, X: np.ndarray) -> np.ndarray:
    z = X @ w
    #sigma = _expit(z)
    sigma = special.expit(z)
    return (1-p)*sigma + p



def _predict(params: np.ndarray, X: np.ndarray) -> np.ndarray:
    p = params[0]
    weights = params[1:]
    return _predict_nparc(p, weights, X)



def _loss(params: np.ndarray, X: np.ndarray, y: np.ndarray, alpha=0.0) -> float:
    p = params[0]
    w = np.array(params[1:])
    z = X @ w
    exp = np.exp(-z)
    term1 = (-y + 1) * (z - np.log1p(-p))
    term2 = special.xlog1py(y, p*exp)
    term3 = np.log1p(exp)
    reg = alpha*(w**2).sum()
    s = term1 + term2 + term3 + reg
    return s.mean()
    

def _w_deriv_nparc(p: float, w: np.ndarray, X: np.ndarray, y: np.ndarray) -> np.ndarray:
    w_derivs = np.empty_like(w)
    for i in range(len(w_derivs)-1):
        x_iplus = X[:, i+1]
        z = x_iplus * w[i]
        exp = np.exp(-z)
        pexp = p * exp
        denom = 1 + pexp
        ratio = (1 + pexp) / (1 + exp)
        num = (-y + ratio) * x_iplus
        w_derivs_i = num / denom
        w_derivs[i] = w_derivs_i.mean()
    return w_derivs

def _p_deriv_nparc(p: float, w: np.ndarray, X: np.ndarray, y: np.ndarray) -> float:
    z = X @ w
    exp = np.exp(-z)
    pexp = p*exp
    one_minp_exp = (1.0-p)*exp
    ratio1 = (1+pexp)/(1+exp)
    left = (-y/ratio1) + 1
    right = exp/(one_minp_exp)
    p_derivs = left*right
    p_deriv = p_derivs.mean()
    return p_deriv


def _jac_nparc(p: float, w: np.ndarray, X: np.ndarray, y: np.ndarray, alpha=0.0) -> np.ndarray:
    jac = np.empty((w.shape[0]+1,))
    p_deriv = _p_deriv_nparc(p, w, X, y)
    w_derivs = _w_deriv_nparc(p, w, X, y)
    jac[0] = p_deriv + 2*alpha*p
    jac[1:] = w_derivs + 2*alpha*w
    return jac


def _gradient(params: np.ndarray, X: np.ndarray, y: np.ndarray, alpha=0.0) -> np.ndarray:
    p = params[0]
    w = np.array(params[1:])
    return _jac_nparc(p, w, X, y, alpha)


def _fit_nparc(X: np.ndarray, y: np.ndarray, alpha=0.0) -> optimize.OptimizeResult:
    pshape = (X.shape[1] + 1,)
    init_params = np.zeros(pshape) - 0.01
    # may be able to avoid pathology in p-derivative
    init_params[0] = 0.0001 # p should be almost zero in most cases
    
    bounds = [(-np.inf, np.inf)] * len(init_params)
    # may be able to avoid pathology in p-derivative
    bounds[0] = (0, 1)
    
    bfgs_res = optimize.minimize(_loss,
                                   init_params,
                                   (X, y, alpha),
                                   'L-BFGS-B',
                                   bounds=bounds)#,
                                   #jac=_gradient)
                                  # jac=_gradient)
    return bfgs_res




class ModelBase(BaseEstimator):
    
    def __init__(self, alpha :float = 0.0):
        self.alpha = alpha
    
    def fit(self, X :np.ndarray, y :np.ndarray):
        if y is None:
            raise ValueError("requires y to be passed, but the target y is None")
        validation.check_array(X, accept_sparse=False)
        X = np.array(X)
        y = np.array(y)
        if issubclass(X.dtype.type, np.complexfloating) or issubclass(y.dtype.type, np.complexfloating):
            raise ValueError("Complex data not supported")
        X = X.astype(float)
        y = y.astype(float)
        if np.isnan(X).any() or np.isnan(y).any():
            raise ValueError("NaN")
        if np.isinf(X).any() or np.isinf(y).any():
            raise ValueError("inf")
        if len(y.shape) > 1:
            warnings.warn(
                'A column-vector y was passed when a 1d array was expected',
                DataConversionWarning)
    
    def predict(self, X :np.ndarray) -> np.ndarray:
        if not self.__sklearn_is_fitted__():
            raise NotFittedError("Not fitted yet")
        if np.isnan(X).any():
            raise ValueError("NaN")
        if np.isinf(X).any():
            raise ValueError("inf")
        if len(X.shape) < 2:
            raise ValueError("Reshape your data")
    
    def __sklearn_is_fitted__(self):
        return hasattr(self, "params_")


class NPARCModel(ModelBase):
    
    def __init__(self, alpha :float = 0.0):
        super().__init__(alpha)
        self.__resid_response = None
    
    def fit(self, X :np.ndarray, y :np.ndarray):
        super().fit(X, y)
        fit_res = _fit_nparc(X, y, self.alpha)
        self.params_ = fit_res.x
        self.fit_success_ = fit_res.success
        self.fit_message_ = fit_res.message
        self.X_ = X
        self.y_ = y
        self.n_features_in_ = X.shape[1]
        hess_inv_op = fit_res.hess_inv
        identity = np.identity(hess_inv_op.shape[0])
        self.inv_hess_ = hess_inv_op @ identity
        return self
    
    def predict(self, X: np.ndarray) -> np.ndarray:
        super().predict(X)
        return _predict(self.params_, X)
    
    def __str__(self):
        cname = type(self).__name__
        
        if self.__sklearn_is_fitted__():
            paramstr = str(self.params_)
            return f"{cname}(w={paramstr};α={self.alpha})".replace('\n', ' ')
        else:
            return f"{cname}(not fitted)"
    
    @property
    def treatment1_decreasing(self):
        if self.__sklearn_is_fitted__():
            sumparam = self.params_[2] + self.params_[5]
            return sumparam < 0
        else:
            raise NotFittedError("tried to check if treatment 1 decreasing before fit")
    
    @property
    def treatment2_decreasing(self):
        if self.__sklearn_is_fitted__():
            sumparam = self.params_[2] + self.params_[6]
            return sumparam < 0
        else:
            raise NotFittedError("tried to check if treatment 2 decreasing before fit")



class ScaledNPARCModel(NPARCModel):
    
    def fit(self, X :np.ndarray, y :np.ndarray):
        self.y_max_ = np.max(y)
        self.small_y_ = y / self.y_max_
        return super().fit(X, self.small_y_)
    
    def predict(self, X: np.ndarray) -> np.ndarray:
        small_preds = super().predict(X)
        return small_preds * self.y_max_
    

        


if __name__ == '__main__':
    
    import permutation_test as ptest
    
    import load_monocyte_cetsa_data as load
    from data_prepper import DataPreparer
    import time
    
    data, can = load.prepare_data()
    groups = data.groupby(by=['PG.ProteinAccessions', 'PG.Genes'],
                                     sort=False)

        
    start = time.time_ns()  
    i = 0
    for ident, table in groups:
        
        if i >=  400:
            break
        dprep = DataPreparer(table)
        try:
            ins, outs, treats, _protids = dprep.transform(table, 'Quercetin', 'DMSO')
        except KeyError:
            i -= 1
            continue
        if np.max(outs) > 1:
            #breakpoint()
            pass
            
        #model = NPARCModel()
        model = ScaledNPARCModel(alpha=1e-6)
        model.fit(ins, outs)
        pds = model.predict(ins)
        #rng = np.random.default_rng((start + i % 97))
        i += 1
        perm_test = ptest.PermutationTest(model,
                                          ins,
                                          outs,
                                          treats,
                                          'Quercetin',
                                          'DMSO',
                                          ident)
        test = perm_test.permutation_test()
        
        if i > 7:
            break
        
        print(i)
        print(test)

