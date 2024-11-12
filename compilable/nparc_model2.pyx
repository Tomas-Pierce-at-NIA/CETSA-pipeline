# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:15:23 2024

@author: piercetf
"""

import warnings
import numpy as np
from scipy import optimize
from scipy import special
from sklearn.base import BaseEstimator
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils import validation
from sklearn.exceptions import DataConversionWarning, NotFittedError

cimport numpy as np

def _pos_expit(np.ndarray[np.float64_t, ndim=1] z) -> np.ndarray:
    return 1 / (1 + np.exp(-z))


def _neg_expit(np.ndarray[np.float64_t, ndim=1] z) -> np.ndarray:
    exp = np.exp(z)
    return exp / (exp + 1)


def _expit(np.ndarray[np.float64_t, ndim=1] z) -> np.ndarray:
    pos = z >=0
    neg = ~pos
    result = np.empty_like(z, dtype=float)
    result[pos] = _pos_expit(z[pos])
    result[neg] = _neg_expit(z[neg])
    return result


def _predict_nparc(p: float, 
                   np.ndarray[np.float64_t, ndim=1] w, 
                   np.ndarray[np.float64_t, ndim=2] X) -> np.ndarray:
    z = X @ w
    sigma = _expit(z)
    return (1-p)*sigma + p


def _predict(np.ndarray[np.float64_t, ndim=1] params, 
             np.ndarray[np.float64_t, ndim=2] X: np.ndarray) -> np.ndarray:
    p = params[0]
    weights = params[1:]
    return _predict_nparc(p, weights, X)



def _loss(np.ndarray[np.float64_t, ndim=1] params, 
          np.ndarray[np.float64_t, ndim=2] X, 
          np.ndarray[np.float64_t, ndim=1] y, 
          alpha=0.0) -> float:
    p = params[0]
    w = np.array(params[1:])
    z = X @ w
    term1 = special.xlog1py(y,-p)
    term2 = -special.xlogy(y, p + np.exp(z))
    term3 = -np.log1p(-p)
    term4 = np.log1p(np.exp(z))
    reg = alpha*(w**2).sum()
    s = term1 + term2 + term3 + term4 + reg
    return s.mean()
    


def _fit_nparc(np.ndarray[np.float64_t, ndim=2] X, 
               np.ndarray[np.float64_t, ndim=1] y, 
               alpha=0.0) -> tuple:
    pshape = (X.shape[1] + 1,)
    init_params = np.zeros(pshape) - 0.01
    # may be able to avoid pathology in p-derivative
    init_params[0] = 0.0001 # p should be almost zero in most cases
    
    bounds = [(-np.inf, np.inf)] * len(init_params)
    # may be able to avoid pathology in p-derivative
    bounds[0] = (1e-16, 1 - 1e-16)
    
    nelder_res = optimize.minimize(_loss,
                                   init_params,
                                   (X, y, alpha),
                                   'Nelder-Mead',
                                   bounds=bounds)
    return nelder_res
    
    


def _mse(np.ndarray[np.float64_t, ndim=1] params, 
         np.ndarray[np.float64_t, ndim=2] X, 
         np.ndarray[np.float64_t, ndim=1] y) -> float:
    preds = _predict(params, X)
    return ((y - preds)**2).mean()



class NPARCModel(BaseEstimator):
    """Model which fits a sigmoid of the form called for by Childs' NPARC
    """
    
    def __init__(self, alpha :float =0.0):
        self.alpha=alpha
    
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
        if np.isnan(X).any():
            raise ValueError("NaN")
        if np.isinf(X).any():
            raise ValueError("inf")
        if np.isinf(y).any():
            raise ValueError("inf")
        if np.isnan(y).any():
            raise ValueError("NaN")
        if len(y.shape) > 1:
            warnings.warn(
                'A column-vector y was passed when a 1d array was expected',
                DataConversionWarning)
        fit_res = _fit_nparc(X, y, self.alpha)
        self.params_ = fit_res.x
        self.fit_success_ = fit_res.success
        self.fit_message_ = fit_res.message
        self.X_ = X
        self.y_ = y
        self.n_features_in_ = X.shape[1]
        self.resid_response_ = y - self.predict(X)
        return self
    
    def predict(self, X: np.ndarray) -> np.ndarray:
        if not self.__sklearn_is_fitted__():
            raise NotFittedError("Not fitted yet")
        
        X = np.array(X)
        X = X.astype(float)
        
        if len(X.shape) < 2 or X.shape[1] != self.n_features_in_:
            raise ValueError("Reshape your data")
        
        if np.isnan(X).any():
            raise ValueError("NaN")
        if np.isinf(X).any():
            raise ValueError("inf")
        X = X.astype(float)
        return _predict(self.params_, X)
    
    def __sklearn_is_fitted__(self):
        return hasattr(self, 'params_')
    
    @property
    def converged(self) -> bool:
        if self.__sklearn_is_fitted__():
            return self.fit_success_
        else:
            return False
    
    @property
    def resid_response(self) -> bool:
        if self.__sklearn_is_fitted__():
            return self.resid_response_
        else:
            raise ValueError("Not fitted yet")

def test1():
    return check_estimator(NPARCModel())

