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
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils import validation
from sklearn.exceptions import DataConversionWarning, NotFittedError


def _pos_expit(z :np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(-z))

def _neg_expit(z :np.ndarray) -> np.ndarray:
    exp = np.exp(z)
    return exp / (exp + 1)

def _expit(z :np.ndarray) -> np.ndarray:
    pos = z >=0
    neg = ~pos
    result = np.empty_like(z, dtype=float)
    result[pos] = _pos_expit(z[pos])
    result[neg] = _neg_expit(z[neg])
    return result


def _predict_nparc(p: float, w: np.ndarray, X: np.ndarray) -> np.ndarray:
    z = X @ w
    sigma = _expit(z)
    return (1-p)*sigma + p



def _predict(params: np.ndarray, X: np.ndarray) -> np.ndarray:
    p = params[0]
    weights = params[1:]
    return _predict_nparc(p, weights, X)



def _loss(params: np.ndarray, X: np.ndarray, y: np.ndarray, alpha=0.0) -> float:
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
    



def _fit_nparc(X: np.ndarray, y: np.ndarray, alpha=0.0) -> optimize.OptimizeResult:
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
    
    


# def _mse(params: np.ndarray, X: np.ndarray, y:np.ndarray) -> float:
#     preds = _predict(params, X)
#     return ((y - preds)**2).mean()



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


if __name__ == '__main__':
    
    from scipy import stats
    
    def mse(predicteds, actuals):
        diffs = actuals - predicteds
        sqdif = diffs**2
        sse = sqdif.sum()
        return sse / len(predicteds)
    
    def ymax_logit(y_act_max, y_values):
        return np.log(y_values) - np.log((y_act_max + 1) - y_values)
    
    
    def ymax_logit_mse(y_predicted, y_actuals):
        y_act_max = max(y_actuals)
        y_pred_logit = ymax_logit(y_act_max, y_predicted)
        y_act_logit = ymax_logit(y_act_max, y_actuals)
        return mse(y_pred_logit, y_act_logit)
    
    
    
    def _permutation_test(model, indata, outdata, treatdata, cat1, cat2, ident, rng, n=50_000):
        
        # need to denote tables where we don't have the raw data
        # required to even do permutation testing
        # and skip them to further cut down on execution time
        if indata[:,2].sum() == 0 or indata[:,3].sum() == 0:
            pvalue = 1.0
            interval = (0.0, 1.0)
            return pvalue, *interval, 'data-insufficient', -1, ident, cat1, cat2
    
        base_pred = model.predict(indata)
        neg_base_mse = -ymax_logit_mse(base_pred, outdata)
        permuted_indata = indata.copy()
    
        fast_better_count = 0
        fast_mmts = np.empty((500,))
        interact_permuted = np.empty_like(permuted_indata)
        
        for j in range(500):
            cat1_permuted = rng.permutation(indata[:,2])
            cat2_permuted = 1 - cat1_permuted
            permuted_indata[:,2] = cat1_permuted
            permuted_indata[:,3] = cat2_permuted
    
            interact_permuted[:,0:4] = permuted_indata[:,0:4]
            interact_permuted[:,4] = interact_permuted[:,1] * interact_permuted[:,2]
            interact_permuted[:,5] = interact_permuted[:,1] * interact_permuted[:,3]
            
    
            permuted_pred = model.predict(interact_permuted)
            neg_perm_mse = -ymax_logit_mse(permuted_pred, outdata)
            
            if neg_perm_mse >= neg_base_mse:
                fast_better_count += 1
            fast_mmts[j] = neg_perm_mse
        
        fast_pvalue = (fast_better_count + 1) / (501)
        if fast_pvalue > 0.1:
            p_variance = fast_pvalue * (1 - fast_pvalue) / 500
            interval = stats.norm.interval(0.95, fast_pvalue, p_variance)
            return fast_pvalue, *interval, 'ecdf-faststop', -1, ident, cat1, cat2
        
        better_count = 0
        mmts = np.empty((n,))
        for i in range(n):
            cat1_permuted = rng.permutation(indata[:,2])
            cat2_permuted = 1 - cat1_permuted
            permuted_indata[:,2] = cat1_permuted
            permuted_indata[:,3] = cat2_permuted
            
            interact_permuted[:,0:4] = permuted_indata[:,0:4]
            interact_permuted[:,4] = interact_permuted[:,1] * interact_permuted[:,2]
            interact_permuted[:,5] = interact_permuted[:,1] * interact_permuted[:,3]
            
            permuted_pred = model.predict(interact_permuted)
            neg_perm_mse = -ymax_logit_mse(permuted_pred, outdata)
            if neg_perm_mse >= neg_base_mse:
                better_count += 1
            mmts[i] = neg_perm_mse
    
        
        if better_count >= 10:
            p_value = (better_count + 1) / (n + 1)
            p_variance = p_value * (1 - p_value) / n
            interval = stats.norm.interval(0.95, p_value, p_variance)
            return p_value, *interval, 'ecdf', -1, ident, cat1, cat2
        
        else:
            sorted_measures = np.sort(mmts)
            thresh_idx = -250
            largest = sorted_measures[thresh_idx:]
            thresh = (sorted_measures[thresh_idx - 1] + sorted_measures[thresh_idx]) / 2
            exceeds = largest - thresh
            gpd_params = stats.genpareto.fit(exceeds)
            # use the 1-sample Kolmogorov-Smirnov test 
            # to check whether the exceedances are from a GPD
            test = stats.ks_1samp(exceeds, stats.genpareto.cdf, args=gpd_params)
            iterations = 0
            
            # we want to have a GPD, which is KS null hypothesis
            while test.pvalue < 0.05:
                iterations += 1
                thresh_idx += 10
                if thresh_idx >= 0: # if we cannot find a GPD, fallback to basic method
                    p_value = (better_count + 1) / (n + 1)
                    p_variance = p_value * (1 - p_value) / n
                    interval = stats.norm.interval(0.95, p_value, p_variance)
                    return (better_count + 1) / (n + 1), *interval, 'ecdf-fallback', iterations, ident, cat1, cat2
                largest = sorted_measures[thresh_idx:]
                thresh = (sorted_measures[thresh_idx - 1] + sorted_measures[thresh_idx]) / 2
                exceeds = largest - thresh
                gpd_params = stats.genpareto.fit(exceeds)
                test = stats.ks_1samp(exceeds, stats.genpareto.cdf, args=gpd_params)
            
            sf = stats.genpareto.sf(neg_base_mse - thresh, *gpd_params)
            p_value = (len(exceeds) / n) * sf
            
            # calculate confidence interval using binomial distribution
            count_low, count_high = stats.binom.interval(0.95, n, p_value)
            low_bound = count_low / n
            high_bound = count_high / n
            
            return p_value, low_bound, high_bound, 'gcd', iterations, ident, cat1, cat2
    


    import cProfile
    import load_monocyte_cetsa_data as load
    from data_prepper import DataPreparer
    import time
    
    profile = cProfile.Profile()
    profile.enable()
    
    try:
        data, can = load.prepare_data()
        start = time.time_ns()  
        i = 0
        for ident, table in data.groupby(by=['PG.ProteinAccessions', 'PG.Genes'],
                                         sort=False):
            
            if i >=  200:
                break
            dprep = DataPreparer(table)
            try:
                ins, outs, treats, _protids = dprep.transform(table, 'Quercetin', 'DMSO')
            except KeyError:
                i -= 1
                continue
            model = NPARCModel()
            model.fit(ins, outs)
            pds = model.predict(ins)
            rng = np.random.default_rng((start + i % 97))
            #profile.enable()
            # p = profile.runcall(_permutation_test, model, ins, outs, treats, 'Quercetin', 'DMSO', ident, rng, 10_000)
            #p = _permutation_test(model, ins, outs, treats, 'Quercetin', 'DMSO', ident, rng, 10_000)
            #profile.disable()
            i += 1
            p = _permutation_test(model, ins, outs, treats, 'Quercetin', 'DMSO', ident, rng, 10_000)
            print(i)
            print(p)
    finally:
        profile.disable()
        profile.dump_stats(r"C:\Users\piercetf\Documents\modeltest.profile")
    