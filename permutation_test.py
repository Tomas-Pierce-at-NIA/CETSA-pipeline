# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 22:02:59 2025

@author: piercetf
"""

import numpy as np
from scipy import special, stats
from sklearn import metrics
import time
import heapq

class PermutationTest:
    
    def __init__(self, model, indata, outdata, treatdata, cat1, cat2, ident):
        
        self.model = model
        self.indata = indata
        self.outdata = outdata
        self.treatdata = treatdata
        self.cat1 = cat1
        self.cat2 = cat2
        self.ident = ident
        
        self.seed = time.time_ns()
        bitgen = np.random.SFC64(self.seed)
        self.rng = np.random.default_rng(bitgen)
        
        self.working = np.copy(indata)
        
        self.y_max = np.max(outdata)
        
        
        self.__base_mse = None
        
    
    def permute(self):
        # randomly permute the category assignments in place in the working copy
        self.rng.shuffle(self.working[:,2])
        self.working[:,3] = 1 - self.working[:,2]
        self.working[:,4] = self.working[:,1] * self.working[:,2]
        self.working[:,5] = self.working[:,1] * self.working[:,3]
    
    
    @property 
    def data_insufficient(self):
        totals = self.indata.sum(axis=0)
        return totals[2] == 0 or totals[3] == 0
    
    @property
    def base_mse(self):
        if self.__base_mse is None:
            base_pred = self.model.predict(self.indata)
            self.__base_mse = self.scaled_logit_mse(base_pred, self.outdata)
        
        return self.__base_mse
    
    
    def data_insufficient_row(self):
        return (1.0, 0.0, 1.0, 'data-insufficient',
                -1, self.ident, self.cat1, self.cat2)
    
    
    def scaled_logit(self, data):
        scaled = data / self.y_max
        return special.logit(scaled)
    
    
    def scaled_logit_mse(self, y_pred, y_act):
        scaled_logit_pred = self.scaled_logit(y_pred)
        scaled_logit_act = self.scaled_logit(y_act)
        
        pred_inf = np.isinf(scaled_logit_pred)
        act_inf = np.isinf(scaled_logit_act)
        
        pred_nan = np.isnan(scaled_logit_pred)
        act_nan = np.isnan(scaled_logit_act)
        
        problem = pred_inf | act_inf | pred_nan | act_nan
        okay = ~problem
        
        okay_act = scaled_logit_act[okay]
        okay_pred = scaled_logit_pred[okay]
        
        #breakpoint()
        return metrics.mean_squared_error(okay_act, okay_pred)
    
    
    def fast_ecdf(self):
        neg_base_mse = -self.base_mse
        better_count = 0
        for i in range(500):
            self.permute()
            perm_pred = self.model.predict(self.working)
            neg_mse = -self.scaled_logit_mse(perm_pred, self.outdata)
            if neg_mse >= neg_base_mse:
                better_count += 1
        pval = (better_count + 1) / 501
        if pval > 0.1:
            pvar = pval * (1 - pval) / 500
            interval = stats.norm.interval(0.95, pval, pvar)
            return (pval, *interval, 'ecdf-faststop', 
                    -1, self.ident, self.cat1, self.cat2)
        else:
            return None
    
    def permutation_test(self, n: int=50_000):
        fast = self.fast_ecdf()
        if fast is not None:
            return fast
        
        neg_base_mse = -self.base_mse
        better_count = 0
        mmts = np.empty((n,))
        #mmts = queue.PriorityQueue(n)
        for i in range(n):
            self.permute()
            perm_pred = self.model.predict(self.working)
            neg_mse = -self.scaled_logit_mse(perm_pred, self.outdata)
            if neg_mse >= neg_base_mse:
                better_count += 1
            mmts[i] = neg_mse
        
        if better_count >= 10:
            pval = (better_count + 1) / (n + 1)
            pvar = pval * (1 - pval) / n
            interval = stats.norm.interval(0.95, pval, pvar)
            return (pval, *interval, 'ecdf', -1, 
                    self.ident, self.cat1, self.cat2)
        
        
        thresh_idx = -250
        largest_idcs = mmts.argpartition(thresh_idx-1)[thresh_idx-1:]
        largest = mmts[largest_idcs]
        sorted_largest = np.sort(largest)
        
        thresh = (sorted_largest[thresh_idx - 1] + sorted_largest[thresh_idx]) / 2
        
        exceeds = sorted_largest - thresh
        
        gpd_params = stats.genpareto.fit(exceeds)
        
        test = stats.ks_1samp(exceeds, stats.genpareto.cdf, args=gpd_params)
        
        iterations = 0
        
        
        while test.pvalue < 0.05:
            iterations += 1
            thresh_idx += 10
            if thresh_idx >= 0:
                p_value = (better_count + 1) / (n + 1)
                pvar = pval * (1 - p_value) / n
                interval = stats.norm.interval(0.95, p_value, pvar)
                return (p_value, *interval, 'ecdf-fallback',
                        iterations, self.ident, self.cat1, self.cat2)
            retained = sorted_largest[thresh_idx:]
            thresh = (sorted_largest[thresh_idx - 1] + sorted_largest[thresh_idx]) / 2
            exceeds = retained - thresh
            gpd_params = stats.genpareto.fit(exceeds)
            test = stats.ks_1samp(exceeds, stats.genpareto.cdf, args=gpd_params)
        
        sf = stats.genpareto.sf(neg_base_mse - thresh, *gpd_params)
        pval = (len(exceeds) / n) * sf
        
        count_low, count_high = stats.binom.interval(0.95, n, pval)
        low_bound = count_low / n
        high_bound = count_high / n
        
        return (pval, low_bound, high_bound, 'gcd', iterations,
                self.ident, self.cat1, self.cat2)
