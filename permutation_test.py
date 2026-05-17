# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 22:02:59 2025

@author: piercetf
"""

# why not avoid cybersecurity implications
#import os
import numpy as np
from scipy import special, stats

# def _get_rand_seed():
#     b = os.urandom(8)
#     return int.from_bytes(b, 'little', signed=False)


class PermutationTest:
    
    def __init__(self, model, indata, outdata, treatdata, cat1, cat2, ident):
        
        self.model = model
        self.indata = indata
        self.outdata = outdata
        self.treatdata = treatdata
        self.cat1 = cat1
        self.cat2 = cat2
        self.ident = ident
        
        self.rng = np.random.default_rng()
        
        self.working = indata.to_numpy()
        
        self.y_max = np.max(outdata)
        
        
        self.__base_mse = None
        
    
    
    def permute(self):
        # randomly permute the category assignments in place in the working copy
        #self.rng.shuffle(self.working[:,2])
        # have to change because changed other aspects of the implementation
        working_np = np.copy(self.working)
        # what new order should be set
        permutation = self.rng.permutation(len(working_np))
        # which columns should be re-ordered
        permutable = working_np[:, [1,2]]
        # treatments are randomly reassinged by permuting
        permuted_treatments = permutable[permutation]
        # temperature is not reassigned
        temp_var = working_np[:, [0]]
        # temperature*treatment products must now be recomputed
        interact_vars = temp_var * permuted_treatments
        # recombine
        permuted_treatment_data = np.concat([temp_var,
                                             permuted_treatments,
                                             interact_vars],
                                            axis=1)
        self.working = permuted_treatment_data
        return permuted_treatment_data
        
    
    
    @property 
    def data_insufficient(self):
        totals = self.indata.sum()
        g1_cnt = totals.item(row=0, column=2)
        g2_cnt = totals.item(row=0, column=3)
        return g1_cnt == 0 or g2_cnt == 0

    
    @property
    def base_mse(self):
        if self.__base_mse is None:
            base_pred = self.model.predict(self.indata)
            self.__base_mse = PermutationTest.mse(base_pred, self.outdata)
            #self.__base_mse = self.scaled_logit_mse(base_pred, self.outdata)
        
        return self.__base_mse
    
    
    def data_insufficient_row(self):
        return (1.0, 0.0, 1.0, 'data-insufficient',
                -1, self.ident, self.cat1, self.cat2)
    
    
    def scaled_logit(self, data):
        scaled = data / (self.y_max * 1.1)
        return special.logit(scaled)
    
    @staticmethod
    def mse(predicted, actual):
        diffs = predicted - actual
        sqdif = diffs**2
        return sqdif.mean()
    
    def adj_logit(self, data):
        return np.log(data) - np.log1p(-data + self.y_max)
    
    def scaled_logit_mse(self, y_pred, y_act):
        scaled_logit_pred = self.adj_logit(y_pred)
        scaled_logit_act = self.adj_logit(y_act)
        
        pred_inf = np.isinf(scaled_logit_pred)
        act_inf = np.isinf(scaled_logit_act)
        
        pred_nan = np.isnan(scaled_logit_pred)
        act_nan = np.isnan(scaled_logit_act)
        
        problem = pred_inf | act_inf | pred_nan | act_nan
        okay = ~problem
        
        okay_act = scaled_logit_act[okay]
        okay_pred = scaled_logit_pred[okay]
        
        mse = np.mean((okay_act - okay_pred)**2)
        #breakpoint()
        return mse
    
    
    def permutation_test(self, n: int=50_000):
        #breakpoint()
        if self.data_insufficient:
            return self.data_insufficient_row()
        
        neg_base_mse = -self.base_mse
        better_count = 0
        mmts = np.empty((n,))
        #mmts = queue.PriorityQueue(n)
        for i in range(n):
            self.permute()
            perm_pred = self.model.predict(self.working)
            neg_mse = -self.scaled_logit_mse(perm_pred, self.outdata)
            # if mse <= base_mse:
            if neg_mse >= neg_base_mse:
                # count # of times model performance is not worse
                # under permutation of drug ID
                better_count += 1
            mmts[i] = neg_mse
        
        # if under treatment permutation, model perf is usually as good or better
        # then treatment is probably not important, so null hypothesis
        # so p-value should be high
        # if under treatment permutation, model perf is usually worse
        # then treatment is probably important, so alt hypothesis
        # so p-value should be low
        
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


if __name__ == '__main__':
    
    import cetsa_paths
    import load_monocyte_cetsa_data as load
    from data_prepper import DataPreparer
    import nparc_model as nparc
    
    data_fname = cetsa_paths.data_filename()
    can_fname = cetsa_paths.candidate_filename()
    lz_dat, lz_can = load.prep_data2(data_fname, can_fname)
    data = lz_dat.collect()
    
    dprep = DataPreparer(data)
    
    ins, outs, treats, _protids = dprep.transform(data, 'Quercetin', 'Myricetin')
    
    mask = (_protids['PG.Genes'] == 'CD38')
    
    ins_m = ins.filter(mask)
    outs_m = outs[mask]
    treats_m = treats.filter(mask)
    
    model = nparc.ScaledNPARCModel()
    model.fit(ins_m, outs_m)
    
    ptest = PermutationTest(model,
                            ins_m,
                            outs_m,
                            treats_m,
                            'Quercetin',
                            'Myricetin',
                            'CD38')
    
    res = ptest.permutation_test()
    
    print(res)

