# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 22:02:59 2025

@author: piercetf
"""

import numpy as np
from scipy import special, stats
from sklearn import metrics
import unittest


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
    
    
    def baseline_perf(self):
        return self.current_perf(self.indata, self.outdata)
    
    def current_perf(self, inputs, target):
        preds = self.model.predict(inputs)
        return metrics.mean_squared_error(target, preds)
    
    def permute(self):
        working_np = np.copy(self.indata)
        treatment_idcols = working_np[:, [2,3]]
        # permute in place
        self.rng.shuffle(treatment_idcols, axis=0)
        # bias term
        bias = np.ones(shape=[working_np.shape[0], 1])
        # temperature term
        temp = working_np[:, [1]]
        # temperature treatment interaction product must be recomputed
        interact = temp * treatment_idcols
        # recombine items
        permuted_data = np.concat([bias,
                                   temp,
                                   treatment_idcols,
                                   interact],
                                  axis=1)
        return permuted_data
    
    
    def data_insufficient(self):
        #totals = self.indata.sum()
        input_arr = np.array(self.indata)
        totals = np.sum(input_arr, axis=0)
        g1_cnt = totals[2]
        g2_cnt = totals[3]
        #g1_cnt = totals.item(row=0, column=2)
        #g2_cnt = totals.item(row=0, column=3)
        return g1_cnt == 0 or g2_cnt == 0 or len(self.indata) == 0
    
    def data_insufficient_row(self):
        return (1.0, 
                0.0, 
                1.0, 
                'data-insufficient',
                -1, 
                self.ident, 
                self.cat1, 
                self.cat2
                )
    
    def permutation_test(self, n_iterations=50_000):
        if self.data_insufficient():
            return self.data_insufficient_row()
        
        baseline = self.baseline_perf()
        
        better_count = 0
        
        measurements = np.empty(shape=(n_iterations,))
        
        for i in range(n_iterations):
            permuted_inputs = self.permute()
            #permuted_preds = self.model.predict(permuted_inputs)
            permuted_perf = self.current_perf(permuted_inputs, self.outdata)
            
            measurements[i] = -permuted_perf
            
            #smaller is better
            if permuted_perf <= baseline:
                better_count += 1
            
            # this should be really rare
            if permuted_perf == baseline:
                pass
                #breakpoint()
            
        if better_count >= 10:
            pval = (better_count + 1) / (n_iterations + 1)
            pvar = pval * (1 - pval) / n_iterations
            p_low, p_high = stats.norm.interval(0.95, pval, pvar)
            return (pval, p_low, p_high, 'ecdf', -1, self.ident, self.cat1, self.cat2)
        
        
        thresh_idx = -250
        largest_idcs = measurements.argpartition(thresh_idx-1)[thresh_idx-1:]
        largest = measurements[largest_idcs]
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
                p_value = (better_count + 1) / (n_iterations + 1)
                pvar = p_value * (1 - p_value) / n_iterations
                interval = stats.norm.interval(0.95, p_value, pvar)
                return (p_value, *interval, 'ecdf-fallback',
                        iterations, self.ident, self.cat1, self.cat2)
            retained = sorted_largest[thresh_idx:]
            thresh = (sorted_largest[thresh_idx - 1] + sorted_largest[thresh_idx]) / 2
            exceeds = retained - thresh
            gpd_params = stats.genpareto.fit(exceeds)
            test = stats.ks_1samp(exceeds, stats.genpareto.cdf, args=gpd_params)
        
        sf = stats.genpareto.sf(-baseline - thresh, *gpd_params)
        pval = (len(exceeds) / n_iterations) * sf
        
        count_low, count_high = stats.binom.interval(0.95, n_iterations, pval)
        low_bound = count_low / n_iterations
        high_bound = count_high / n_iterations
        
        return (pval, 
                low_bound, 
                high_bound, 
                'gcd', 
                iterations,
                self.ident, 
                self.cat1, 
                self.cat2)


class PermutationTestUnitTests(unittest.TestCase):
    
    from sklearn import linear_model
    import nparc_model
    
    # generate linear data where the two groups have different intercepts
    @staticmethod
    def _gen_lin_data_helper_different():
        rng = np.random.default_rng()
        x = rng.normal(loc=2.0, scale=3.0, size=(50,1))
        groups = rng.binomial(1, 0.5, size=(50,1))
        inputs = np.concat((x, groups, (1 - groups)), axis=1)
        
        slope = 5
        intercept1 = -3
        intercept2 = 7
        params = np.array([slope, intercept1, intercept2])
        
        y = (params * inputs).sum(axis=1)
        
        return (inputs, y)
    
    # generate linear data where the two groups have the same intercepts
    @staticmethod
    def _gen_lin_data_helper_same():
        rng = np.random.default_rng()
        x = rng.normal(loc=2.0, scale=3.0, size=(50,1))
        groups = rng.binomial(1, 0.5, size=(50,1))
        inputs = np.concat((x, groups, (1 - groups)), axis=1)
        
        slope = 5
        
        intercept1 = -3
        intercept2 = -3
        params = np.array([slope, intercept1, intercept2])
        
        y = (params * inputs).sum(axis=1)
        
        return (inputs, y)
    
    # generate data suitable for checking whether permutation testing
    # is giving sensible p-values for NPARC data
    # case: identical effects for all groups for all parameters
    @classmethod
    def _gen_nparc_data_helper_same(cls):
        gen_params = np.array([0.01, -0.6, -0.3, 0.01, 0.01, 0.01, 0.01])
        rng = np.random.default_rng()
        biases = np.ones(shape=(50,1))
        temps = rng.uniform(size=(50,1))
        group1 = rng.binomial(1, 0.5, size=(50,1))
        group2 = 1 - group1
        columns = [biases,
                   temps,
                   group1,
                   group2,
                   group1*temps,
                   group2*temps]
        
        X = np.concat(columns, axis=1)
        y = cls.nparc_model._predict(gen_params, X)
        return (X, y)
    
    # generate data suitable for checking whether permutation testing
    # is giving sensible p-values for NPARC data
    # case: different offsets per group inside parameter position
    @classmethod
    def _gen_nparc_data_helper_dif_offset(cls):
        gen_params = np.array([0.01, -0.6, -0.3, -0.5, 7, 0.01, 0.01])
        rng = np.random.default_rng()
        biases = np.ones(shape=(50,1))
        temps = rng.uniform(size=(50,1))
        group1 = rng.binomial(1, 0.5, size=(50,1))
        group2 = 1 - group1
        columns = [biases,
                   temps,
                   group1,
                   group2,
                   group1*temps,
                   group2*temps]
        
        X = np.concat(columns, axis=1)
        y = cls.nparc_model._predict(gen_params, X)
        return (X, y)
    
    # generate data suitable for checking whether permutation testing
    # is giving sensible p-values for NPARC data
    # case: different temperature interactions per group inside parameter position
    @classmethod
    def _gen_nparc_data_helper_dif_effect(cls):
        gen_params = np.array([0.01, -0.6, -0.3, 0.0, 0.0, -2.0, 3.5])
        rng = np.random.default_rng()
        biases = np.ones(shape=(50,1))
        temps = rng.uniform(size=(50,1))
        group1 = rng.binomial(1, 0.5, size=(50,1))
        group2 = 1 - group1
        columns = [biases,
                   temps,
                   group1,
                   group2,
                   group1*temps,
                   group2*temps]
        
        X = np.concat(columns, axis=1)
        y = cls.nparc_model._predict(gen_params, X)
        return (X, y)
    
    
    # generate data suitable for checking whether permutation testing
    # is giving sensible p-values for NPARC data
    # case: different offsets and group temperature interactions
    @classmethod
    def _gen_nparc_data_helper_diff_all(cls):
        gen_params = np.array([0.01, -0.6, -0.3, -2.0, 4.0, -2.0, 3.5])
        rng = np.random.default_rng()
        biases = np.ones(shape=(50,1))
        temps = rng.uniform(size=(50,1))
        group1 = rng.binomial(1, 0.5, size=(50,1))
        group2 = 1 - group1
        columns = [biases,
                   temps,
                   group1,
                   group2,
                   group1*temps,
                   group2*temps]
        
        X = np.concat(columns, axis=1)
        y = cls.nparc_model._predict(gen_params, X)
        return (X, y)
    
    
    def test_nparc_same_not_significant(self):
        X, y = self._gen_nparc_data_helper_same()
        model = self.nparc_model.ScaledNPARCModel()
        model.fit(X,y)
        ptest = PermutationTest(model,
                                X,
                                y,
                                None,
                                "Treat1",
                                "Treat2",
                                "ident")
        res = ptest.permutation_test()
        p_val = res[0]
        self.assertGreater(p_val, 0.05)
    
    def test_nparc_diff_offset_significant(self):
        X, y = self._gen_nparc_data_helper_dif_offset()
        model = self.nparc_model.ScaledNPARCModel()
        model.fit(X,y)
        ptest = PermutationTest(model,
                                X,
                                y,
                                None,
                                "Treat1",
                                "Treat2",
                                "ident")
        res = ptest.permutation_test()
        p_val = res[0]
        self.assertLess(p_val, 0.05)
    
    def test_nparc_diff_slope_significant(self):
        X, y = self._gen_nparc_data_helper_dif_effect()
        model = self.nparc_model.ScaledNPARCModel()
        model.fit(X,y)
        ptest = PermutationTest(model,
                                X,
                                y,
                                None,
                                "Treat1",
                                "Treat2",
                                "ident")
        res = ptest.permutation_test()
        p_val = res[0]
        self.assertLess(p_val, 0.05)
    
    def test_nparc_diff_all_significant(self):
        X, y = self._gen_nparc_data_helper_diff_all()
        model = self.nparc_model.ScaledNPARCModel()
        model.fit(X,y)
        ptest = PermutationTest(model,
                                X,
                                y,
                                None,
                                "Treat1",
                                "Treat2",
                                "ident")
        res = ptest.permutation_test()
        p_val = res[0]
        self.assertLess(p_val, 0.05)
    
    
    def test_diff_intercept_significant(self):
        X, y = self._gen_lin_data_helper_different()
        model = self.linear_model.LinearRegression(fit_intercept=False)
        model.fit(X,y)
        rng = np.random.default_rng()
        
        base_score = model.score(X, y)
        
        x_group = X[:,[1,2]]
        
        better_count = 0
        
        for i in range(1_000):
            rng.shuffle(x_group, axis=0)
            X[:,[1,2]] = x_group
            perm_score = model.score(X, y)
            # larger is better in this case
            if perm_score >= base_score:
                better_count += 1
        
        p_val = (better_count + 1) / (1_000 + 1)
        print(p_val)
        self.assertLess(p_val, 0.05)
        
    
    def test_same_intercept_not_significant(self):
        X, y = self._gen_lin_data_helper_same()
        model = self.linear_model.LinearRegression(fit_intercept=False)
        model.fit(X,y)
        rng = np.random.default_rng()
        
        base_score = model.score(X, y)
        
        x_group = X[:,[1,2]]
        
        better_count = 0
        
        for i in range(1_000):
            rng.shuffle(x_group, axis=0)
            X[:,[1,2]] = x_group
            perm_score = model.score(X, y)
            # larger is better in this case
            if perm_score >= base_score:
                better_count += 1
        
        p_val = (better_count + 1) / (1_000 + 1)
        print(p_val)
        self.assertGreater(p_val, 0.05)
    
    
        
        
    


if __name__ == '__main__':
    
    unittest.main()
    
    import cetsa_paths
    import load_monocyte_cetsa_data as load
    from data_prepper import DataPreparer
    import nparc_model as nparc
    
    TEST_PROT = 'PHIP'
    CONTROL = 'DMSO'
    
    data_fname = cetsa_paths.data_filename()
    can_fname = cetsa_paths.candidate_filename()
    lz_dat, lz_can = load.prep_data2(data_fname, can_fname)
    data = lz_dat.collect()
    
    dprep = DataPreparer(data)
    
    ins, outs, treats, _protids = dprep.transform(data, 'Fisetin', CONTROL)
    
    mask = (_protids['PG.Genes'] == TEST_PROT)
    
    ins_m = ins.filter(mask)
    outs_m = outs[mask]
    treats_m = treats.filter(mask)
    
    model = nparc.ScaledNPARCModel(1e-8)
    model.fit(ins_m, outs_m)
    
    ptest = PermutationTest(model,
                            ins_m,
                            outs_m,
                            treats_m,
                            'Fisetin',
                            CONTROL,
                            TEST_PROT)
    
    permuted = ptest.permute()
    #print(permuted == ins_m)
    
    res = ptest.permutation_test(10_000)
    
    print(res)
    
    import polars as pl
    import seaborn
    from matplotlib import pyplot
    subdata = data.filter(pl.col('PG.Genes').eq(TEST_PROT))
    seaborn.scatterplot(
        subdata.filter(
            pl.col('Treatment').is_in(pl.lit([CONTROL, 'Fisetin']))
            ).to_pandas(), 
        x='Temperature', 
        y='Normalized_FG_Quant', 
        hue='Treatment')
    pyplot.show()
