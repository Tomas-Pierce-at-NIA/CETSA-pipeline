# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 15:40:56 2024

@author: piercetf
"""

from sklearn.preprocessing import OneHotEncoder, MinMaxScaler, PolynomialFeatures
import seaborn
NORMPROT = 'Normalized_FG_Quantity'


class DataPreparer:
    
    __slots__ = ['onehot',
                 'category_names',
                 'scaler',
                 'aware_poly']
    
    def __init__(self, focused_data):
        self.onehot = OneHotEncoder()
        self.onehot.fit(focused_data.loc[:, ['Treatment']])
        self.category_names = self.onehot.categories_[0]
        self.scaler = MinMaxScaler()
        self.scaler.fit(focused_data[['Temperature']])
        self.aware_poly = PolynomialFeatures(interaction_only=True,
                                             include_bias=True)

    
    def transform(self, focused_data, category1, category2):
        categories = self.onehot.transform(focused_data.loc[:, ['Treatment']]).toarray()
        focused_data.loc[:, self.category_names] = categories
        scaled_temps = self.scaler.transform(focused_data[['Temperature']])
        focused_data.loc[:, ['ScaledTemp']] = scaled_temps
        focal_data = focused_data[(focused_data[category1] == 1) | (focused_data[category2] == 1)]
        treatments = focal_data['Treatment']
        prot_ids = focal_data[['PG.ProteinAccessions', 'PG.Genes']]
        interact_vars = self.aware_poly.fit_transform(focal_data[['ScaledTemp',
                                                                  category1,
                                                                  category2]])
        
        col_count = 6 # bias + temp + cat1 + cat2 + temp*cat1 + temp*cat2
        #idx = np.argwhere(np.all(interact_vars[...,:] == 0, axis=0))
        #interact_data = np.delete(interact_vars, idx, axis=1)
        interact_data = interact_vars[:, :col_count]
        outputs = focal_data[NORMPROT].to_numpy()
        return interact_data, outputs, treatments, prot_ids
    
    def palette(self):
        colors = seaborn.color_palette('hls', len(self.category_names))
        pal = dict(zip(self.category_names, colors))
        return pal
    
    
    def null_model_cols_transform(self, focused_data, category1, category2):
        categories = self.onehot.transform(focused_data.loc[:, ['Treatment']]).toarray()
        focused_data.loc[:, self.category_names] = categories
        scaled_temps = self.scaler.transform(focused_data[['Temperature']])
        focused_data.loc[:, ['ScaledTemp']] = scaled_temps
        focal_data = focused_data[(focused_data[category1] == 1) | (focused_data[category2] == 1)]
        treatments = focal_data['Treatment']
        prot_ids = focal_data[['PG.ProteinAccessions', 'PG.Genes']]
        interact_vars = self.aware_poly.fit_transform(focal_data[['ScaledTemp',
                                                                  category1,
                                                                  category2]])
        col_count = 2 # bias + temp
        nh_data = interact_vars[:, :col_count]
        outputs = focal_data[NORMPROT].to_numpy()
        return nh_data, outputs, treatments, prot_ids