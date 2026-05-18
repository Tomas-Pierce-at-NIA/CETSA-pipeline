# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 15:40:56 2024

@author: piercetf
"""

from sklearn.preprocessing import OneHotEncoder, MinMaxScaler
import seaborn
import polars as pl
NORMPROT = 'Normalized_FG_Quant'


class DataPreparer:
    
    __slots__ = ['onehot',
                 'category_names',
                 'scaler',
                 'aware_poly']
    
    def __init__(self, focused_data):
        self.onehot = OneHotEncoder(sparse_output=False).set_output(transform='polars')
        self.onehot.fit(focused_data.select(pl.col('Treatment')))
        self.category_names = self.onehot.categories_[0]
        self.scaler = MinMaxScaler().set_output(transform='polars')
        self.scaler.fit(focused_data[['Temperature']])

    
    def transform(self, focused_data, category1, category2):
        my_data = focused_data.filter(
            pl.col('Treatment').eq(category1) |
            pl.col('Treatment').eq(category2)
            )
        treatment_only = my_data.select(pl.col('Treatment'))
        categories = self.onehot.transform(treatment_only)
        
        temp_only = my_data.select(pl.col('Temperature'))
        scaled_temps = self.scaler.transform(temp_only)
        
        base_inputs = pl.concat([scaled_temps, categories], how='horizontal')
        
        interact_inputs = base_inputs.select(
            # including a shared bias competes with treatment category biases
            # we need it anyway for permutation testing reasons
            pl.lit(1).alias('bias'),
            pl.col('Temperature'),
            pl.col(f"Treatment_{category1}"),
            pl.col(f"Treatment_{category2}"),
            (pl.col('Temperature')*pl.col(f"Treatment_{category1}")).alias(f"Temperature Treatment_{category1}"),
            (pl.col('Temperature')*pl.col(f"Treatment_{category2}")).alias(f"Temperature Treatment_{category2}")
            )
        outputs = my_data[NORMPROT].to_numpy()
        prot_ids = my_data.select(pl.col('PG.ProteinAccessions'),
                                  pl.col('PG.Genes')
                                  )
        return interact_inputs, outputs, my_data['Treatment'], prot_ids

    
    def palette(self):
        colors = seaborn.color_palette('hls', len(self.category_names))
        pal = dict(zip(self.category_names, colors))
        return pal
    
    
    def null_model_cols_transform(self, focused_data, category1, category2):
        my_data = focused_data.filter(
            pl.col('Treatment').eq(category1) |
            pl.col('Treatment').eq(category2)
            )
        treatment_only = my_data.select(pl.col('Treatment'))
        categories = self.onehot.transform(treatment_only)
        
        temp_only = my_data.select(pl.col('Temperature'))
        scaled_temps = self.scaler.transform(temp_only)
        
        base_inputs = pl.concat([scaled_temps, categories], how='horizontal')
        
        interact_inputs = base_inputs.select(
            pl.lit(1),
            pl.col('Temperature')
            )
        prot_ids = my_data.select(pl.col('PG.ProteinAccessions'),
                                  pl.col('PG.Genes')
                                  )
        outputs = my_data[NORMPROT].to_numpy()
        return interact_inputs, outputs, my_data['Treatment'], prot_ids
    

if __name__ == '__main__':
    import cetsa_paths
    import load_monocyte_cetsa_data as load
    
    dpath = cetsa_paths.get_data_filepath()
    cpath = cetsa_paths.get_candidates_filepath()
    
    lzdat, lzcan = load.prep_data2(dpath, cpath)
    
    data = lzdat.collect()
    
    subdata = data.filter(
        pl.col('PG.Genes').eq(pl.lit('MTREX'))
        )
    
    dprep = DataPreparer(subdata)
    
    transformed = dprep.transform(subdata, 'Fisetin', 'DMSO')
    null_transformed = dprep.null_model_cols_transform(subdata, 'Fisetin', 'DMSO')
    
    