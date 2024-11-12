# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 11:22:12 2024

@author: piercetf
"""

from load_monocyte_cetsa_data import prepare_data

d, c = prepare_data()

aifm1 = d[d['PG.Genes'] == 'AIFM1']

prot = aifm1.assign(idx = aifm1.groupby('R.Replicate').cumcount())

fis = prot[prot['Treatment'] == 'Fisetin']
qur = prot[prot['Treatment'] == 'Quercetin']
dmso = prot[prot['Treatment'] == 'DMSO']
myr = prot[prot['Treatment'] == 'Myricetin']

fis_tot_pivot = fis.pivot(index='Temperature', columns='R.Replicate',values='Total_FG_Quantity').reset_index()

quer_tot_pivot = qur.pivot(index='Temperature', columns='R.Replicate',values='Total_FG_Quantity').reset_index()

dmso_tot_pivot = dmso.pivot(index='Temperature', columns='R.Replicate',values='Total_FG_Quantity').reset_index()

myr_tot_pivot = myr.pivot(index='Temperature', columns='R.Replicate',values='Total_FG_Quantity').reset_index()
