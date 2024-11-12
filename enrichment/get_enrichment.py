# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 13:14:18 2024

@author: piercetf
"""

import pandas
import requests

import load_monocyte_cetsa_data as load
import cetsa_paths

def get_background():
    rawdata, rawcan = load.prepare_data()
    all_genes = set(rawcan['Genes'])
    return all_genes

mydir = cetsa_paths.get_outdir()

nparc_fis_name = mydir / "nparc_unshared_fisetin_Oct2024.csv"
nparc_quer_name = mydir / "nparc_unshared_quercetin_Oct2024.csv"

ita_fis_name = mydir / "ITA_fisetin_Student_Oct2024.csv"
ita_quer_name = mydir / "ITA_quercetin_Student_Oct2024.csv"

nparc_fisetin = pandas.read_csv(nparc_fis_name)
nparc_quercetin = pandas.read_csv(nparc_quer_name)

ita_fisetin = pandas.read_csv(ita_fis_name)
ita_quercetin = pandas.read_csv(ita_quer_name)

bg_genes = get_background()

nparc_fis_genes = set(nparc_fisetin['PG.Genes'])
nparc_quer_genes = set(nparc_quercetin['PG.Genes'])

ita_fis_genes = set(ita_fisetin['PG.Genes'])
ita_quer_genes = set(ita_quercetin['PG.Genes'])

# See
# https://maayanlab.cloud/Enrichr/help#api


BASE_URL = "https://maayanlab.cloud/speedrichr"

description = "cetsa data"

res = requests.post(
    BASE_URL + "/api/addList",
    files=dict(list=(None, "\n".join(nparc_fis_genes)))
    )

print(res.ok)

