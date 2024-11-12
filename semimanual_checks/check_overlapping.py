# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 09:19:00 2024

@author: piercetf
"""

import pandas

FIS_ITA = "C:\\Users\\piercetf\\Documents\\ITA_fisetin_Student_Oct2024.csv"

QUER_ITA = "C:\\Users\\piercetf\\Documents\\ITA_quercetin_Student_Oct2024.csv"

FIS_NPARC = "C:\\Users\\piercetf\\Documents\\nparc_unshared_fisetin_Oct2024.csv"

QUER_NPARC = "C:\\Users\\piercetf\\Documents\\nparc_unshared_quercetin_Oct2024.csv"

fis_ita = pandas.read_csv(FIS_ITA)

quer_ita = pandas.read_csv(QUER_ITA)

fis_nparc = pandas.read_csv(FIS_NPARC)

quer_nparc = pandas.read_csv(QUER_NPARC)

fis_ita_genes = set(fis_ita['PG.Genes'])

quer_ita_genes = set(quer_ita['PG.Genes'])

fis_nparc_genes = set(fis_nparc['PG.Genes'])

quer_nparc_genes = set(quer_nparc['PG.Genes'])

two_method_fis = fis_ita_genes.intersection(fis_nparc_genes)

two_method_quer = quer_ita_genes.intersection(quer_nparc_genes)

any_method_fis = fis_ita_genes.union(fis_nparc_genes)

any_method_quer = quer_ita_genes.union(quer_nparc_genes)

print("Detected by both methods for Fisetin:")
print(two_method_fis)
print("IoU: {}".format(len(two_method_fis) / len(any_method_fis)))

print("Detected by two methods for quercetin:")
print(two_method_quer)
print("IoU: {}".format(len(two_method_quer) / len(any_method_quer)))

print("Total unique detections")
print("Fisetin: {}".format(len(any_method_fis)))
print("Quercetin: {}".format(len(any_method_quer)))

with open("C:\\Users\\piercetf\\Documents\\summary.txt", "w") as handle:
    handle.write("Detected by both methods for Fisetin:\n")
    handle.write(str(two_method_fis) + "\n")
    handle.write("IoU: {}\n".format(len(two_method_fis) / len(any_method_fis)))
    handle.write("Detected by two methods for quercetin:\n")
    handle.write(str(two_method_quer) + "\n")
    handle.write("IoU: {}\n".format(len(two_method_quer) / len(any_method_quer)))
    handle.write("Total unique detections\n")
    handle.write("Fisetin: {}\n".format(len(any_method_fis)))
    handle.write("Quercetin: {}\n".format(len(any_method_quer)))
    

print("yeet")

