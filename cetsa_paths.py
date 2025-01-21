# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:02:19 2024

@author: piercetf
"""

from pathlib import Path
import platform
import os

def get_userdir() -> Path:
    if platform.system() == 'Windows':
        userdir = os.environ['USERPROFILE']
        return Path(userdir)
    elif platform.system() == 'Linux':
        userdir = os.environ['HOME']
        datadir = userdir.replace('home', 'data')
        return Path(datadir)
    else:
        failsys()


def get_outdir() -> Path:
    if platform.system() == 'Windows':
        userdir = get_userdir()
        documents = userdir / 'Documents'
        return documents
    elif platform.system() == 'Linux':
        datadir = get_userdir() / '2024_CETSA_MS' / 'outdata'
        return datadir
    else:
        failsys()


def failsys():
    system = platform.system()
    raise NotImplementedError("{} system is not supported".format(system))

def get_candidates_filepath(cached=False) -> Path:
    
    if platform.system() == 'Windows':
        if cached:
            cachepath = r"C:\Users\piercetf\Projects\CachedCETSAData\Candidates.tsv"
            return Path(cachepath)
        else:
            canonpath = r'T:\TGB\LSS\TGU\Users\Tomas\2024_CETSA_MS\Monocyte_CETSA_Statistical_Analysis\CETSA_ind_temp_analysis_starting_files\Candidates.tsv'
            return Path(canonpath)
    
    elif platform.system() == 'Linux':
        datadir = get_userdir()
        candidate = datadir / "2024_CETSA_MS" / "indata" / "Candidates.tsv"
        return candidate
    else:
        failsys()
    

def get_data_filepath(cached=False) -> Path:
    if platform.system() == 'Windows':
        if cached:
            cachedpath = r"C:\Users\piercetf\Projects\CachedCETSAData\Complete CETSA analysis w F-37-4_Report_Delaney_Default (Normal).tsv"
            return Path(cachedpath)
        else:
            canonpath = r"T:\TGB\LSS\TGU\Users\Tomas\2024_CETSA_MS\Monocyte_CETSA_Statistical_Analysis\CETSA_ind_temp_analysis_starting_files\Complete CETSA analysis w F-37-4_Report_Delaney_Default (Normal).tsv"
            return Path(canonpath)
    elif platform.system() == 'Linux':
        datadir = get_userdir()
        indata_dir = datadir / '2024_CETSA_MS' / 'indata'
        indata_name = indata_dir / 'Complete_CETSA_analysis_w_F_37_4_Report_Delaney_Default_Normal.tsv'
        return indata_name
    else:
        failsys()


def get_logging_path(logname) -> Path:
    userdir = get_userdir()
    if platform.system() == 'Windows':
        logdir = userdir / 'Documents'
        logfilename = logdir / logname
        return logfilename
    elif platform.system() == 'Linux':
        return userdir / logname
    else:
        failsys()



