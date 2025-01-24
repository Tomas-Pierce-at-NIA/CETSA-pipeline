# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:02:19 2024

@author: piercetf
"""

import pathlib
import toml

def loadparams() -> dict:
    configpath = pathlib.Path(__file__).with_name("cetsa_config.toml")
    with configpath.open("r") as config:
        params = toml.load(config)
    return params

def get_outdir() -> pathlib.Path:
    params = loadparams()
    return pathlib.Path(params['outputs']['outdir'])


def get_candidates_filepath(cached=False) -> pathlib.Path:
    params = loadparams()
    if cached:
        return pathlib.Path(params['inputs']['cached']['candidates'])
    else:
        return pathlib.Path(params['inputs']['uncached']['candidates'])
    

def get_data_filepath(cached=False) -> pathlib.Path:
    params = loadparams()
    if cached:
        return pathlib.Path(params['inputs']['cached']['data'])
    else:
        return pathlib.Path(params['inputs']['uncached']['data'])


def get_logging_path(logname) -> pathlib.Path:
    params = loadparams()
    logdir = pathlib.Path(params['outputs']['logdir'])
    return logdir / logname



if __name__ == '__main__':
    params = loadparams()
    outdir = get_outdir()
    can_cache = get_candidates_filepath(True)
    can = get_candidates_filepath(False)
    data_cache = get_data_filepath(True)
    data = get_data_filepath(False)
    logpath = get_logging_path('mylog.log')
    