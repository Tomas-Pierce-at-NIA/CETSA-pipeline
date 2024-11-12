# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:54:00 2024

@author: piercetf
"""

import requests

import asyncio
import aiohttp
import aiofiles

import json
import io

import pandas

BASEURL = "https://rest.uniprot.org"

RUN = '/idmapping/run'

STATUS = '/idmapping/status'

RESULTS = '/idmapping/results'


SEARCH = '/uniprotkb/search'      

FISETIN = r'T:\TGB\LSS\TGU\Users\Tomas\2024_CETSA_MS\Monocyte_CETSA_Statistical_Analysis\outputs\nparc_based_Aug2024\sigcomps_DMSO_Fisetin.csv'
QUERCETIN = r'T:\TGB\LSS\TGU\Users\Tomas\2024_CETSA_MS\Monocyte_CETSA_Statistical_Analysis\outputs\nparc_based_Aug2024\sigcomps_DMSO_Quercetin.csv'

async def start_job(session :aiohttp.ClientSession, job_params :dict) -> dict:
    async with session.post(RUN, params=job_params) as resp:
        bodystr = await resp.text()
        data = json.loads(bodystr)
    return data

async def checkout_results(session :aiohttp.ClientSession, jobid :str):
    async with session.get(STATUS + '/' + jobid, allow_redirects=True) as resp:
        bodystr = await resp.text()
        data = json.loads(bodystr)
    return data


async def domain_search(session : aiohttp.ClientSession, ident: str):
    params = {'query': ident,
              'fields' : ['accession', 
                          'id', 
                          'gene_names',
                          'organism_name',
                          'xref_interpro_full'],
              'format' : 'tsv'
              }
    async with session.get(SEARCH, params=params) as resp:
        bodystr = await resp.text()
    
    return bodystr

async def domain_searchs(session :aiohttp.ClientSession, idlist: list):
    #domaintxt = {}
    checkall = " OR ".join(idlist)
    checkall_human = f"({checkall}) AND organism_name:human"
    alltext = await domain_search(session, checkall_human)
            
    return alltext




async def main():
    
    genes = []
    accessions = []
    
    async with aiofiles.open(FISETIN) as fisetinhandle:
        await fisetinhandle.readline()
        async for line in fisetinhandle:
            components = line.split(',')
            genes.append(components[2])
            accessions.append(components[1])
    
    async with aiohttp.ClientSession(BASEURL) as session:
        dbody = await domain_searchs(session, accessions)
    
    text = io.StringIO(dbody)
    table = pandas.read_csv(text, sep='\t')
    unquoted = table['InterPro'].str.replace('"', '')
    table.loc[:, 'Domains'] = unquoted.str.split(';')
    
    return table




if __name__ == '__main__':
    asyncio.run(main())