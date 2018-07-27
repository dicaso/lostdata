# -*- coding: utf-8 -*-
"""ensembl

Reference: https://www.ensembl.org/
"""
from bidali import LSD
from bidali.LSD import retrieveSources, cacheableTable, processedDataStorage
import os, gzip, pandas as pd
from io import StringIO
import requests

# KEGG API
KEGG_BASE_URL = 'http://rest.kegg.jp/'

def list_KEGG_pathways(species='hsa'):
    """List KEGG pathways for specified species

    Args:
        species (str): KEGG species code. Default `hsa` == human.
    """
    r = requests.get(f'{KEGG_BASE_URL}list/pathway/{species}')
    data = StringIO(r.text)
    table = pd.read_table(data,names=['pathway_id','description'])
    table['pathway_name'] = table.description.apply(lambda x: x[:x.rindex('-')].strip())
    return table

## pathway regex
def get_KEGG_pathway(pathway_id):
    """Get a KEGG pathway

    Args:
        pathway_id (str): Pathway ID, e.g. `'path:hsa00010'`
    """
    p = requests.get(f'{KEGG_BASE_URL}get/{pathway_id}')
    try: genes = p.text[p.text.index('\nGENE')+5:p.text.index('\nCOMPOUND')]
    except ValueError: genes = p.text[p.text.index('\nGENE')+5:p.text.index('\nREFERENCE')]
    genes = [g.strip().split()[1][:-1] for g in genes.split('\n')]
    return genes
