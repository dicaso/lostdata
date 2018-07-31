# -*- coding: utf-8 -*-
"""lostdata: LOcalised STructured DAta TAbles
Module that preprocesses commonly used datasets
and makes them available for use in python.
"""
import gzip, pickle, os
from zipfile import ZipFile
from io import TextIOWrapper, StringIO
import pandas as pd, numpy as np, xarray as xr
from os.path import expanduser, exists

## Defaults
from .config import config

# Check file locations
if not (
        os.path.exists(config['LSD']['privatedir']) and
        os.path.exists(config['LSD']['cachedir'])
        ):
    import warnings
    warnings.warn(
        'Create directory %s or %s for lostdata to function properly, or configure cachedir/privatedir'
        % (config['LSD']['privatedir'], config['LSD']['cachedir'])
    )

## Utility functions
def getLSDataset(name,**kwargs):
    try: return name(**kwargs)
    except TypeError:
        return globals()[name](**kwargs)

def listLSDatasets(subpackages=True):
    print(*[i for i in sorted(globals()) if i.startswith('get_')],sep='\n')
    if subpackages:
        print('\nWithin subpackages:')
        import pkgutil as pk
        from inspect import getmembers
        for mod in pk.walk_packages(__path__,__name__+'.'):
            if not mod.ispkg:
                mol = pk.importlib.import_module(mod[1])
                print('>',mod[1])
                print(*[i for i in sorted(dict(getmembers(mol)))
                        if i.startswith('get_')],sep='\n',end='\n\n')

# Structured dataset classes
from .formats import Dataset, IntegratedDataset, DatasetRepo

# Processing
from .processing import processedDataStorage, retrieveSources, cacheable, \
     cacheableTable, storeDatasetLocally

