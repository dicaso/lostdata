# -*- coding: utf-8 -*-
"""cancer related data
"""
import lostdata as LSD
import pandas as pd, numpy as np, gzip, os
from itertools import count
from lostdata import storeDatasetLocally, Dataset

@LSD.retrieveSources
def get_cancer_incidence_in_5_continents():
    """Cancer incidence data
    
    Info: http://gco.iarc.fr/
    Source: http://ci5.iarc.fr/CI5-X/CI5-Xd.zip
    """
    import zipfile
    from io import TextIOWrapper
    zip = zipfile.ZipFile(
        os.path.join(
            LSD.config['LSD']['cachedir'],
            'CI5-Xd.zip'
        )
    )
    registry = pd.read_table(
        TextIOWrapper(
            zip.open('registry.txt'),
            encoding='Windows-1252'
        ),
        header=None, names=['code','region']
    )
    registry = registry.set_index('code').region
    cancers = zip.open('cancer.TXT').readlines()
    cfiles = zip.infolist()
    regional_data = {}
    for cf in cfiles:
        if cf.filename.endswith('.csv'):
            regional_data[
                registry[int(cf.filename[:-4])].strip()
            ] = pd.read_csv(
                TextIOWrapper(zip.open(cf)),
                header = None,
                names = ['sex', 'cancer_site', 'age_group',
                             'number_of_cases', 'person_years_at_risk'
                             ]
            )
    zip.close()
    return {
        'cancers': cancers,
        'registry': registry,
        'regional_data': regional_data
    }
