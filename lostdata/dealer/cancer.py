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
    cancers = [c.decode().rstrip() for c in cancers]
    cancers = [(int(c[:3]),c[3:]) for c in cancers if c.strip()]
    for i,c in enumerate(cancers):
        if c[1].lstrip(' ').startswith('\t'):
            cancers[i] = (c[0],cancers[lvi][1]+c[1])
        else: lvi = i
    cancers = pd.DataFrame(cancers).set_index(0)[1]
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

def process_cancerregion(cancerregiondata,region):
    pediatrics = cancerregiondata.age_group <= 5
    yadults = (cancerregiondata.age_group > 5)&(cancerregiondata.age_group<=10)
    madults = (cancerregiondata.age_group > 10)&(cancerregiondata.age_group<=15)
    pedcases = cancerregiondata[pediatrics].number_of_cases.sum()
    yadcases = cancerregiondata[yadults].number_of_cases.sum()
    pediatricVSadolescentRatio = pedcases/yadcases
    pedagcases = cancerregiondata[pediatrics][['cancer_site','number_of_cases']].groupby('cancer_site').sum().loc[210].number_of_cases # 210     Adrenal gland (C74)
    #pedagVSotheryad = pedagcases/cancerregiondata[yadults][['cancer_site','number_of_cases']].groupby('cancer_site').sum().number_of_cases
    pedagVSotheryad = cancerregiondata[yadults][
        ['cancer_site','number_of_cases']].groupby('cancer_site').sum(
            ).number_of_cases/cancerregiondata[madults][['cancer_site','number_of_cases']].groupby('cancer_site').sum().number_of_cases
    pedagVSotheryad.name = region
    return (pedcases,yadcases,pediatricVSadolescentRatio,pedagVSotheryad)
