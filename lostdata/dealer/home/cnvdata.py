#!/usr/bin/env python
from bidali import LSD
import pandas as pd, numpy as np

@LSD.storeDatasetLocally
def get_UHRprofiles():
    """
    PDP meta cnv dataset from Ultra High Risk project. 
    Samples are all high risk (not ultra in se, that was project aim)

    Source: ~/Dropbiz/Basecamp/BRIP1/2017/2017_003_Data_analysis_BRIP1_paper/profiles_UHR/
    """
    datacn = LSD.expanduser("~/Dropbiz/Basecamp/BRIP1/2017/2017_003_Data_analysis_BRIP1_paper/profiles_UHR/")
    centromereshg38 = LSD.get_centromeres()
    lo = LSD.get_lift19to38()
    
    # All samples => TODO filter patients that only have whole chromosome gains
    samples = pd.read_table(datacn+'samples_UHR.txt')
    profiles = pd.read_table(datacn+'profiles_UHR.txt')
    profiles.chromosome = profiles.chromosome.apply(str).apply(lambda x: 'X' if x=='23' else x)
    profiles['size'] = profiles['max']-profiles['min']
    profiles['min38'] = profiles.apply(lambda x: lo.convert_coordinate('chr{}'.format(x['chromosome']),x['min']),axis=1).apply(lambda x: x[0][1] if x else np.nan)
    profiles['max38'] = profiles.apply(lambda x: lo.convert_coordinate('chr{}'.format(x['chromosome']),x['max']),axis=1).apply(lambda x: x[0][1] if x else np.nan)
    print('alterations lost per',profiles[profiles[['min38','max38']].isnull().any(axis=1)].groupby('chromosome').size())
    profiles = profiles[~profiles[['min38','max38']].isnull().any(axis=1)]
    del profiles['min']
    del profiles['max']

    #Study only non-normal gains
    profiles = profiles[profiles.annotation != 'normal']

    #Annotate chr arm
    profiles['chrarm'] = profiles.apply(lambda x: 'chr{}{}'.format(x['chromosome'],'p' if x['max38'] < centromereshg38.ix['chr'+x['chromosome']]['left_base'] else
                                                                   ('q' if x['min38'] > centromereshg38.ix['chr'+x['chromosome']]['right_base'] else 'p+q')), axis=1)
    #Annotate in profiles mycn status
    samples.index = samples.Name
    profiles['MYCN'] = profiles.profile_id.apply(lambda x: samples.ix[x].MYCN)
    #profsperMYCNstatus = profiles[profiles.chromosome != '2'].groupby(('MYCN','profile_id')).size()
    #profsperMYCNstatus[0].mean()
    #(profsperMYCNstatus[0].var())**.5

    return LSD.Dataset(profiles=profiles,metadata=samples)
    
