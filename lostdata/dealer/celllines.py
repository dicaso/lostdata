#!/usr/bin/env python
from bidali import LSD
import pandas as pd, numpy as np, gzip
from os.path import expanduser, exists
from itertools import count
from bidali.LSD import storeDatasetLocally, datadir, Dataset

def get_NB39():
    """
    39 neuroblastoma cell lines + RPE1 and HU.FETAL.BRAIN

    Reference: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89413
    Source: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE89413&format=file&\
file=GSE89413%5F2016%2D10%2D30%2DNBL%2Dcell%2Dline%2DSTAR%2Dfpkm%2Etxt%2Egz
    """
    exprdata = pd.read_table(
        gzip.open(
            datadir+'GEO/NB39_celllines_Maris/GSE89413_2016-10-30-NBL-cell-line-STAR-fpkm.txt.gz',
            'rt',encoding='UTF-8'
        ),
        index_col='GeneID'
    )

    return Dataset(exprdata=exprdata)

def get_CCLE():
    """
    Reference: https://portals.broadinstitute.org/ccle
    """
    metadata = pd.read_table(datadir+'CCLE/CCLE_sample_info_file_2012-10-18.txt',index_col='CCLE name')
    exprdata = pd.read_table(datadir+'CCLE/CCLE_Expression_Entrez_2012-09-29.gct',skiprows=2,index_col='Description')
    return Dataset(exprdata=exprdata,metadata=metadata)
