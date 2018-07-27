# -*- coding: utf-8 -*-
"""Data from genenames.org

Reference: https://www.genenames.org/
"""
from bidali import LSD
from bidali.LSD import retrieveSources,cacheableTable,processedDataStorage,datadir
import os, gzip, pandas as pd
from io import TextIOWrapper, StringIO

@retrieveSources
def get_genenames():
    """
    Source: genenames.tsv ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt
    """
    return pd.read_table(
        os.path.join(processedDataStorage, 'genenames.tsv'), low_memory = False
    )#, index_col='GeneID')

@retrieveSources
def get_genefamilies():
    """
    Source: genefamilies.tsv https://www.genenames.org/cgi-bin/genefamilies/download-all/tsv
    """
    return pd.read_table(
        os.path.join(processedDataStorage, 'genefamilies.tsv'), low_memory = False
    )
