# -*- coding: utf-8 -*-
"""Wellcome Sanger Institute datasets

Reference: http://www.sanger.ac.uk/
"""
from bidali import LSD
import os, gzip, pandas as pd

@LSD.retrieveSources
def get_census():
    """Cancer census genes

    Info: https://cancer.sanger.ac.uk/census
    Locked source: https://cancer.sanger.ac.uk/cosmic/s3download?data=GRCh38%2Fcosmic%2Fv84%2Fcancer_gene_census.csv
    """
    return pd.read_csv(os.path.join(LSD.processedDataStorage,'cancer_gene_census.csv'),index_col=0)
