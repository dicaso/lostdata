# -*- coding: utf-8 -*-
"""entrez

References:
- https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html 
- https://www.ncbi.nlm.nih.gov/home/develop/api/
- ftp://ftp.ncbi.nih.gov/gene/
"""
from bidali import LSD
from bidali.LSD import retrieveSources,cacheableTable,processedDataStorage,datadir
import os, gzip, pandas as pd
from io import TextIOWrapper, StringIO
from urllib.parse import parse_qsl

@retrieveSources
def get_refseq():
    """
    Source: ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene
    """
    entrez = pd.read_table(processedDataStorage+'gene_RefSeqGene', index_col='GeneID')
    return entrez

@retrieveSources
def get_gene_info():
    """
    Source: ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

    All species -> ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
    """
    gene_info = pd.read_table(
        TextIOWrapper(gzip.open(processedDataStorage+'Homo_sapiens.gene_info.gz')),
        index_col='GeneID'
    )
    return gene_info
