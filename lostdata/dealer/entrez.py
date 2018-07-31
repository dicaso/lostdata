# -*- coding: utf-8 -*-
"""entrez

References:
- https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html 
- https://www.ncbi.nlm.nih.gov/home/develop/api/
- ftp://ftp.ncbi.nih.gov/gene/
"""
from lostdata import retrieveSources,cacheableTable,processedDataStorage
import os, gzip, pandas as pd
from io import TextIOWrapper, StringIO
from urllib.parse import parse_qsl

@retrieveSources
def get_liftover(frm=19,to=38):
    """
    Info: http://hgdownload.cse.ucsc.edu/downloads.html
    """
    from pyliftover import LiftOver
    liftoverfile = 'hg{}ToHg{}.over.chain.gz'.format(frm,to)
    try: return LiftOver(processedDataStorage+liftoverfile)
    except FileNotFoundError:
        raise FileNotFoundError('Source: http://hgdownload.cse.ucsc.edu/gbdb/hg{}/liftOver/{}'.format(frm,liftoverfile))

def get_lift19to38():
    """
    Source: http://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
    """
    return get_liftover(frm=19,to=38)

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

#@retrieveSources -> not working login required
def get_msigdb6():
    """
    Source: http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.2/msigdb_v6.0.xml
    """
    import xml.etree.ElementTree as ET
    import pickle
    parser = ET.parse(processedDataStorage+'msigdb_v6.0.xml')
    root = parser.getroot()
    genesetsCollections = {} 
    for geneset in root:
        if (geneset.attrib['CATEGORY_CODE'] == 'ARCHIVED' or
            geneset.attrib['ORGANISM'] != 'Homo sapiens'): continue
        try:
            genesetsCollections[geneset.attrib['CATEGORY_CODE']
            ][geneset.attrib['STANDARD_NAME']] = geneset.attrib['MEMBERS_SYMBOLIZED'].split(',')
        except KeyError as e:
            genesetsCollections[geneset.attrib['CATEGORY_CODE']] = {}
            genesetsCollections[geneset.attrib['CATEGORY_CODE']
            ][geneset.attrib['STANDARD_NAME']] = geneset.attrib['MEMBERS_SYMBOLIZED'].split(',')
    genesetsCollections = {'version':'6','MSigDB':genesetsCollections}        
    print('MSigDB {}'.format(genesetsCollections['version']))
    mdb = genesetsCollections['MSigDB']
    return mdb
