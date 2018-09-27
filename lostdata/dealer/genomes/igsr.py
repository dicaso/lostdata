# -*- coding: utf-8 -*-
"""1000 genomes data
International Genome Sample Resource (IGSR)

Info:
    - http://www.internationalgenome.org/
    - http://www.internationalgenome.org/about
"""
import lostdata as LSD
from lostdata.formats import SamplesDataset
from lostdata.processing import processedDataStorage, retrieveSources, download_ftp_resource
import os, re, glob, gzip, pandas as pd

def get_vcf_1000g(includeIndices=True):
    targetdir = os.path.join(processedDataStorage,'1000g_vcf')
    if not os.path.exists(targetdir): os.mkdir(targetdir)
    vcf_re = re.compile('\.vcf\.gz' if includeIndices else '\.vcf\.gz$')
    download_ftp_resource(
        server = 'ftp.1000genomes.ebi.ac.uk',
        directory = 'vol1/ftp/release/20130502',
        filter = vcf_re,
        targetdir = targetdir
    )
    return glob.glob(os.path.join(targetdir,'*.vcf.gz'))

@retrieveSources
def get_500g_RNAseq():
    """Get RNA seq expression data
    for the 462 subset (mostly European) of the 1000genome project

    Source: https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
    """
    return pd.read_table(
        gzip.open(
            os.path.join(
                processedDataStorage,
                'GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz'
            ),'rt',encoding='UTF-8'
        )
    )
