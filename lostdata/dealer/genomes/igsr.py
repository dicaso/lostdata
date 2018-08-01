# -*- coding: utf-8 -*-
"""1000 genomes data
International Genome Sample Resource (IGSR)

Info:
    - http://www.internationalgenome.org/
    - http://www.internationalgenome.org/about
"""
import lostdata as LSD
from lostdata.formats import SamplesDataset
from lostdata.processing import processedDataStorage, download_ftp_resource
import os, re, glob

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
