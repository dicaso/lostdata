# -*- coding: utf-8 -*-
"""uniprot

Reference: https://www.uniprot.org/
"""
from lostdata import retrieveSources, cacheableTable, processedDataStorage
import os, gzip, pandas as pd
from io import TextIOWrapper, StringIO
from urllib.parse import parse_qsl

@retrieveSources
def get_uniprot():
    """
    Info: https://www.uniprot.org/downloads
    Source: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz
    """
    import xml.etree.ElementTree as ET
    with gzip.open(processedDataStorage+'uniprot_sprot.xml.gz') as gf:
        tree = ET.parse(gf)
        root = tree.getroot()
    return root
