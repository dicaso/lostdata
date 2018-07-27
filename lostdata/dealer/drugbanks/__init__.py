# -*- coding: utf-8 -*-
"""Drugging databases

Drug names, medication and procedure links.
"""
from bidali import LSD
from bidali.LSD import retrieveSources,cacheableTable,processedDataStorage,datadir
import os, gzip, zipfile, pandas as pd
from io import TextIOWrapper, StringIO
from urllib.parse import parse_qsl

class DrugBank:
    """Convenience class to work with
    https://www.drugbank.ca drugbank xml.
    
    Mainly methods to access specific data
    """
    def __init__(self,drugtree):
        self.tree = drugtree
        self.root = self.tree.getroot()
        self.namespace_symbol = 'ns0'
        self.tree_namespace = {self.namespace_symbol:'http://www.drugbank.ca'}

    def __getitem__(self,xpath):
        #prepending namespace to all path elements
        xpath = '/'.join(
            ':'.join((self.namespace_symbol,d)) for d in xpath.split('/')
        )
        return self.root.findall(xpath, self.tree_namespace)

    def get_tags(self):
        return {
            e.tag.replace(
                '{'+self.tree_namespace[self.namespace_symbol]+'}',
                ''
            )
            for e in self.tree.iter()
        }

    def get_drugnameset(self):
        return {
            d.text.strip() for d in self['drug/name']
        }
        
@retrieveSources
def get_drugbank():
    """
    Info: https://www.drugbank.ca/docs/drugbank.xsd

    Source: drugbank_all_full_database.xml.zip $locked$https://www.drugbank.ca/releases/5-1-1/downloads/all-full-database
    Source: https://www.drugbank.ca/docs/drugbank.xsd
    """
    import xml.etree.ElementTree as ET:
    with zipfile.ZipFile(os.path.join(processedDataStorage,'drugbank_all_full_database.xml.zip')) as z:
        with z.open('full database.xml') as f:
            tree = ET.parse(f)
    return DrugBank(tree)
