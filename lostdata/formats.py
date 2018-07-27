# -*- coding: utf-8 -*-
"""Format types for localised data structures.
"""
from bidali import LSD
from bidali.LSD import retrieveSources,cacheableTable,processedDataStorage,datadir
import os, gzip, zipfile, pandas as pd

class XML(object):
    """Class to work with xml files.
    
    Mainly methods to access specific data easier, than directry to Python standard lib.

    Args:
        xmlfile (str): Can either be a filename or an already parsed tree (but not the root element).
    """
    def __init__(self,xmlfile,namespace_symbol=None,namespace_url=None,schema=None):
        import xml.etree.ElementTree as ET
        if isintance(xmlfile,str):
            with open(xmlfile) as f:
                self.tree = ET.parse(f)
        else:
            self.tree = xmlfile
        self.root = self.tree.getroot()
        self.namespace_symbol = namespace_symbol
        if namespace_symbol:
            self.tree_namespace = {self.namespace_symbol:namespace_url}
        self.schema = schema

    def __getitem__(self,xpath):
        #prepending namespace to all path elements
        if self.namespace_symbol: 
            xpath = '/'.join(
                ':'.join((self.namespace_symbol,d)) for d in xpath.split('/')
            )
        return self.root.findall(xpath, self.tree_namespace)

    def get_tags(self):
        return {
            e.tag.replace(
                '{'+self.tree_namespace[self.namespace_symbol]+'}',
                ''
            ) if self.namespace_symbol else e
            for e in self.tree.iter()
        }

    def validate(self):
        import xmlschema
        xs = xmlschema.XMLSchema(self.schema, converter=xmlschema.BadgerFishConverter)
        xs.validate(self.tree)
