# -*- coding: utf-8 -*-
"""Biological network datasets
"""
from lostdata import retrieveSources,cacheableTable,processedDataStorage,datadir
import os, gzip, pandas as pd
from io import TextIOWrapper, StringIO

@storeDatasetLocally
def get_proteinNetworks():
    """
    Source: https://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.147/BIOGRID-ALL-3.4.147.tab2.zip
    Source: http://string-db.org/download/protein.links.v10/9606.protein.links.v10.txt.gz
    Source: http://string-db.org/mapping_files/entrez_mappings/entrez_gene_id.vs.string.v10.28042015.tsv
    """
    import networkx as nx
    from lostdata.dealer import entrez
    
    #Biogrid
    with ZipFile(datadir+'ProteinNetworks/BIOGRID-ALL-3.4.147.tab2.zip') as biogridzip:
        ds = pd.read_table(TextIOWrapper(biogridzip.open('BIOGRID-ALL-3.4.147.tab2.txt','r')),low_memory=False)
    ds = ds[ds['Organism Interactor A'] == 9606]
    Gbio = nx.Graph()
    ds.T.apply(lambda x: Gbio.add_edge(x['Official Symbol Interactor A'],x['Official Symbol Interactor B']))

    #String-DB
    stringdb = pd.read_table(gzip.open(datadir+'ProteinNetworks/9606.protein.links.v10.txt.gz','rt'),sep=' ')
    stringids = pd.read_table(datadir+'ProteinNetworks/entrez_gene_id.vs.string.v10.28042015.tsv',index_col='STRING_Locus_ID')
    entrez = entrez.get_refseq()
    stringids = stringids[stringids['#Entrez_Gene_ID'].isin(entrez.index)]
    stringdb = stringdb[stringdb.combined_score > 400] #study stringdb.combined_score.hist(bins='auto') to set threshold
    stringdb = stringdb[stringdb.protein1.isin(stringids.index) & stringdb.protein2.isin(stringids.index)]
    stringdb.protein1 = stringdb.protein1.apply(lambda x: stringids.loc[x]['#Entrez_Gene_ID'])
    stringdb.protein2 = stringdb.protein2.apply(lambda x: stringids.loc[x]['#Entrez_Gene_ID'])
    stringdb.protein1 = stringdb.protein1.apply(lambda x: entrez.loc[x].Symbol)
    stringdb.protein2 = stringdb.protein2.apply(lambda x: entrez.loc[x].Symbol)
    Gstring = nx.Graph()
    stringdb.T.apply(lambda x: Gstring.add_edge(x.protein1,x.protein2))
    
    return Dataset(biogridnx = Gbio, biogrid = ds,
                   stringnx = Gstring, string = stringdb)
