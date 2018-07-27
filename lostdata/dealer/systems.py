#!/usr/bin/env python
from bidali import LSD
import gzip, tarfile, biomart, pandas as pd
from io import TextIOWrapper, StringIO
from .ensembl import get_biomart

def get_proteinHalfLifes(datadir=LSD.datadir+'PMID_published_tables/21593866_protein_half_lifes/'):
    """
    Mouse protein half lifes

    >>> phl = get_proteinHalfLifes()
    """
    phl = pd.read_excel(datadir+"nature10098-s5.xls",sheetname=0)
    return phl

def get_cyclebase(datadir=LSD.datadir+'PMID_published_tables/25378319_cyclebase/'):
    """
    Info: https://cyclebase.org/
    Info: http://genome-www.stanford.edu/Human-CellCycle/HeLa/data/dataPlusScores_all5.txt

    cyb = get_cyclebase()
    """
    tf = tarfile.open(datadir+'human_timecourses.tar')
    #metadata = pd.read_table(TextIOWrapper(tf.extractfile('human_metadata.tsv')))
    hexp = pd.read_table(TextIOWrapper(tf.extractfile('human_experiments.tsv')))
    #whitfieldData = pd.read_table(datadir+'dataPlusScores_all5.txt',index_col='UID')
    #whitfieldData['gene_name'] = whitfieldData.NAME.fillna('').apply(lambda x: x.split()[1] if len(x.split())>1 else '')
    #hexp['gene_name'] = hexp.identifier.apply(lambda x: whitfieldData.ix[x].gene_name)
    mapping = pd.read_table(datadir+'cyb_probeset_mapping.txt',header=None,names=('species','protein','code_id'))

    # Periodic genes
    tf = tarfile.open(datadir+'human_periodic.tar')
    periodicGenes = pd.read_table(TextIOWrapper(tf.extractfile('human_periodic.tsv')))
    bmrt = get_biomart()
    protgenenames = bmrt[['Gene name','Protein stable ID']].dropna().set_index('Protein stable ID')
    periodicGenes = periodicGenes.join(other=protgenenames,on='gene').set_index('Gene name')

    return LSD.Dataset(
        periodicGenes=periodicGenes,
        hexp=hexp,
        mapping=mapping
    )
