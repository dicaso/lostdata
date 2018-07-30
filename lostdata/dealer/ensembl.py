# -*- coding: utf-8 -*-
"""ensembl

Reference: https://www.ensembl.org/
"""
from lostdata import retrieveSources,cacheableTable,processedDataStorage,datadir
import os, gzip, pandas as pd
from io import TextIOWrapper, StringIO
from urllib.parse import parse_qsl

## Biomart
@cacheableTable
def get_biomart(atts=None,dataset='hsapiens_gene_ensembl'):
    """
    Get biomart id mappings

    Info:
        For possible attributes, see attribute table in
        https://www.bioconductor.org/packages/2.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf

    >>> bmrt = get_biomart()
    """
    import biomart
    server = biomart.BiomartServer("http://www.ensembl.org/biomart")
    # default atts
    if not atts:
        atts = ['external_gene_name','external_gene_source','ensembl_gene_id',
                'ensembl_transcript_id','ensembl_peptide_id']
    seda = server.datasets[dataset]
    s = seda.search({'attributes': atts}, header=1)
    data = pd.read_table(StringIO(s.content.decode()))
    return data

## Human resources
def get_ensembl(onlyGeneLabeled=True,onlyInChromosomes=None):
    import warnings
    warnings.warn("deprecated, use get_biomart instead", DeprecationWarning)
    ensembl = pd.read_table(datadir+'Genomes/Ensembl/Biomart/idmapping_extended.txt',
                        header=None,index_col=0,names=('egid','etid','start','stop','gcC','chr','strand','TSS','typeg','typet','gene_label','entrez'))
    ensembl = ensembl[~ensembl.index.duplicated()]
    if onlyGeneLabeled: ensembl = ensembl[ensembl.gene_label.isnull().apply(lambda x: not(x))]
    ensembl.chr = ensembl.chr.apply(lambda x: 'chr'+x.replace('MT','M'))
    if onlyInChromosomes: ensembl = ensembl[ensembl.chr.isin(onlyInChromosomes)]
    return ensembl

@retrieveSources
def get_ensemblGeneannot():
    """
    Info: http://www.ensembl.org/info/data/ftp/index.html
    Source: ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz
    """
    import gffutils
    try: db = gffutils.FeatureDB(processedDataStorage+'Homo_sapiens.GRCh38.88.sqlite3')
    except ValueError:
        if not os.path.exists(processedDataStorage+'Homo_sapiens.GRCh38.88.gtf.gz'):
            raise FileNotFoundError
        db = gffutils.create_db(processedDataStorage+'Homo_sapiens.GRCh38.88.gtf.gz',
                                processedDataStorage+'Homo_sapiens.GRCh38.88.sqlite3',
                                disable_infer_genes=True,disable_infer_transcripts=True)
    return db

## Cross species resources
### Mouse
@retrieveSources
def get_mouseEnsemblSet():
    """
    Source: ftp://ftp.ensembl.org/pub/release-90/gff3/mus_musculus/Mus_musculus.GRCm38.90.gff3.gz
    """
    db = pd.read_table(TextIOWrapper(gzip.open(processedDataStorage+'Mus_musculus.GRCm38.90.gff3.gz')),comment='#',
                       names='seqid,source,type,start,end,score,strand,phase,attribute'.split(','),low_memory=False)
    genes = db[db.type == 'gene'].copy()
    transcripts = db[db.type == 'mRNA'].copy()
    del db
    transcripts['TRANSCRIPT_ID'] = transcripts.attribute.apply(lambda x: parse_qsl(x)[0][1].split(':')[1])
    transcripts['GENE_ID'] = transcripts.attribute.apply(lambda x: parse_qsl(x)[1][1].split(':')[1])
    transcripts.set_index('TRANSCRIPT_ID',inplace=True)
    return (genes,transcripts)
