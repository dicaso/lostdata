#!/usr/bin/env python
from bidali import LSD
import gzip, pandas as pd
from io import TextIOWrapper

def get_cancerrxgene(datadir=LSD.datadir+'PMID_published_tables/27397505_pharmacogenomics_cancer/', gene_cn=False, exprdata=False):
    """
    Reference: http://www.cancerrxgene.org/downloads

    cn data:
    Four comma seperated peices of data are presented for each gene cell line combination (n1,n2,n3,n4), hyphen (-) is used where the value is unknown:-

    n1 - Maximum copy number of any genomic segment containing coding sequence of the gene (-1 indicates a value could not be assigned).
    n2 - Minimum copy number of any genomic segment containing coding sequence of the gene (-1 indicates a value could not be assigned).
    n3 - Zygosity - (H) if all segments containing gene sequence are heterozygous, (L) if any segment containing coding sequence has LOH, 
         (0) if the complete coding sequence of the gene falls within a homozygous deletion.
    n4 - Disruption (D) if the gene spans more than 1 genomic segment (-) if no disruption occures.

    >>> cdx = get_cancerrxgene()
    """
    # Metadata
    cell_lines = pd.read_excel(datadir+'Cell_Lines_Details.xlsx',sheetname='Cell line details')
    cell_tissues = pd.read_excel(datadir+'Cell_Lines_Details.xlsx',
                                 sheetname='COSMIC tissue classification',index_col='Line')
    cell_lines = cell_lines.join(cell_tissues,on='Sample Name')
    del cell_lines['COSMIC_ID']
    del cell_tissues

    # Compounds
    compounds = pd.read_excel(datadir+'Screened_Compounds.xlsx',index_col='Drug Name')

    # Compound responses
    compound_responses = pd.read_excel(datadir+'v17_fitted_dose_response.xlsx',index_col='COSMIC_ID')

    # Expression data
    if exprdata:
        exprdata = pd.read_table(TextIOWrapper(gzip.open(datadir+'sanger1018_brainarray_ensemblgene_rma.txt.gz')),index_col='ensembl_gene')
    
    # CN data
    if gene_cn:
        gene_cn = pd.read_excel(datadir+'Gene_level_CN.xlsx',sheetname='Gene_level_CN',index_col='gene')

    return LSD.Dataset(**locals())
