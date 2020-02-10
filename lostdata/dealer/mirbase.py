"""miRBase resources

- mature.fa
"""
from lostdata import processedDataStorage, retrieveSources

@retrieveSources
def get_mature_miRNA_sequences(filterSpecies=''):
    """Get the mature miRNA sequences

    Source: ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

    Args:
        filterSpecies (str): E.g. 'Homo sapiens'
    """
    import gzip
    import os
    with gzip.open(
        os.path.join(processedDataStorage, 'mature.fa.gz'),
        mode='rt'
    ) as micrornafile:
        micrornas = {}
        for line in micrornafile:
            if line.startswith('>'):
                hs_micro = line.strip() if filterSpecies in line else False
            else:
                if hs_micro: micrornas[hs_micro] = line.strip()
    return micrornas
