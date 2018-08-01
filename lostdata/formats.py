# -*- coding: utf-8 -*-
"""Format types for localised data structures.
"""
from .processing import retrieveSources, cacheableTable, processedDataStorage
import os, gzip, zipfile, pandas as pd
from collections import OrderedDict

## General
class Dataset(object):
    """A Dataset object is a collection of data tables, accessible as
    attributes from the Dataset object.

    to_R pushes the sub datasets to R making them available in the
    global namespace
    """
    def __init__(self,**kwargs):
        self.__datasets__ = set(kwargs)
        for kw in kwargs:
            self.__setattr__(kw,kwargs[kw])

## Integrated
class IntegratedDataset(object):
    """IntegratedDataset checks and assures that for a set of samples, 
    features and datatypes all data is internally cross-referencable and
    then stores it as an xarray.
    """
    def __init__(self, xarray_kwargs = {}, **kwargs):
        self.data = xr.Dataset(data_vars = kwargs, **xarray_kwargs)

    def save(self,filename,csvsep=',',csvdec='.'):
        """Save IntegratedDataset to disk

        Iterates over the individual data variables and saves them
        in a zipfolder

        Args:
            filename (str): Save location. Extension '.zip' will be added.
        """
        import zipfile, io
        if filename.endswith('.zip'): filename = filename[:-4]
        with zipfile.ZipFile('{}.zip'.format(filename), mode='w') as zdir:
            for datatablename in self.data.data_vars:
                b = io.StringIO()
                self.data.data_vars[datatablename].to_pandas().to_csv(b,sep=csvsep,decimal=csvdec)
                b.seek(0)
                zdir.writestr('{}.csv'.format(datatablename),b.read())
                
    @classmethod
    def load(cls,filename,csvsep=',',csvdec='.'):
        """Load IntegratedDataset from disk zip folder.

        Args:
            filename (str): Save location. Extension '.zip' will be added.
        """
        import zipfile, io
        if filename.endswith('.zip'): filename = filename[:-4]
        with zipfile.ZipFile('{}.zip'.format(filename), mode='r') as zdir:
            tables = {
                zf.filename[:-4]:pd.read_csv(
                    io.TextIOWrapper(zdir.open(zf)),
                    sep=csvsep, decimal=csvdec, index_col=0
                )
                for zf in zdir.infolist()
            }
        return cls(**tables)

    # def save_netcdf(self, filename):
    #     """Save IntegratedDataset to disk

    #     Uses xr.Dataset netCDF method

    #     Args:
    #         filename (str): Save location. Suggested extension '.nc'.
    #     """
    #     raise NotImplementedError('Issues with netcdf output')
    #     self.data.to_netcdf(filename)

    # @staticmethod
    # def load_from_netcdf(filename):
    #     """Load IntegratedDataset from disk

    #     Uses xr.Dataset netCDF method

    #     Args:
    #         filename (str): Filename of dataset to load.
    #     """
    #     raise NotImplementedError('serializing and loading not yet working')
    #     return xr.open_dataset(filename)
    

### Specific integrated datasets
class CoexDataset(IntegratedDataset):
    """IntegratedDataset that contains data on
    copy number alteration status of genes (optionally
    including the regions), and gene expression data.

    *conureal* is optional, but if provided and *conugeal* is None,
    *conugeal* is calculated based on *conureal* and the genes in
    *expression*.

    Args:
        conugeal (pd.DataFrame): Copy number gene alterations, samples should be columns and genes rows.
        expression (pd.DataFrame): Expression data, samples should be columns and genes rows.
        metadata (pd.DataFrame): Annotation data should be columns and samples rows.
            Recommended columns are `DOD` (bool, "dead of disease"), `POD` (bool, "progress of disease"), `overallSurvival` (int), `progressionFreeSurvival` (int)
        metadataColMapping (dict): dict if metadata column names need to be changed to follow recommendations, e.g. {'survived':'DOD'}.
        conureal (pd.DataFrame): Copy number region alterations. 
            Required columns are `chrom`, `chromStart`, `chromEnd` and `samples`.
        conurealColMapping (dict): dict if conureal column names need to be changed to follow requirements , e.g. {'chr':'chrom'}
        filterUnknownExpressed (bool): If True, genes for which there is no expression data, will be filtered from conugeal.
    """
    def __init__(self,conugeal,expression,metadata,metadataColMapping={},conureal=None,conurealColMapping={},filterUnknownExpressed=True):
        # Applying col mapping
        if metadataColMapping: metadata = metadata.rename(metadataColMapping,axis=1)
        if conurealColMapping: conureal = conureal.rename(conurealColMapping,axis=1)

        # Checking data
        import warnings
        from pandas.api import types
        typefns = {bool:types.is_bool_dtype,int:types.is_integer_dtype}
        for n,t in (
                ('DOD',bool),
                ('POD',bool),
                ('overallSurvival',int),
                ('progressionFreeSurvival', int)
        ):
            if n in metadata and not typefns[t](metadata[n]):
                warnings.warn('%s in metadata is not of correct type %s' % (n,t))
            elif n not in metadata:
                warnings.warn('%s not in metadata' % n)
        
        # If only regions is provided, and conugeal is None, calculate gene copy number values for all expression genes
        if conureal and not conugeal:
            conugeal = self.__class__.calculate_conugeal(conureal,expression.index)
        self.conureal = conureal

        # Filter genes for which there is no expression data
        if filterUnknownExpressed:
            conugeal = conugeal.filter(items=expression.index,axis=0)
            
        #Setting dimension names
        conugeal.columns.name = expression.columns.name = metadata.index.name = 'sample'
        conugeal.index.name = expression.index.name = 'gene'
        metadata.columns.name = 'annotation'
        if conureal:
            conureal.columns.name = 'sample'
            conureal.index.name = 'region'
            
        super().__init__(
            conugeal=conugeal,
            expression=expression,
            metadata=metadata
        )

    @staticmethod
    def calculate_conugeal(conureal,genes):
        # Assign genes to regions
        genannot = get_ensemblGeneannot()
        segments['genes'] = segments.T.apply(
            lambda x: {f.attributes['gene_name'][0] for f in genannot.region('{}:{}-{}'.format(
                x.chromosome[3:],int(x.Start38),int(x.End38)),featuretype='gene')}
        )
        segments['nrGenes'] = segments.genes.apply(len)
        del genannot#, get_genes
        conugeal = {}
        for sname,sample in segments.groupby('sample'):
            conugeal[sname] = {}
            for i in sample.T:
                for g in sample.ix[i].genes: conugeal[sname][g] = sample.ix[i].value
        conugeal = pd.DataFrame(conugeal)
        return conugeal

## Dataset repository
class DatasetRepo:
    """
    Object containing:
      - a Dataset object
      - a report on how the dataset was generated
      - file location
      - a history of the code that generated the current
        dataset and earlier versions
      - a hash of the code that generated the current dataset
      - archived earlier dataset objects

    Method to keep only the most recent version of the
    dataset to free up hard disk space
    """
    def __init__(self,dataset,code,report,filename):
        import hashlib, pickle
        self.dataset = dataset
        self.report = report
        self.currentHash = hashlib.md5(code.encode()).hexdigest()
        self.code = OrderedDict([(self.currentHash,code)])
        self.archive = OrderedDict()
        self.filename = filename
        pickle.dump(self,open(filename,'wb'))

    def update(self,dataset,code,report):
        import hashlib, pickle
        # Put previous object in archive
        self.archive[self.currentHash] = self.dataset

        # Update
        self.dataset = dataset
        self.report = report
        self.currentHash = hashlib.md5(code.encode()).hexdigest()
        self.code[self.currentHash] = code
        pickle.dump(self,open(self.filename,'wb'))

    def wipeArchive(self):
        self.archive = OrderedDict()
        pickle.dump(self,open(self.filename,'wb'))

## Dataset class for working with samples collections
class SamplesDataset(object):
    """Samples dataset

    Class to create a collection of samples. Intantiating an object does
    not start downloading data, but only configures the object. This
    design choice has been made because the lostdata.dealer module
    containes several SamplesDataset resources.

    Args:
        name (str): Dataset name.
        source (str): Where the original dataset is located.
        protocol (str): Protocol for accessing/downloading the data.
        target (str): Directory where data will be stored. If not provided,
          will be stored in 'name' subdir of `processedDataStorage`.
          
    """
    class Samples(object):
        def __init__(self):
            self.__list = []

        def append(self,obj):
            self.__list.append(obj)

        def __getitem__(self,key):
            return self.__list[key]
    
    def __init__(self,name,source,protocol='requests',target=None):
        self.name = name
        self.source = source
        self.protocol = protocol
        self.samples = Samples()
        self.target = target if target else os.path.join(processedDataStorage,name)
        self.targetsamples = os.path.join(self.target,'samples')
        
    def setuptarget(self):
        """Setup target directories"""
        if not os.path.exists(self.target):
            os.mkdir(self.target)
            os.mkdir(self.targetsamples)
    

## XML dataset
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
