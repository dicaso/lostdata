# -*- coding: utf-8 -*-
"""lostdata: LOcalised STructured DAta TAbles
Module that preprocesses commonly used datasets
and makes them available for use in python.
"""
import gzip, pickle, time, os
from zipfile import ZipFile
from io import TextIOWrapper, StringIO
import pandas as pd, numpy as np, xarray as xr
from os.path import expanduser, exists
from collections import OrderedDict
from contextlib import redirect_stdout, redirect_stderr

## Defaults
from .config import config
processedDataStorage = config['LSD']['cachedir']
datadir = config['LSD']['privatedir']

# Check file locations
if not (
        os.path.exists(datadir) and
        os.path.exists(processedDataStorage)
        ):
    import warnings
    warnings.warn(
        'Create directory %s or %s for lostdata to function properly, or configure cachedir/privatedir'
        % (datadir, processedDataStorage)
    )

## Utility functions
def getLSDataset(name,**kwargs):
    try: return name(**kwargs)
    except TypeError:
        return globals()[name](**kwargs)

def listLSDatasets(subpackages=True):
    print(*[i for i in sorted(globals()) if i.startswith('get_')],sep='\n')
    if subpackages:
        print('\nWithin subpackages:')
        import pkgutil as pk
        from inspect import getmembers
        for mod in pk.walk_packages(__path__,__name__+'.'):
            if not mod.ispkg:
                mol = pk.importlib.import_module(mod[1])
                print('>',mod[1])
                print(*[i for i in sorted(dict(getmembers(mol)))
                        if i.startswith('get_')],sep='\n',end='\n\n')

# Structured dataset classes
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
    

## Specific integrated datasets
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
        
def retrieveSources(dataset_getfunction):
    """
    A dataset_getfunction function that contains 'Source:' lines
    in the docstring, can be decorated with this function.
    If a source is not locally available, it will be downloaded
    and added to the processedDataStorage location.

    A source line has to be formatted accordingly:
    Source: [filename] url

    If filename is not provided, the last part of the url (after last '/')
    is taken as filename.

    Source lines can also be provided as arguments to the FileNotFoundError
    that the dataset_getfunction throws.
    """
    import inspect, requests
    from urllib.request import urlopen, urlretrieve

    def wrapper(*args, **kwargs):
        try:
            return dataset_getfunction(*args, **kwargs)
        except FileNotFoundError as fnf:
            for docline in inspect.getdoc(dataset_getfunction).split('\n') + list(fnf.args):
                if docline.startswith('Source:'):
                    docline = docline.split()
                    if len(docline) == 2:
                        url = docline[1]
                        filename = url[url.rindex('/')+1:]
                    elif len(docline) == 3:
                        url = docline[2]
                        filename = docline[1]
                    if not exists(os.path.join(processedDataStorage,filename)):
                        print('Downloading {}:'.format(filename))
                        if url.startswith('ftp://'):
                            urlretrieve(url,os.path.join(processedDataStorage,filename),
                                        lambda x,y,z: print("\r[{}{}]".format('=' * int(50*x*y/z),
                                                                              ' ' * (50-int(50*x*y/z))),
                                                            end='',flush=True))
                        else:
                            r = requests.get(url,stream=True)
                            total_length = r.headers.get('content-length')
                            with open(processedDataStorage+filename,'wb') as f:
                                if total_length is None: f.write(r.content)
                                else:
                                    dl = 0
                                    total_length = int(total_length)
                                    for data in r.iter_content(chunk_size=4096):
                                        dl += len(data)
                                        f.write(data)
                                        done = int(50 * dl / total_length)
                                        print("\r[{}{}]".format('=' * done, ' ' * (50-done)),end='',flush=True)
            try: return dataset_getfunction(*args, **kwargs)
            except FileNotFoundError:
                print('Either not all source files are documented correctly in docstring,',
                      'or there is a source file unrelated issue')
                raise

    return wrapper

def cacheable(importer,exporter,extension='.cache'):
    """
    Produces decorators with a specified importer and exporter function.
    For example, if a function produces a pandas DataFrame, the importer,
    and exporter could be respectively a lambda for DataFrame.read_csv 
    and DataFrame.to_csv

    extension specifies the file extension used by the caching functions
    """
    def cachedecorator(function):
        def wrapper(*args,cache=True,cache_name='',**kwargs):
            #produce function call hash -> function name, args and kwargs should be merged, stringified and hashed
            import inspect, hashlib, time, datetime
            signature = inspect.signature(function)
            boundArgs = signature.bind(*args,**kwargs)
            boundArgs.apply_defaults()
            callhash = '{}_{}'.format(
                function.__name__ if not cache_name else cache_name,
                hashlib.md5(boundArgs.__repr__().encode()).hexdigest()
            )
            cachedir = config['LSD']['cachedir']
            cachefile = os.path.join(cachedir,'{}{}'.format(callhash,extension))
            # If cache, check if cache exists and how long it is allowed to exist in config
            timeMap = {'h': 'hours', 'd': 'days', 'w': 'weeks'}
            cacheAllowedTime = config['LSD']['cachetime']
            cacheAllowedTime = datetime.timedelta(
                **{timeMap[t]:int(cacheAllowedTime[:-1]) for t in timeMap if cacheAllowedTime.endswith(t)}
            )
            tTBM = time.time() - cacheAllowedTime.total_seconds() # time To Be Modified
            #time.strftime('%c',time.gmtime(tTBM))
            if cache and os.path.exists(cachefile) and os.path.getmtime(cachefile) > tTBM: #within allowed cache time
                return importer(cachefile)
            elif cache:
                # Check if cachedir exists
                if not os.path.exists(cachedir):
                    raise FileNotFoundError(
                        "LSD cache dir ({}) does not exist. Create, change in config or run function with cache=False".format(cachedir)
                    )
                # Redirect stdout and stderr
                stouterr_redirect = StringIO()
                with redirect_stdout(stouterr_redirect), redirect_stderr(stouterr_redirect):
                    functionData = function(*args,**kwargs)
                stouterr_function = stouterr_redirect.getvalue().strip()
                if stouterr_function: print(stouterr_function)
                exporter(functionData,cachefile,stouterr_function)
                return functionData
            else:
                return function(*args,**kwargs)
            
        return wrapper
    return cachedecorator

#cacheableTable reads from and writes to csv, does not log function output
cacheableTable = cacheable(
    importer = lambda x: pd.read_csv(x,index_col=0),
    exporter = lambda x,y,z: x.to_csv(y),
    extension = '.csv'
)

def storeDatasetLocally(dataset_getfunction):
    """
    Can be used as a decorator for 'get_dataset' functions.
    It will check if a processed dataset is locally available,
    and if so, load that one instead of processing from the source
    files.

    Should only be used for functions that do not process the data
    differently depending on the 'get_dataset' function arguments.
    The wrapper raises a warning if there are arguments to pass to 
    the 'get_dataset' function.

    As opposed to cacheable decorated functions, storeDatasetLocally
    is intended for volatile source code, such as e.g. your own scripts
    for which it is nonetheless useful being able to cache results. 
    """
    import inspect, hashlib, re
    from plumbum import colors
    dependency = re.compile(r'\W*Dependenc(y|ies): (.+)')
    
    def wrapper(*args, verbose=True, **kwargs):
        if args or kwargs:
            import warnings
            warnings.warn(
                'This decorated function is not designed to use with arguments. Be warned!'
            )
        # Check if data was already processed
        ## Prepare hash
        functionSource = inspect.getsource(dataset_getfunction).encode()
        ### Check if dependencies
        docstr = inspect.getdoc(dataset_getfunction)
        if docstr:
            dependencies = [e for d in (d.groups()[1].split()
                                        for d in (dependency.search(l)
                                                  for l in docstr.split('\n')) if d)
                            for e in d]
            for d in dependencies:
                d = d.split('.')
                d = getattr(globals()[d[0]],d[-1],globals()[d[0]]) #hack to get d with globals and getattr
                try:
                    functionSource+=inspect.getsource(d).encode()
                except TypeError:
                    functionSource+=pickle.dumps(d)
        hashvalue = hashlib.md5(functionSource).hexdigest()
        datastorage = '{}{}.pickle'.format(processedDataStorage,
                                              dataset_getfunction.__name__.replace('get_','')
        )
        
        if exists(datastorage):
            with open(datastorage,'rb') as openedDatastorage:
                datasetrepo = pickle.load(openedDatastorage)
            if datasetrepo.currentHash == hashvalue:
                print(datasetrepo.report)
                if verbose:
                    print(colors.green & colors.bold | 'Repo size {:.1f}MB, archive contains {} other versions'.format(
                        os.stat(datastorage).st_size/1024**2,
                        len(datasetrepo.archive)
                    ))
                return datasetrepo.dataset
            else:
                print(colors.cyan | 'Dataset content out of date, updating:')
                updateRepo = True
        else:
            print(colors.cyan | 'Dataset not locally available, generating:')
            updateRepo = False

        # Redirect stdout and stderr
        stouterr_redirect = StringIO()
        with redirect_stdout(stouterr_redirect), redirect_stderr(stouterr_redirect):
            start = time.process_time()
            # Run the function that generates the dataset
            dataset = dataset_getfunction(*args, **kwargs)
            duration = time.process_time() - start
            print('Dataset',datastorage,'generated',time.strftime('%c'))
            print('Processing time',time.strftime('%H:%M:%S', time.gmtime(duration)))
            
        print(stouterr_redirect.getvalue())

        if updateRepo:
            datasetrepo.update(dataset,functionSource.decode(),stouterr_redirect.getvalue().strip())
        else:
            try: DatasetRepo(dataset,functionSource.decode(),stouterr_redirect.getvalue().strip(),datastorage)
            except FileNotFoundError:
                print('Not possible to store dataset locally. Create',processedDataStorage,
                      'if you want to avoid reprocessing dataset on every call.')
            
        return dataset

    return wrapper

## Datasets
### References/annotations
from .dealer.ensembl import get_ensembl, get_ensemblGeneannot

@retrieveSources
def get_liftover(frm=19,to=38):
    """
    Info: http://hgdownload.cse.ucsc.edu/downloads.html
    """
    from pyliftover import LiftOver
    liftoverfile = 'hg{}ToHg{}.over.chain.gz'.format(frm,to)
    try: return LiftOver(processedDataStorage+liftoverfile)
    except FileNotFoundError:
        raise FileNotFoundError('Source: http://hgdownload.cse.ucsc.edu/gbdb/hg{}/liftOver/{}'.format(frm,liftoverfile))

def get_lift19to38():
    """
    Source: http://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
    """
    return get_liftover(frm=19,to=38)

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
    
#@retrieveSources -> not working login required
def get_msigdb6():
    """
    Source: http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.2/msigdb_v6.0.xml
    """
    import xml.etree.ElementTree as ET
    import pickle
    parser = ET.parse(processedDataStorage+'msigdb_v6.0.xml')
    root = parser.getroot()
    genesetsCollections = {} 
    for geneset in root:
        if (geneset.attrib['CATEGORY_CODE'] == 'ARCHIVED' or
            geneset.attrib['ORGANISM'] != 'Homo sapiens'): continue
        try:
            genesetsCollections[geneset.attrib['CATEGORY_CODE']
            ][geneset.attrib['STANDARD_NAME']] = geneset.attrib['MEMBERS_SYMBOLIZED'].split(',')
        except KeyError as e:
            genesetsCollections[geneset.attrib['CATEGORY_CODE']] = {}
            genesetsCollections[geneset.attrib['CATEGORY_CODE']
            ][geneset.attrib['STANDARD_NAME']] = geneset.attrib['MEMBERS_SYMBOLIZED'].split(',')
    genesetsCollections = {'version':'6','MSigDB':genesetsCollections}        
    print('MSigDB {}'.format(genesetsCollections['version']))
    mdb = genesetsCollections['MSigDB']
    return mdb

