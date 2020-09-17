# -*- coding: utf-8 -*-
"""Functions for processing data and lostdata Datasets
"""
from .config import config
import os, time
from os.path import expanduser, exists
from io import TextIOWrapper, StringIO
from contextlib import redirect_stdout, redirect_stderr
from collections import OrderedDict

processedDataStorage = config['LSD']['cachedir']

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
            dataset_getfunction_result = dataset_getfunction(*args, **kwargs)
        except FileNotFoundError as fnf:
            for docline in inspect.getdoc(dataset_getfunction).split('\n') + list(fnf.args):
                # Test if docline is str, as sometimes non-str arguments get thrown
                if isinstance(docline, str) and docline.startswith('Source:'):
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
                        print('Downloaded', filename, exists(os.path.join(processedDataStorage,filename)))
            try: dataset_getfunction_result = dataset_getfunction(*args, **kwargs)
            except FileNotFoundError:
                print('Either not all source files are documented correctly in docstring,',
                      'or there is a source file unrelated issue')
                raise

        # Returning the result produced by the wrapped function
        return dataset_getfunction_result

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
    from .formats import DatasetRepo
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

## Communication protocols
def download_ftp_resource(server,directory='.',filelist=None,filter=None,targetdir='.'):
    """Download FTP resource. Files that are already downloaded are skipped.

    In case of an interruption this does imply those files need to be manually removed.

    TODO: - instead of simply skipping downloaded file, check if it was
    downloaded correctly (e.g. +- same size)

    Args:
      server (str): FTP server location. Should not start with 'fpt://'
      directory (str): Directory on the remote ftp server.
      filelist (list of str): List of files to download.
      filter (re): Compiled regex for filtering files in `directory`.
      targetdir (str): Destination where files will be saved.
    """
    from ftplib import FTP
    ftp = FTP(server)
    ftp.login()
    ftp.cwd(directory)
    if not filelist:
        filelist = []
        ftp.retrlines('LIST', lambda x: filelist.append(x.split()[-1]))
    for f in filelist:
        if not filter or filter.search(f):
            foutpath = os.path.join(targetdir,f)
            if not os.path.exists(foutpath):
                print('downloading',f)
                ftp.retrbinary('RETR %s' % f, open(foutpath, 'wb').write)
    ftp.quit()
