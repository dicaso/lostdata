if __name__ == '__main__':
    import argparse, inspect, LSD
    from rpy2.rinterface import RRuntimeError
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    #Activate automatic pandas/r conversion
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    parser = argparse.ArgumentParser(description=
            '''
Generate LSD dataset and store as RData.
datasets have to be provided with their full path e.g. LSD.get_centromeres
            '''
    )
    parser.add_argument('datasets', type=str, nargs='*', help='datasets to generate')
    parser.add_argument('--list', action='store_true', help='list available datasets')

    args = parser.parse_args()

    if args.list:
        LSD.listLSDatasets()

    for ds in args.datasets:
        if not ds.startswith('LSD.'):
            raise Exception('Not an LSD dataset. `'+ds+'`should start with `LSD.`')
        datasetname = ds.split('.')[-1].replace('get_','')
        filename_r = '{}{}.RData'.format(LSD.processedDataStorage,datasetname)

        # Generate/load dataset
        m = inspect.importlib.import_module(ds[:ds.rindex('.')])
        f = m.__dict__[ds.split('.')[-1]]
        ds = f()
        ds_r = ro.ListVector(
            {d:ds.__getattribute__(d) for d in ds.__datasets__}
        )
        ro.globalenv[datasetname] = ds_r
        ro.r.save(datasetname,file=filename_r)
        print('{} now available for use in R with `load("{}")`'.format(datasetname,filename_r))
