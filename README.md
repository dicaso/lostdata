# LOcalized STructured DAta TAbles

Package for retrieving useful biological datasets. Datasets are cached
on a local directory, and served in a structured format, most
frequently a pandas DataFrame.

## Install

    pip3 install lostdata

## Configuration

To make use of lostdata, it expects a few directories where datasets
can be stored. You either have to create those directories, or
redirect the lostdata config file to existing ones.

Create default directories:

    mkdir -p ~/LSData/{cache,private}

Create config file `~/.lsd.cfg` with custom chosen directories:

    [LSD]
    cachetime=4w
    cachedir=~/LSData/cache/
    privatedir=~/LSData/private/


