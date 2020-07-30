# -*- coding: utf-8 -*-
"""LSD configuration

This module organizes the configuration for the lostdata package

Todo:
    * Give example symlinking privatedir ln -s ~/Dropbiz/Lab/z_archive/Datasets ~/LSData/private
"""

import configparser, os, appdirs

configdir = appdirs.AppDirs('lostdata').user_config_dir
appdatadir = appdirs.AppDirs('lostdata').user_data_dir

# Check existance of app dirs
if not os.path.exists(configdir):
    os.makedirs(configdir)
    print('created lostdata configdir', configdir)
if not os.path.exists(appdatadir):
    os.makedirs(appdatadir)
    print('created lostdata data dir', appdatadir)

configFileOptions = [
    'lsd.cfg', # in current working dir
    os.path.join(configdir, 'lsd.cfg')
]

# Default configuration
config = configparser.ConfigParser()
config['LSD'] = {
    'cachetime': '4w', #supports w[eeks], d[ays] or h[ours]
    'cachedir': os.path.join(appdatadir, 'cache'),
    'privatedir': os.path.join(appdatadir, 'private')
}

# Read configuration file
config.read(configFileOptions)

# Check cache and privatedir
if not os.path.exists(config['LSD']['cachedir']):
    os.makedirs(config['LSD']['cachedir'])
    print('created lostdata cache dir', config['LSD']['cachedir'])
if not os.path.exists(config['LSD']['privatedir']):
    os.makedirs(config['LSD']['privatedir'])
    print('created lostdata private dir', config['LSD']['privatedir'])

# Secrets: config for storing user API keys and other sensitive/personal information
from kindi import Secrets
secrets = Secrets(default_section=__package__)

def makeLSDataLinks():
    if not os.path.exists(os.path.expanduser('~/LSData')):
        os.makedirs(os.path.expanduser('~/LSData'))
    os.symlink(
        config['LSD']['cachedir'],
        os.path.expanduser('~/LSData/cache/')
    )
    os.symlink(
        config['LSD']['privatedir'],
        os.path.expanduser('~/LSData/private/')
    )
