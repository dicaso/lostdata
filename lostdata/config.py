# -*- coding: utf-8 -*-
"""LSD configuration

This module organizes the configuration for the lostdata package

Todo:
    * Give example symlinking privatedir ln -s ~/Dropbiz/Lab/z_archive/Datasets ~/LSData/private
"""

import configparser, os
configFileOptions = [
    'lsd.cfg', # in current working dir
    os.path.expanduser('~/.lsd.cfg'),
    '/usr/local/etc/lsd.cfg'
]

# Default configuration
config = configparser.ConfigParser()
config['LSD'] = {
    'cachetime': '4w', #supports w[eeks], d[ays] or h[ours]
    'cachedir': os.path.expanduser('~/LSData/cache/'),
    'privatedir': os.path.expanduser('~/LSData/private/')
}

# Read configuration file
config.read(configFileOptions)

# Secrets: config for storing user API keys and other sensitive/personal information
from kindi import Secrets
secrets = Secrets(default_section=__package__)
