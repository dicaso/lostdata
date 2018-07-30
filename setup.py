from setuptools import setup, find_packages

setup(name = 'lostdata',
      version = '0.0.2',
      description = 'LOcalized STructured DAta TAbles',
      url = 'https://github.com/dicaso/lostdata',
      author = 'Christophe Van Neste',
      author_email = 'christophe.vanneste@ugent.be',
      license = 'MIT',
      packages = find_packages(),
      python_requires='>3.6',
      install_requires = [
          #Generated with `pipreqs .` and then moved here from requirements.txt
          'numpy',
          'pandas',
          'xarray',
          'kindi',
          'networkx',
          'biomart',
          'gffutils',
          'pyliftover',
          'requests',
          'plumbum',
      ],
      extras_require = {
          #for netCDF4 output, depends on: libnetcdf-dev libhdf5-dev
          'netcdf': ['netCDF4'],
          'documentation': ['Sphinx']
      },
      zip_safe = False,
      #entry_points = {
      #    'console_scripts': ['getLSDataset=LSD.command_line:main'],
      #},
      test_suite = 'nose.collector',
      tests_require = ['nose']
)

#To install with symlink, so that changes are immediately available:
#pip install -e . 
