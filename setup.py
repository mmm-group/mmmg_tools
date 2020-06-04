from os.path import join

from setuptools import setup
from setuptools.config import read_configuration

py_config = read_configuration('setup.cfg')

setup(**{
    **py_config['metadata'],
    **py_config['options'],
})
