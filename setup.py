import setuptools
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
from os.path import join
from setuptools.config import read_configuration

py_config = read_configuration('setup.cfg')
fort_config = Configuration('mmmg_tools')

#-- pyBader compilation setup -----------------------------------------#
lib = [
    join('src','pybader','kind_mod.f90'),
    join('src','pybader','matrix_mod.f90'),
    join('src','pybader','charge_mod.f90'),
    join('src','pybader','bader_mod.f90'),
    join('src','pybader','weight_mod.f90'),
    join('src','pybader','interface_mod.f90'),
]
src = [
        join('src','pybader','pybader.pyf'),
        join('src','pybader','pybader.f90'),
]
fort_config.add_library('bader',
        sources = lib,
        extra_f90_compile_args = ['-O2']
)
fort_config.add_extension('pybader',
        sources = src,
        libraries = ['bader'],
        depends = lib,
        include_dirs = [join('build','temp')],
)
#-- End pyBader -------------------------------------------------------#

setup(**{
        **py_config['metadata'],
        **py_config['options'],
        **fort_config.todict(),
})
