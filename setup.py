import setuptools
from distutils.core import setup, Extension

module1 = Extension('libtimetag',
                    sources = ['./src/algos.cpp', './src/getline.cpp', './src/python_bindings.cpp', './src/sstt_file.cpp', './src/sstt_file2.cpp'], 
                    extra_compile_args=['-std=c++11','-stdlib=libc++'],
                    include_dirs = ['.','./include'],
					define_macros=[('LIBTIMETAG_COMPILE_PYTHON', None), ('BUILDING_LIBTIMETAG', None)])

setup (name = 'libtimetag',
       version = '0.8',
       author='Stijn Hinterding',
       author_email='git@stijnhinterding.nl',
       description = 'Library for storage and processing of time-correlated single-photon counting (TCSPC) data.',
       ext_modules = [module1])
