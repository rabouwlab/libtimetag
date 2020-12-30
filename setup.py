import setuptools
from setuptools import setup, Extension, find_packages

module1 = Extension('_libtimetag',
                    sources = ['./src/algos.cpp', './src/getline.cpp', './src/python_bindings.cpp', './src/sstt_file.cpp', './src/sstt_file2.cpp'], 
                    extra_compile_args=['-std=c++11','-stdlib=libc++'],
                    include_dirs = ['.','./include'],
					define_macros=[('LIBTIMETAG_COMPILE_PYTHON', None), ('BUILDING_LIBTIMETAG', None)])

with open("README.md","r", encoding="utf-8") as fh:
    long_description = fh.read()

install_requires = [
    'numpy',
    'python-dateutil',
    ]
    
setup (name = 'libtimetag',
       version = '0.9',
       author='Stijn Hinterding',
       author_email='git@stijnhinterding.nl',
       description = 'Library for storage and processing of time-correlated single-photon counting (TCSPC) data.',
       long_description=long_description,
       packages=['libtimetag'],
       package_dir={'libtimetag':'src'},
       long_description_content_type="text/markdown",
       url="https://www.github.com/rabouwlab/libtimetag",
       classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            ],
        python_requires='>=3.6',
        install_requires=install_requires,
       ext_modules = [module1])
