# About
Libtimetag is a library for opening time-correlated single-photon counting data stored in the "small simple time-tagged" (SSTT) file format, and also provides algorithms to process these data, such as fast cross-correlation functions. These functions are fast enough that anti-bunching and fluorescence-correlation spectroscopy (FCS) cross-correlation curves can be visualized live; no hardware correlator is required.

The library was made for two use cases: (i) as a supporting library for other (C/C++/...) applications (e.g., [PHAST](https://www.github.com/rabouwlab/phast)), and (ii) as Python library for data processing and analysis. It is used on a daily basis for data acquisition, processing, and analysis at the lab of [Freddy Rabouw](https://www.uu.nl/medewerkers/FTRabouw) at Utrecht University. The library was originally developed by Stijn Hinterding.


# Installing as Python module
Simply install from pypi using pip:

    pip install libtimetag

	
# Usage as Python module
Please see the [python notebook](https://github.com/rabouwlab/libtimetag/blob/master/examples/libtimetag_example-notebook.ipynb).


# Installing as C++ library
Libtimetag was developed on Windows using [MSYS2](https://www.msys2.org) (more specifically, the MSYS2 MinGW 64-bit toolchain). The following should work under Windows (MSYS2), and Linux-like OS's. Inside a shell, situated in the project root directory, execute:

	mkdir build
	cd build
	cmake ..
	make
	make install

	
# Acknowledgements
This work was supported by The Netherlands Center for Multiscale Catalytic Energy Conversion (MCEC), an NWO Gravitation Programme funded by the Ministry of Education, Culture and Science of the government of The Netherlands.


# License
This project is licensed under the MIT license. See the LICENSE file in the project root.

