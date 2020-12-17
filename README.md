# About
Libtimetag is a library for opening time-correlated single-photon counting data stored in the "small simple time-tagged" (SSTT) file format, and also provides algorithms to process these data, such as fast cross-correlation functions. These functions are fast enough that anti-bunching and fluorescence-correlation spectroscopy (FCS) cross-correlation curves can be visualized live; no hardware correlator is required.

The library was made for two use cases: (i) as a supporting library for other (C/C++/...) applications, and (ii) as Python library for data processing and analysis. It is used on a daily basis for data acquisition, processing, and analysis at the lab of [Freddy Rabouw](https://www.uu.nl/medewerkers/FTRabouw) at Utrecht University. The library was originally developed by Stijn Hinterding.


# Installing as Python module
## Windows
Make sure you have a C++ compiler installed. If you have Visual Studio installed you should be good. Don't have a clue? Install the [Visual C++ 2015 Build Tools](http://go.microsoft.com/fwlink/?LinkId=691126&fixForIE=.exe). You are likely to run into plenty of errors along the way -- this document should be updated to give a more detailed explanation of the installation process.

While situated inside the project root directory, execute the following command (note: if you are using Anaconda, you should execute this command from inside the [Anaconda Prompt](https://docs.anaconda.com/anaconda/user-guide/getting-started/):

	python setup.py install

## Linux/MacOS/...
This project was developed on Windows, but known to work on MacOS as well.
While situated inside the project root directory, execute:

	python setup.py install

# Installing as C++ library
Libtimetag was developed on Windows using [MSYS2](www.msys2.org) (more specifically, the MSYS2 MinGW 64-bit toolchain). The following should work under Windows (MSYS2), and Linux-like OS's. Inside a shell, situated in the project root directory, execute:

	mkdir build
	cd build
	cmake ..
	make
	make install

# Basic usage as Python module
(To be written)

# License
This project is licensed under the MIT license. See the LICENSE file in the project root.

