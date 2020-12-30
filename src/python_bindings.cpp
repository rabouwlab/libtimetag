/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

/**
 * \file    python_bindings.cpp
 * \brief   Python bindings for this project
 * \author  Stijn Hinterding
*/

#ifdef LIBTIMETAG_COMPILE_PYTHON

#include <iostream>
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

#include "sstt_file.h"
#include "sstt_file2.h"
#include "algos.h"

namespace py = pybind11;

typedef void destr(void);

PYBIND11_MODULE(_libtimetag, m) {

    m.doc() = "Library for opening and processing Time-Correlated Single-Photon Counting data\n"
			  "\n"
			  "This module can open and process small simple time-tagged (SSTT) time-correlated\n"
			  "single-photon counting (TCSPC) datasets, and also provides algorithms to process\n"
			  "these data, such as fast cross-correlation functions.\n"
			  "\n"
			  "Data is most easily imported using the import_data() function, as this imports\n"
			  "the data as well as the header information, and does basic pre-processing.\n"
			  "More advanced use-cases may benefit from the read_sstt_data() and gen_micro_times()\n"
			  "functions.";

    m.def("gen_micro_times", [](const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& pulse_times,
          const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& data_times,
          uint64_t total_sync_divider) -> py::array {

        uint64_t ret_size = data_times.size();
        int64_t* ret = new int64_t[ret_size]{0};


        int success = gen_microtimes(pulse_times.data(0),
                                       pulse_times.size(),
                                       data_times.data(0),
                                       data_times.size(),
                                       ret,
                                       ret_size,
                                       total_sync_divider);
        if (success != 0) {
            delete[] ret;
        }

        if (success == 1) {
            throw std::runtime_error("Invalid input");
        }

        if (success == 2) {
            throw std::runtime_error("Internal error :-(");
        }

        auto capsule = py::capsule(ret, [](void *v) { delete[] (int64_t*)v; });


        return py::array(ret_size, ret, capsule);
    },"Generate micro timestamps (timestamps relative to a reference channel)\n"
      "\n"
	  "Note: if you use import_data() you will most likely not need this function\n"
	  "\n"
      "Parameters\n"
      "----------\n"
      "ref_timestamps : array_like\n"
      "     Array containing timestamps of the reference channel, e.g. the laser sync channel.\n"
      "data_timestamps : array_like\n"
      "     Array containing timestamps of the channel for which the micro timestamps should be calculated\n"
      "total_sync_divider : positive integer\n"
      "     The total sync divider that was applied to the reference channel during data acquisition. The sync divider\n"
      "     determines how many of the recorded events are discarded; where the ratio total_num_events/total_sync_divider\n"
      "     gives the number of events which are not discarded. For example: with a sync divider of unity, all events\n"
      "     are recorded; with a sync divider of 2, every other event is recorded; with a sync divider of 3, every\n"
      "     third event is recorded, et cetera.\n"
      "\n"
      "Returns\n"
      "-------\n"
      "microtimestamps : array_type\n"
      "     Returns an array containing the micro timestamps, corresponding to the combination of the supplied reference\n"
      "     and data channel.",
        py::arg("ref_timestamps"), py::arg("data_timestamps"), py::arg("total_sync_divider"));

    m.def("read_sstt_data", [](std::string filepath,uint64_t n_photons_to_skip, uint64_t n_overflow_events) -> py::tuple{
        std::vector<int64_t>* macrotimes = new std::vector<int64_t>();
        std::vector<int64_t>* microtimes = new std::vector<int64_t>();
        uint64_t n_overflows = 0;
        int success = 0;

        if (test_is_sstt2_file(filepath)) {
            // py::print("Reading SSTT V2 file");

            success = read_data_file_sstt2(filepath, macrotimes, n_photons_to_skip, n_overflow_events, &n_overflows);
        } else {
            // py::print("Reading SSTT V1 file");

            success = read_data_file(filepath, macrotimes, microtimes);
        }

        if (success != 0) {
            delete macrotimes;
            delete microtimes;
        }

        if (success == 1) {
            throw std::runtime_error("Failed to open file '" + std::string(filepath) + "'");
        }

        if (success == 3) {
            throw std::runtime_error("Did not recognize file format as either SSTT v1 or v2!");
        }

        if (success != 0) {
            throw std::runtime_error("Unknown error");
        }

        auto capsule_macro = py::capsule(macrotimes, [](void *v) { delete reinterpret_cast<std::vector<uint64_t>*>(v); });
        py::array py_macrotimes(macrotimes->size(), macrotimes->data(), capsule_macro);

        auto capsule_micro = py::capsule(microtimes, [](void *v) { delete reinterpret_cast<std::vector<uint64_t>*>(v); });
        py::array py_microtimes(microtimes->size(), microtimes->data(), capsule_micro);

        return py::make_tuple(py_macrotimes, py_microtimes,n_overflows);
    },"Reads in data from a single *.sstt.c* data file\n"
    "\n"
	"Note: the import_data() function is generally more convenient to use.\n"
	"\n"
    "Parameters\n"
    "----------\n"
    "filepath : string\n"
    "     Path to the *.sstt.c* data file to open.\n"
    "n_photons_to_skip : uint64 (optional)\n"
    "     The number of photon events to skip when\n"
    "     reading this file. Useful when reading in\n"
    "     a file which is still being updated. Be\n"
    "     careful to also specify the correct number\n"
    "     of overflow events.\n"
    "n_overflow_events : uint64_t (optional)\n"
    "     The number of overflow events already\n"
    "     encountered in this file. Should be used\n"
    "     in combination with n_photons_to_skip\n"
    "\n"
    "Returns\n"
    "-------\n"
	"py_macrotimes : list\n"
	"		List of macro timestamps stored in the data file.\n"
	"py_microtimes : list\n"
	"		List of micro timestamps stored in the data file.\n"
	"		Only legacy SSTT files (v1) save the microtimes\n"
	"		explicitely, so this list is likely to be empty.\n"
	"		Generate the microtimes using the\n"
	"		gen_micro_times() function.\n"
	"n_overflows : integer\n"
	"		Number of overflow events encountered while\n"
	"		reading the data file. Use this information,\n"
	" 		if desired, in subsequent calls to\n"
	"		read_sstt_data(), to read in only a portion\n"
	"		of the data file.",
    py::arg("filepath"),py::arg("n_photons_to_skip")=0,py::arg("n_overflow_events")=0);

    m.def("correlate_fcs", [](const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& bin_edges,
                                const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& left_list,
                                    const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& right_list) -> py::array {
        if (bin_edges.size() <= 1) {
            throw std::runtime_error("bin_edges should have a minimum length of two");
        }

        uint64_t ret_size = bin_edges.size() - 1;
        int64_t* ret = new int64_t[ret_size]{0};
        auto capsule = py::capsule(ret, [](void *v) { delete[] (int64_t*)v; });

        if (left_list.size() == 0 || right_list.size() == 0) {
            return py::array(ret_size, ret, capsule);
        }

        int success = correlate_many_per_bin(bin_edges.data(0),
                                        bin_edges.size(),
                                        left_list.data(0),
                                        left_list.size(),
                                        right_list.data(0),
                                        right_list.size(),
                                        ret,
                                        bin_edges.size() - 1);

        if (success == 1) {
            throw std::runtime_error("Internal error #1");
        } else if (success == 2) {
            throw std::runtime_error("bin_edges should have a minimum length of two");
        } else if (success == 3) {
            throw std::runtime_error("Internal error #3");
        } else if (success != 0) {
            throw std::runtime_error("Unknown error");
        }

        return py::array(ret_size, ret, capsule);
    }, "Cross-correlates two arrays. Optimized for cases with many photons per bin.\n"
    "\n"
	"Correlates two arrays containing timestamped data. Optimized for cases.\n"
	"where there will be many photons per bin in the correlation histogram,\n"
	"or when correlating datasets over large lag times (e.g., up to seconds).\n"
	"Useful for Fluorescence Correlation Spectroscopy (FCS) data sets.\n"
	"Note that the output of this function should be normalized using the\n"
	"norm_corr() function."
	"\n"
    "Parameters\n"
    "----------\n"
    "bin_edges : list\n"
    "     Edges for the correlation histogram. Size of bins is allowed to vary\n"
	"  	  within the histogram.\n"
    "left_array : list\n"
    "     List containing the timestamps of the 'left' dataset\n"
    "right_array : uint64_t\n"
    "     List containing the timestamps of the 'right' dataset\n"
    "\n"
    "Returns\n"
    "-------\n"
    "data : list\n"
    "     List containing the non-normalized cross-correlation histogram.",
    py::arg("bin_edges"), py::arg("left_array"), py::arg("right_array"));

    m.def("correlate_lin", [](const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& bin_edges,
                                    const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& left_list,
                                    const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& right_list) {
        if (bin_edges.size() <= 1) {
            throw std::runtime_error("bin_edges should have a minimum length of two");
        }

        uint64_t ret_size = bin_edges.size() - 1;
        std::vector<int64_t> ret(ret_size, 0);

        if (left_list.size() == 0 || right_list.size() == 0) {
            return ret;
        }

        int success = correlate_unit_bins(bin_edges.data(0),
                                            bin_edges.size(),
                                            left_list.data(0),
                                            left_list.size(),
                                            right_list.data(0),
                                            right_list.size(),
                                            ret.data(),
                                            ret_size);

        if (success == 1) {
            throw std::runtime_error("Internal error #1");
        } else if (success == 2) {
            throw std::runtime_error("bin_edges should have a minimum length of two");
        } else if (success == 3) {
            throw std::runtime_error("Internal error #3");
        } else if (success == 4) {
            throw std::runtime_error("Bins should have a size of unity");
        } else if (success != 0) {
            throw std::runtime_error("Unknown error");
        }

        return ret;
    }, "Cross-correlates two arrays. Optimized for small lag times/few photons per bin.\n"
    "\n"
	"Cross-correlates two arrays containing timestamped data. Optimized for cases.\n"
	"where there will be few photons per bin in the correlation histogram,\n"
	"or when correlating datasets over short lag times.\n"
	"Useful for generating cross-correlation functions to check for anti-bunching.\n"
	"Note that the output of this function is not normalized. If normalization is\n"
	"desired, use the norm_corr() function."
	"\n"
    "Parameters\n"
    "----------\n"
    "bin_edges : list\n"
    "     Edges for the correlation histogram. Bin width must be unity.\n"
    "left_array : list\n"
    "     List containing the timestamps of the 'left' dataset\n"
    "right_array : uint64_t\n"
    "     List containing the timestamps of the 'right' dataset\n"
    "\n"
    "Returns\n"
    "-------\n"
    "data : list\n"
    "     List containing the non-normalized cross-correlation histogram.",
    py::arg("bin_edges"), py::arg("left_array"), py::arg("right_array"));

    m.def("norm_corr", [](const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& data,
                          const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& bin_edges,
                          uint64_t T_min,
                          uint64_t T_max,
                          uint64_t n_photons_left_channel,
                          uint64_t n_photons_right_channel) {
        if (bin_edges.size() <= 1) {
            throw std::runtime_error("bin_edges should have a minimum length of two");
        }

        if (data.size() != bin_edges.size() - 1) {
            throw std::runtime_error("histogram should be exactly one element shorter than bin_edges");
        }

        double* ret = new double[data.size()]{0};
        auto capsule = py::capsule(ret, [](void *v) { delete[] (double*)v; });

        int success = normalize_correlation(data.data(0),
                                            data.size(),
                                            bin_edges.data(0),
                                            bin_edges.size(),
                                            T_min,
                                            T_max,
                                            n_photons_left_channel,
                                            n_photons_right_channel,
                                            ret);

        if (success == 1) {
            throw std::runtime_error("Internal error #1");
        } else if (success != 0) {
            throw std::runtime_error("Unknown error");
        }

        return py::array(data.size(), ret, capsule);
    }, "Normalizes a cross-correlation histogram\n"
    "\n"
	"Cross-correlation functions generated using\n"
	"correlate_fcs() and correlate_lin() are not normalized.\n"
	"This means that the correlation amplitudes are effectively\n"
	"arbitrary. This function normalizes cross-correlation\n"
	"functions, so that for two completely non-correlated\n"
	"signals the correlation amplitude is unity, and for lag\n"
	"times where there are no photon counts, the correlation\n"
	"amplitude is zero. Correlation curves will then tend to\n"
	"zero in the case of anti-bunching, and to unity for\n"
	"long lag times. Some fields use a different definition\n"
	"of the cross-correlation function, where it tends to\n"
	"zero for non-correlated signals, and to -1 for lag\n"
	"times where there is no signal. If this convention\n"
	"is desired, simply subtract 1 from the values returned\n"
	"by this function.\n"
	"\n"
    "Parameters\n"
    "----------\n"
    "data : list\n"
    "	Non-normalized cross-correlation histogram.\n"
    "bin_edges : list\n"
    " 	The bin edges of the cross-correlation histogram\n"
	"T_min : positive integer\n"
	"	Time of experiment start, or minimum timestamp\n"
	"	value within the entire dataset.\n"
	"T_max : positive integer\n"
	"	Time of experiment end, or maximum timestamp\n"
	"	value within the entire dataset.\n"
	"n_photons_left_channel : positive integer\n"
	"	Number of timestamps in the 'left' data set."
	"n_photons_right_channel : positive integer\n"
	"	Number of timestamps in the 'right' data set."
    "\n"
    "Returns\n"
    "-------\n"
    "data : list\n"
    "     List containing the normalized cross-correlation histogram.",
       py::arg("data"), py::arg("bin_edges"), py::arg("T_min"), py::arg("T_max"), py::arg("n_photons_left_chan"), py::arg("n_photons_right_chan"));

    m.def("rebin", [](const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& data,
                      uint64_t new_bin_size) -> py::array {
        if ((int64_t)new_bin_size > data.size()) {
            throw std::runtime_error("n cannot be larger than the total number of bins");
        }

        uint64_t remainder = data.size() % new_bin_size;
        uint64_t ret_size = (data.size() - remainder) / new_bin_size;

        if (ret_size < 1) {
            throw std::runtime_error("Invalid n: the resulting histogram would have not even have a single bin");
        }

        int64_t* ret = new int64_t[ret_size]{0};
        auto capsule = py::capsule(ret, [](void *v) { delete[] (uint64_t*)v; });

        int success = rebin(data.data(0),
                    data.size(),
                    new_bin_size,
                    ret,
                    ret_size);

        if (success != 0) {
            throw std::runtime_error("Internal error");
        }


        return py::array(ret_size, ret, capsule);
    }, "Rebins a histogram\n"
    "\n"
	"Returns a new histogram in which the value of each new bin\n"
	"is the sum of n original bins (n >= 1). Any original bins\n"
	"that together do not make up an entire new bin will be\n"
	"discarded."
	"\n"
    "Parameters\n"
    "----------\n"
    "data : list\n"
    "	Original histogram\n"
    "new_bin_size : positive integer\n"
    " 	New bin width, i.e., how many original bins are\n"
	"	combined to make a bin in the new histogram\n"
    "\n"
    "Returns\n"
    "-------\n"
    "ret : list\n"
    "     Re-binned histogram.",
    py::arg("histogram"), py::arg("n"));

    m.def("rebin_bin_edges", [](const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& bin_edges, uint64_t new_bin_size) -> py::array {
        if (bin_edges.size() <= 1) {
            throw std::runtime_error("bin_edges should have a minimum length of two");
        }

        if ((int64_t)new_bin_size > bin_edges.size() - 1) {
            throw std::runtime_error("n cannot be larger than the total number of bins");
        }

        uint64_t remainder = (bin_edges.size() - 1) % new_bin_size;
        uint64_t ret_size = (bin_edges.size() - 1 - remainder) / new_bin_size + 1;

        if (ret_size <= 1) {
            throw std::runtime_error("Invalid n: the resulting histogram would have not even have a single bin");
        }

        int64_t* ret = new int64_t[ret_size]{0};
        auto capsule = py::capsule(ret, [](void *v) { delete[] (uint64_t*)v; });

        int success = rebin_bin_edges(bin_edges.data(0),
                        bin_edges.size(),
                        new_bin_size,
                        ret,
                        ret_size);

        if (success == 1) {
            throw std::runtime_error("Internal error #1");
        } else if (success == 2) {
            throw std::runtime_error("bin_edges should have a minimum length of two");
        } else if (success == 3) {
            throw std::runtime_error("Internal error #3");
        } else if (success != 0) {
            throw std::runtime_error("Unknown error");
        }

        return py::array(ret_size, ret, capsule);
    }, "Rebins histogram bin edges\n"
    "\n"
	"Returns the bin edges of a histogram\n"
	"for which the value of each new bin\n"
	"is the sum of n original bins (n >= 1). Any original bins\n"
	"that together do not make up an entire new bin will be\n"
	"discarded."
	"\n"
    "Parameters\n"
    "----------\n"
    "bin_edges : list\n"
    "	Original bin edges\n"
    "new_bin_size : positive integer\n"
    " 	New bin width, i.e., how many original bins are\n"
	"	combined to make a bin in the new histogram\n"
    "\n"
    "Returns\n"
    "-------\n"
    "ret : list\n"
    "     Re-binned bin edges.",
    py::arg("bin_edges"), py::arg("n"));
}

#endif
