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

PYBIND11_MODULE(libtimetag, m) {

    m.doc() = "Contains many useful functions related to Time-Correlated Single-Photon Counting experiments"; // optional

    py::class_<channel_info>(m, "channel_info")
            .def(py::init())
            .def_readwrite("ID", &channel_info::ID)
			.def_readwrite("n_photons", &channel_info::n_photons)
            .def_readwrite("filename", &channel_info::filename);

    m.def("get_sstt_info", [](const char *filename) {
        int error = 0;
        std::vector<channel_info> infos = get_sstt_info(filename, &error);

        if (error == 1) {
            throw std::runtime_error("Failed to open file '" + std::string(filename) + "'");
        } else if (error == 2) {
            throw std::runtime_error("Channel data in sstt file appears malformed");
        } else if (error == 3) {
            throw std::runtime_error("Could not find channel data in sstt file");
        }

        py::list ret;

        for (uint64_t i = 0; i < infos.size(); i++) {
            ret.append(infos[i]);
        }

        return ret;

    }, "Obtain information on an *.sstt header file\n"
       "\n"
       "Parameters\n"
       "----------\n"
       "filepath : string\n"
       "    path to the *.sstt file to open\n"
       "\n"
       "Returns\n"
       "-------\n"
       "channel_infos : array_like\n"
       "    Return an array containing channel_info structs", py::arg("filepath"));

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
    },"Generate micro timestamps (timestamps relative to a reference channel, e.g. laser pulses)\n"
      "\n"
      "Parameters\n"
      "----------\n"
      "ref_timestamps : array_like\n"
      "     Array containing timestamps of the reference channel\n"
      "data_timestamps : array_like\n"
      "     Array containing timestamps of the channel of which the micro timestamps should be calculated\n"
      "total_sync_divider : positive integer\n"
      "     The total sync divider, applied to the reference channel, during data acquisition. The sync divider\n"
      "     determines how many of the recorded events are discarded; where the ratio total_num_events/total_sync_divider\n"
      "     gives the number of events which are not discarded. For example: with a sync divider of unity, all events\n"
      "     are recorded; with a sync divider of 2, every other event is recorded; with a sync divider of 3, every\n"
      "     third event is recorded, et cetera.\n"
      "\n"
      "Returns\n"
      "-------\n"
      "microtimestamps : array_type\n"
      "     Return an array containing the micro timestamps, corresponding to the combination of the supplied reference\n"
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
    },"Reads in data from a *.sstt.c* data file\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "filepath : string\n"
    "     Path to the *.sstt.c* data file to open.\n"
    "n_photons_to_skip : uint64\n"
    "     The number of photon events to skip when\n"
    "     reading this file. Usefule when reading in\n"
    "     a file which is still being updated. Be\n"
    "     careful to also specify the correct number\n"
    "     of overflow events!\n"
    "n_overflow_events : uint64_t\n"
    "     The number of overflow events already\n"
    "     encountered in this file. Should be used\n"
    "     in combination with n_photons_to_skip\n"
    "\n"
    "Returns\n"
    "-------\n"
    "data : tuple\n"
    "       Return a tuple consisting of two arrays. The first array corresponds to the macro timestamps\n"
    "       stored in the data file. The second array corresponds to the micro timestamps stored in the\n"
    "       data file. Depending on the SSTT file version, the data file may or may not contain\n"
    "       micro timestamps. If the data file does not contain any micro timestamps, the returned\n"
    "       micro timestamp array is empty. In such cases, where the macro timestamps of the desired\n"
    "       reference channel are available, the gen_micro_times() function can be used to generate\n"
    "       the micro timestamps",
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
    }, "Correlates two arrays. This function is optimised for the case where there are many photons per bin, e.g. in Fluorescence Correlation Spectroscopy",
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
    }, "Correlates two arrays. This function is optimised for the case where there are few photons per bin (e.g. in an anti-bunching curve). The size of each bin should be unity (otherwise an exception will be thrown).",
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
    }, "Normalises a photon correlation histogram, returns the normalised histogram.",
       py::arg("data"), py::arg("bin_edges"), py::arg("T_min"), py::arg("T_max"), py::arg("n_photons_left_chan"), py::arg("n_photons_right_chan"));

    m.def("bindata_interp_seq", [](const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& bin_edges,
                                    const py::array_t<int64_t,py::array::c_style|py::array::forcecast>& data) {
        if (bin_edges.size() <= 1) {
            throw std::runtime_error("bin_edges should have a minimum length of two");
        }

        int64_t* ret = new int64_t[bin_edges.size() - 1]{0};
        auto capsule = py::capsule(ret, [](void *v) { delete[] (int64_t*)v; });

        int success = bindata_interp_seq(bin_edges.data(0),
                                        bin_edges.size(),
                                        data.data(0),
                                        data.size(),
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

        return py::array(bin_edges.size() - 1, ret, capsule);
    }, "Bins the supplied data into the supplied histogram. This function is optimized for linear bins (e.g. with constant size) and might perform poorly on bins with variable size (e.g. logarithmic).",
    py::arg("bin_edges"), py::arg("data"));

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
    }, "Returns a new histogram wherein each new bin corresponds to n original bins (n >= 1). Any original bins which together do not make up an entire new bin will be dropped",
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
    }, "Returns new bin edges, wherein each new bin corresponds to n original bins (n >= 1). Any original bins which together do not make up an entire new bin will be dropped.",
    py::arg("bin_edges"), py::arg("n"));
}

#endif
