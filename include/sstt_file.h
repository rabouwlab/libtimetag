/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

/**
 * \file    sstt_file.h
 * \brief   Defines and implements the "small simple time-tagged" (SSTT) file format, version 1 (deprecated!).
 * \author  Stijn Hinterding
*/

/*
	This is the legacy SSTT file format (version 1), and is not 
	supported any more. 
*/

#ifndef SSTT_FILE_H
#define SSTT_FILE_H

#define SSTT_N_BITS_TOT         64
#define SSTT_N_BITS_SIGNAL      2
#define SSTT_N_BITS_MICRO       34
#define SSTT_N_BITS_MACRO       28
#define SSTT_N_BITS_OVERFLOW    62

#include <stdint.h>
#include <string>
#include <vector>

#ifdef _WIN32
#ifdef BUILDING_LIBTIMETAG
#define LIBTIMETAG_DLL __declspec(dllexport)
#else
#define LIBTIMETAG_DLL __declspec(dllimport)
#endif
#else
#define LIBTIMETAG_DLL
#endif

#define SSTT_MASK_SIGNAL        (((uint64_t)1 << SSTT_N_BITS_SIGNAL) - 1)
#define SSTT_MASK_MICRO         (((uint64_t)1 << SSTT_N_BITS_MICRO) - 1)
#define SSTT_MASK_MACRO         (((uint64_t)1 << SSTT_N_BITS_MACRO) - 1)
#define SSTT_MASK_OVERFLOW      (((uint64_t)1 << SSTT_N_BITS_OVERFLOW) - 1)

#define SSTT_OVERFLOW_VAL       ((uint64_t)1 << SSTT_N_BITS_MACRO)

struct channel_info
{
public:
    uint64_t ID;
    uint64_t n_photons;
    std::string filename;

    channel_info() :
        ID(0),
        n_photons(0),
        filename("")
    {
    }
};

int LIBTIMETAG_DLL read_data_file(const std::string &filepath,
                   std::vector<int64_t> *macrotimes,
                   std::vector<int64_t> *microtimes);

int LIBTIMETAG_DLL n_photons_in_datafile(const char* directory,
                          const char* filename,
                          uint64_t* ret);

std::vector<channel_info> LIBTIMETAG_DLL get_sstt_info(const char* filename, int *error_code);

#endif // SSTT_FILE_H
