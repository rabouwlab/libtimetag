/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

/**
 * \file    sstt_file2.h
 * \brief   Defines and implements the "small simple time-tagged" (SSTT) file format, read-only, version 2.
 * \author  Stijn Hinterding
*/

/*
	This file format uses 6 bytes to represent a time-tag 'event'.
	In each event, the first two bits signify the type of event:
		bit #0:		1: overflow event, 0: not an overflow event
		bit #1:		<reserved, not used in current implementation>
		
	The next 46 bit store the event time, in units of the intrinsic
	time unit of the time-to-digital converter used
	(e.g., for QuTools quTAG: 1 ps; for QuTools quTAU: 81 ps).
	
	Since there are only 46 bits available to represent an
	event time, at some point during the experiment the time
	counter will likely overflow. In very rare cases, the
	time counter may overflow multiple times in between time-tag
	events. Any overflow events are communicated by setting
	the first bit to 1. The 46 data bits then hold the number
	of overflows that have occurred since the last event.
*/


#ifndef SSTT_FILE2_H
#define SSTT_FILE2_H

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

#define SSTT2_N_BYTES_TOT         6
#define SSTT2_N_BYTES_HEADER     (SSTT2_N_BYTES_TOT*3)
#define SSTT2_MAGIC_INFO         "Simple Small Time Tagged (V2)\n"
#define SSTT2_MAGIC              "SSTT2\0"
#define SSTT2_N_BITS_TOT         (SSTT2_N_BYTES_TOT*8)
#define SSTT2_N_BITS_SIGNAL      2
#define SSTT2_N_BITS_MACRO       (SSTT2_N_BITS_TOT-SSTT2_N_BITS_SIGNAL)
#define SSTT2_N_BITS_OVERFLOW    (SSTT2_N_BITS_TOT-SSTT2_N_BITS_SIGNAL)
#define SSTT2_MASK_SIGNAL        (((uint64_t)1 << SSTT2_N_BITS_SIGNAL) - 1)
#define SSTT2_MASK_MACRO         (((uint64_t)1 << SSTT2_N_BITS_MACRO) - 1)
#define SSTT2_MASK_OVERFLOW      (((uint64_t)1 << SSTT2_N_BITS_OVERFLOW) - 1)
#define SSTT2_OVERFLOW_VAL       ((uint64_t)1 << SSTT2_N_BITS_MACRO)

struct channel_info_sstt2
{
public:
    uint64_t ID;
    uint64_t n_photons;
    std::string filename;

    bool channel_has_microtime;
    bool is_pulses_channel;
    bool has_pulses_channel;
    uint64_t corresponding_pulses_channel;
    uint64_t sync_divider;
    uint64_t additional_sync_divider;

    channel_info_sstt2() :
        ID(0),
        n_photons(0),
        filename(),
        channel_has_microtime(false),
        is_pulses_channel(false),
        has_pulses_channel(false),
        corresponding_pulses_channel(0),
        sync_divider(1),
        additional_sync_divider(1)
    {
    }
};

struct exp_info_sstt2
{
public:
    double time_unit_seconds;
    std::string device_type;
};

int LIBTIMETAG_DLL test_is_sstt2_info_file(const std::string& filepath);

int LIBTIMETAG_DLL test_is_sstt2_file(const std::string& filepath);

int LIBTIMETAG_DLL read_data_file_sstt2(const std::string &filepath,
                   std::vector<int64_t> *macrotimes, uint64_t n_events_to_skip,
                                        uint64_t n_overflows_had, uint64_t *n_overflows_in_file);

int LIBTIMETAG_DLL n_photons_in_datafile_sstt2(const char* directory,
                          const char* filename,
                          uint64_t* ret);

std::vector<channel_info_sstt2> LIBTIMETAG_DLL get_sstt2_info(const char* filename, int *error_code, exp_info_sstt2* exp_info);

#endif // SSTT_FILE2_H
