/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

/**
 * \file    sstt_file2.cpp
 * \brief   Defines and implements the "small simple time-tagged" (SSTT) file format, read-only, version 2.
 * \author  Stijn Hinterding
*/

#include "sstt_file2.h"
#include "getline.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#define SSTT2_CHAN_HEADER_TEXT "CHANNEL_HEADER\n"

#define SSTT2_HEADER_DELIMITER       "\t\n"
#define SSTT2_HEADER_CHANID          "ChannelID"
#define SSTT2_HEADER_FILENAME        "Filename"
#define SSTT2_HEADER_NUMPHOTONS      "NumPhotons"
#define SSTT2_HEADER_SYNCDIV         "HardwareSyncDivider"
#define SSTT2_HEADER_ADDI_SYNCDIV    "AdditionalSyncDivider"
#define SSTT2_HEADER_IS_PULSES       "IsPulsesChannel"
#define SSTT2_HEADER_HAS_PULSES      "HasPulsesChannel"
#define SSTT2_HEADER_CORR_PULSECHAN  "CorrespondingPulsesChannel"

#define SSTT2_EXP_HEADER_TEXT "EXPERIMENT_HEADER\n"

#define SSTT2_HEADER_TIMEUNIT       "Time_unit_seconds"
#define SSTT2_HEADER_DEV_TYPE       "device_type"

struct processed_event_sstt2
{
    int64_t macrotime;
    uint64_t n_overflows;
};

struct processed_event_sstt2 handle_event_sstt2(int64_t e, uint64_t num_overflows)
{
    // Prepare the return value
    struct processed_event_sstt2 ret;
    ret.macrotime = 0;
    ret.n_overflows = 0;

    if ( (e & 1) && !( (e >> 1) & 1 ) ) {
        // Overflow event
        ret.n_overflows = (e >> SSTT2_N_BITS_SIGNAL) & SSTT2_MASK_OVERFLOW;
    } else if ( !(e & 1) && !( (e >> 1) & 1)) {
        // Photon event
        ret.macrotime = (e >> SSTT2_N_BITS_SIGNAL) & SSTT2_MASK_MACRO;

        // Also account for the overflows
        ret.macrotime += num_overflows * SSTT2_OVERFLOW_VAL;
    }

    return ret;
}

int LIBTIMETAG_DLL n_photons_in_datafile_sstt2(const char* directory,
                          const char* filename,
                          uint64_t* ret)
{
    if (ret == NULL) {
        return 1;
    }

    FILE* f = fopen((std::string(directory) + std::string(filename)).c_str(), "rb");

    if (f == NULL) {
        // Error opening the file
        return 2;
    }

    // Read in the header
    std::vector<char> header(18,'\0');

    size_t n_read = fread(header.data(), SSTT2_N_BYTES_HEADER, 1, f);

    uint64_t n_overflows = 0;
    uint64_t event = 0;
   n_read = fread(&event, SSTT2_N_BYTES_TOT, 1, f);

    while (n_read == 1) {
        struct processed_event_sstt2 pe = handle_event_sstt2(event, n_overflows);

        if (pe.n_overflows == 0)
            (*ret)++;

        n_read = fread(&event, SSTT2_N_BYTES_TOT, 1, f);
    }

    fclose(f);

    return 0;
}

int LIBTIMETAG_DLL test_is_sstt2_file(const std::string &filepath)
{
    FILE* f = fopen(filepath.c_str(), "rb");

    if (f == NULL) {
        // Error opening the file
        return 0;
    }

    // Read in the header
    std::vector<char> header(18,'\0');

    fread(header.data(), SSTT2_N_BYTES_HEADER, 1, f);

    if (std::string(header.data()) == std::string(SSTT2_MAGIC)) {
        fclose(f);
        return 1;
    }

    fclose(f);
    return 0;
}

int LIBTIMETAG_DLL read_data_file_sstt2(const std::string& filepath,
                   std::vector<int64_t>* macrotimes,
                         uint64_t n_events_to_skip,
                         uint64_t n_overflows_had,
                         uint64_t* n_overflows_in_file)
{
    if (macrotimes == nullptr) {
        return 2;
    }

    if (!test_is_sstt2_file(filepath)) {
        return 3;
    }

    FILE* f = fopen(filepath.c_str(), "rb");

    if (f == NULL) {
        // Error opening the file
        return 1;
    }

    // Read in the header
    std::vector<char> header(18,'\0');

    size_t n_read = fread(header.data(), SSTT2_N_BYTES_HEADER, 1, f);
    uint64_t n_overflows = 0;

    if (n_events_to_skip != 0) {
        int success = fseek(f, SSTT2_N_BYTES_TOT * (n_events_to_skip + n_overflows_had), SEEK_CUR);

        if (success != 0) {
            printf("ERROR: cound not skip required number (%lld photons, %lld overflows) of events! \n", n_events_to_skip, n_overflows_had);
            return 4;
        }

        n_overflows = n_overflows_had;
    }

    int64_t event = 0;

    n_read = fread(&event, SSTT2_N_BYTES_TOT, 1, f);

    while (n_read == 1) {
        struct processed_event_sstt2 pe = handle_event_sstt2(event, n_overflows);

        if (pe.n_overflows != 0) {
            // Update the number of overflows
            n_overflows += pe.n_overflows;
        } else {
            // Store the photon
            macrotimes->push_back(pe.macrotime);
        }

        n_read = fread(&event, SSTT2_N_BYTES_TOT, 1, f);
    }

    fclose(f);

    if (n_overflows_in_file != nullptr) {
        *n_overflows_in_file = n_overflows;
    }

    return 0;
}

int LIBTIMETAG_DLL test_is_sstt2_info_file(const std::string &filepath)
{
    FILE* f = fopen(filepath.c_str(), "r");

    if (f == nullptr) {
        return 0;
    }

    char* line = nullptr;
    size_t len = 0;

    getline(&line, &len, f);

    if (std::string(line) == std::string(SSTT2_MAGIC_INFO)) {
        fclose(f);
        return 1;
    }

    if (line != nullptr)
        free(line);

    fclose(f);
    return 0;

}
std::vector<channel_info_sstt2> LIBTIMETAG_DLL get_sstt2_info(const char* filename, int* error_code, exp_info_sstt2 *exp_info)
{
    std::vector<channel_info_sstt2> ret;

    if (error_code == nullptr) {
        return ret; // Because screw you. You best take notice of errors.
    }

    FILE* f = fopen(filename, "r");

    if (f == nullptr) {
        *error_code = 1;
        return ret;
    }

    long long read = 0;
    char* line = nullptr;
    size_t len = 0;

    int start_exp_header = 0;
    int index_timeunit = 0;
    int index_dev_type = 0;

    int start_chan_header = 0;

    int start_exp_data = 0;
    int start_chan_data = 0;

    int index_chan_id = -1;
    int index_filename = -1;
    int index_num_photons = -1;
    int index_sync_div = -1;
    int index_add_sync_div = -1;
    int index_is_pulsechan = -1;
    int index_has_pulsechan = -1;
    int index_corr_pulsechan = -1;

    long start_chan_data_seek_pos = 0;
    int64_t n_channels = 0;

    double time_unit_seconds = 0.0;

    std::string dev_type = "";

    // Find out how many channels there are,
    // how the header columns are distributed
    while ((read = getline(&line, &len, f)) != -1) {

        if (strcmp(line, SSTT2_EXP_HEADER_TEXT) == 0) {
            // Channel header starts
            start_exp_header = 1;
            continue;
        }

        if (start_exp_header) {
            start_exp_header = 0;
            char* substring = strtok(line, SSTT2_HEADER_DELIMITER);

            unsigned int index = 0;

            // Loop through all column titles. If we find a hit,
            // store the index
            while (substring != NULL) {
                if (strcmp(substring, SSTT2_HEADER_TIMEUNIT) == 0) {
                    index_timeunit = index;
                } else if (strcmp(substring, SSTT2_HEADER_DEV_TYPE) == 0) {
                    index_dev_type = index;
                }

                index++;
                substring = strtok (NULL, SSTT2_HEADER_DELIMITER);
            }

            start_exp_data = 1;
            continue;
        }

        if (start_exp_data) {
            start_exp_data = 0;

            char* substring = strtok(line, SSTT2_HEADER_DELIMITER);

            int index = 0;

            // Loop through the info for this channel
            while (substring != NULL) {
                if (index == index_timeunit) {
                    time_unit_seconds = atof(substring);
                } else if (index == index_dev_type) {
                    dev_type = std::string(substring);
                }

                index++;
                substring = strtok (NULL, SSTT2_HEADER_DELIMITER);
            }
            continue;
        }

        if (strcmp(line, SSTT2_CHAN_HEADER_TEXT) == 0) {
            // Channel header starts
            start_chan_header = 1;
            continue;
        }

        if (start_chan_header) {
            start_chan_header = 0;
            char* substring = strtok(line, SSTT2_HEADER_DELIMITER);

            unsigned int index = 0;

            // Loop through all column titles. If we find a hit,
            // store the index
            while (substring != NULL) {
                if (strcmp(substring, SSTT2_HEADER_CHANID) == 0) {
                    index_chan_id = index;
                } else if (strcmp(substring, SSTT2_HEADER_FILENAME) == 0) {
                    index_filename = index;
                } else if (strcmp(substring, SSTT2_HEADER_NUMPHOTONS) == 0) {
                    index_num_photons = index;
                } else if (strcmp(substring, SSTT2_HEADER_SYNCDIV) == 0) {
                    index_sync_div = index;
                } else if (strcmp(substring, SSTT2_HEADER_ADDI_SYNCDIV) == 0) {
                    index_add_sync_div = index;
                } else if (strcmp(substring, SSTT2_HEADER_IS_PULSES) == 0) {
                    index_is_pulsechan = index;
                } else if (strcmp(substring, SSTT2_HEADER_HAS_PULSES) == 0) {
                    index_has_pulsechan= index;
                } else if (strcmp(substring, SSTT2_HEADER_CORR_PULSECHAN) == 0) {
                    index_corr_pulsechan = index;
                }

                index++;
                substring = strtok (NULL, SSTT2_HEADER_DELIMITER);
            }

            start_chan_data = 1;

            start_chan_data_seek_pos = ftell(f);
            continue;
        }

        if (start_chan_data) {
            if (index_chan_id == -1 ||
                    index_filename == -1 ||
                    index_num_photons == -1) {
                // The channel data is malformed
                *error_code = 2;
                return ret;
            }

            // Count the number of channels
            if (strlen(line) > 1) {
                n_channels++;
            } else {
                // Empty line signals the end of the table
                break;
            }
        }
    }

    if (!start_chan_data) {
        *error_code = 3; // Could not read channel data

        return ret;
    }


    // Now read in the channel info
    fseek(f, start_chan_data_seek_pos, 0);

    int64_t chan_counter = 0;

    while ((read = getline(&line, &len, f)) != -1 && chan_counter < n_channels) {
        char* substring = strtok(line, SSTT2_HEADER_DELIMITER);

        int index = 0;

        channel_info_sstt2 ci;
        ci.filename = "";
        ci.n_photons = 0;
        ci.ID = 0;
        ci.additional_sync_divider = 1;
        ci.channel_has_microtime = false;
        ci.corresponding_pulses_channel = 0;
        ci.has_pulses_channel = false;
        ci.is_pulses_channel = false;
        ci.sync_divider = 1;

        // Loop through the info for this channel
        while (substring != NULL) {
            if (index == index_chan_id) {
                ci.ID = atoi(substring);
            } else if (index == index_filename) {
                ci.filename = std::string(substring);
            } else if (index == index_num_photons) {
                ci.n_photons = atoi(substring);
            } else if (index == index_is_pulsechan) {
                ci.is_pulses_channel = (bool)atoi(substring);
            } else if (index == index_sync_div) {
                ci.sync_divider = atoi(substring);
            } else if (index == index_add_sync_div) {
                ci.additional_sync_divider = atoi(substring);
            } else if (index == index_has_pulsechan) {
                ci.has_pulses_channel= (bool)atoi(substring);
            } else if (index == index_corr_pulsechan) {
                ci.corresponding_pulses_channel = atoi(substring);
            }

            index++;
            substring = strtok (NULL, SSTT2_HEADER_DELIMITER);
        }

        if (ci.filename.size() >= 2) {
            ci.filename = ci.filename.substr(1, ci.filename.size() - 2);
        }
        ret.push_back(ci);
        chan_counter++;
    }

    fclose(f);

    *error_code = 0;

    exp_info->time_unit_seconds = time_unit_seconds;
    exp_info->device_type = dev_type;

    return ret;
}
