/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

/**
 * \file    sstt_file.cpp
 * \brief   Defines and implements the "small simple time-tagged" (SSTT) file format, version 1 (deprecated!).
 * \author  Stijn Hinterding
*/

#include "sstt_file.h"
#include "getline.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#define SSTT_CHAN_HEADER_TEXT "CHANNEL_HEADER\n"

#define SSTT_HEADER_DELIMITER       "\t"
#define SSTT_HEADER_CHANID          "ChannelID"
#define SSTT_HEADER_FILENAME        "Filename"
#define SSTT_HEADER_NUMPHOTONS      "NumPhotons"

struct processed_event
{
    int64_t macrotime;
    int64_t microtime;
    uint64_t n_overflows;
};

struct processed_event handle_event(int64_t e, uint64_t num_overflows)
{
    // Prepare the return value
    struct processed_event ret;
    ret.macrotime = 0;
    ret.microtime = 0;
    ret.n_overflows = 0;

    if ( (e & 1) && !( (e >> 1) & 1 ) ) {
        // Overflow event
        ret.n_overflows = (e >> SSTT_N_BITS_SIGNAL) & SSTT_MASK_OVERFLOW;
    } else if ( !(e & 1) && !( (e >> 1) & 1)) {
        // Photon event
        ret.microtime = (e >> SSTT_N_BITS_SIGNAL) & SSTT_MASK_MICRO;
        ret.macrotime = (e >> (SSTT_N_BITS_SIGNAL + SSTT_N_BITS_MICRO)) & SSTT_MASK_MACRO;

        // Also account for the overflows
        ret.macrotime += num_overflows * SSTT_OVERFLOW_VAL;
    }

    return ret;
}

int LIBTIMETAG_DLL n_photons_in_datafile(const char* directory,
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

    uint64_t n_overflows = 0;
    uint64_t event = 0;
    uint64_t n_read = fread(&event, 8, 1, f);

    while (n_read == 1) {
        struct processed_event pe = handle_event(event, n_overflows);

        if (pe.n_overflows == 0)
            (*ret)++;

        n_read = fread(&event, 8, 1, f);
    }

    fclose(f);

    return 0;
}

int LIBTIMETAG_DLL read_data_file(const std::string& filepath,
                   std::vector<int64_t>* macrotimes,
                   std::vector<int64_t>* microtimes)
{
    bool save_macro = false;
    bool save_micro = false;

    // Check if input is valid
    if (macrotimes != NULL) {
        save_macro = true;
    }

    if (microtimes != NULL) {
        save_micro = true;
    }

    if (!save_macro && !save_micro) {
        return 0;
    }

    FILE* f = fopen(filepath.c_str(), "rb");

    if (f == NULL) {
        // Error opening the file
        return 1;
    }

    int64_t event = 0;

    size_t n_read = fread(&event, 8, 1, f);
    uint64_t n_overflows = 0;

    while (n_read == 1) {
        struct processed_event pe = handle_event(event, n_overflows);

        if (pe.n_overflows != 0) {
            // Update the number of overflows
            n_overflows += pe.n_overflows;
        } else {
            // Store the photon
            macrotimes->push_back(pe.macrotime);
            microtimes->push_back(pe.microtime);
        }

        n_read = fread(&event, 8, 1, f);
    }

    fclose(f);

    return 0;
}

std::vector<channel_info> LIBTIMETAG_DLL get_sstt_info(const char* filename, int* error_code)
{
    std::vector<channel_info> ret;

    if (error_code == NULL) {
        return ret; // Because screw you. You best take notice of errors.
    }

    FILE* f = fopen(filename, "r");

    if (f == NULL) {
        *error_code = 1;
        return ret;
    }

    long long read = 0;
    char* line = NULL;
    size_t len = 0;

    int start_chan_header = 0;
    int start_chan_data = 0;

    int index_chan_id = -1;
    int index_filename = -1;
    int index_num_photons = -1;
    long start_chan_data_seek_pos = 0;
    int64_t n_channels = 0;

    // Find out how many channels there are,
    // how the header columns are distributed
    while ((read = getline(&line, &len, f)) != -1) {

        if (strcmp(line, SSTT_CHAN_HEADER_TEXT) == 0) {
            // Channel header starts
            start_chan_header = 1;
            continue;
        }

        if (start_chan_header) {
            start_chan_header = 0;
            char* substring = strtok(line, SSTT_HEADER_DELIMITER);

            unsigned int index = 0;

            // Loop through all column titles. If we find a hit,
            // store the index
            while (substring != NULL) {
                if (strcmp(substring, SSTT_HEADER_CHANID) == 0) {
                    index_chan_id = index;
                } else if (strcmp(substring, SSTT_HEADER_FILENAME) == 0) {
                    index_filename = index;
                } else if (strcmp(substring, SSTT_HEADER_NUMPHOTONS) == 0) {
                    index_num_photons = index;
                }
                index++;
                substring = strtok (NULL, SSTT_HEADER_DELIMITER);
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
        char* substring = strtok(line, SSTT_HEADER_DELIMITER);

        int index = 0;

        channel_info ci;
        ci.filename = "";
        ci.n_photons = 0;
        ci.ID = 0;

        // Loop through the info for this channel
        while (substring != NULL) {
            if (index == index_chan_id) {
                ci.ID = atoi(substring);
            } else if (index == index_filename) {
                ci.filename = std::string(substring);
            } else if (index == index_num_photons) {
                ci.n_photons = atoi(substring);
            }

            index++;
            substring = strtok (NULL, SSTT_HEADER_DELIMITER);
        }

        if (ci.filename.size() >= 2) {
            ci.filename = ci.filename.substr(1, ci.filename.size() - 2);
        }
        ret.push_back(ci);
        chan_counter++;
    }

    fclose(f);

    *error_code = 0;
    return ret;
}
