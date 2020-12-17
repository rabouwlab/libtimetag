/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

/**
 * \file    algos.cpp
 * \brief   Contains algorithms useful in Time-Correlated Single-Photon counting experiments
 * \author  Stijn Hinterding
*/

#include "algos.h"

#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <vector>

template <typename T, typename U>
int _normalize_correlation(const int64_t* corr_hist,
                           uint64_t hist_len,
                           const T* bin_edges,
                           uint64_t n_bin_edges,
                           U T_min,
                           U T_max,
                           uint64_t n_photons_left,
                           uint64_t n_photons_right,
                           double* ret)
{
    if (hist_len != n_bin_edges - 1) {
        return 1;
    }

    double n_photons_squared = (double)(n_photons_left * n_photons_right);
    double mult = n_photons_squared / ((double)pow(T_max - T_min, 2.0));

    for (uint64_t i = 0; i < hist_len; i++) {
        double A = (double)(bin_edges[i + 1] - bin_edges[i]) * (T_max - T_min + 0.5 - 0.5*(bin_edges[i] + bin_edges[i + 1]));

        double divider = A * mult;
        double val = (divider == 0) ? 0.0 : (corr_hist[i] / divider);

        ret[i] = val;
    }

    return 0;
}

template <typename T>
uint64_t _seq_search_left(const T* a, T value, uint64_t guess_i, uint64_t len_a)
{
    if (value < a[0]) {
        return 0;
    }

    if (value > a[len_a - 1]) {
        return len_a;
    }

    if (guess_i >= len_a) {
        guess_i = len_a - 1;
    }

    if (a[guess_i] >= value || guess_i == len_a - 1) {
        for (int64_t j = guess_i; j >= 0; j--) {
            if (a[j] < value) {
                return j + 1;
            }
        }

        return 0;
    } else {
        for (uint64_t j = guess_i; j < len_a; j++) {
            if (a[j] >= value) {
                return j;
            }
        }

        return len_a - 1;
    }

    return len_a;
}

template <typename T>
uint64_t _seq_search(const T* a, T value, uint64_t guess_i, uint64_t len_a, int64_t side)
{
    if (value < a[0]) {
        return 0;
    }

    if (value > a[len_a - 1]) {
        return len_a;
    }

    if (guess_i >= len_a) {
        guess_i = len_a - 1;
    }

    if (a[guess_i] > value) {
        for (int64_t j = guess_i - 1; j >= 0; j--) {
            if (a[j] <= value) {
                return j + side;
            }
        }
    } else {
        for (uint64_t j = guess_i + 1; j < len_a; j++) {
            if (a[j] > value) {
                return j + side - 1;
            }
        }
    }

    return len_a;
}

template <typename T>
uint64_t _interp_seq_search_left(const T* a, T value, uint64_t len_a)
{
    double guess_rel = (double)(value - a[0])/(double)(a[len_a - 1]-a[0]);

    if (guess_rel < 0) {
        return 0;
    }

    if (guess_rel > 1) {
        return len_a;
    }

    uint64_t guess_i = (uint64_t)(guess_rel * (len_a - 1));

    return _seq_search_left(a, value, guess_i, len_a);
}

template <typename T>
uint64_t _interp_seq_search(const T* a, T value, uint64_t len_a, int side)
{
    double guess_rel = (double)(value - a[0])/(double)(a[len_a - 1]-a[0]);

    if (guess_rel < 0) {
        return 0;
    }

    if (guess_rel > 1) {
        return len_a;
    }

    uint64_t guess_i = (uint64_t)(guess_rel * (len_a - 1));

    return _seq_search(a, value, guess_i, len_a, side);
}

template <typename T>
int _correlate_many_per_bin(const T* bin_edges,
                    uint64_t n_bin_edges,
                    const T* left_list,
                    uint64_t left_list_len,
                    const T* right_list,
                    uint64_t right_list_len,
                    int64_t* histogram_ret,
                    uint64_t histogram_ret_len)
{
    if (bin_edges == nullptr || left_list == nullptr || right_list == nullptr || histogram_ret == nullptr)
        return 1; // Input is invalid

    if (n_bin_edges <= 1)   // We should have at least one bin
        return 2;

    if (histogram_ret_len != n_bin_edges - 1)   // The return histogram and the bin edges should match
        return 3;

    if (left_list_len == 0 || right_list_len == 0) // We are finished
        return 0;

    uint64_t* prev_indices = (uint64_t*)malloc(n_bin_edges * sizeof(uint64_t));
    memset(prev_indices, 0, n_bin_edges * sizeof(uint64_t));

    for (uint64_t i = 0; i < n_bin_edges; i++) {
        prev_indices[i] = _interp_seq_search_left(right_list, bin_edges[i] + left_list[0], right_list_len);
    }

    for (uint64_t i = 0; i < left_list_len; i++) {
        uint64_t prev_index = _seq_search_left(right_list, left_list[i] + bin_edges[0], prev_indices[0], right_list_len);

        prev_indices[0] = prev_index;

        for (uint64_t j = 1; j < n_bin_edges; j++) {
            uint64_t found_index = _seq_search_left(right_list, left_list[i] + bin_edges[j], prev_indices[j], right_list_len);

            prev_indices[j] = found_index;

            histogram_ret[j - 1] += found_index - prev_index;
            prev_index = found_index;
        }
    }

    free(prev_indices);
    return 0;
}

#ifdef __cplusplus
extern "C" {
#endif

void LIBTIMETAG_DLL logspace(double start, double stop, uint64_t num, double base, double* ret)
{
    double real_start = pow(base, start);
    double real_base = pow(base, (stop - start)/(double)num);

    double cur_value = real_start;

    for (uint64_t i = 0; i < num ; i++) {
        ret[i] = cur_value;
        cur_value *= real_base;
    }
}

int64_t LIBTIMETAG_DLL linspace_len(int64_t start,
                 int64_t stop,
                 int64_t step_size,
                 int right_inclusive,
                 int list_must_contain_stop)
{
    if (start > stop) {
        return -1;
    }

    if (start == stop && step_size != 1) {
        return -2;
    }

    if (step_size < 0) {
        return -3;
    }

    if (step_size == 0) {
        return 0;
    }

    int64_t n_entire_bins = 0;

    if (list_must_contain_stop) {
        right_inclusive = 1;
        // We need to add the leftovers, so that we can also have the stop value
        int64_t leftover = (stop - start) % step_size;

        if (leftover != 0 && leftover < step_size) {
            n_entire_bins += 1;
        } else {
            n_entire_bins += leftover;
        }
    }

    if (right_inclusive) {
        // This truncates the division because we are not using floating point numbers
        n_entire_bins += (int64_t)((stop - start) / step_size) + 1;
    } else {
        n_entire_bins += (int64_t)((stop - 1 - start) / step_size) + 1;
    }

    return n_entire_bins;
}

int64_t LIBTIMETAG_DLL linspace(int64_t start,
                 int64_t stop,
                 int64_t step_size,
                 int right_inclusive,
                 int list_must_contain_stop,
                 int64_t* result,
                 int64_t result_len)
{
    int64_t len = linspace_len(start, stop, step_size, right_inclusive, list_must_contain_stop);

    if (len <= 0) {
        return len;
    }

    if (result_len != len) {
        return -1337;
    }

    for (int64_t i = 0; i < result_len; i++) {
        result[i] = start + i * step_size;
    }

    return len;
}


uint64_t LIBTIMETAG_DLL seq_search(const int64_t* a, int64_t value, uint64_t guess_i, uint64_t len_a, int64_t side)
{
    return _seq_search(a, value, guess_i, len_a, side);
}

uint64_t LIBTIMETAG_DLL interp_seq_search(const int64_t* a, int64_t value, uint64_t len_a, int side)
{
    return _interp_seq_search(a, value, len_a, side);
}


int LIBTIMETAG_DLL correlate_many_per_bin(const int64_t* bin_edges,
                    uint64_t n_bin_edges,
                    const int64_t* left_list,
                    uint64_t left_list_len,
                    const int64_t* right_list,
                    uint64_t right_list_len,
                    int64_t* histogram_ret,
                    uint64_t histogram_ret_len)
{
    return _correlate_many_per_bin(bin_edges, n_bin_edges, left_list, left_list_len, right_list, right_list_len, histogram_ret, histogram_ret_len);
}

int LIBTIMETAG_DLL correlate_many_per_bin_double(const double* bin_edges,
                    uint64_t n_bin_edges,
                    const double* left_list,
                    uint64_t left_list_len,
                    const double* right_list,
                    uint64_t right_list_len,
                    int64_t* histogram_ret,
                    uint64_t histogram_ret_len)
{
    return _correlate_many_per_bin(bin_edges, n_bin_edges, left_list, left_list_len, right_list, right_list_len, histogram_ret, histogram_ret_len);
}

int LIBTIMETAG_DLL correlate_unit_bins(const int64_t* bin_edges,
                        uint64_t n_bin_edges,
                        const int64_t* left_list,
                        uint64_t left_list_len,
                        const int64_t* right_list,
                        uint64_t right_list_len,
                        int64_t* histogram_ret,
                        uint64_t histogram_ret_len)
{
    if (bin_edges == NULL || left_list == NULL || right_list == NULL || histogram_ret == NULL)
        return 1; // Input is invalid

    if (n_bin_edges <= 1)   // We should have at least one bin
        return 2;

    if (histogram_ret_len != n_bin_edges - 1)   // The return histogram and the bin edges should match
        return 3;

    if (bin_edges[1] - bin_edges[0] != 1) {   // Imperfect check to see if the input bins are OK
        return 4;
    }

    if (left_list_len == 0 || right_list_len == 0) // We are finished
        return 0;

    // Now do the correlation
    uint64_t next_photon_to_check = 0;

    // We loop through the left photon list.
    for (uint64_t i = 0; i < left_list_len; i++) {

        // We now take the current photon in the left photon list as
        // our 'origin', for building up the current correlation histogram.
        //
        // We loop through the right photon list. Note that we do not start
        // at the very first photon: we start at a point we previously
        // determined.
        for (uint64_t j = next_photon_to_check; j < right_list_len; j++) {

            // Here we determine the time difference between
            // the current photon (in the right list) and the reference
            // photon (in the left list). We use the very first bin edge
            // here as offset.
            int64_t index2 = right_list[j] - (left_list[i] + bin_edges[0]);

            // Now see what we need to do
            // if index2 < 0:
            //      This photon is not in our histogram (because it is too small)
            //      However, future photon times may be large enough, so we continue
            // if index2 >= n_bins:
            //      This photon is not in our histogram (because it is too big)
            //      Future photons will be too large, so we are completely done
            //      with the current histogram. We break.
            // else:
            //      The photon is within our histogram. We update the histogram.
            if (index2 < 0) {
                next_photon_to_check = j;
                continue;
            } else if ((uint64_t)index2 >= n_bin_edges) {
                break;
            } else {
                histogram_ret[index2]++;
            }
        }
    }

    return 0;
}

int LIBTIMETAG_DLL bindata_interp_seq(const int64_t* bin_edges,
                        uint64_t n_bin_edges,
                        const int64_t* data,
                        uint64_t data_len,
                        int64_t* histogram_ret,
                        uint64_t histogram_ret_len)
{
    if (bin_edges == NULL || data == NULL || histogram_ret == NULL)
        return 1;

    if (n_bin_edges <= 1)
        return 2;

    if (histogram_ret_len != n_bin_edges - 1)
        return 3;

    if (data_len == 0)
        return 0;

    int64_t leftmost_bin_edge = bin_edges[0];
    int64_t rightmost_bin_edge = bin_edges[n_bin_edges - 1];

    for (uint64_t i = 0; i < data_len; i++) {
        int64_t d = data[i];

        if (d < leftmost_bin_edge || d > rightmost_bin_edge)
            continue;

        uint64_t index = interp_seq_search(bin_edges, data[i], n_bin_edges, 0);

        if (index >= histogram_ret_len)
            continue;

        histogram_ret[index]++;
    }

    return 0;
}

uint64_t LIBTIMETAG_DLL rebin_bin_edges_len(uint64_t n_org_bin_edges, uint64_t new_bin_size)
{
    uint64_t remainder = (n_org_bin_edges - 1) % new_bin_size;
    return (n_org_bin_edges - 1 - remainder) / new_bin_size + 1;
}

int LIBTIMETAG_DLL rebin_bin_edges(const int64_t* org_bin_edges,
                    uint64_t n_org_bin_edges,
                    uint64_t new_bin_size,
                    int64_t* new_bin_edges,
                    uint64_t n_new_bin_edges)
{
    if (org_bin_edges == NULL || new_bin_edges == NULL)
        return 1;

    if (n_org_bin_edges <= 1)
        return 2;

    uint64_t ret_size = rebin_bin_edges_len(n_org_bin_edges, new_bin_size);

    if (n_new_bin_edges != ret_size)
        return 3;

    uint64_t counter = 0;

    for (uint64_t i = 0; i < n_org_bin_edges; i++) {
        if (i % new_bin_size == 0) {
            new_bin_edges[counter] = org_bin_edges[i];
            counter++;
        }
    }

    return 0;
}

uint64_t LIBTIMETAG_DLL rebin_len(uint64_t binned_data_len, uint64_t new_bin_size)
{
    uint64_t remainder = binned_data_len % new_bin_size;
    return (binned_data_len - remainder) / new_bin_size;
}

int LIBTIMETAG_DLL rebin( const int64_t* binned_data,
           uint64_t binned_data_len,
           uint64_t new_bin_size,
           int64_t* ret_hist,
           uint64_t ret_hist_len)
{
    if (binned_data == NULL || ret_hist == NULL)
        return 1;

    uint64_t ret_size = rebin_len(binned_data_len, new_bin_size);

    if (ret_hist_len != ret_size)
        return 2;

    if (binned_data_len == 0)
        return 0;

    if (new_bin_size == 1) {
        memcpy(ret_hist, binned_data, sizeof(uint64_t) * ret_hist_len);
        return 0;
    }

    int64_t temp_val = 0;
    uint64_t counter = 0;
    for (uint64_t i = 0; i < binned_data_len; i++) {
        temp_val += binned_data[i];

        if ((i+1) % new_bin_size == 0) {
            ret_hist[counter] += temp_val;
            temp_val = 0;
            counter++;
        }
    }

    return 0;
}

int LIBTIMETAG_DLL normalize_correlation(const int64_t* corr_hist,
                           uint64_t hist_len,
                           const int64_t* bin_edges,
                           uint64_t n_bin_edges,
                           uint64_t T_min,
                           uint64_t T_max,
                           uint64_t n_photons_left,
                           uint64_t n_photons_right,
                           double* ret)
{
    return _normalize_correlation(corr_hist, hist_len, bin_edges, n_bin_edges, T_min, T_max, n_photons_left, n_photons_right, ret);
}

int LIBTIMETAG_DLL normalize_correlation_double(const int64_t* corr_hist,
                           uint64_t hist_len,
                           const double* bin_edges,
                           uint64_t n_bin_edges,
                           double T_min,
                           double T_max,
                           uint64_t n_photons_left,
                           uint64_t n_photons_right,
                           double* ret)
{
    return _normalize_correlation(corr_hist, hist_len, bin_edges, n_bin_edges, T_min, T_max, n_photons_left, n_photons_right, ret);
}

int LIBTIMETAG_DLL gen_microtimes(const int64_t *pulses_macrotimes,
                   uint64_t pulses_macrotimes_len,
                   const int64_t *data_macrotimes,
                   uint64_t data_macrotimes_len,
                   int64_t *results_buffer,
                   uint64_t results_buffer_len,
                   uint64_t total_sync_divider)
{
    if (pulses_macrotimes_len == 0 ||
            data_macrotimes_len == 0 ||
            pulses_macrotimes == nullptr ||
            data_macrotimes == nullptr ||
            results_buffer == nullptr ||
            data_macrotimes_len != results_buffer_len) {
        return 1; // Input invalid.
    }

    // See if we need to generate extra pulse times
    std::vector<int64_t> extra_pulses;

    double avg_pulse_duration = (double)(pulses_macrotimes[pulses_macrotimes_len - 1] - pulses_macrotimes[0])/((double)pulses_macrotimes_len - 1);
    int64_t pulse_duration = (int64_t)round(avg_pulse_duration);

    int64_t latest_gen = pulses_macrotimes[0];

    while (latest_gen > data_macrotimes[0]) {
        latest_gen -= pulse_duration;
        extra_pulses.push_back(latest_gen);
    }

    latest_gen = pulses_macrotimes[pulses_macrotimes_len-1];

    while (latest_gen <= data_macrotimes[data_macrotimes_len - 1]) {
        latest_gen += pulse_duration;
        extra_pulses.push_back(latest_gen);
    }

    // Now copy the existing pulses to the extra pulses vector
    uint64_t cur_vector_len = extra_pulses.size();
    extra_pulses.resize(cur_vector_len + pulses_macrotimes_len);

    memcpy(&extra_pulses[cur_vector_len], pulses_macrotimes, pulses_macrotimes_len * sizeof(int64_t));

    // Sort the pulses vector
    std::sort(extra_pulses.begin(), extra_pulses.end());

    int64_t prev_found_pulse_index = 0;

    // Now actually generate the microtimes
    for (uint64_t i = 0; i < data_macrotimes_len; i++) {
        int64_t macro_t = data_macrotimes[i];

        // Now find the pulse we are most closely related to
        auto found_index = std::upper_bound(extra_pulses.begin() + prev_found_pulse_index, extra_pulses.end(), macro_t);
        prev_found_pulse_index =  found_index - extra_pulses.begin() - 1;

        // See if we found anything
        if (found_index == extra_pulses.end()) {
            // We failed to find anything.
            // This means we will not find anything in the future... Quit.
            // (this should not happen)
            return 2;
        }

        // We found something :-)
        int64_t found_pulse_t = *(found_index-1);
        int64_t dt = macro_t - found_pulse_t;

        // TODO: maybe do not use the average pulse duration here?
        double div = avg_pulse_duration / (double)total_sync_divider;
        double rem = fmod((double)dt, div) ;

        results_buffer[i] = (int64_t) rem;
    }

    return 0;
}

#ifdef __cplusplus
}
#endif
