/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

/**
 * \file    algos.h
 * \brief   Contains algorithms useful in Time-Correlated Single-Photon counting experiments
 * \author  Stijn Hinterding
*/

#ifndef ALGOS_H
#define ALGOS_H
#include <stdint.h>

#ifdef _WIN32
#ifdef BUILDING_LIBTIMETAG
#define LIBTIMETAG_DLL __declspec(dllexport)
#else
#define LIBTIMETAG_DLL __declspec(dllimport)
#endif
#else
#define LIBTIMETAG_DLL
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief   Finds the index where a new element should be inserted to maintain order
 *
 * Finds the index into a sorted array \p a such that, if the \p value were inserted before the index, the order of \p a would be preserved.
 * This function iterates through \p a sequentially, as such it is not efficient for very long arrays.
 *
 * \param   a         Array to search in
 * \param   value     The value to search the index for
 * \param   guess_i   Index to start the search at
 * \param   len_a     The number of elements in array a
 * \param   side      0: return found index, 1: return found index plus one
 * \returns The found index. If the value is smaller than any value in the array, zero is returned. If the value is larger than any value in the array, the length of the array is returned.
 * \note    This function is similar to numpy's searchsorted() function.
*/
uint64_t LIBTIMETAG_DLL seq_search(const int64_t *a, int64_t value, uint64_t guess_i, uint64_t len_a, int64_t side);

/**
 * \brief   Finds the index where a new element should be inserted to maintain order
 *
 * Finds the index into a sorted array \p a such that, if the \p value were inserted before the index, the order of \p a would be preserved.
 * This function guesses an initial index by linear interpolation, then uses the seq_search() function to find the exact value.
 * As such, this function is relatively efficient for large arrays, wherein the values are spaced more-or-less equidistantly.
 *
 * \param   a         Array to search in
 * \param   value     The value to search the index for
 * \param   len_a     The number of elements in array a
 * \param   side      0: return found index, 1: return found index plus one
 * \returns The found index. If the value is smaller than any value in the array, zero is returned. If the value is larger than any value in the array, the length of the array is returned.
 * \note    This function is similar to numpy's searchsorted() function.
*/
uint64_t LIBTIMETAG_DLL interp_seq_search(const int64_t *a, int64_t value, uint64_t len_a, int side);

/**
 * \brief   Correlates two arrays with each other
 *
 * Correlated the sorted array \p left_list with the sorted array \p right_list, at an interval determined by the \p bin_edges.
 * This function is optimised to work for situations wherein there are many occurences per bin (e.g. in a Fluorescence Correlation Spectroscopy curve).
 * For data sets in which there are only few occurrences per bin (e.g. fluorescence intensity decay curves), use the correlate_unit_bins() function.
 * To normalise the resulting histogram, use the normalize_correlation() function.
 *
 * \param   bin_edges       An array containing the edges of the bins
 * \param   n_bin_edges     The number of bin edges
 * \param   left_list       An array containing the first data set
 * \param   left_list_len   The number of data points in the first data set
 * \param   right_list      An array containing the second data set
 * \param   right_list_len  The number of data points in the second data set
 * \param   histogram_ret   The array to store the correlation data in. Each new value will be added to the corresponding existing element.
 * \param   histogram_ret_len   The number of histogram bins, should be one smaller than \p n_bin_edges
 * \returns On success: 0. Else: 1: NULL pointer supplied as input; 2: \p n_bin_edges <= 1; 3: \p histogram_ret_len != \p n_bin_edges - 1.
 * \note    The \p histogram_ret array does not need to consist of zeroes. This may be useful in cases where you need to sum multiple histograms.
*/
int LIBTIMETAG_DLL correlate_many_per_bin(const int64_t *bin_edges,
                    uint64_t n_bin_edges,
                    const int64_t *left_list,
                    uint64_t left_list_len,
                    const int64_t *right_list,
                    uint64_t right_list_len,
                    int64_t *histogram_ret,
                    uint64_t histogram_ret_len);

int LIBTIMETAG_DLL correlate_many_per_bin_double(const double *bin_edges,
                    uint64_t n_bin_edges,
                    const double *left_list,
                    uint64_t left_list_len,
                    const double *right_list,
                    uint64_t right_list_len,
                    int64_t *histogram_ret,
                    uint64_t histogram_ret_len);
/**
 * \brief   Correlates two arrays with each other
 *
 * Correlated the sorted array \p left_list with the sorted array \p right_list, at an interval determined by the \p bin_edges. The bins must have a size of unity.
 * This function is optimised to work for situations wherein there are few occurences per bin (e.g. in fluorescence intensity decay curves).
 * For data sets in which there are many occurrences per bin (e.g. Fluorescence Correlation Spectroscopy curves), use the correlate_many_per_bin() function.
 * To normalise the resulting histogram, use the normalize_correlation() function.
 *
 * \param   bin_edges       An array containing the edges of the bins
 * \param   n_bin_edges     The number of bin edges
 * \param   left_list       An array containing the first data set
 * \param   left_list_len   The number of data points in the first data set
 * \param   right_list      An array containing the second data set
 * \param   right_list_len  The number of data points in the second data set
 * \param   histogram_ret   The array to store the correlation data in. Each new value will be added to the corresponding existing element.
 * \param   histogram_ret_len   The number of histogram bins, should be one smaller than \p n_bin_edges
 * \returns On success: 0. Else: 1: NULL pointer supplied as input; 2: \p n_bin_edges <= 1; 3: \p histogram_ret_len != \p n_bin_edges - 1; 4: bins are not unity-sized.
 * * \note    The \p histogram_ret array does not need to consist of zeroes. This may be useful in cases where you need to sum multiple histograms.
*/
int LIBTIMETAG_DLL correlate_unit_bins(const int64_t *bin_edges,
                        uint64_t n_bin_edges,
                        const int64_t *left_list,
                        uint64_t left_list_len,
                        const int64_t *right_list,
                        uint64_t right_list_len,
                        int64_t *histogram_ret,
                        uint64_t histogram_ret_len);

/**
 * \brief   Finds the index of the bins corresponding to the supplied data values
 *
 * Bins the supplied data values into the supplied bins.
 * This function makes an initial guess using linear interpolation, as such, it is most efficient for linear bins (bins all having the same size).
 *
 * \param   bin_edges       An array containing the edges of the bins
 * \param   n_bin_edges     The number of bin edges
 * \param   data            An array containing the data values
 * \param   data_len        The number of data points in the data set
 * \param   histogram_ret   The array to store the correlation data in. Each new value will be added to the corresponding existing element.
 * \param   histogram_ret_len   The number of histogram bins, should be one smaller than \p n_bin_edges
 * \returns On success: 0. Else: 1: NULL pointer supplied as input; 2: \p n_bin_edges <= 1; 3: \p histogram_ret_len != \p n_bin_edges - 1; 4: bins are not unity-sized.
 * \note    The \p histogram_ret array does not need to consist of zeroes. This may be useful in cases where you need to sum multiple histograms.
*/
int LIBTIMETAG_DLL bindata_interp_seq(const int64_t *bin_edges,
                        uint64_t n_bin_edges,
                        const int64_t *data,
                        uint64_t data_len,
                        int64_t *histogram_ret,
                        uint64_t histogram_ret_len);

/**
 * \brief   Rebins a histogram according to a new bin size
 *
 * Takes already binned data and computes a new histogram, based on a new bin size, which is a multiple of the original bin size.
 * If not all original bins fit in the new histogram (i.e., if there are some bins left over, which together cannot form a new bin),
 * these leftover bins are discarded.
 * The number of new bins is calculated as: (\p binned_data_len - (\p binned_data_len % \p new_bin_size) ) / \p new_bin_size
 *
 * \param   binned_data     An array containing the existing histogram
 * \param   binned_data_len The number of bins in the existing histogram
 * \param   new_bin_size    The size of the bins in the new histogram, expressed in units of the number of original bins.
 * \param   ret_hist        The array to store the new histogram in
 * \param   ret_hist_len    The number of elements in (capacity of) the \p ret_hist array.
 * \returns On success: 0. Else: 1: NULL pointer supplied as input; 2: the value of \p ret_hist_len is not correct in combination with the supplied \p new_bin_size.
 * \note    The \p ret_hist array does not need to consist of zeroes. This may be useful in cases where you need to sum multiple histograms.
*/
int LIBTIMETAG_DLL rebin(const int64_t *binned_data,
           uint64_t binned_data_len,
           uint64_t new_bin_size,
           int64_t *ret_hist,
           uint64_t ret_hist_len);

uint64_t LIBTIMETAG_DLL rebin_len(uint64_t binned_data_len,
                                  uint64_t new_bin_size);

uint64_t LIBTIMETAG_DLL rebin_bin_edges_len(uint64_t n_org_bin_edges,
                                            uint64_t new_bin_size);
/**
 * \brief   Determines the bins corresponding to a rebinned histogram
 *
 * Takes the bin edges of an original histogram and computes new bin edges, based on a new bin size.
 * The number of new bin edges is calculated as: (\p n_org_bin_edges - 1 - ((\p n_org_bin_edges - 1) % \p new_bin_size) ) / (\p new_bin_size + 1)
 *
 * \param   org_bin_edges   An array containing the existing bin edges
 * \param   n_org_bin_edges The number of existing bin edges
 * \param   new_bin_size    The size of the bins in the new histogram, expressed in units of the number of original bins
 * \param   new_bin_edges   The array to store the new bin edges in
 * \param   n_new_bin_edges The number of elements in (capacity of) the \p ne_bin_edges array.
 * \returns On success: 0. Else: 1: NULL pointer supplied as input; 2: \p n_org_bin_edges <= 1; 3: the value of \p n_new_bin_edges is not correct in combination with the supplied \p new_bin_size.
*/
int LIBTIMETAG_DLL rebin_bin_edges(const int64_t *org_bin_edges,
                    uint64_t n_org_bin_edges,
                    uint64_t new_bin_size,
                    int64_t *new_bin_edges,
                    uint64_t n_new_bin_edges);

void LIBTIMETAG_DLL logspace(double start, double stop, uint64_t num, double base, double* ret);

int64_t LIBTIMETAG_DLL linspace_len(int64_t start,
                                     int64_t stop,
                                     int64_t step_size,
                                     int right_inclusive,
                                     int list_must_contain_stop);

int64_t LIBTIMETAG_DLL linspace(int64_t start,
                             int64_t stop,
                             int64_t step_size,
                             int right_inclusive,
                             int list_must_contain_stop,
                             int64_t* result,
                             int64_t result_len);

int LIBTIMETAG_DLL normalize_correlation(const int64_t *corr_hist,
                           uint64_t hist_len,
                           const int64_t *bin_edges,
                           uint64_t n_bin_edges,
                           uint64_t T_min,
                           uint64_t T_max,
                           uint64_t n_photons_left,
                           uint64_t n_photons_right,
                           double* ret);

int LIBTIMETAG_DLL normalize_correlation_double(const int64_t *corr_hist,
                           uint64_t hist_len,
                           const double *bin_edges,
                           uint64_t n_bin_edges,
                           double T_min,
                           double T_max,
                           uint64_t n_photons_left,
                           uint64_t n_photons_right,
                           double* ret);

int LIBTIMETAG_DLL gen_microtimes(const int64_t* pulses_macrotimes,
                                  uint64_t pulses_macrotimes_len,
                                  const int64_t* data_macrotimes,
                                  uint64_t data_macrotimes_len,
                                  int64_t* results_buffer,
                                  uint64_t results_buffer_len, uint64_t total_sync_divider);
#ifdef __cplusplus
}
#endif

#endif // ALGOS_H
