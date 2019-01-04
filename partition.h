/*
 * Implements variants of parallel partition and serial partition.
 */

#include <cstdint>
#include "params.h"

#ifndef PARTITION_H
#define PARTITION_H
int64_t serial_partition(int64_t *array, uint64_t n, int64_t pivot);
void partition(int64_t *array, uint64_t n, int64_t pivot);
int64_t in_place_partition(int64_t *array, uint64_t n, int64_t pivot);
void small_prefix_partition(int64_t *array, uint64_t n, int64_t pivot);
void test_prefix_sum();
void parallel_quicksort (int64_t* array, uint64_t num_elts);
void parallel_quicksort_high_span (int64_t* array, uint64_t num_elts);
void high_span_partition(int64_t* array, uint64_t n, int64_t pivot);
#endif
