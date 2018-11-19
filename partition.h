/*
 * Implements variants of parallel partition and serial partition.
 */

#include <cstdint>
#include "params.h"

#ifndef SORT_H
#define SORT_H
void serial_partition(int64_t *array, uint64_t n, int64_t pivot);
void serial_partition2(int64_t *array, uint64_t n, int64_t pivot);
void partition(int64_t *array, uint64_t n, int64_t pivot);
void in_place_partition(int64_t *array, uint64_t n, int64_t pivot);
void small_prefix_partition(int64_t *array, uint64_t n, int64_t pivot);
void test_prefix_sum();
#endif
