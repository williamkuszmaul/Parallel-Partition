#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <cstdlib>
#include <math.h>
#include <assert.h>
#include <unordered_set>
#include <cstdint>
#include <bits/stdc++.h> 
#include <sys/time.h>
#include <cilk/cilk.h>
#include <string>
#include <cilk/cilk_api.h>

using namespace std;


struct timeval tp;

#define BLOCK_SIZE 16


#define NUM_THREADS 2
//#define USE_OPENMP
#define USE_CILK
//#define RUN_SERIAL
#ifdef USE_OPENMP
#define parallel_for for
#define USE_OPENMP_FOR omp parallel for num_threads(NUM_THREADS)
#endif
#ifdef USE_CILK
#define parallel_for cilk_for
#define USE_OPENMP_FOR 
#endif
#ifdef RUN_SERIAL
#define parallel_for for
#define USE_OPENMP_FOR 
#endif


 // often when an algorithm that is not naturally in-place is made in
 // place it pays significantly in practical-performance. The
 // simplicity of our algorithmic techniques instead allows for us to
 // build an (almost) in-place partition that achieves modest
 // \emph{improvements} in runtime, running a factor of roughly a
 // third faster in serial than an optimized non-in-place
 // implementation with the same span, and running only a factor of
 // two slower than the textbook serial algorithm. In fact, it closes
 // more than half the gap. (Note span is log n for our implementations)


 // When discussing why in-place algorithmsa re important: For
 // algorithms that are unable to fit their entire input in RAM,
 // portions of the problem are often handled in RAM. The size of the
 // portion that can be handled is determined by the memory footprint
 // of the algorithm. Also in systems where many processes or threads
 // may be running in parallel, smaller memory footprint reduces the
 // likelyhood of there being a moment in which the system runs out of
 // memory and is forced to page and improves cache performance. 

/*
 * Replaces array with a prefix sum of itself. Implements parallel
 * prefix sum algorithm but is currently written in serial. Uses extra
 * memory to be cache friendly.
 */
void prefix_sum_unoptimized_cache_friendly(int64_t *array, uint64_t n) {
  assert((n & (n - 1)) == 0); // Assert n is a pwer of two
  if (n <= 1) return;
  int64_t *small_array = (int64_t*)malloc(sizeof(int64_t) * n);
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n / 2; i++) {
    small_array[i] = array[2 * i] + array[2 * i + 1];
  }
  prefix_sum_unoptimized_cache_friendly(small_array, n / 2);
  array[1] = array[1] + array[0];
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 1; i < n / 2; i++) {
    array[2 * i] = small_array[i - 1] + array[2 * i];
    array[2 * i + 1] = small_array[i];
  }
  free(small_array);
}

/*
 * First reduces the size of the prefix sum by BLOCK_SIZE, and then
 * outsources it to prefix_sum_unoptimized2.
 */
void prefix_sum_optimized_cache_friendly(int64_t *array, uint64_t n) {
  int64_t *small_array = (int64_t*)malloc(sizeof(int64_t) * n / BLOCK_SIZE);
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n / BLOCK_SIZE; i++) {
    int64_t sum = 0;
    for (int64_t j = i * BLOCK_SIZE; j < i * BLOCK_SIZE + BLOCK_SIZE; j++) {
      sum += array[j];
    }
    small_array[i] = sum;
  }
  prefix_sum_unoptimized_cache_friendly(small_array, n / BLOCK_SIZE);
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n / BLOCK_SIZE; i++) {
    int64_t prefix_sum = 0;
    if (i > 0) prefix_sum = small_array[i - 1];
    for (int64_t j = i * BLOCK_SIZE; j < i * BLOCK_SIZE + BLOCK_SIZE; j++) {
      prefix_sum += array[j];
      array[j] = prefix_sum;
    }
  }
}


/*
 * Replaces array with a prefix sum of itself. Implements parallel
 * prefix sum algorithm but is currently written in serial.
 */
void prefix_sum_unoptimized(int64_t *array, uint64_t n) {
  prefix_sum_unoptimized_cache_friendly(array, n);
  return;
  
  // The above code also works and does the prefix sum in-place. The
  // function call above runs slightly faster, however.
  
  assert((n & (n - 1)) == 0); // Assert n is a pwer of two
  int64_t jump;
  for (jump = 2; jump <= n; jump *= 2) {
    #pragma USE_OPEN_MP_FOR
    parallel_for (int64_t i = n; i >= jump; i -= jump) {
      array[i - 1] += array[i - jump / 2 - 1];
    }
  }
  for (; jump >= 2; jump /= 2) {
    #pragma USE_OPEN_MP_FOR
    parallel_for (int64_t i = jump / 2 + jump; i <= n; i+= jump) {
      array[i - 1] = array[i - jump / 2 - 1] + array[i - 1];
    }
  }
}

/*
 * Replaces array with a prefix sum of itself. Implements parallel
 * prefix sum algorithm but is currently written in serial.
 */
void prefix_sum_optimized(int64_t *array, uint64_t n) {
  prefix_sum_optimized_cache_friendly(array, n);
  return;

  // The above code also works and does the prefix sum in-place. The
  // function call above runs slightly faster, however.

  assert((n & (n - 1)) == 0); // Assert n is a pwer of two
  if (n < BLOCK_SIZE) {
    prefix_sum_unoptimized(array, n);
    return;
  }
  int64_t jump;
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n / BLOCK_SIZE; i++) {
    int64_t sum = 0;
    for (int64_t j = i * BLOCK_SIZE; j < i * BLOCK_SIZE + BLOCK_SIZE; j++) {
      sum += array[j];
      array[j] = sum;
    }
  }
  for (jump = BLOCK_SIZE * 2; jump <= n; jump *= 2) {
    #pragma USE_OPEN_MP_FOR
    parallel_for (int64_t i = n; i >= jump; i -= jump) {
      array[i - 1] += array[i - jump / 2 - 1];
    }
  }
  for (; jump >= BLOCK_SIZE * 2; jump /= 2) {
    #pragma USE_OPEN_MP_FOR
    parallel_for (int64_t i = jump / 2 + jump; i <= n; i+= jump) {
      array[i - 1] = array[i - jump / 2 - 1] + array[i - 1];
    }
  }
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 1; i < n / BLOCK_SIZE; i++) {
    for (int64_t j = i * BLOCK_SIZE; j < i * BLOCK_SIZE + BLOCK_SIZE - 1; j++) {
      array[j] += array[i * BLOCK_SIZE - 1];
    }
  }
}



/*
 * A simple test for prefix_sum.
 */
void test_prefix_sum() {
  int64_t array[16] = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
  prefix_sum_optimized(array, 16);
  for (int64_t i = 0; i < 16; i++) assert(array[i] = i / 2 + 1);

  int64_t array2[16] = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
  prefix_sum_unoptimized(array2, 16);
  for (int64_t i = 0; i < 16; i++) assert(array2[i] = i / 2 + 1);

}

/*
 * Partitions array so that elements <= pivot appear before other
 * elts.
 */
void partition(int64_t *array, uint64_t n, int64_t pivot) {
  int64_t *prefixes = (int64_t*) malloc(sizeof(int64_t) * n);
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n; i++) {
    prefixes[i] = (array[i] <= pivot);
  }
  prefix_sum_optimized(prefixes, n);
  int64_t *output_tmp = (int64_t*) malloc(sizeof(int64_t) * n);
  int64_t num_preds = prefixes[n - 1];
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n; i++) {
    if (array[i] <= pivot) {
      output_tmp[prefixes[i] - 1] = array[i];
    } else {
      output_tmp[num_preds + (i - prefixes[i])] = array[i]; 
    }
  }
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n; i++) {
    array[i] = output_tmp[i];
  }
  free(output_tmp);
  free(prefixes);
}

/*
 * Performs the reordering step of the partitioning algorithm
 * in-place; requires the preprocessing step that ensures all prefixes
 * are successor-heavy.
 */
void in_place_reorder_with_blocks(int64_t* array, int64_t* prefixes, uint64_t n, int64_t pivot) {
  assert ((n & (BLOCK_SIZE - 1)) == 0);
  if (n <= BLOCK_SIZE * 10) {
    // We just take the easy way out here
    sort(array, array + n);
    return;
  }
  int64_t m = (n * 3 / 4 + BLOCK_SIZE) & (~(BLOCK_SIZE - 1));
  in_place_reorder_with_blocks(array, prefixes, m, pivot);

  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t j = m / BLOCK_SIZE; j < n / BLOCK_SIZE; j++) {
    int64_t new_pos = prefixes[j - 1];
    for (int64_t i = j * BLOCK_SIZE; i < j * BLOCK_SIZE + BLOCK_SIZE; i++) {
      if (array[i] <= pivot) {
        int64_t tmp = array[i];
        array[i] = array[new_pos];
        array[new_pos] = tmp;
        new_pos++; 
      }
    }
  }
}

/*
 * Performs the reordering step of the partitioning algorithm
 * in-place; requires the preprocessing step that ensures all prefixes
 * are successor-heavy.
 */
void in_place_reorder(int64_t* array, int64_t* prefixes, uint64_t n, int64_t pivot) {
  if (n <= 10) {
    // We just take the easy way out here
    sort(array, array + n);
    return;
  }
  int64_t m = n * 3 / 4 + 1;
  in_place_reorder(array, prefixes, m, pivot);
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = m; i < n; i++) {
    int64_t new_pos = prefixes[i] - 1;
    if (array[i] <= pivot) {
      int64_t tmp = array[i];
      array[i] = array[new_pos];
      array[new_pos] = tmp;
    }
  }
}


/*
 * Partitions array so that elements <= pivot appear before other
 * elts. Uses algorithmic techniques to reduce memory footprint
 * for prefix sums
 */
void small_prefix_partition(int64_t *array, uint64_t n, int64_t pivot) {
  int64_t *prefixes = (int64_t*) malloc(sizeof(int64_t) * n / BLOCK_SIZE);
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n / BLOCK_SIZE; i++) {
    int64_t sum = 0;
    for (int64_t j = i * BLOCK_SIZE; j < i * BLOCK_SIZE + BLOCK_SIZE; j++) {
      sum += (array[j] <= pivot);
    }
    prefixes[i] = sum;
  }
  // By using the unoptimized prefix sum, we make it so the advantage
  // we're getting really is from our use of less memory.
  prefix_sum_unoptimized(prefixes, n / BLOCK_SIZE);
  int64_t *output_tmp = (int64_t*) malloc(sizeof(int64_t) * n);
  int64_t num_preds = prefixes[n / BLOCK_SIZE - 1];
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n / BLOCK_SIZE; i++) {
    int64_t prefix_sum = 0;
    if (i > 0) prefix_sum = prefixes[i - 1];
    for (int64_t j = i * BLOCK_SIZE; j < i * BLOCK_SIZE + BLOCK_SIZE; j++) {
      if (array[j] <= pivot) {
        prefix_sum++;
        output_tmp[prefix_sum - 1] = array[j];
      } else {
        output_tmp[num_preds + (j - prefix_sum)] = array[j]; 
      }
    }
  }
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n; i++) {
    array[i] = output_tmp[i];
  }
  free(output_tmp);
  free(prefixes);
}

/*
 * Partitions array so that elements <= pivot appear before other
 * elts. Uses algorithmic techniques to reduce memory footprint
 * for the reordering step. For simplicity, asserts case
 * that number of predecessors is at least as large as number of
 * successors.
 */
void small_reorder_partition(int64_t *array, uint64_t n, int64_t pivot) {
  int64_t *prefixes = (int64_t*) malloc(sizeof(int64_t) * n);
  int64_t num_pred = 0;
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n; i++) {
    num_pred += (array[i] <= pivot);
  }
  assert(num_pred <= n / 2);
  for (int64_t m = n; m > 1; m /= 2) {
    #pragma USE_OPEN_MP_FOR
    parallel_for (int64_t i = 0; i < m / 2; i++) {
      if (array[i] <= pivot) {
        int64_t tmp = array[i];
        array[i] = array[m / 2 + i];
        array[m / 2 + i] = tmp;
      }
    }
  }
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n; i++) {
    prefixes[i] = (array[i] <= pivot);
  }
  prefix_sum_optimized(prefixes, n);
  in_place_reorder(array, prefixes, n, pivot);
  free(prefixes);
}



/*
 * Partitions array so that elements <= pivot appear before other
 * elts. Uses algorithmic techniques to reduce memory footprint
 * (though not all the way in place). For simplicity, asserts case
 * that number of predecessors is at least as large as number of
 * successors.
 */
void in_place_partition(int64_t *array, uint64_t n, int64_t pivot) {
  int64_t *prefixes = (int64_t*) malloc(sizeof(int64_t) * n / BLOCK_SIZE);
  int64_t num_pred = 0;
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n; i++) {
    num_pred += (array[i] <= pivot);
  }
  assert(num_pred <= n / 2);
  for (int64_t m = n; m > 1; m /= 2) {
    #pragma USE_OPEN_MP_FOR
    parallel_for (int64_t i = 0; i < m / 2; i++) {
      if (array[i] <= pivot) {
          int64_t tmp = array[i];
          array[i] = array[m / 2 + i];
          array[m / 2 + i] = tmp;
      }
    }
  }
  #pragma USE_OPEN_MP_FOR
  parallel_for (int64_t i = 0; i < n / BLOCK_SIZE; i++) {
    int64_t sum = 0;
    for (int64_t j = i * BLOCK_SIZE; j < i * BLOCK_SIZE + BLOCK_SIZE; j++) {
      sum += (array[j] <= pivot);
    }
    prefixes[i] = sum;
  }
  // By using the unoptimized prefix sum, we make it so the advantage
  // we're getting really is from our use of less memory.
  prefix_sum_unoptimized(prefixes, n / BLOCK_SIZE);
  in_place_reorder_with_blocks(array, prefixes, n, pivot);
  free(prefixes);
}

void serial_partition(int64_t *array, uint64_t n, int64_t pivot) {
  int64_t swap_index = n - 1;
  while (array[swap_index] > pivot) swap_index--;
  for (int64_t i = 0; i <= swap_index; i++) {
    if (array[i] > pivot) {
      int64_t tmp = array[i];
      array[i] = array[swap_index];
      array[swap_index] = tmp;
      while (array[swap_index] > pivot) swap_index--;
    }
  }
}


/*
 * A simple test for partition and in_place_partition.
 */
void test_partition(int64_t in_place, uint64_t n) {
  int64_t *array = (int64_t*) malloc(sizeof(int64_t) * n);
  uint64_t num_zeros;
  do {
    num_zeros = 0;
    for (int64_t i = 0; i < n; i++) {
      array[i] = rand() % 2;
      if (array[i] == 0) num_zeros++;
    }
  } while (num_zeros > n / 2); // For simplicity, we make number of predecessors at most n / 2
  if (in_place == -1) {
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    serial_partition(array, n, 0);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    cout<<(ms2 - ms1)<<" milliseconds for serial implementation test"<<endl;
  }
  if (in_place == 0) {
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    in_place_partition(array, n, 0);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    cout<<(ms2 - ms1)<<" milliseconds for in-place test"<<endl;
  }
  if (in_place == 1) {
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    small_prefix_partition(array, n, 0);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    cout<<(ms2 - ms1)<<" milliseconds for partially in-place prefix sum test"<<endl;
  }
  if (in_place == 2) {
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    small_reorder_partition(array, n, 0);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    cout<<(ms2 - ms1)<<" milliseconds for partially in-place reorder test"<<endl;
  }
  if (in_place == 3) {
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    partition(array, n, 0);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    cout<<(ms2 - ms1)<<" milliseconds for non-in-place test"<<endl;
  }
  for (int64_t i = 0; i < num_zeros; i++) {
    assert(array[i] == 0);
  }
  for (int64_t i = num_zeros; i < n; i++) {
    assert(array[i] == 1);
  }
  free(array);
}

int main() {
  string cilk_var = "CILK_NWORKERS=" + to_string(NUM_THREADS);
  assert(putenv((char*)cilk_var.c_str()) == 0);
  assert(__cilkrts_get_nworkers() == NUM_THREADS);
  srand (time(NULL)); // initialize random seed
  test_prefix_sum();
  for (uint64_t n = 64; n < 80000000; n *= 2) {
    cout<<"-----------"<<endl;
    cout<<n<<endl;
    test_partition(-1, n);
    test_partition(0, n);
    test_partition(1, n);
    test_partition(2, n);
    test_partition(3, n);
  }
  return 0;
}
