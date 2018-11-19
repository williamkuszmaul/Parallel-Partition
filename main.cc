/*
 * Runs timing (and correctness) tests on the partition functions in partition.cc/h
 */

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
#include <cilk/cilk_api.h>
#include "omp.h"
#include <string>
#include "partition.h"

using namespace std;

struct timeval tp;

void set_num_threads(int num_threads) {
  __cilkrts_end_cilk();
  string str1 = "nworkers";
  string str2 = to_string(num_threads);
  assert(__cilkrts_set_param((char*)str1.c_str(), (char*)str2.c_str()) == 0);
  assert(__cilkrts_get_nworkers() == num_threads);
}


/*
 * A simple test for variants of partition. The argument in_place
 * determines which version of the algorithm we run. Returns number of
 * milliseconds that test took. Depending on the value of more_succ,
 * we either guarantee that the input has more predecessors than
 * successors or more successors than predecessors. For the in-place
 * algorithm this results in two different cases; our implementation
 * for the case where more_succ == false is slightly slower, so we
 * use it for all the experiments.
 *
 */
int64_t test_partition(int64_t in_place, uint64_t n, bool more_succ) {
  int64_t answer;
  int64_t *array = (int64_t*) malloc(sizeof(int64_t) * n);
  uint64_t num_zeros;
  do {
    num_zeros = 0;
    // We fill the array with zeros and 100s.
    for (int64_t i = 0; i < n; i++) {
      array[i] = rand() % 2 * 100;
      if (array[i] == 0) num_zeros++;
    }
  } while ((((num_zeros > n / 2) && more_succ) || ((num_zeros < n / 2) && !more_succ)));
  if (in_place == -2) { 
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    serial_partition2(array, n, 50);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    answer = (ms2 - ms1);
  }
  if (in_place == -1) { 
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    serial_partition(array, n, 50);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    answer = (ms2 - ms1);
  }
  if (in_place == 0) {
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    in_place_partition(array, n, 50); // We use pivot 50
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    answer = (ms2 - ms1);
  }
  if (in_place == 1) {
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    small_prefix_partition(array, n, 50);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    answer = ms2 - ms1;
  }
  if (in_place == 2) {
    gettimeofday(&tp, NULL);
    long int ms1 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    partition(array, n, 50);
    gettimeofday(&tp, NULL);
    long int ms2 = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    answer = ms2 - ms1;
  }
  // We verify correctness of output:
  for (int64_t i = 0; i < num_zeros; i++) {
    assert(array[i] == 0);
  }
  for (int64_t i = num_zeros; i < n; i++) {
    assert(array[i] == 100);
  }
  free(array);
  return answer;
}


/*
 * Averages answer of test_partition over multiple trials.
 */
double test_partition_multiple_trials(int64_t in_place, uint64_t n, uint64_t num_trials, bool more_succ) {
  int64_t answer = 0;
  for (int i = 0; i < num_trials; i++) {
    answer += test_partition(in_place, n, more_succ);
  }
  if (answer == 0) answer++;
  return (double)answer / (double)num_trials;
}

/* 
 *  Builds the table CILKtable in the paper.
 */ 
void parallel_tests(string tablename) {
  double serial_baseline = test_partition_multiple_trials(-1, MAX_INPUT_SIZE, NUM_TRIALS, false);
  cout<<"\\def \\"<<tablename<<"serialbaseline {"<<serial_baseline<<"}"<<endl;
  cout<<"\\def \\"<<tablename<<"blocksize {"<<BLOCK_SIZE<<"}"<<endl;
  cout<<"\\def \\"<<tablename<<"numtrials {"<<NUM_TRIALS<<"}"<<endl;
  cout<<"\\def \\"<<tablename<<"inputsize {"<<MAX_INPUT_SIZE<<"}"<<endl;
  cout<<"\\def \\"<<tablename<<"table {"<<endl;
  cout<<"\\begin{tikzpicture}[scale = .8]"
      <<endl<<"\\begin{axis}["
      <<endl<<"title={Speedup Versus Number of Threads},"
      <<endl<<"xtick pos=left,"
      <<endl<<"ytick pos=left,"
      <<endl<<"legend style={draw=none},"
      <<endl<<"axis line style = { draw = none },"
      <<endl<<"legend pos= north west,"
      <<endl<<"xtick = data,"
      <<endl<<"xlabel={Number of Threads},"
      <<endl<<"ylabel={Speedup Over Serial Partition},"
      <<endl<<"ymax = 4,"
      <<endl<<"legend columns = 2,"
      <<endl<<"scatter/classes=%"
      <<endl<<"{a={mark=o,draw=blue}}]"<<endl;
  
  
  for (int i = 0; i <= 2; i++) {
    if (i == 0) cout<<"%% In-Place"<<endl;
    if (i == 1) cout<<"%% In-Place Prefix-Sum"<<endl;
    if (i == 2) cout<<"%% Out-of-Place"<<endl;
    cout<<"\\addplot coordinates {";
    for (int num_cores = 1; num_cores <= NUM_THREADS_DEFAULT; num_cores++) {
      set_num_threads(num_cores);
      double millisecs = test_partition_multiple_trials(i, MAX_INPUT_SIZE, NUM_TRIALS, false);
      double speedup = serial_baseline / millisecs;
      cout << "( " << num_cores  << ", " << speedup << ") ";
    }
    cout<<"};"<<endl;
  }
  cout<<"\\legend{Low-Space, Med-Space, High-Space}"<<endl;
  cout<<"\\end{axis}"<<endl<<"\\end{tikzpicture}"<<endl<<"}"<<endl;
}


/* 
 *  Builds the table CILKtabletwo in the paper.
 */ 
void parallel_tests2(string tablename) {
  cout<<"\\def \\"<<tablename<<"blocksizetwo {"<<BLOCK_SIZE<<"}"<<endl;
  cout<<"\\def \\"<<tablename<<"numtrialstwo {"<<NUM_TRIALS<<"}"<<endl;
  cout<<"\\def \\"<<tablename<<"numcorestwo {"<<NUM_THREADS_DEFAULT<<"}"<<endl;
  cout<<"\\def \\"<<tablename<<"inputsizetwo {"<<MAX_INPUT_SIZE<<"}"<<endl;
  cout<<"\\def \\"<<tablename<<"tabletwo {"<<endl;
  cout<<"\\begin{tikzpicture}[scale = .8]"
      <<endl<<"\\begin{axis}["
      <<endl<<"title={Speedup Versus Input Size},"
      <<endl<<"xtick pos=left,"
      <<endl<<"ytick pos=left,"
      <<endl<<"legend style={draw=none},"
      <<endl<<"axis line style = { draw = none },"
      <<endl<<"legend pos= north west,"
      <<endl<<"xtick = data,"
      <<endl<<"xlabel={Number of Threads},"
      <<endl<<"ylabel={Speedup Over Serial Partition},"
      <<endl<<"ymax = 5,"
      <<endl<<"ymin = 0,"
      <<endl<<"legend columns = 2,"
      <<endl<<"scatter/classes=%"
      <<endl<<"{a={mark=o,draw=blue}}]"<<endl;
  vector<double> baselines(100, 0);
  set_num_threads(NUM_THREADS_DEFAULT);  
  for (int i = -1; i <= 2; i++) {
    if (i == -1) cout<<"%% Serial Baseline"<<endl<<"%% baselines in ms: ";
    if (i == 0) cout<<"%% In-Place"<<endl;
    if (i == 1) cout<<"%% In-Place Prefix-Sum"<<endl;
    if (i == 2) cout<<"%% Out-of-Place"<<endl;
    cout<<"\\addplot coordinates {";
    for (int input_size_log = 23; ((uint64_t)1 << input_size_log) <= MAX_INPUT_SIZE; input_size_log++) {
      int input_size = (1 << input_size_log);
      double millisecs = test_partition_multiple_trials(i, input_size, NUM_TRIALS, false);
      if (i == -1) {
        baselines[input_size_log] = millisecs;
        assert (input_size_log < 100); // Because I was lazy in allocating space to baselines
        cout << "( " << input_size_log  << ", " << millisecs << " ) ";
      } else {
        double speedup = baselines[input_size_log] / millisecs;
        cout << "( " << input_size_log  << ", " << speedup << ") ";
      }
    }
    cout<<"};"<<endl;
  }
  cout<<"\\legend{Low-Space, Med-Space, High-Space}"<<endl;
  cout<<"\\end{axis}"<<endl<<"\\end{tikzpicture}"<<endl<<"}"<<endl;
}


/*
 * Builds the table serialtable in the paper.
 */
void tests_in_serial() {
  cout<<"\\def \\serialnumtrials {"<<NUM_TRIALS<<"}"<<endl;
  cout<<"\\def \\serialblocksize {"<<BLOCK_SIZE<<"}"<<endl;
  cout<<"\\def \\serialtable {"<<endl;
  cout<<"\\begin{tikzpicture}[scale = .8]"
      <<endl<<"\\begin{axis}["
      <<endl<<"title={Slowdown Versus Input Size in Serial},"
      <<endl<<"xtick pos=left,"
      <<endl<<"ytick pos=left,"
      <<endl<<"ymax = 4,"
      <<endl<<"ymin = 0,"
      <<endl<<"legend style={draw=none},"
      <<endl<<"axis line style = { draw = none },"
      <<endl<<"legend pos= north west,"
      <<endl<<"xtick = data,"
      <<endl<<"xlabel={Log Input Size},"
      <<endl<<"ylabel={Slowdown Over Serial Partition},"
      <<endl<<"legend columns = 2,"
      <<endl<<"scatter/classes=%"
      <<endl<<"{a={mark=o,draw=blue}}]"<<endl;

  vector<double> baselines(100, 0);
  for (int i = -1; i <= 2; i++) {
    if (i == -1) cout<<"%% Serial Baseline"<<endl<<"%% baselines in ms: ";
    if (i == 0) cout<<"%% In-Place"<<endl;
    if (i == 1) cout<<"%% In-Place Prefix-Sum"<<endl;
    if (i == 2) cout<<"%% Out-of-Place"<<endl;
    cout<<"\\addplot coordinates {";
    for (int input_size_log = 23; ((uint64_t)1 << input_size_log) <= MAX_INPUT_SIZE; input_size_log++) {
      int input_size = (1 << input_size_log);
      double millisecs = test_partition_multiple_trials(i, input_size, NUM_TRIALS, false);
      if (i == -1) {
        baselines[input_size_log] = millisecs;
        assert (input_size_log < 100); // Because I was lazy in allocating space to baselines
        cout << "( " << input_size_log  << ", " << millisecs << " ) ";
      } else {
        double slowdown = millisecs / baselines[input_size_log];
        cout << "( " << input_size_log  << ", " << slowdown << ") ";
      }
    }
    cout<<"};"<<endl;
  }
  cout<<"\\legend{Low-Space, Med-Space, High-Space}"<<endl;
  cout<<"\\end{axis}"<<endl<<"\\end{tikzpicture}"<<endl<<"}"<<endl;
}

int main() {
  srand (time(NULL)); // initialize random seed
  // BLOCK_SIZE needs to be a power of two:
  assert(((BLOCK_SIZE - 1) & BLOCK_SIZE) == 0);

  // We begin by running tests of correctness on the functions
  test_prefix_sum();
  test_partition(-2, 141123, false); // Make sure to test on a non-power of two here for completeness.
  test_partition(-1, 141123, false); 
  test_partition(0, 141123, false);
  test_partition(0, 141123, true); // This is the case where the value of more_succ changes which case of the code is run
  test_partition(1, 141123, false);
  test_partition(2, 141123, false);
#ifdef USE_CILK
  parallel_tests("CILK"); // Run the parallel tests!
  parallel_tests2("CILK"); // Run the parallel tests!
#endif
#ifdef RUN_SERIAL
  tests_in_serial(); // Run the serial tests!
#endif
  // This is just some code I like to run when messing around with things:
  // set_num_threads(2);
  //int64_t n = (1 << 25);
  //cout<<test_partition(-1, n, false)<<endl;
  //cout<<test_partition(0, n, false)<<endl;
  // cout<<test_partition(0, n, true)<<endl;
  //cout<<test_partition(1, n, false)<<endl;
  //cout<<test_partition(2, n, false)<<endl;
  return 0;
}
