/**
 * Serial implementation of radix sort.
 *
 * Author: Yourii Martiak
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "timing.h"

#define b 32           // number of bits for integer, usually four words
#define g 8            // group of bits for each scan, one word
#define B (1 << g)     // number of buckets, 2^g

void usage() {
  fprintf(stderr, "Incorrect usage!\n");
  fprintf(stderr, "Usage: radix_sort [f] [n]\n");
  fprintf(stderr, "  [f] - input file to be sorted\n");
  fprintf(stderr, "  [n] - number of elements in the file\n");
  fprintf(stderr, "  [r] - print sorted results 0/1, 0 by default\n");
}

void print_array(int *a, const int n) {
  for (int i = 0; i < n; i++) {
    printf("%d\n", a[i]);
  } 
}

int init_array(char* file, const int begin, const int n, int *a) {

  // open file in read-only mode and check for errors
  FILE *file_ptr;
  file_ptr = fopen(file, "r");
  if (file_ptr == NULL) {
    return EXIT_FAILURE;
  }

  // read n numbers from a file into array a starting at begin position
  int skip;

  // first skip to the begin position
  for (int i = 0; i < begin; i++) {
    int s = fscanf(file_ptr, "%d", &skip); 
  }
  // then read numbers into array a
  for (int i = 0; i < n; i++) {
    int s = fscanf(file_ptr, "%d", &a[i]);
  }

  return EXIT_SUCCESS;
}

// Compute j bits which appear k bits from the right in x
// Ex. to obtain rightmost bit of x call bits(x, 0, 1)
unsigned bits(unsigned x, int k, int j) {
  return (x >> k) & ~(~0 << j);
}

void radix_sort(int *a, const int n) {
  int* t = malloc(n*sizeof(int));     // temp array used for sorting
  int count[B];                       // array of counts per bucket

  for (int pass = 0; pass < b/g; pass++) {       // each pass
    for (int j = 0; j < B; j++) count[j] = 0;    // init counts array
    for (int i = 0; i < n; i++) {
      count[bits(a[i], pass*g, g)]++;            // count keys per bucket
    }
    for (int j = 1; j < B; j++) {
      count[j] = count[j-1] + count[j];          // compute prefix sum
    }
    for (int i = n-1; i >= 0; i--) {
      int idx = --count[bits(a[i], pass*g, g)];
      t[idx] = a[i];                             // transpose to temp array
    }
    for (int i = 0; i < n; i++) a[i] = t[i];     // copy back to master
  }
  free(t);
}

int main(int argc, char** argv)
{
  int print_results = 0;

  // check for correct number of arguments
  if (argc < 3) {
    usage();
    return EXIT_FAILURE;
  } else if (argc > 3) {
    print_results = atoi(argv[3]);
  }

  // initialize vars and allocate memory
  const int n = atoi(argv[2]);
  int* a = malloc(sizeof(int) * n);

  // initialize local array
  if (init_array(argv[1], 0, n, &a[0]) != EXIT_SUCCESS) {
    printf("File %s could not be opened!\n", argv[1]);
    return EXIT_FAILURE;
  }

  // take a timestamp before the sort starts
  timestamp_type time1, time2;
  get_timestamp(&time1);

  // sort elements
  radix_sort(&a[0], n);
  
  // take a timestamp after the process finished sorting
  get_timestamp(&time2);

  // calculate fish updates per second
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("%f s\n", elapsed);
  printf("%d elements sorted\n", n);
  printf("%f elements/s\n", n / elapsed);

  // print sorted resutls
  if (print_results) {
    print_array(&a[0], n);
  }

  // release resources no longer used
  free(a);

  return 0;
}
