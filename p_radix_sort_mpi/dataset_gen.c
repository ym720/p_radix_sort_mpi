#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/**
 * This program is used to generate a dataset of random numbers by providing
 * number of elements, min-max range and optionally, number of digits if we
 * need to generate elements of fixed size.
 */
void usage() {
  fprintf(stderr, "Incorrect usage!\n");
  fprintf(stderr, "Usage: dataset_gen [n] [min] [max]\n");
  fprintf(stderr, "  [n] - number of elements in the dataset\n");
  fprintf(stderr, "  [min] - min value of element\n");
  fprintf(stderr, "  [max] - max value of element\n");
  fprintf(stderr, "  <digits> - optional, number of digits in each element\n");
}

// generate random number within a given range
int random_in_range(int min, int max) {
  return (rand() % (max + 1 - min)) + min;
}

int main(int argc, char** argv) {

  // check for correct number of arguments
  if (argc < 4) {
    usage();
    return -1;
  }

  // initialize seed
  srand(time(0));

  int n = atoi(argv[1]);
  int min = atoi(argv[2]);
  int max = atoi(argv[3]);
  int digits = -1;
  if (argc > 4) {
    digits = atoi(argv[4]);
  }

  // if fixed digits length supplied, generate numbers pad with zeros
  // otherwise generate random numbers of variable length
  if (digits > 0) {
    for (int i = 0; i < n; i++) {
      printf("%d\n", random_in_range(min, max));
    }
  } else {
    for (int i = 0; i < n; i++) {
      printf("%0*d\n", digits, random_in_range(min, max));
    }
  }

  return 0;
}
