#include <bsp.hpp>
#include <cstdio>
#include <cmath>
#include <unistd.h>

long P; // number of processors requested
long n;

bool isPrime(unsigned int n) {
  bool *list = new bool[n - 1]();
  // list[i] == false means i+2 is not crossed out.
  // new bool[n-1]() automatically initializes the array with false.
  double maxIter = sqrt(n) - 2;
  for (int i = 0; i <= maxIter; i++) {
    if (list[i] == 0) { // cross out all larger multiples of i+2
      for (int j = (i + 2) * (i + 2) - 2; j < n - 1; j += (i + 2)) {
        list[j] = true;
      }
    }
  }
  bool result = !list[n - 2];
  delete[] list;
  return result;
}

void bspsieve_cyc() {
  bsp_begin(P);
  double time0 = bsp_time();

  long p = bsp_nprocs(); // p = number of processors
  long s = bsp_pid();    // s = processor number

  long arrayLength = (long) ceil(((double) (n - 1)) / p);

  bool *list = new bool[arrayLength]();
  bsp_push_reg(list, (arrayLength) * sizeof(bool));
  bsp_sync();
  bool putValue = true;

  double maxIter = sqrt(n) - 2;
  for (int iter = 0; iter <= maxIter / p; iter++) {
    int i = s + p * iter;
    if (i <= maxIter && !list[iter]) { // cross out all larger multiples of i+2
      for (int j = (i + 2) * (i + 2) - 2; j < n - 1; j += (i + 2)) {
        int dest = j % p;
        int offset = j / p;
        if (dest == s)
          list[offset] = true;
        else
          bsp_put((long) dest, &putValue, list, offset * sizeof(bool), sizeof(bool));
      }
    }
    bsp_sync();
  }

  double time1 = bsp_time();
  if (s == 0)
    printf("cyc needed %f ms. \n", time1 - time0);


  /*for (unsigned int i = 0; i < n-1; i++) {
    if (i % p == s && !list[i / p]) {
      printf("%i is prime on proc %li.\n", i + 2, s);
    }
    bsp_sync();
  }*/

  bsp_pop_reg(list);
  delete[] list;
  bsp_end();
}

void bspsieve_block() {
  bsp_begin(P);
  double time0 = bsp_time();

  long p = bsp_nprocs(); // p = number of processors
  long s = bsp_pid();    // s = processor number

  unsigned int arrayLength = (unsigned int) ceil(((double) (n - 1)) / p);

  bool *list = new bool[arrayLength]();
  bsp_push_reg(list, (arrayLength) * sizeof(bool));
  bsp_sync();
  bool putValue = true;

  double maxIter = sqrt(n) - 2;
  for (int iter = 0; iter <= maxIter / p; iter++) {
    int i = s * arrayLength + iter;
    if (i <= maxIter && !list[iter]) { // cross out all larger multiples of i+2
      for (int j = (i + 2) * (i + 2) - 2; j < n - 1; j += (i + 2)) {
        int dest = j / arrayLength;
        int offset = j % arrayLength;
        if (dest == s)
          list[offset] = true;
        else
          bsp_put((long) dest, &putValue, list, offset * sizeof(bool), sizeof(bool));
      }
    }
    bsp_sync();
  }

  double time1 = bsp_time();
  if (s == 0)
    printf("block needed %f ms. \n", time1 - time0);


  /*for (unsigned int i = 0; i < n-1; i++) {
    if (i / arrayLength == s && !list[i % arrayLength]) {
      printf("%i is prime on proc %li.\n", i + 2, s);
    }
    bsp_sync();
  }*/

  bsp_pop_reg(list);
  delete[] list;
  bsp_end();
}

void bspsieve_block_cyc() {
  bsp_begin(P);
  double time0 = bsp_time();

  long p = bsp_nprocs(); // p = number of processors
  long s = bsp_pid();    // s = processor number

  long arrayLength = (long) ceil(((double) (n - 1)) / p);

  bool *list = new bool[arrayLength]();
  bsp_push_reg(list, (arrayLength) * sizeof(bool));
  bsp_sync();
  bool putValue = true;

  double maxIter = sqrt(n) - 2;
  long blockSize = (long) ceil(maxIter / p);

  for (long iter = 0; iter < blockSize; iter++) {
    long i = s * blockSize + iter;
    if (i <= maxIter && !list[iter]) { // cross out all larger multiples of i+2
      for (long j = (i + 2) * (i + 2) - 2; j < n - 1; j += (i + 2)) {
        int dest, offset;
        if (j < blockSize * p) {
          dest = j / blockSize;
          offset = j % blockSize;
        } else {
          dest = (j - blockSize * p) % p;
          offset = blockSize + (j - blockSize * p) / p;
        }
        if (dest == s)
          list[offset] = true;
        else
          bsp_put((long) dest, &putValue, list, offset * sizeof(bool), sizeof(bool));
      }
    }
    bsp_sync();
  }

  double time1 = bsp_time();
  if (s == 0)
    printf("block_cyc needed %f ms. \n", time1 - time0);

  /*for (unsigned int i = 0; i < n-1; i++) {
    if (i < p*blockSize && s == i / blockSize && !list[i % blockSize]) {
      printf("%i is prime on proc %li under blockSize.\n", i + 2, s);
    } else if (i >= p*blockSize && s == (i - p*blockSize) % p && !list[blockSize + (i - p*blockSize) / p]) {
      printf("%i is prime on proc %li over blockSize.\n", i + 2, s);
    }
    bsp_sync();
  }*/

  bsp_pop_reg(list);
  delete[] list;
  bsp_end();
}

void bspsieve_block_block() {
  bsp_begin(P);
  double time0 = bsp_time();

  long p = bsp_nprocs(); // p = number of processors
  long s = bsp_pid();    // s = processor number

  long arrayLength = (long) ceil(((double) (n - 1)) / p);

  bool *list = new bool[arrayLength]();
  bsp_push_reg(list, (arrayLength) * sizeof(bool));
  bsp_sync();
  bool putValue = true;

  double maxIter = sqrt(n) - 2;
  long block1Size = (long) ceil(((double) maxIter) / p);
  long block2Size = (long) ceil(((double) n - 1 - p * block1Size) / p);

  for (long iter = 0; iter < block1Size; iter++) {
    long i = s * block1Size + iter;
    if (i <= maxIter && !list[iter]) { // cross out all larger multiples of i+2
      for (long j = (i + 2) * (i + 2) - 2; j < n - 1; j += (i + 2)) {
        int dest, offset;
        if (j < block1Size * p) {
          dest = j / block1Size;
          offset = j % block1Size;
        } else {
          dest = (j - block1Size * p) / block2Size;
          offset = block1Size + (j - block1Size * p) % block2Size;
        }
        if (dest == s)
          list[offset] = true;
        else
          bsp_put(dest, &putValue, list, offset * sizeof(bool), sizeof(bool));
      }
    }
    bsp_sync();
  }

  double time1 = bsp_time();
  if (s == 0)
    printf("block_block needed %f ms.\n", time1 - time0);

  /*
  for (unsigned int i = 0; i < n-1; i++) {
    if (i < p*block1Size && s == i / block1Size && !list[i % block1Size]) {
      printf("%i is prime on proc %li under blockSize.\n", i + 2, s);
    } else if (i >= p*block1Size && s == (i - p*block1Size) / block2Size && !list[block1Size + (i - p*block1Size) % block2Size]) {
      printf("%i is prime on proc %li over blockSize.\n", i + 2, s);
    }
    bsp_sync();
  }*/

  bsp_pop_reg(list);
  delete[] list;
  bsp_end();
}

int main(int argc, char **argv) {
  printf("How many processors do you want to use?\n");
  fflush(stdout);
  scanf("%ld", &P);
  if (P > bsp_nprocs()) {
    printf("Sorry, only %u processors available.\n",
           bsp_nprocs());
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  n = 100000000;

  bsp_init(bspsieve_block, argc, argv);
  bspsieve_block();
  bsp_init(bspsieve_block_cyc, argc, argv);
  bspsieve_block_cyc();
  bsp_init(bspsieve_cyc, argc, argv);
  bspsieve_cyc();
  bsp_init(bspsieve_block_block, argc, argv);
  bspsieve_block_block();
  exit(EXIT_SUCCESS);

}
