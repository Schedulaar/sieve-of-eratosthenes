#include <bsp.hpp>
#include <cstdio>
#include <cmath>
#include <cstring>

bool PRINT_PRIMES = false;
bool COUNT_PRIMES = false;

long P; // number of processors requested
long n;
double timeTaken;
long totalNumberPrimes;

void bspsieve_coord() {
  bsp_begin(P);
  double time0 = bsp_time();

  long p = bsp_nprocs(); // p = number of processors
  long s = bsp_pid();    // s = processor number

  long arrayLength = (long) ceil(((double) (n - 1)) / p);
  long blockStart = s * arrayLength + 2;

  bool *list = new bool[arrayLength]();
  long currPrimeAndCoord[2] = {2, 0};
  long lastCoordinator = 0;
  bsp_push_reg(currPrimeAndCoord, 2 * sizeof(long));
  bsp_sync();

  double maxSievePrime = sqrt(n);

  while (currPrimeAndCoord[0] < maxSievePrime) {
    // Sieve, if the coordinator didn't change.
    if (lastCoordinator == currPrimeAndCoord[1]) {
      long p = currPrimeAndCoord[0];
      long j = p * p - blockStart;
      if (j < 0) {
        long mod = blockStart % p;
        j = (mod == 0) ? 0 : p - mod;
      }
      while (j < arrayLength) {
        list[j] = true;
        j += p;
      }
    }

    // Change current prime or coordinator
    if (s == currPrimeAndCoord[1]) {
      // Search for the next prime to sieve
      long j = currPrimeAndCoord[0] - blockStart + (lastCoordinator == currPrimeAndCoord[1]);
      while (j < arrayLength && list[j]) j++;
      currPrimeAndCoord[0] = blockStart + j;

      lastCoordinator = currPrimeAndCoord[1];
      currPrimeAndCoord[1] = j < arrayLength ? s : s + 1;

      long k = s + 1;
      while (k < p) {
        bsp_put(k, currPrimeAndCoord, currPrimeAndCoord, 0, 2 * sizeof(long));
        k++;
      }
      bsp_sync();
    } else {
      lastCoordinator = currPrimeAndCoord[1];
      bsp_sync();
    }
  }

  double time1 = bsp_time();
  if (s == 0) {
    timeTaken = time1 - time0;
    printf("coord needed %fs. \n", timeTaken);
  }

  if (PRINT_PRIMES) {
    for (long i = 0; i < n - 1; i++) {
      if (i / arrayLength == s && !list[i % arrayLength]) {
        printf("%li is prime on proc %li.\n", i + 2, s);
      }
      bsp_sync();
    }
  }

  bsp_pop_reg(currPrimeAndCoord);
  delete[] list;
  bsp_end();
}

void bspsieve_coord_ignore_even() {
  bsp_begin(P);
  double time0 = bsp_time();

  long p = bsp_nprocs(); // p = number of processors
  long s = bsp_pid();    // s = processor number
  long globalArrayLength = (n - 1) / 2;

  long arrayLength = (long) ceil(((double) globalArrayLength) / p);
  long blockStart = 2 * s * arrayLength + 3;

  bool *list = new bool[arrayLength]();
  long currPrimeAndCoord[2] = {3, 0};
  long lastCoordinator = 0;
  bsp_push_reg(currPrimeAndCoord, 2 * sizeof(long));
  bsp_sync();

  double maxSievePrime = sqrt(n);

  while (currPrimeAndCoord[0] <= maxSievePrime) {
    // Sieve, if the coordinator didn't change.
    if (lastCoordinator == currPrimeAndCoord[1]) {
      long q = currPrimeAndCoord[0];
      long j = (q * q - blockStart) / 2;
      if (j < 0) {
        long mod = blockStart % q;
        if (mod == 0) j = 0;
        else if (mod % 2 == 1) j = (q - mod) / 2;
        else j = (2 * q - mod) / 2;
      }
      while (j < arrayLength) {
        list[j] = true;
        j += q;
      }
    }

    // Let coordinator change current prime or coordinator
    if (s == currPrimeAndCoord[1]) {
      // Search for the next prime to sieve
      long j = (currPrimeAndCoord[0] - blockStart) / 2 + (lastCoordinator == currPrimeAndCoord[1]);
      while (j < arrayLength && list[j]) j++;
      currPrimeAndCoord[0] = blockStart + j * 2;

      lastCoordinator = currPrimeAndCoord[1];
      currPrimeAndCoord[1] = j < arrayLength ? s : s + 1;

      long k = s + 1;
      while (k < p) {
        bsp_put(k, currPrimeAndCoord, currPrimeAndCoord, 0, 2 * sizeof(long));
        k++;
      }
      bsp_sync();
    } else {
      lastCoordinator = currPrimeAndCoord[1];
      bsp_sync();
    }
  }

  bsp_pop_reg(currPrimeAndCoord);

  double time1 = bsp_time();
  if (s == 0) {
    timeTaken = time1 - time0;
  }

  if (COUNT_PRIMES) {
    long numberPrimes[p];
    memset(numberPrimes, 0, sizeof(numberPrimes));
    bsp_push_reg(numberPrimes, p * sizeof(long));
    bsp_sync();

    for (long i = 0; i < arrayLength && blockStart + 2 * i <= n; i++) {
      numberPrimes[s] += (list[i] == 0);
    }
    bsp_put(0, &numberPrimes[s], numberPrimes, s * sizeof(long), sizeof(long));
    bsp_sync();

    if (s == 0) {
      totalNumberPrimes = 1;
      for (long j = 0; j < p; j++) {
        totalNumberPrimes += numberPrimes[j];
      }
    }
    bsp_pop_reg(numberPrimes);
  }

  if (PRINT_PRIMES) {
    if (s == 0)
      printf("2 is prime by definition.\n");
    for (long i = 0; i < (n - 1) / 2; i++) {
      if (i / arrayLength == s && !list[i % arrayLength]) {
        printf("%li is prime on proc %li.\n", 2 * i + 3, s);
      }
      bsp_sync();
    }
  }
  delete[] list;
  bsp_end();
}

int main(int argc, char **argv) {
  long NITERS = 100;
  long maxP = bsp_nprocs();
  for (n = 10; n <= 10e8; n*= 10) {
    printf("----------------------\n");
    printf("       n=%li\n", n);
    printf("----------------------\n");
    for (P = 1; P <= maxP; P = (2*P <= maxP || P == maxP) ? 2*P : maxP ) {
      double averageTime = 0;
      for (int i = 0; i < NITERS; i++) {
        bsp_init(bspsieve_coord_ignore_even, argc, argv);
        bspsieve_coord_ignore_even();
        averageTime += timeTaken;
      }
      printf("p=%li: t=%f\n", P, averageTime / NITERS);
    }
  }
  return 0;

  /*printf("How many processors do you want to use?\n");
  fflush(stdout);
  scanf("%ld", &P);
  if (P > bsp_nprocs()) {
    printf("Sorry, only %u processors available.\n", bsp_nprocs());
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  n = 1000000000;

  bsp_init(bspsieve_coord_ignore_even, argc, argv);
  bspsieve_coord_ignore_even();
  exit(EXIT_SUCCESS);*/
}