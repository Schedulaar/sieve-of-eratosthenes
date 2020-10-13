#include <bsp.hpp>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>

int PRINT_PRIMES = 0;
bool COUNT_PRIMES = false;
bool GENERATE_TWINS = false;

long P; // number of processors requested
long n;
double timeTaken;
long totalNumberPrimes;

void bsp_sieve() {
  bsp_begin(P);
  double time0 = bsp_time();

  long p = bsp_nprocs(); // p = number of processors
  long s = bsp_pid();    // s = processor number

  long arrayLength = (long) ceil(((double) (n - 1)) / p);
  long blockStart = s * arrayLength + 2;

  bool *crosses = new bool[arrayLength]();
  long currPrime = 2, currCoordinator = 0, lastCoordinator = 0;
  bsp_push_reg(&currPrime, sizeof(long));
  bsp_push_reg(&currCoordinator, sizeof(long));
  bsp_sync();

  double maxSievePrime = (long) sqrt(n);
  while (currPrime < maxSievePrime) {
    // Sieve, if the coordinator didn't change.
    if (lastCoordinator == currCoordinator) {
      long j = currPrime * currPrime - blockStart;
      if (j < 0) {
        long mod = blockStart % currPrime;
        j = (mod == 0) ? 0 : currPrime - mod;
      }
      while (j < arrayLength) {
        crosses[j] = true;
        j += currPrime;
      }
    }

    // Change current prime or coordinator
    if (s == currCoordinator) {
      // Search for the next prime to sieve
      long j = currPrime - blockStart + (lastCoordinator == currCoordinator);
      while (j < arrayLength && crosses[j]) j++;
      currPrime = blockStart + j;

      lastCoordinator = currCoordinator;
      currCoordinator = j < arrayLength ? s : s + 1;

    } else lastCoordinator = currCoordinator;
    bsp_sync();

    if (s == lastCoordinator) {
      for (long k = 0; k < p; k++) {
        bsp_put(k, &currCoordinator, &currCoordinator, 0, sizeof(long));
        bsp_put(k, &currPrime, &currPrime, 0, sizeof(long));
      }
    }
    bsp_sync();
  }

  double time1 = bsp_time();
  if (s == 0) {
    timeTaken = time1 - time0;
  }

  // Accumulate all local primes in list
  std::vector<long> primes;
  primes.reserve(n / log(n)); // n/ln(n) is average number of primes <= n
  for (long j = 0; j < arrayLength; j++) {
    if (!crosses[j] && blockStart+j <= n) // We might have few numbers >n
      primes.push_back(blockStart+j);
  }
  delete[] crosses;
  bsp_pop_reg(&currPrime);
  bsp_pop_reg(&currCoordinator);
  bsp_sync();

  if (PRINT_PRIMES) {
    for (long i = 0; i < n - 1; i++) {
      if (i / arrayLength == s && !crosses[i % arrayLength]) {
        printf("%li is prime on proc %li.\n", i + 2, s);
      }
      bsp_sync();
    }
  }

  primes.clear();
  bsp_end();
}

void sieve_optimized(long prime, bool * crosses,
                     long blockStart, long arrayLength) {
  // Search for first multiple of prime in local array
  // Start with prime^2 and check if it is prior to the local array
  long j = (prime*prime - blockStart) / 2;
  if (j < 0) {
    // If so, take the first multiple of prime in local array:
    long mod = blockStart % prime;
    // If the remainder `mod` is ...
    if (mod == 0) j = 0;// zero, j=0 represents the first multiple
    else if (mod%2 == 1) j = (prime - mod)/2; // odd, advance by prime-mod.
    else j = (2 * prime - mod) / 2; // even, advance by 2*prime-mod
  }
  // Now cross out all multiples of the prime starting with j
  while (j < arrayLength) {
    crosses[j] = true;
    j += prime; // Step size: 2*prime/2; ignore even multiples.
  }
}

/* Returns list of twins represented by their mean */
std::vector<long> bsp_twins(std::vector<long> primes) {
  long p = bsp_nprocs(), s = bsp_pid();
  // Get largest prime of previous processor
  long prevPrime = s == 0 ? 2 : 0;
  bsp_push_reg(&prevPrime, sizeof(long));
  bsp_sync();

  if (s < p - 1 && primes.size() > 0)
    bsp_put(s + 1, &primes[primes.size() - 1],
            &prevPrime, 0, sizeof(long));
  bsp_sync();

  std::vector<long> twins;
  if (primes.size() == 0)
    return twins;

  if (prevPrime + 2 == primes[0])
    twins.push_back(prevPrime + 1);
  for (long i = 0; i < primes.size() - 1; i++) {
    if (primes[i] + 2 == primes[i + 1])
      twins.push_back(primes[i] + 1);
  }

  bsp_pop_reg(&prevPrime);
  bsp_sync();
  return twins;
}

void bsp_sieve_optimized() {
  bsp_begin(P);
  double time0 = bsp_time();
  long p = bsp_nprocs(), s = bsp_pid();

  // We only consider odd numbers to save on memory and enhance performance
  long globalArrayLength = (n - 1) / 2;
  long arrayLength = (long) ceil(((double) globalArrayLength) / p);
  long blockStart = 2 * s * arrayLength + 3;
  bool *crosses = new bool[arrayLength]();
  // We start sieving with prime 3; the first coordinator is processor 0.
  long currPrime = 3, currCoordinator = 0, lastCoordinator = 0;
  bsp_push_reg(&currPrime, sizeof(long));
  bsp_push_reg(&currCoordinator, sizeof(long));
  bsp_sync();

  long maxSievePrime = (long) sqrt(n);
  while (currPrime <= maxSievePrime) {
    // Sieve, if the coordinator did not change
    if (lastCoordinator == currCoordinator)
      sieve_optimized(currPrime, crosses, blockStart, arrayLength);

    // Let coordinator change current prime or coordinator
    if (s == currCoordinator) {
      // Search for the next prime in local array
      long j = (currPrime - blockStart) / 2
               + (lastCoordinator == currCoordinator);
      while (j < arrayLength && crosses[j]) j++;
      // Translate back to normal representation
      currPrime = blockStart + j * 2;
      lastCoordinator = currCoordinator;
      // If j exceeds the local array, transfer coordination to next proc.
      currCoordinator = j < arrayLength ? s : s + 1;
    } else lastCoordinator = currCoordinator;
    bsp_sync();

    // Distribute information
    if (s == lastCoordinator) {
      for (long k = 0; k < p; k++) {
        bsp_put(k, &currCoordinator, &currCoordinator, 0, sizeof(long));
        bsp_put(k, &currPrime, &currPrime, 0, sizeof(long));
      }
    }
    bsp_sync();
  }

  // Accumulate all local primes in list
  std::vector<long> primes;
  primes.reserve(n / log(n)); // n/ln(n) is average number of primes <= n
  for (long j = 0; j < arrayLength; j++) {
    if (!crosses[j] && blockStart + 2 * j <= n) // We might have few numbers >n
      primes.push_back(blockStart + 2 * j);
  }
  delete[] crosses;
  bsp_pop_reg(&currPrime);
  bsp_pop_reg(&currCoordinator);
  bsp_sync();

  double time1 = bsp_time();
  if (s == 0) {
    timeTaken = time1 - time0;
  }

  if (GENERATE_TWINS) {
    std::vector<long> twins = bsp_twins(primes);
    for (long k = 0; k < p; k++) {
      if (s == k) {
        for (long i = 0; i < twins.size(); i++)
          printf("%li±1 is twin-prime on proc %li.\n", twins[i], s);
      }
      bsp_sync();
    }
    twins.clear();
  }


  if (COUNT_PRIMES) {
    long numberPrimes[p];
    memset(numberPrimes, 0, sizeof(numberPrimes));
    bsp_push_reg(numberPrimes, p * sizeof(long));
    bsp_sync();

    numberPrimes[s] = primes.size();
    bsp_put(0, &numberPrimes[s], numberPrimes, s*sizeof(long), sizeof(long));
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
    for (long k = 0; k < p; k++) {
      if (s == k) {
        for (long i = 0; i < primes.size(); i++)
          printf("%li is prime on proc %li.\n", primes[i], s);
      }
      bsp_sync();
    }
  }

  primes.clear();
  bsp_end();
}

int main(int argc, char **argv) {
  printf("Please enter your config (or d for defaults): minP,maxP,minN,maxN,nIters,print,opt\n");
  fflush(stdout);
  long maxP = bsp_nprocs();
  long minP = 1;
  long minN = 1;
  long maxN = 10e8;
  long nIters = 100;
  int opt = 1;
  scanf("%ld,%ld,%ld,%ld,%ld,%i, %i", &minP, &maxP, &minN, &maxN, &nIters, &PRINT_PRIMES, &opt);
  if (P > bsp_nprocs()) {
    printf("Sorry, only %u processors available.\n", bsp_nprocs());
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  printf("[\n");
  for (n = minN; n <= maxN; n *= 10) {
    for (P = minP; P <= maxP; P = (2*P <= maxP || P == maxP) ? 2*P : maxP ) {
      for (int i = 0; i < nIters; i++) {
        if (opt) {
          bsp_init(bsp_sieve_optimized, argc, argv);
          bsp_sieve_optimized();
        } else {
          bsp_init(bsp_sieve, argc, argv);
          bsp_sieve();
        }
        printf("  { \"p\": %li, \"n\": %li, \"t\": %f, \"iter\": %i },\n", P, n, timeTaken, i);
      }
    }
  }
  printf("]\n");
  return 0;
}
