#include <cstdio>
#include <cmath>
#include <vector>

long numberOps = 0;

std::vector<long> sieve(long n) {
  bool *crosses = new bool[n - 1]();
  // list[i] == true means i+2 is crossed out.
  // new bool[n-1]() automatically initializes the array with false.
  double maxIter = sqrt(n) - 2;

  numberOps += n + 1;
  for (long i = 0; i <= maxIter; i++) {
    if (!crosses[i]) { // cross out all larger multiples of i+2
      for (int j = (i+2)*(i+2) - 2; j < n - 1; j += (i+2)) {
        crosses[j] = true;
        numberOps += 5;
      }
    }
    numberOps += 3;
  }

  std::vector<long> primes;
  primes.reserve(n / log(n)); // n/ln(n) is average number of primes <= n
  for (long i = 0; i < n - 1; i++) {
    if (!crosses[i]) primes.push_back(i+2);
    numberOps += 4;
  }

  delete [] crosses;
  return primes;
}


int main() {
  long lastN = 0;
  long lastOps = 0;
  for (long n = 2; n <= 10e10; n *= 2) {
    std::vector<long> primes = sieve(n);
    double differenceQuotient = ((double) numberOps - lastOps) / (n - lastN);
    printf("(%f,%f)\n", ((double) n) / 2 + 0.001, numberOps, differenceQuotient);
    printf("(%li,%f)\n", n, numberOps, differenceQuotient);
    primes.clear();
    numberOps = 0;
  }
  return 0;
}
