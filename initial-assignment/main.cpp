#include <cstdio>
#include <cmath>

bool isPrime(unsigned int n) {
  bool *list = new bool[n - 1]();
  // list[i] == false means i+2 is not crossed out.
  // new bool[n-1]() automatically initializes the array with false.
  double maxIter = sqrt(n) - 2;
  for (int i = 0; i <= maxIter; i++) {
    if (!list[i]) { // cross out all larger multiples of i+2
      for (int j = (i + 2) * (i + 2) - 2; j < n; j += (i + 2)) {
        list[j] = true;
      }
    }
  }

  for (int i = 0; i <= n-1; i++) {
    if (!list[i]) {
      printf("%i is prime. \n", i + 2);
    }
  }
  bool result = !list[n - 2];
  delete [] list;
  return result;
}


int main() {
  isPrime(1000000000);

  printf("Type in some number to check for primality!\n");
  unsigned int n;
  if (!scanf("%i", &n)) {
    printf("Something went wrong.");
    return 1;
  }
  if (isPrime(n)) {
    printf("The number you entered is prime.");
  } else {
    printf("The number you entered is not prime.");
  }
  return 0;
}
