#include <cstdio>
#include <cmath>

bool isPrime(unsigned int n) {
  bool *list = new bool[n - 1]();
  // list[i] == false means i+2 is not crossed out.
  // new bool[n-1]() automatically initializes the array with false.
  double maxIter = sqrt(n) - 2;
  for (int i = 0; i <= maxIter; i++) {
    if (list[i] == 0) { // cross out all larger multiples of i+2
      for (int j = (i + 2) * (i + 2) - 2; j < n; j += (i + 2)) {
        list[j] = true;
      }
    }
  }
  bool result = !list[n - 2];
  delete [] list;
  return result;
}


int main() {
  for (int i = 2; i < 101; i++) {
    if (isPrime(i)) {
      printf("%i is a prime.\n", i);
    }
  }

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
