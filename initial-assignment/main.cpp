#include <stdio.h>
#include <stdbool.h>
#include <malloc.h>
#include <math.h>

bool isPrime(unsigned int n) {
    bool *list = calloc(n, sizeof(bool)); // list[i] == false means i+2 is not crossed out.
    // calloc automatically initializes the array with false.
    double maxIter = sqrt(n) - 2;
    for (int i = 0; i <= maxIter; i++) {
        if (list[i] == 0) { // cross out all larger multiples of i+2
            for (int j = i + (i + 2); j < n; j += (i + 2)) {
                list[j] = true;
            }
        }
    }
    bool result = list[n - 2] == false;
    free(list);
    return result;
}


int main() {
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
