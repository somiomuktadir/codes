/*
 * Recursion Example
 * Calculates factorial using a recursive function.
 */

#include <stdio.h>

long long int factorial(int n) {
    if (n == 0 || n == 1) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}

int main() {
    int n;

    printf("Enter a positive integer: ");
    if (scanf("%d", &n) != 1 || n < 0) {
        printf("Invalid input. Please enter a non-negative integer.\n");
        return 1;
    }

    printf("Factorial of %d is %lld\n", n, factorial(n));

    return 0;
}
