/*
 * Fibonacci Sequence Generator
 * Generates the first n numbers in the Fibonacci sequence
 * Uses iteration for efficiency
 */

#include <stdio.h>

int main() {
    int n, i;
    long long int first = 0, second = 1, next;

    // Get number of terms from user
    printf("Enter the number of Fibonacci terms to generate: ");
    if (scanf("%d", &n) != 1) {
        printf("Error: Invalid input. Please enter a positive integer.\n");
        return 1;
    }

    // Validate input
    if (n <= 0) {
        printf("Error: Please enter a positive integer.\n");
        return 1;
    }

    if (n > 93) {
        printf("Warning: Values beyond 93 terms may overflow.\n");
    }

    printf("\nFibonacci sequence (first %d terms):\n", n);
    
    for (i = 0; i < n; i++) {
        if (i <= 1) {
            next = i;
        } else {
            next = first + second;
            first = second;
            second = next;
        }
        printf("%d: %lld\n", i + 1, next);
    }
    
    return 0;
}
