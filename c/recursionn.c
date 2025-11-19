/*
 * Recursive Fibonacci Calculator
 * Calculates the nth Fibonacci number using recursion
 * 
 * Note: This recursive implementation is inefficient for large n
 * due to repeated calculations. For production use, consider
 * iterative approach or memoization.
 * 
 * Key difference: Iteration uses repetition; Recursion uses selection
 */

#include <stdio.h>

// Function prototype
unsigned long long int fibonacci(unsigned int n);

int main() {
    unsigned long long int result;
    unsigned int number;

    // Get input from user
    printf("Enter a number to find its Fibonacci value: ");
    if (scanf("%u", &number) != 1) {
        printf("Error: Invalid input. Please enter a positive integer.\n");
        return 1;
    }

    // Warn about performance for large numbers
    if (number > 40) {
        printf("Warning: Recursive calculation for n > 40 may take a long time.\n");
        printf("Consider using an iterative approach for better performance.\n");
    }

    // Calculate and display result
    result = fibonacci(number);
    printf("Fibonacci(%u) = %llu\n", number, result);

    return 0;
}

/*
 * Function: fibonacci
 * Recursively calculates the nth Fibonacci number
 * Base cases: F(0) = 0, F(1) = 1
 * Recursive case: F(n) = F(n-1) + F(n-2)
 */
unsigned long long int fibonacci(unsigned int n) {
    if (n == 0 || n == 1) {
        return n;
    } else {
        return fibonacci(n - 1) + fibonacci(n - 2);
    }
}