/*
 * Arrays 1D Example
 * Demonstrates basic 1D array operations: input, display, sum, average, min/max.
 */

#include <stdio.h>

#define MAX_SIZE 100

int main() {
    int arr[MAX_SIZE];
    int n, i;
    int sum = 0;
    float average;
    int min, max;

    // Input size
    printf("Enter size of array (1-%d): ", MAX_SIZE);
    if (scanf("%d", &n) != 1 || n <= 0 || n > MAX_SIZE) {
        printf("Invalid size.\n");
        return 1;
    }

    // Input elements
    printf("Enter %d elements:\n", n);
    for (i = 0; i < n; i++) {
        printf("Element %d: ", i + 1);
        scanf("%d", &arr[i]);
    }

    // Calculate sum, min, max
    min = arr[0];
    max = arr[0];
    
    for (i = 0; i < n; i++) {
        sum += arr[i];
        
        if (arr[i] < min) {
            min = arr[i];
        }
        if (arr[i] > max) {
            max = arr[i];
        }
    }

    average = (float)sum / n;

    // Display results
    printf("\nArray elements: ");
    for (i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    printf("Sum: %d\n", sum);
    printf("Average: %.2f\n", average);
    printf("Minimum: %d\n", min);
    printf("Maximum: %d\n", max);

    return 0;
}
