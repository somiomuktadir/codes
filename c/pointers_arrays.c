/*
 * Pointers and Arrays
 * Demonstrates pointer arithmetic and traversing arrays using pointers.
 */

#include <stdio.h>

int main() {
    int arr[] = {10, 20, 30, 40, 50};
    int n = sizeof(arr) / sizeof(arr[0]);
    int *ptr = arr; // Points to the first element
    int i;

    printf("Array elements using pointer arithmetic:\n");
    for (i = 0; i < n; i++) {
        printf("*(ptr + %d) = %d\n", i, *(ptr + i));
    }

    printf("\nTraversing array by incrementing pointer:\n");
    ptr = arr; // Reset pointer
    for (i = 0; i < n; i++) {
        printf("Address: %p, Value: %d\n", (void*)ptr, *ptr);
        ptr++; // Move to next integer
    }

    return 0;
}
