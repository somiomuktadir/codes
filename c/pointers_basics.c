/*
 * Pointers Basics
 * Demonstrates introduction to pointers and swapping variables using call-by-reference.
 */

#include <stdio.h>

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

int main() {
    int x = 10, y = 20;
    int *ptrX = &x;

    printf("Initial values:\n");
    printf("x = %d, y = %d\n", x, y);
    
    printf("\nPointer demonstration:\n");
    printf("Address of x: %p\n", (void*)&x);
    printf("Value of ptrX: %p\n", (void*)ptrX);
    printf("Value pointed to by ptrX: %d\n", *ptrX);

    printf("\nSwapping using pointers...\n");
    swap(&x, &y);

    printf("Values after swap:\n");
    printf("x = %d, y = %d\n", x, y);

    return 0;
}
