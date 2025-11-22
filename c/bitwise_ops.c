/*
 * Bitwise Operations Example
 * Demonstrates bitwise AND, OR, XOR, NOT, Left Shift, Right Shift.
 */

#include <stdio.h>

void printBinary(unsigned int n) {
    for (int i = 7; i >= 0; i--) {
        if ((n >> i) & 1)
            printf("1");
        else
            printf("0");
    }
    printf("\n");
}

int main() {
    unsigned int a = 12; // 0000 1100
    unsigned int b = 25; // 0001 1001

    printf("a = %u, b = %u\n", a, b);
    printf("Binary a: "); printBinary(a);
    printf("Binary b: "); printBinary(b);

    printf("\nBitwise AND (a & b): %u\n", a & b);
    printf("Binary: "); printBinary(a & b);

    printf("\nBitwise OR (a | b): %u\n", a | b);
    printf("Binary: "); printBinary(a | b);

    printf("\nBitwise XOR (a ^ b): %u\n", a ^ b);
    printf("Binary: "); printBinary(a ^ b);

    printf("\nBitwise NOT (~a): %u (Result depends on integer size)\n", ~a);
    
    printf("\nLeft Shift (a << 2): %u\n", a << 2);
    printf("Binary: "); printBinary(a << 2);

    printf("\nRight Shift (b >> 2): %u\n", b >> 2);
    printf("Binary: "); printBinary(b >> 2);

    return 0;
}
