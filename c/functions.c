/*
 * Maximum of Three Numbers
 * Demonstrates function usage to find the largest of three integers
 */

#include <stdio.h>

// Function prototype
int maximum(int x, int y, int z);

int main() {
    int a, b, c;
    
    // Get three numbers from user
    printf("Enter three integers separated by spaces: ");
    if (scanf("%d%d%d", &a, &b, &c) != 3) {
        printf("Error: Invalid input. Please enter three integers.\n");
        return 1;
    }
    
    // Display result
    printf("The maximum of %d, %d, and %d is: %d\n", a, b, c, maximum(a, b, c));
    
    return 0;
}

/*
 * Function: maximum
 * Parameters: x, y, z - three integers to compare
 * Returns: the largest of the three integers
 */
int maximum(int x, int y, int z) {
    int max = x;
    
    if (y > max) {
        max = y;
    }
    if (z > max) {
        max = z;
    }
    
    return max;
}