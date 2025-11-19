#include <stdio.h>

int main(void) {
    int passes = 0;
    int failures = 0;
    int student_counter = 1;
    int result;

    while (student_counter <= 10) {
        printf("Enter result (1 = pass, 2 = fail): ");
        
        // Check if scanf successfully read an integer
        if (scanf("%d", &result) != 1) {
            // Clear the input buffer to prevent infinite loop on invalid input
            while (getchar() != '\n');
            printf("Invalid input. Please enter a number.\n");
            continue; 
        }

        if (result == 1) {
            passes = passes + 1;
        } else {
            failures = failures + 1;
        }

        student_counter = student_counter + 1;
    }

    printf("Passed: %d\n", passes);
    printf("Failed: %d\n", failures);

    if (passes >= 9) {
        printf("Congo!\n");
    } else {
        printf("Meow!\n");
    }

    return 0;
}