/*
 * Student Grade Tracker
 * Tracks pass/fail results for 10 students
 * Displays summary and congratulatory message if 9+ students pass
 */

#include <stdio.h>

int main(void) {
    int passes = 0;
    int failures = 0;
    int student_counter = 1;
    int result;
    const int TOTAL_STUDENTS = 10;
    const int PASS_THRESHOLD = 9;

    printf("=== Student Grade Tracker ===\n");
    printf("Enter results for %d students\n\n", TOTAL_STUDENTS);

    while (student_counter <= TOTAL_STUDENTS) {
        printf("Student %d - Enter result (1 = pass, 2 = fail): ", student_counter);
        
        // Check if scanf successfully read an integer
        if (scanf("%d", &result) != 1) {
            // Clear the input buffer to prevent infinite loop on invalid input
            while (getchar() != '\n');
            printf("Invalid input. Please enter 1 or 2.\n");
            continue; 
        }

        // Validate input is either 1 or 2
        if (result == 1) {
            passes++;
        } else if (result == 2) {
            failures++;
        } else {
            printf("Invalid result. Please enter 1 (pass) or 2 (fail).\n");
            continue;
        }

        student_counter++;
    }

    // Display results summary
    printf("\n=== Results Summary ===\n");
    printf("Passed: %d\n", passes);
    printf("Failed: %d\n", failures);
    printf("Pass Rate: %.1f%%\n", (passes * 100.0) / TOTAL_STUDENTS);

    // Display performance message
    if (passes >= PASS_THRESHOLD) {
        printf("\nCongratulations! Excellent class performance!\n");
    } else {
        printf("\nClass needs improvement. More support may be needed.\n");
    }

    return 0;
}