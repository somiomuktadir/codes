/*
 * Student Pass/Fail Analyzer
 * Analyzes student results and determines if pass rate meets threshold
 * Pass threshold: 80% or higher
 */

#include <stdio.h>

int main() {
    int passes = 0, failures = 0, student = 1;
    int result, num_students;
    float pass_rate;
    const float PASS_THRESHOLD = 0.8; // 80%

    // Get number of students
    printf("Enter the number of students: ");
    if (scanf("%d", &num_students) != 1 || num_students <= 0) {
        printf("Error: Invalid number of students.\n");
        return 1;
    }

    printf("\nEnter results for each student (1 = pass, 0 = fail):\n");

    // Collect results for all students
    while (student <= num_students) {
        printf("Student %d: ", student);
        if (scanf("%d", &result) != 1) {
            printf("Error: Invalid input.\n");
            return 1;
        }

        if (result == 1) {
            passes++;
        } else if (result == 0) {
            failures++;
        } else {
            printf("Warning: Invalid result. Please enter 1 or 0. Treating as fail.\n");
            failures++;
        }
        student++;
    }

    // Calculate pass rate
    pass_rate = (float)passes / num_students;

    // Display results
    printf("\n=== Results Summary ===\n");
    printf("Total students: %d\n", num_students);
    printf("Passed: %d\n", passes);
    printf("Failed: %d\n", failures);
    printf("Pass rate: %.1f%%\n", pass_rate * 100);

    // Check if pass rate meets threshold
    if (pass_rate >= PASS_THRESHOLD) {
        printf("\nExcellent! The class meets the %.0f%% pass threshold.\n", PASS_THRESHOLD * 100);
    } else {
        printf("\nThe class does not meet the %.0f%% pass threshold.\n", PASS_THRESHOLD * 100);
    }

    return 0;
}