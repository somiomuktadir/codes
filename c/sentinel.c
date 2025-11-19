/*
 * Grade Average Calculator (Sentinel-Controlled Loop)
 * Calculates the average of grades entered by the user
 * Uses -1 as the sentinel value to stop input
 */

#include <stdio.h>

int main() {
    int grade_count = 0;
    int total_sum = 0;
    int grade;
    float average;

    printf("=== Grade Average Calculator ===\n");
    printf("Enter grades one by one (enter -1 to finish):\n\n");

    // Get first grade
    printf("Enter grade: ");
    if (scanf("%d", &grade) != 1) {
        printf("Error: Invalid input.\n");
        return 1;
    }

    // Continue reading grades until sentinel value (-1) is entered
    while (grade != -1) {
        // Validate grade is non-negative (except sentinel)
        if (grade < 0) {
            printf("Warning: Negative grades are not valid. Skipping.\n");
        } else {
            total_sum += grade;
            grade_count++;
        }

        printf("Enter grade: ");
        if (scanf("%d", &grade) != 1) {
            printf("Error: Invalid input.\n");
            return 1;
        }
    }

    // Calculate and display results
    if (grade_count > 0) {
        average = (float)total_sum / grade_count;
        printf("\n=== Results ===\n");
        printf("Total grades entered: %d\n", grade_count);
        printf("Sum of grades: %d\n", total_sum);
        printf("Average grade: %.2f\n", average);
    } else {
        printf("\nNo valid grades were entered.\n");
    }

    return 0;
}