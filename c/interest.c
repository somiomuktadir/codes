/*
 * Compound Interest Calculator
 * Calculates compound interest over a period of years
 * Formula: A = P(1 + r)^t
 */

#include <stdio.h>
#include <math.h>

int main() {
    float amount, principal, rate;
    int year, num_years;

    // Get principal amount
    // Get principal amount
    while (1) {
        printf("Enter principal amount: $");
        if (scanf("%f", &principal) == 1 && principal > 0) {
            break;
        } else {
            printf("Error: Invalid principal amount. Please enter a positive number.\n");
            while (getchar() != '\n'); // Clear input buffer
        }
    }

    // Get interest rate
    while (1) {
        printf("Enter annual interest rate (as decimal, e.g., 0.05 for 5%%): ");
        if (scanf("%f", &rate) == 1 && rate >= 0 && rate <= 1) {
            break;
        } else {
            printf("Error: Invalid interest rate. Please enter a value between 0 and 1.\n");
            while (getchar() != '\n'); // Clear input buffer
        }
    }

    // Get number of years
    while (1) {
        printf("Enter number of years: ");
        if (scanf("%d", &num_years) == 1 && num_years > 0) {
            break;
        } else {
            printf("Error: Invalid number of years. Please enter a positive integer.\n");
            while (getchar() != '\n'); // Clear input buffer
        }
    }

    // Display results in table format
    printf("\n=== Compound Interest Calculation ===\n");
    printf("Principal: $%.2f\n", principal);
    printf("Annual Rate: %.2f%%\n\n", rate * 100);
    printf("%-10s %-15s\n", "Year", "Amount");
    printf("----------------------------------------\n");

    for (year = 1; year <= num_years; year++) {
        amount = principal * pow(1.0 + rate, year);
        printf("%-10d $%-14.2f\n", year, amount);
    }

    return 0;
}