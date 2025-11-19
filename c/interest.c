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
    printf("Enter principal amount: $");
    if (scanf("%f", &principal) != 1 || principal <= 0) {
        printf("Error: Invalid principal amount.\n");
        return 1;
    }

    // Get interest rate
    printf("Enter annual interest rate (as decimal, e.g., 0.05 for 5%%): ");
    if (scanf("%f", &rate) != 1 || rate < 0 || rate > 1) {
        printf("Error: Invalid interest rate.\n");
        return 1;
    }

    // Get number of years
    printf("Enter number of years: ");
    if (scanf("%d", &num_years) != 1 || num_years <= 0) {
        printf("Error: Invalid number of years.\n");
        return 1;
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