/*
 * Simple Calculator Program
 * Performs basic arithmetic operations: +, -, *, /
 * Author: Code Review Improvement
 */

#include <stdio.h>

int main() {
    double num1, num2, result;
    char operation;

    // Get first number
    printf("Enter first number: ");
    if (scanf("%lf", &num1) != 1) {
        printf("Error: Invalid input for first number.\n");
        return 1;
    }

    // Get second number
    printf("Enter second number: ");
    if (scanf("%lf", &num2) != 1) {
        printf("Error: Invalid input for second number.\n");
        return 1;
    }

    // Get operation
    printf("Enter operation (+, -, *, /): ");
    if (scanf(" %c", &operation) != 1) {
        printf("Error: Invalid input for operation.\n");
        return 1;
    }

    // Perform calculation based on operation
    switch (operation) {
        case '+':
            result = num1 + num2;
            printf("%.2lf + %.2lf = %.2lf\n", num1, num2, result);
            break;

        case '-':
            result = num1 - num2;
            printf("%.2lf - %.2lf = %.2lf\n", num1, num2, result);
            break;

        case '*':
            result = num1 * num2;
            printf("%.2lf * %.2lf = %.2lf\n", num1, num2, result);
            break;

        case '/':
            if (num2 != 0) {
                result = num1 / num2;
                printf("%.2lf / %.2lf = %.2lf\n", num1, num2, result);
            } else {
                printf("Error: Division by zero is not allowed.\n");
                return 1;
            }
            break;

        default:
            printf("Error: Invalid operator '%c'. Use +, -, *, or /\n", operation);
            return 1;
    }
    
    return 0;
}