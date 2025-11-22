/*
 * Simple Calculator Program
 * Performs basic arithmetic operations: +, -, *, /
 * Author: Code Review Improvement
 */

#include <stdio.h>

int main() {
    double num1, num2, result;
    char operation;
    int valid_input;

    // Get first number
    do {
        printf("Enter first number: ");
        if (scanf("%lf", &num1) == 1) {
            valid_input = 1;
        } else {
            valid_input = 0;
            printf("Error: Invalid input. Please enter a number.\n");
            while (getchar() != '\n'); // Clear input buffer
        }
    } while (!valid_input);

    // Get second number
    do {
        printf("Enter second number: ");
        if (scanf("%lf", &num2) == 1) {
            valid_input = 1;
        } else {
            valid_input = 0;
            printf("Error: Invalid input. Please enter a number.\n");
            while (getchar() != '\n'); // Clear input buffer
        }
    } while (!valid_input);

    // Get operation
    do {
        printf("Enter operation (+, -, *, /): ");
        if (scanf(" %c", &operation) == 1) {
            if (operation == '+' || operation == '-' || operation == '*' || operation == '/') {
                valid_input = 1;
            } else {
                valid_input = 0;
                printf("Error: Invalid operator. Use +, -, *, or /.\n");
            }
        } else {
            valid_input = 0;
            while (getchar() != '\n'); // Clear input buffer
        }
    } while (!valid_input);

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
            }
            break;
    }
    
    return 0;
}