// This program calculates the average mark of a class test and comments on it.

#include <stdio.h>

#define CLASS_SIZE 10

int main() {
    int counter = 1;
    int sum = 0;
    int grade;
    float average;

    printf("Enter grades for %d students:\n", CLASS_SIZE);

    while (counter <= CLASS_SIZE) {
        printf("Student %d: ", counter);
        if (scanf("%d", &grade) != 1) {
            printf("Invalid input. Please enter a number.\n");
            return 1;
        }
        sum += grade;
        counter++;
    }

    average = (float)sum / CLASS_SIZE;
    printf("Average is %.2f\n", average);

    if (average >= 90) {
        printf("Excellent!\n");
    }
    else if (average >= 50) {
        printf("Average!\n");
    }
    else {
        printf("Meow!\n");
    }

    return 0;
}