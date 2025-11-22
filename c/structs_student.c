/*
 * Structures Example
 * Demonstrates defining a Student structure, array of structures, and input/output.
 */

#include <stdio.h>

#define MAX_STUDENTS 50

struct Student {
    int id;
    char name[50];
    float marks;
};

int main() {
    struct Student students[MAX_STUDENTS];
    int n, i;

    printf("Enter number of students (1-%d): ", MAX_STUDENTS);
    if (scanf("%d", &n) != 1 || n <= 0 || n > MAX_STUDENTS) {
        printf("Invalid number.\n");
        return 1;
    }

    // Input
    for (i = 0; i < n; i++) {
        printf("\nStudent %d:\n", i + 1);
        printf("ID: ");
        scanf("%d", &students[i].id);
        getchar(); // Consume newline
        printf("Name: ");
        scanf("%[^\n]", students[i].name);
        printf("Marks: ");
        scanf("%f", &students[i].marks);
    }

    // Output
    printf("\n--- Student Details ---\n");
    printf("%-5s %-20s %-10s\n", "ID", "Name", "Marks");
    for (i = 0; i < n; i++) {
        printf("%-5d %-20s %-10.2f\n", students[i].id, students[i].name, students[i].marks);
    }

    return 0;
}
