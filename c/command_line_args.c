/*
 * Command Line Arguments Example
 * Demonstrates argc and argv usage.
 * Run this program from terminal: ./command_line_args arg1 arg2
 */

#include <stdio.h>

int main(int argc, char *argv[]) {
    int i;

    printf("Number of arguments passed: %d\n", argc);
    printf("Program name: %s\n", argv[0]);

    if (argc > 1) {
        printf("\nArguments passed:\n");
        for (i = 1; i < argc; i++) {
            printf("argv[%d]: %s\n", i, argv[i]);
        }
    } else {
        printf("\nNo extra arguments passed.\n");
    }

    return 0;
}
