/*
 * File I/O Example
 * Demonstrates writing to a file and reading back its content.
 */

#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *fp;
    char filename[] = "example_file.txt";
    char buffer[100];

    // Writing to file
    fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file for writing.\n");
        return 1;
    }

    fprintf(fp, "Hello, File I/O!\n");
    fprintf(fp, "This is a test line.\n");
    fprintf(fp, "C programming is fun.\n");
    
    printf("Data written to %s successfully.\n", filename);
    fclose(fp);

    // Reading from file
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error opening file for reading.\n");
        return 1;
    }

    printf("\nReading content from %s:\n", filename);
    while (fgets(buffer, sizeof(buffer), fp) != NULL) {
        printf("%s", buffer);
    }

    fclose(fp);
    
    // Clean up (optional, keeping it for user to see)
    // remove(filename); 

    return 0;
}
