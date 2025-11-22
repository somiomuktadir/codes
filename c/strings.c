/*
 * Strings Example
 * Demonstrates string manipulation (length, copy, concat) without string.h library functions.
 */

#include <stdio.h>

int stringLength(char str[]) {
    int len = 0;
    while (str[len] != '\0') {
        len++;
    }
    return len;
}

void stringCopy(char source[], char dest[]) {
    int i = 0;
    while (source[i] != '\0') {
        dest[i] = source[i];
        i++;
    }
    dest[i] = '\0';
}

void stringConcat(char str1[], char str2[]) {
    int i = 0, j = 0;
    
    // Move to end of str1
    while (str1[i] != '\0') {
        i++;
    }

    // Copy str2 to end of str1
    while (str2[j] != '\0') {
        str1[i] = str2[j];
        i++;
        j++;
    }
    str1[i] = '\0';
}

int main() {
    char str1[100], str2[100], copy[100];

    printf("Enter first string: ");
    scanf("%[^\n]", str1); // Read line including spaces
    getchar(); // Consume newline

    printf("Enter second string: ");
    scanf("%[^\n]", str2);

    printf("\nLength of String 1: %d\n", stringLength(str1));
    printf("Length of String 2: %d\n", stringLength(str2));

    stringCopy(str1, copy);
    printf("Copy of String 1: %s\n", copy);

    stringConcat(str1, str2);
    printf("Concatenation of String 1 and 2: %s\n", str1);

    return 0;
}
