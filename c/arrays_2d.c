/*
 * Arrays 2D Example
 * Demonstrates matrix addition and display.
 */

#include <stdio.h>

#define ROWS 3
#define COLS 3

void inputMatrix(int matrix[ROWS][COLS], char *name) {
    int i, j;
    printf("Enter elements for Matrix %s (%dx%d):\n", name, ROWS, COLS);
    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            printf("Enter element [%d][%d]: ", i, j);
            scanf("%d", &matrix[i][j]);
        }
    }
}

void printMatrix(int matrix[ROWS][COLS]) {
    int i, j;
    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            printf("%4d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    int mat1[ROWS][COLS], mat2[ROWS][COLS], sum[ROWS][COLS];
    int i, j;

    inputMatrix(mat1, "A");
    inputMatrix(mat2, "B");

    // Addition
    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            sum[i][j] = mat1[i][j] + mat2[i][j];
        }
    }

    printf("\nMatrix A:\n");
    printMatrix(mat1);

    printf("\nMatrix B:\n");
    printMatrix(mat2);

    printf("\nSum of Matrix A and B:\n");
    printMatrix(sum);

    return 0;
}
