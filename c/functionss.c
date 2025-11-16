// Stack is an important data structure in C programming
// Stacks are known as LIFO (last in, first out) structures
// A stack frame contains all the return address

// The following program demonstrates the function call stack
// And stack frame uses a function square

// This program is simple

#include <stdio.h>

int square (int); // prototype for function square

int main(){
    int a = 10; // value to square
    printf("%d squared is %d\n", a, square(a)); // display a squared
} // end main

int square(int x) // x is a local variable
{
    return x * x;
} //end fucntion square