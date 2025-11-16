#include <stdio.h>

int square (int y); // function prototype

int main(){
    int x; //counter

    for (x=1; x<=10; x=x+1){
        printf("%d\n", square(x)); //function call
    }
} // end main

int square (int y)// y is a copy of the argument to the function
{
    return y*y; // returns the square of y as an int
}