// Let us recursively calculate factorials

#include <stdio.h>

unsigned long long int factorial (unsigned int number);

int main(){
    unsigned int i; //counter
    // during each iteration, calculate
    // factorial(i) and display result

    for (i=0; i <=20; i=i+1){
        printf("%d ! = %llu\n", i, factorial(i));
    }
}// end main

unsigned long long int factorial (unsigned int number){
    if (number <= 1){
        return 1;
    }
    else{
        return (number * factorial(number-1));
    }
}