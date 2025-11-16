// Let's get the fibonacci numbers recursively

#include <stdio.h>

unsigned long long int fibonacci (unsigned int n);

int main(){
    unsigned long long int result;
    unsigned int number;
    scanf("%u", &number);

    result = fibonacci(number);
    printf("%u = %llu\n", number, result);
}

unsigned long long int fibonacci (unsigned int n){
    if (0==n ||  1 == n){
        return n;
    }
    else {
        return fibonacci (n-1)+fibonacci(n-2);
    }
} 

// The core difference between iteration and recursion is :
// Iteration - repetition structure; Recursion - selection structure