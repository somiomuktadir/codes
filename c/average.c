// This program calculates the average mark of a class test and comments on it.

#include <stdio.h>

int main(){
    int counter=1, sum=0, grade, average;

    while (counter <= 10){
        scanf("%d", &grade);
        sum = sum+grade;
        counter+=1;
    }
    average=sum/10;
    printf("Average is %d\n", average);

    if (average>=90){
        printf("Excellent!\n");
    }
    else if (average >= 50 && average <90){
        printf("Average!\n");
    }
    else {
        printf("Meow!\n");
    }
    return 0;
}