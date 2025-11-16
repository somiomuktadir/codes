#include <stdio.h>

int main(){
    int counter =0, sum= 0, grade;
    float average;
    scanf("%d", &grade);

    while (grade!=-1){
        sum=sum+grade;
        counter+=1;
        scanf("%d", &grade);
    }

    if (counter!=0){
        average=(float)sum/counter;
        printf("Average = %f\n", average);
    }
    else {
        printf("Meow!\n");
    }

    return 0;
}