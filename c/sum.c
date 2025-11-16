#include <stdio.h>

int main(){
    int sum=0, count;
    for (count = 2; count <=100; count=count+2){
        sum = sum+count;
    }
    printf("%d\n", sum);
    return 0;
}