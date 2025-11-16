#include <stdio.h>

int main(){
    int n, sum = 0, counter;
    scanf("%d", &n);
    for (counter = 1; counter <=n; counter=counter+2){
        sum=sum+counter;
    }
    printf("%d\n", sum);
    return 0;
}