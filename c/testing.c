#include <stdio.h>

int main(){
    int passes=0, failures=0, student=1;
    int result,n;
    scanf("%d", &n);

    while (student <=n){
        scanf("%d", &result);
        if (result == 1){
            passes=passes+1;
        }
        else {
            failures=failures+1;
        }
        student=student+1;
    }
    if (passes >= .8*n){
        printf("Meow!\n");
    }
    return 0;
}