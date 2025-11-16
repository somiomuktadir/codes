#include <stdio.h>

int main(){
    int p=0,f=0,count=1;
    int result;

    while (count <= 10){
        scanf("%d", &result);
        if (result == 1){
            p=p+1;
        }
        else{
            f=f+1;
        }
        count=count+1;
    }
    printf("%d\n",p);
    printf("%d\n",f);

    if (p>=9){
        printf("Congo!\n");
    }
    else{
        printf("Meow!\n");
    }
    return 0;
}