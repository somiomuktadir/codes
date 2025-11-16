#include <stdio.h>

int main(){
    int counter=1;
    do {
        printf("%d\n", counter);
        counter=counter+1;
    }
    while (counter<=10);
    return 0;
}