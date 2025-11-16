// Lets write a program that determines the largest of the three numbers

#include <stdio.h>

int maximum(int x, int y, int z); // function prototype

int main(){
    int a,b,c;
    scanf("%d%d%d", &a, &b, &c);
    printf("%d\n", maximum(a,b,c));//calling the function
}

//now lets define the function
//x,y,z are parameters

int maximum(int x, int y, int z){
    int max = x;
    if (y>max){
        max=y;
    }
    if (z>max){
        max=z;
    }
    return max;
}