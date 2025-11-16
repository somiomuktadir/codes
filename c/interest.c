#include <stdio.h>
#include <math.h>

int main(){
    float amount,principal=1000,rate=0.05;
    int year;

    for (year=1; year<=10; year=year+1){
        amount = principal * pow (1.0+ rate, year);
        printf("%d    %f\n", year, amount);
    }
    return 0;
}