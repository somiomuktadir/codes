#include <stdio.h>

int main() {
    int n, i;
    long long int first = 0, second = 1, next;
    scanf("%d", &n);
    for (i = 0; i < n; i=i+1) {
        if (i <= 1)
            next = i;
        else {
            next = first + second;
            first = second;
            second = next;
        }
        printf("%lld\n", next);
    }
    return 0;
}
