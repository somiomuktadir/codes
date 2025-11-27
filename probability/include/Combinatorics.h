#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include <vector>
#include <string>

namespace Combinatorics {
    // Basic
    unsigned long long factorial(int n);
    unsigned long long permutations(int n, int r);
    unsigned long long combinations(int n, int r);

    // Advanced
    unsigned long long permutationsWithRepetition(int n, const std::vector<int>& repetitions);
    unsigned long long circularPermutations(int n);
    unsigned long long derangements(int n);
    unsigned long long catalanNumber(int n);
    unsigned long long stirlingNumberSecondKind(int n, int k);
    unsigned long long bellNumber(int n);
    
    // Extended Features
    unsigned long long fibonacci(int n);
    unsigned long long lucasNumber(int n);
    unsigned long long integerPartitions(int n);
    unsigned long long multinomial(int n, const std::vector<int>& k_counts);

    // Explanations
    void printFormulas();
}

#endif // COMBINATORICS_H
