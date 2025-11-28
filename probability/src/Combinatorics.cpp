#include "Combinatorics.h"
#include "Utils.h"
#include <iostream>
#include <numeric>
#include <cmath>

namespace Combinatorics {

    unsigned long long factorial(int n) {
        if (n < 0) throw std::invalid_argument("Factorial undefined for negative numbers");
        Utils::log("Calculating Factorial of " + std::to_string(n));
        unsigned long long result = 1;
        std::string step = "";
        for (int i = 1; i <= n; ++i) {
            result *= i;
            if (Utils::isVerbose()) {
                if (i == 1) {
                    step = std::to_string(i);
                } else {
                    step += " × " + std::to_string(i);
                }
            }
        }
        if (Utils::isVerbose() && n > 0) {
            Utils::logStep(std::to_string(n) + "! = " + step + " = " + std::to_string(result));
        } else if (Utils::isVerbose()) {
            Utils::logStep("0! = 1");
        }
        return result;
    }

    unsigned long long permutations(int n, int r) {
        if (n < 0 || r < 0 || r > n) throw std::invalid_argument("Invalid n or r for permutations");
        Utils::log("Calculating Permutations nPr: n=" + std::to_string(n) + ", r=" + std::to_string(r));
        Utils::logStep("Formula: n! / (n-r)!");
        Utils::logStep("Substitution: " + std::to_string(n) + "! / (" + std::to_string(n) + "-" + std::to_string(r) + ")!");
        unsigned long long result = Utils::nPr(n, r);
        Utils::logStep("Result: " + std::to_string(result));
        return result;
    }

    unsigned long long combinations(int n, int r) {
        if (n < 0 || r < 0 || r > n) throw std::invalid_argument("Invalid n or r for combinations");
        Utils::log("Calculating Combinations nCr: n=" + std::to_string(n) + ", r=" + std::to_string(r));
        Utils::logStep("Formula: n! / (r! × (n-r)!)");
        Utils::logStep("Substitution: " + std::to_string(n) + "! / (" + std::to_string(r) + "! × " + std::to_string(n-r) + "!)");
        unsigned long long result = Utils::nCr(n, r);
        Utils::logStep("Result: " + std::to_string(result));
        return result;
    }

    unsigned long long permutationsWithRepetition(int n, const std::vector<int>& repetitions) {
        unsigned long long denom = 1;
        int sum_r = 0;
        for (int r : repetitions) {
            denom *= Utils::factorial(r);
            sum_r += r;
        }
        if (sum_r != n) {
            // If sum of repetitions doesn't equal n, it assumes the rest are 1s (unique)
            // or it could be an error depending on interpretation. 
            // Standard formula is n! / (n1! * n2! * ... * nk!) where sum(ni) = n.
             if (sum_r > n) throw std::invalid_argument("Sum of repetitions cannot exceed n.");
        }
        return Utils::factorial(n) / denom;
    }

    unsigned long long circularPermutations(int n) {
        if (n <= 0) return 0;
        return Utils::factorial(n - 1);
    }

    unsigned long long derangements(int n) {
        if (n == 0) return 1;
        if (n == 1) return 0;
        // !n = (n-1) * (!(n-1) + !(n-2))
        unsigned long long d_prev2 = 1; // !0
        unsigned long long d_prev1 = 0; // !1
        unsigned long long d_curr = 0;
        
        for (int i = 2; i <= n; ++i) {
            d_curr = (i - 1) * (d_prev1 + d_prev2);
            d_prev2 = d_prev1;
            d_prev1 = d_curr;
        }
        return d_curr;
    }

    unsigned long long catalanNumber(int n) {
        // C_n = (2n)! / ((n+1)!n!)
        return Utils::nCr(2 * n, n) / (n + 1);
    }

    unsigned long long stirlingNumberSecondKind(int n, int k) {
        // S(n, k) = k*S(n-1, k) + S(n-1, k-1)
        if (k < 0 || k > n) return 0;
        if (k == 0) return (n == 0) ? 1 : 0;
        if (k == n) return 1;
        if (k == 1) return 1;

        // Using DP for larger values if needed, but for reasonable n recursion is fine or iterative
        std::vector<std::vector<unsigned long long>> dp(n + 1, std::vector<unsigned long long>(k + 1, 0));
        dp[0][0] = 1;

        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= k; ++j) {
                dp[i][j] = j * dp[i - 1][j] + dp[i - 1][j - 1];
            }
        }
        return dp[n][k];
    }

    unsigned long long bellNumber(int n) {
        // Sum of Stirling numbers of second kind S(n, k) for k=0 to n
        unsigned long long sum = 0;
        for (int k = 0; k <= n; ++k) {
            sum += stirlingNumberSecondKind(n, k);
        }
        return sum;
    }

    unsigned long long fibonacci(int n) {
        if (n < 0) throw std::invalid_argument("Fibonacci undefined for negative indices.");
        if (n == 0) return 0;
        if (n == 1) return 1;
        
        unsigned long long a = 0, b = 1, c;
        for (int i = 2; i <= n; ++i) {
            c = a + b;
            a = b;
            b = c;
        }
        return b;
    }

    unsigned long long lucasNumber(int n) {
        if (n < 0) throw std::invalid_argument("Lucas number undefined for negative indices.");
        if (n == 0) return 2;
        if (n == 1) return 1;
        
        unsigned long long a = 2, b = 1, c;
        for (int i = 2; i <= n; ++i) {
            c = a + b;
            a = b;
            b = c;
        }
        return b;
    }

    unsigned long long integerPartitions(int n) {
        if (n < 0) return 0;
        if (n == 0) return 1;
        
        // p(n) using Euler's pentagonal number theorem
        // p(n) = sum_{k!=0} (-1)^(k-1) * p(n - k(3k-1)/2)
        std::vector<unsigned long long> p(n + 1, 0);
        p[0] = 1;
        
        for (int i = 1; i <= n; ++i) {
            for (int k = 1; ; ++k) {
                int pent1 = k * (3 * k - 1) / 2;
                int pent2 = k * (3 * k + 1) / 2;
                
                if (pent1 > i) break;
                
                if (k % 2 == 1) p[i] += p[i - pent1];
                else p[i] -= p[i - pent1];
                
                if (pent2 > i) continue;
                
                if (k % 2 == 1) p[i] += p[i - pent2];
                else p[i] -= p[i - pent2];
            }
        }
        return p[n];
    }

    unsigned long long multinomial(int n, const std::vector<int>& k_counts) {
        int sum_k = 0;
        unsigned long long denom = 1;
        for (int k : k_counts) {
            if (k < 0) throw std::invalid_argument("Counts must be non-negative.");
            sum_k += k;
            denom *= Utils::factorial(k);
        }
        if (sum_k != n) throw std::invalid_argument("Sum of counts must equal n.");
        
        return Utils::factorial(n) / denom;
    }

    void printFormulas() {
        std::cout << "\n--- Combinatorics Formulas ---\n";
        std::cout << "Factorial (n!): n * (n-1) * ... * 1\n";
        std::cout << "Permutations (nPr): n! / (n-r)!\n";
        std::cout << "Combinations (nCr): n! / (r! * (n-r)!)\n";
        std::cout << "Permutations with Repetition: n! / (n1! * n2! * ...)\n";
        std::cout << "Circular Permutations: (n-1)!\n";
        std::cout << "Derangements (!n): n! * sum((-1)^k / k!) for k=0 to n\n";
        std::cout << "Catalan Numbers (Cn): (2n)! / ((n+1)!n!)\n";
        std::cout << "Stirling Numbers 2nd Kind S(n,k): Ways to partition set of n elements into k non-empty subsets\n";
        std::cout << "Bell Numbers (Bn): Sum of S(n,k) for k=0 to n\n";
        std::cout << "Fibonacci (Fn): F(n) = F(n-1) + F(n-2)\n";
        std::cout << "Lucas (Ln): L(n) = L(n-1) + L(n-2)\n";
        std::cout << "Integer Partitions p(n): Ways to write n as sum of positive integers\n";
        std::cout << "Multinomial: n! / (k1! * k2! * ... * km!)\n";
        std::cout << "------------------------------\n";
    }
}
