#include "Menu.h"
#include "Utils.h"
#include "Combinatorics.h"
#include "Probability.h"
#include "Distributions.h"
#include <iostream>
#include <vector>
#include <cstdlib>

namespace Menu {

    void clearScreen() {
        #ifdef _WIN32
            system("cls");
        #else
            int ret = system("clear");
            (void)ret;
        #endif
    }

    void pause() {
        std::cout << "\nPress Enter to continue...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cin.get();
    }

    void showCombinatoricsMenu() {
        while (true) {
            clearScreen();
            std::cout << "=== Combinatorics Menu ===\n";
            std::cout << "1. Factorial (n!)\n";
            std::cout << "2. Permutations (nPr)\n";
            std::cout << "3. Combinations (nCr)\n";
            std::cout << "4. Permutations with Repetition\n";
            std::cout << "5. Circular Permutations\n";
            std::cout << "6. Derangements (!n)\n";
            std::cout << "7. Catalan Numbers\n";
            std::cout << "8. Stirling Numbers (2nd Kind)\n";
            std::cout << "9. Bell Numbers\n";
            std::cout << "10. Fibonacci Numbers\n";
            std::cout << "11. Lucas Numbers\n";
            std::cout << "12. Integer Partitions\n";
            std::cout << "13. Multinomial Coefficient\n";
            std::cout << "14. View Formulas\n";
            std::cout << "0. Back to Main Menu\n";
            
            int choice = Utils::getIntInput("Select option: ");

            try {
                if (choice == 0) break;
                else if (choice == 1) {
                    int n = Utils::getIntInput("Enter n: ");
                    std::cout << "Result: " << Combinatorics::factorial(n) << std::endl;
                }
                else if (choice == 2) {
                    int n = Utils::getIntInput("Enter n: ");
                    int r = Utils::getIntInput("Enter r: ");
                    std::cout << "Result: " << Combinatorics::permutations(n, r) << std::endl;
                }
                else if (choice == 3) {
                    int n = Utils::getIntInput("Enter n: ");
                    int r = Utils::getIntInput("Enter r: ");
                    std::cout << "Result: " << Combinatorics::combinations(n, r) << std::endl;
                }
                else if (choice == 4) {
                    int n = Utils::getIntInput("Enter total items n: ");
                    std::vector<double> reps = Utils::getVectorInput("Enter repetitions for each type (space separated): ");
                    std::vector<int> intReps;
                    for(double d : reps) intReps.push_back((int)d);
                    std::cout << "Result: " << Combinatorics::permutationsWithRepetition(n, intReps) << std::endl;
                }
                else if (choice == 5) {
                    int n = Utils::getIntInput("Enter n: ");
                    std::cout << "Result: " << Combinatorics::circularPermutations(n) << std::endl;
                }
                else if (choice == 6) {
                    int n = Utils::getIntInput("Enter n: ");
                    std::cout << "Result: " << Combinatorics::derangements(n) << std::endl;
                }
                else if (choice == 7) {
                    int n = Utils::getIntInput("Enter n: ");
                    std::cout << "Result: " << Combinatorics::catalanNumber(n) << std::endl;
                }
                else if (choice == 8) {
                    int n = Utils::getIntInput("Enter n: ");
                    int k = Utils::getIntInput("Enter k: ");
                    std::cout << "Result: " << Combinatorics::stirlingNumberSecondKind(n, k) << std::endl;
                }
                else if (choice == 9) {
                    int n = Utils::getIntInput("Enter n: ");
                    std::cout << "Result: " << Combinatorics::bellNumber(n) << std::endl;
                }
                else if (choice == 10) {
                    int n = Utils::getIntInput("Enter n: ");
                    std::cout << "Result: " << Combinatorics::fibonacci(n) << std::endl;
                }
                else if (choice == 11) {
                    int n = Utils::getIntInput("Enter n: ");
                    std::cout << "Result: " << Combinatorics::lucasNumber(n) << std::endl;
                }
                else if (choice == 12) {
                    int n = Utils::getIntInput("Enter n: ");
                    std::cout << "Result: " << Combinatorics::integerPartitions(n) << std::endl;
                }
                else if (choice == 13) {
                    int n = Utils::getIntInput("Enter total n: ");
                    std::vector<double> ks = Utils::getVectorInput("Enter k counts (space separated): ");
                    std::vector<int> intKs;
                    for(double d : ks) intKs.push_back((int)d);
                    std::cout << "Result: " << Combinatorics::multinomial(n, intKs) << std::endl;
                }
                else if (choice == 14) {
                    Combinatorics::printFormulas();
                }
                else {
                    std::cout << "Invalid option.\n";
                }
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
            pause();
        }
    }

    void showProbabilityMenu() {
        while (true) {
            clearScreen();
            std::cout << "=== Probability Menu ===\n";
            std::cout << "1. Basic Probability\n";
            std::cout << "2. Complementary Probability\n";
            std::cout << "3. Conditional Probability P(A|B)\n";
            std::cout << "4. Bayes' Theorem\n";
            std::cout << "5. Union Probability P(A U B)\n";
            std::cout << "6. Intersection (Independent Events)\n";
            std::cout << "7. Expected Value (Discrete RV)\n";
            std::cout << "8. Variance (Discrete RV)\n";
            std::cout << "9. View Formulas\n";
            std::cout << "0. Back to Main Menu\n";

            int choice = Utils::getIntInput("Select option: ");

            try {
                if (choice == 0) break;
                else if (choice == 1) {
                    double f = Utils::getDoubleInput("Enter favorable outcomes: ");
                    double t = Utils::getDoubleInput("Enter total outcomes: ");
                    std::cout << "Result: " << Probability::probability(f, t) << std::endl;
                }
                else if (choice == 2) {
                    double p = Utils::getDoubleInput("Enter probability P(E): ");
                    std::cout << "Result: " << Probability::complementaryProbability(p) << std::endl;
                }
                else if (choice == 3) {
                    double pi = Utils::getDoubleInput("Enter P(A n B): ");
                    double pg = Utils::getDoubleInput("Enter P(B) (given): ");
                    std::cout << "Result: " << Probability::conditionalProbability(pi, pg) << std::endl;
                }
                else if (choice == 4) {
                    double pba = Utils::getDoubleInput("Enter P(B|A): ");
                    double pa = Utils::getDoubleInput("Enter P(A): ");
                    double pb = Utils::getDoubleInput("Enter P(B): ");
                    std::cout << "Result: " << Probability::bayesTheorem(pba, pa, pb) << std::endl;
                }
                else if (choice == 5) {
                    double pa = Utils::getDoubleInput("Enter P(A): ");
                    double pb = Utils::getDoubleInput("Enter P(B): ");
                    double pi = Utils::getDoubleInput("Enter P(A n B): ");
                    std::cout << "Result: " << Probability::unionProbability(pa, pb, pi) << std::endl;
                }
                else if (choice == 6) {
                    double pa = Utils::getDoubleInput("Enter P(A): ");
                    double pb = Utils::getDoubleInput("Enter P(B): ");
                    std::cout << "Result: " << Probability::intersectionIndependent(pa, pb) << std::endl;
                }
                else if (choice == 7) {
                    std::vector<double> vals = Utils::getVectorInput("Enter values (x): ");
                    std::vector<double> probs = Utils::getVectorInput("Enter probabilities (p): ");
                    std::cout << "Expected Value: " << Probability::expectedValue(vals, probs) << std::endl;
                }
                else if (choice == 8) {
                    std::vector<double> vals = Utils::getVectorInput("Enter values (x): ");
                    std::vector<double> probs = Utils::getVectorInput("Enter probabilities (p): ");
                    std::cout << "Variance: " << Probability::variance(vals, probs) << std::endl;
                }
                else if (choice == 9) {
                    Probability::printFormulas();
                }
                else {
                    std::cout << "Invalid option.\n";
                }
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
            pause();
        }
    }

    void showDistributionsMenu() {
        while (true) {
            clearScreen();
            std::cout << "=== Distributions Menu ===\n";
            std::cout << "--- Discrete ---\n";
            std::cout << "1. Binomial Distribution\n";
            std::cout << "2. Poisson Distribution\n";
            std::cout << "3. Geometric Distribution\n";
            std::cout << "4. Hypergeometric Distribution\n";
            std::cout << "5. Negative Binomial Distribution\n";
            std::cout << "--- Continuous ---\n";
            std::cout << "6. Uniform Distribution\n";
            std::cout << "7. Normal (Gaussian) Distribution\n";
            std::cout << "8. Exponential Distribution\n";
            std::cout << "9. Gamma Distribution\n";
            std::cout << "10. Beta Distribution\n";
            std::cout << "11. Chi-Square Distribution\n";
            std::cout << "12. Student's t-Distribution\n";
            std::cout << "13. F-Distribution\n";
            std::cout << "--- Info ---\n";
            std::cout << "14. View Formulas\n";
            std::cout << "0. Back to Main Menu\n";

            int choice = Utils::getIntInput("Select option: ");

            try {
                if (choice == 0) break;
                else if (choice == 1) {
                    int n = Utils::getIntInput("Enter trials (n): ");
                    double p = Utils::getDoubleInput("Enter probability (p): ");
                    int k = Utils::getIntInput("Enter successes (k): ");
                    std::cout << "PMF P(X=k): " << Distributions::Discrete::binomialPMF(n, p, k) << std::endl;
                    std::cout << "CDF P(X<=k): " << Distributions::Discrete::binomialCDF(n, p, k) << std::endl;
                    std::cout << "Mean: " << Distributions::Discrete::binomialMean(n, p) << std::endl;
                    std::cout << "Variance: " << Distributions::Discrete::binomialVariance(n, p) << std::endl;
                }
                else if (choice == 2) {
                    double l = Utils::getDoubleInput("Enter lambda: ");
                    int k = Utils::getIntInput("Enter occurrences (k): ");
                    std::cout << "PMF P(X=k): " << Distributions::Discrete::poissonPMF(l, k) << std::endl;
                    std::cout << "CDF P(X<=k): " << Distributions::Discrete::poissonCDF(l, k) << std::endl;
                    std::cout << "Mean: " << Distributions::Discrete::poissonMean(l) << std::endl;
                }
                else if (choice == 3) {
                    double p = Utils::getDoubleInput("Enter probability (p): ");
                    int k = Utils::getIntInput("Enter trial number (k): ");
                    std::cout << "PMF P(X=k): " << Distributions::Discrete::geometricPMF(p, k) << std::endl;
                    std::cout << "CDF P(X<=k): " << Distributions::Discrete::geometricCDF(p, k) << std::endl;
                    std::cout << "Mean: " << Distributions::Discrete::geometricMean(p) << std::endl;
                }
                else if (choice == 4) {
                    int N = Utils::getIntInput("Enter population size (N): ");
                    int K = Utils::getIntInput("Enter successes in population (K): ");
                    int n = Utils::getIntInput("Enter sample size (n): ");
                    int k = Utils::getIntInput("Enter successes in sample (k): ");
                    std::cout << "PMF P(X=k): " << Distributions::Discrete::hypergeometricPMF(N, K, n, k) << std::endl;
                    std::cout << "Mean: " << Distributions::Discrete::hypergeometricMean(N, K, n) << std::endl;
                }
                else if (choice == 5) {
                    int r = Utils::getIntInput("Enter successes (r): ");
                    double p = Utils::getDoubleInput("Enter probability (p): ");
                    int k = Utils::getIntInput("Enter failures (k): ");
                    std::cout << "PMF P(X=k): " << Distributions::Discrete::negativeBinomialPMF(r, p, k) << std::endl;
                    std::cout << "Mean: " << Distributions::Discrete::negativeBinomialMean(r, p) << std::endl;
                }
                else if (choice == 6) {
                    double a = Utils::getDoubleInput("Enter min (a): ");
                    double b = Utils::getDoubleInput("Enter max (b): ");
                    double x = Utils::getDoubleInput("Enter value (x): ");
                    std::cout << "PDF f(x): " << Distributions::Continuous::uniformPDF(a, b, x) << std::endl;
                    std::cout << "CDF F(x): " << Distributions::Continuous::uniformCDF(a, b, x) << std::endl;
                    std::cout << "Mean: " << Distributions::Continuous::uniformMean(a, b) << std::endl;
                }
                else if (choice == 7) {
                    double m = Utils::getDoubleInput("Enter mean: ");
                    double s = Utils::getDoubleInput("Enter std dev: ");
                    double x = Utils::getDoubleInput("Enter value (x): ");
                    std::cout << "PDF f(x): " << Distributions::Continuous::normalPDF(m, s, x) << std::endl;
                    std::cout << "CDF F(x): " << Distributions::Continuous::normalCDF(m, s, x) << std::endl;
                }
                else if (choice == 8) {
                    double l = Utils::getDoubleInput("Enter lambda: ");
                    double x = Utils::getDoubleInput("Enter value (x): ");
                    std::cout << "PDF f(x): " << Distributions::Continuous::exponentialPDF(l, x) << std::endl;
                    std::cout << "CDF F(x): " << Distributions::Continuous::exponentialCDF(l, x) << std::endl;
                    std::cout << "Mean: " << Distributions::Continuous::exponentialMean(l) << std::endl;
                }
                else if (choice == 9) {
                    double k = Utils::getDoubleInput("Enter shape (k): ");
                    double theta = Utils::getDoubleInput("Enter scale (theta): ");
                    double x = Utils::getDoubleInput("Enter value (x): ");
                    std::cout << "PDF f(x): " << Distributions::Continuous::gammaPDF(k, theta, x) << std::endl;
                    std::cout << "Mean: " << Distributions::Continuous::gammaMean(k, theta) << std::endl;
                }
                else if (choice == 10) {
                    double alpha = Utils::getDoubleInput("Enter alpha: ");
                    double beta = Utils::getDoubleInput("Enter beta: ");
                    double x = Utils::getDoubleInput("Enter value (x): ");
                    std::cout << "PDF f(x): " << Distributions::Continuous::betaPDF(alpha, beta, x) << std::endl;
                    std::cout << "Mean: " << Distributions::Continuous::betaMean(alpha, beta) << std::endl;
                }
                else if (choice == 11) {
                    int k = Utils::getIntInput("Enter degrees of freedom (k): ");
                    double x = Utils::getDoubleInput("Enter value (x): ");
                    std::cout << "PDF f(x): " << Distributions::Continuous::chiSquarePDF(k, x) << std::endl;
                    std::cout << "Mean: " << Distributions::Continuous::chiSquareMean(k) << std::endl;
                }
                else if (choice == 12) {
                    int v = Utils::getIntInput("Enter degrees of freedom (v): ");
                    double x = Utils::getDoubleInput("Enter value (x): ");
                    std::cout << "PDF f(x): " << Distributions::Continuous::studentTPDF(v, x) << std::endl;
                    std::cout << "Mean: " << Distributions::Continuous::studentTMean(v) << std::endl;
                }
                else if (choice == 13) {
                    int d1 = Utils::getIntInput("Enter d1: ");
                    int d2 = Utils::getIntInput("Enter d2: ");
                    double x = Utils::getDoubleInput("Enter value (x): ");
                    std::cout << "PDF f(x): " << Distributions::Continuous::fDistributionPDF(d1, d2, x) << std::endl;
                    std::cout << "Mean: " << Distributions::Continuous::fDistributionMean(d1, d2) << std::endl;
                }
                else if (choice == 14) {
                    Distributions::printFormulas();
                }
                else {
                    std::cout << "Invalid option.\n";
                }
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
            pause();
        }
    }

    void showMainMenu() {
        while (true) {
            clearScreen();
            std::cout << "==========================================\n";
            std::cout << "   Combinatorics & Probability Calculator \n";
            std::cout << "==========================================\n";
            std::cout << "1. Combinatorics\n";
            std::cout << "2. Probability\n";
            std::cout << "3. Distributions\n";
            std::cout << "0. Exit\n";
            
            int choice = Utils::getIntInput("Select module: ");

            if (choice == 0) {
                std::cout << "Exiting...\n";
                break;
            }
            else if (choice == 1) showCombinatoricsMenu();
            else if (choice == 2) showProbabilityMenu();
            else if (choice == 3) showDistributionsMenu();
            else {
                std::cout << "Invalid option. Please try again.\n";
                pause();
            }
        }
    }
}
