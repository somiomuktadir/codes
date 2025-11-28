#include "Distributions.h"
#include "Utils.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace Distributions {

    namespace Discrete {
        double binomialPMF(int n, double p, int k) {
            if (k < 0 || k > n) return 0.0;
            if (p < 0 || p > 1) throw std::invalid_argument("Probability p must be between 0 and 1.");
            return Utils::nCr(n, k) * std::pow(p, k) * std::pow(1 - p, n - k);
        }

        double binomialCDF(int n, double p, int k) {
            if (k < 0) return 0.0;
            if (k >= n) return 1.0;
            double sum = 0.0;
            for (int i = 0; i <= k; ++i) {
                sum += binomialPMF(n, p, i);
            }
            return sum;
        }

        double binomialMean(int n, double p) {
            return n * p;
        }

        double binomialVariance(int n, double p) {
            return n * p * (1 - p);
        }

        double poissonPMF(double lambda, int k) {
            if (lambda <= 0) throw std::invalid_argument("Lambda must be positive.");
            if (k < 0) return 0.0;
            return (std::pow(lambda, k) * std::exp(-lambda)) / Utils::factorial(k);
        }

        double poissonCDF(double lambda, int k) {
            if (k < 0) return 0.0;
            double sum = 0.0;
            for (int i = 0; i <= k; ++i) {
                sum += poissonPMF(lambda, i);
            }
            return sum;
        }

        double poissonMean(double lambda) {
            return lambda;
        }

        double poissonVariance(double lambda) {
            return lambda;
        }

        double geometricPMF(double p, int k) {
            if (p <= 0 || p > 1) throw std::invalid_argument("Probability p must be in (0, 1].");
            if (k < 1) return 0.0;
            return std::pow(1 - p, k - 1) * p;
        }

        double geometricCDF(double p, int k) {
            if (k < 1) return 0.0;
            return 1 - std::pow(1 - p, k);
        }

        double geometricMean(double p) {
            return 1.0 / p;
        }

        double geometricVariance(double p) {
            return (1 - p) / (p * p);
        }

        double hypergeometricPMF(int N, int K, int n, int k) {
            if (k < 0 || k > n || k > K || n - k > N - K) return 0.0;
            return (double)(Utils::nCr(K, k) * Utils::nCr(N - K, n - k)) / Utils::nCr(N, n);
        }

        double hypergeometricMean(int N, int K, int n) {
            return (double)(n * K) / N;
        }

        double hypergeometricVariance(int N, int K, int n) {
            double p = (double)K / N;
            return n * p * (1 - p) * ((double)(N - n) / (N - 1));
        }

        double negativeBinomialPMF(int r, double p, int k) {
            if (k < 0) return 0.0;
            if (p <= 0 || p > 1) throw std::invalid_argument("Probability p must be in (0, 1].");
            return Utils::nCr(k + r - 1, k) * std::pow(p, r) * std::pow(1 - p, k);
        }

        double negativeBinomialMean(int r, double p) {
            return r * (1 - p) / p;
        }

        double negativeBinomialVariance(int r, double p) {
            return r * (1 - p) / (p * p);
        }
    }

    namespace Continuous {
        // Uniform
        double uniformPDF(double a, double b, double x) {
            if (a >= b) throw std::invalid_argument("a must be less than b");
            Utils::log("Calculating Uniform PDF");
            Utils::logStep("Range: [" + std::to_string(a) + ", " + std::to_string(b) + "], x: " + std::to_string(x));
            if (x < a || x > b) return 0.0;
            return 1.0 / (b - a);
        }
        double uniformCDF(double a, double b, double x) {
            if (a >= b) throw std::invalid_argument("a must be less than b");
            Utils::log("Calculating Uniform CDF");
            if (x < a) return 0.0;
            if (x > b) return 1.0;
            return (x - a) / (b - a);
        }
        double uniformMean(double a, double b) {
            return (a + b) / 2.0;
        }
        double uniformVariance(double a, double b) {
            return std::pow(b - a, 2) / 12.0;
        }
        double uniformStdDev(double a, double b) {
            return std::sqrt(uniformVariance(a, b));
        }

        // Normal
        double normalPDF(double mean, double stdDev, double x) {
            if (stdDev <= 0) throw std::invalid_argument("stdDev must be positive");
            Utils::log("Calculating Normal PDF");
            Utils::logStep("Mean: " + std::to_string(mean) + ", StdDev: " + std::to_string(stdDev) + ", x: " + std::to_string(x));
            double exponent = -0.5 * std::pow((x - mean) / stdDev, 2);
            return (1.0 / (stdDev * std::sqrt(2 * Utils::PI))) * std::exp(exponent);
        }
        double normalCDF(double mean, double stdDev, double x) {
            if (stdDev <= 0) throw std::invalid_argument("stdDev must be positive");
            Utils::log("Calculating Normal CDF");
            return 0.5 * (1 + std::erf((x - mean) / (stdDev * std::sqrt(2))));
        }
        double normalMean(double mean) {
            return mean;
        }
        double normalVariance(double stdDev) {
            return stdDev * stdDev;
        }
        double normalStdDev(double stdDev) {
            return stdDev;
        }

        // Exponential
        double exponentialPDF(double lambda, double x) {
            if (lambda <= 0) throw std::invalid_argument("lambda must be positive");
            Utils::log("Calculating Exponential PDF");
            if (x < 0) return 0.0;
            return lambda * std::exp(-lambda * x);
        }
        double exponentialCDF(double lambda, double x) {
            if (lambda <= 0) throw std::invalid_argument("lambda must be positive");
            Utils::log("Calculating Exponential CDF");
            if (x < 0) return 0.0;
            return 1.0 - std::exp(-lambda * x);
        }
        double exponentialMean(double lambda) {
            return 1.0 / lambda;
        }
        double exponentialVariance(double lambda) {
            return 1.0 / (lambda * lambda);
        }
        double exponentialStdDev(double lambda) {
            return 1.0 / lambda;
        }

        // Gamma
        double gammaPDF(double k, double theta, double x) {
            if (k <= 0 || theta <= 0) throw std::invalid_argument("k and theta must be positive.");
            Utils::log("Calculating Gamma PDF");
            if (x < 0) return 0.0;
            return (std::pow(x, k - 1) * std::exp(-x / theta)) / (std::pow(theta, k) * std::tgamma(k));
        }
        double gammaMean(double k, double theta) {
            return k * theta;
        }
        double gammaVariance(double k, double theta) {
            return k * theta * theta;
        }
        double gammaStdDev(double k, double theta) {
            return std::sqrt(gammaVariance(k, theta));
        }

        // Beta
        double betaPDF(double alpha, double beta, double x) {
            if (alpha <= 0 || beta <= 0) throw std::invalid_argument("alpha and beta must be positive.");
            Utils::log("Calculating Beta PDF");
            if (x < 0 || x > 1) return 0.0;
            return (std::pow(x, alpha - 1) * std::pow(1 - x, beta - 1)) / Utils::betaFunction(alpha, beta);
        }
        double betaMean(double alpha, double beta) {
            return alpha / (alpha + beta);
        }
        double betaVariance(double alpha, double beta) {
            return (alpha * beta) / (std::pow(alpha + beta, 2) * (alpha + beta + 1));
        }
        double betaStdDev(double alpha, double beta) {
            return std::sqrt(betaVariance(alpha, beta));
        }

        // Chi-Square
        double chiSquarePDF(int k, double x) {
            if (k <= 0) throw std::invalid_argument("Degrees of freedom k must be positive.");
            Utils::log("Calculating Chi-Square PDF");
            if (x < 0) return 0.0;
            return (std::pow(x, (k / 2.0) - 1) * std::exp(-x / 2.0)) / (std::pow(2, k / 2.0) * std::tgamma(k / 2.0));
        }
        double chiSquareMean(int k) {
            return k;
        }
        double chiSquareVariance(int k) {
            return 2 * k;
        }
        double chiSquareStdDev(int k) {
            return std::sqrt(2 * k);
        }

        // Student's t
        double studentTPDF(int v, double x) {
            if (v <= 0) throw std::invalid_argument("Degrees of freedom v must be positive.");
            Utils::log("Calculating Student's t PDF");
            double num = std::tgamma((v + 1) / 2.0);
            double den = std::sqrt(v * Utils::PI) * std::tgamma(v / 2.0);
            return (num / den) * std::pow(1 + (x * x) / v, -(v + 1) / 2.0);
        }
        double studentTMean(int v) {
            if (v <= 1) return std::numeric_limits<double>::quiet_NaN(); // Undefined for v <= 1
            return 0.0;
        }
        double studentTVariance(int v) {
            if (v <= 2) return std::numeric_limits<double>::infinity(); // Undefined for v <= 1, infinite for 1 < v <= 2
            return (double)v / (v - 2);
        }
        double studentTStdDev(int v) {
            if (v <= 2) return std::numeric_limits<double>::infinity();
            return std::sqrt(studentTVariance(v));
        }

        // F-Distribution
        double fDistributionPDF(int d1, int d2, double x) {
            if (d1 <= 0 || d2 <= 0) throw std::invalid_argument("Degrees of freedom d1, d2 must be positive.");
            Utils::log("Calculating F-Distribution PDF");
            if (x < 0) return 0.0;
            double num = std::sqrt(std::pow(d1 * x, d1) * std::pow(d2, d2) / std::pow(d1 * x + d2, d1 + d2));
            double den = x * Utils::betaFunction(d1 / 2.0, d2 / 2.0);
            return num / den;
        }
        double fDistributionMean(int d1, int d2) {
            (void)d1; // d1 is not used in the F-distribution mean formula
            if (d2 <= 2) return std::numeric_limits<double>::quiet_NaN(); // Undefined for d2 <= 2
            return (double)d2 / (d2 - 2);
        }
        double fDistributionVariance(int d1, int d2) {
            if (d2 <= 4) return std::numeric_limits<double>::quiet_NaN(); // Undefined for d2 <= 4
            double num = 2.0 * d2 * d2 * (d1 + d2 - 2);
            double den = d1 * std::pow(d2 - 2, 2) * (d2 - 4);
            return num / den;
        }
        double fDistributionStdDev(int d1, int d2) {
            if (d2 <= 4) return std::numeric_limits<double>::quiet_NaN();
            return std::sqrt(fDistributionVariance(d1, d2));
        }
    }

    void printFormulas() {
        std::cout << "\n--- Distribution Formulas ---\n";
        std::cout << "Binomial PMF: nCk * p^k * (1-p)^(n-k)\n";
        std::cout << "Poisson PMF: (lambda^k * e^-lambda) / k!\n";
        std::cout << "Geometric PMF: (1-p)^(k-1) * p\n";
        std::cout << "Negative Binomial PMF: (k+r-1)Ck * p^r * (1-p)^k\n";
        std::cout << "Uniform PDF: 1/(b-a)\n";
        std::cout << "Normal PDF: (1 / (sigma * sqrt(2pi))) * e^(-(x-mu)^2 / (2sigma^2))\n";
        std::cout << "Exponential PDF: lambda * e^(-lambda * x)\n";
        std::cout << "Gamma PDF: (x^(k-1) * e^(-x/theta)) / (theta^k * Gamma(k))\n";
        std::cout << "Beta PDF: (x^(alpha-1) * (1-x)^(beta-1)) / B(alpha, beta)\n";
        std::cout << "Chi-Square PDF: Special case of Gamma(k/2, 2)\n";
        std::cout << "Student's t PDF: Bell-shaped, heavier tails than Normal\n";
        std::cout << "F-Distribution PDF: Ratio of two Chi-Squares\n";
        std::cout << "-----------------------------\n";
    }
}
