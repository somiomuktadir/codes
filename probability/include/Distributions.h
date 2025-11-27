#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <string>

namespace Distributions {

    // Discrete Distributions
    namespace Discrete {
        // Binomial: P(X=k) = nCk * p^k * (1-p)^(n-k)
        double binomialPMF(int n, double p, int k);
        double binomialCDF(int n, double p, int k);
        double binomialMean(int n, double p);
        double binomialVariance(int n, double p);

        // Poisson: P(X=k) = (lambda^k * e^-lambda) / k!
        double poissonPMF(double lambda, int k);
        double poissonCDF(double lambda, int k);
        double poissonMean(double lambda);
        double poissonVariance(double lambda);

        // Geometric: P(X=k) = (1-p)^(k-1) * p (trials until first success)
        double geometricPMF(double p, int k);
        double geometricCDF(double p, int k);
        double geometricMean(double p);
        double geometricVariance(double p);

        // Hypergeometric: P(X=k) = (KCk * (N-K)C(n-k)) / NCn
        // N: Population size, K: Successes in population, n: Sample size, k: Successes in sample
        double hypergeometricPMF(int N, int K, int n, int k);
        double hypergeometricMean(int N, int K, int n);
        double hypergeometricVariance(int N, int K, int n);

        // Negative Binomial: P(X=k) = (k+r-1)C(k) * p^r * (1-p)^k
        double negativeBinomialPMF(int r, double p, int k);
        double negativeBinomialMean(int r, double p);
        double negativeBinomialVariance(int r, double p);
    }

    // Continuous Distributions
    namespace Continuous {
        // Uniform: f(x) = 1/(b-a) for a <= x <= b
        double uniformPDF(double a, double b, double x);
        double uniformCDF(double a, double b, double x);
        double uniformMean(double a, double b);
        double uniformVariance(double a, double b);

        // Normal (Gaussian): f(x) = (1 / (sigma * sqrt(2pi))) * e^(-(x-mu)^2 / (2sigma^2))
        double normalPDF(double mean, double stdDev, double x);
        double normalCDF(double mean, double stdDev, double x); // Uses error function
        double normalMean(double mean);
        double normalVariance(double stdDev);

        // Exponential: f(x) = lambda * e^(-lambda * x) for x >= 0
        double exponentialPDF(double lambda, double x);
        double exponentialCDF(double lambda, double x);
        double exponentialMean(double lambda);
        double exponentialVariance(double lambda);

        // Gamma: f(x) = (x^(k-1) * e^(-x/theta)) / (theta^k * Gamma(k))
        double gammaPDF(double k, double theta, double x);
        double gammaMean(double k, double theta);
        double gammaVariance(double k, double theta);

        // Beta: f(x) = (x^(alpha-1) * (1-x)^(beta-1)) / B(alpha, beta)
        double betaPDF(double alpha, double beta, double x);
        double betaMean(double alpha, double beta);
        double betaVariance(double alpha, double beta);

        // Chi-Square: Special case of Gamma
        double chiSquarePDF(int k, double x);
        double chiSquareMean(int k);
        double chiSquareVariance(int k);

        // Student's t
        double studentTPDF(int v, double x);
        double studentTMean(int v);
        double studentTVariance(int v);

        // F-Distribution
        double fDistributionPDF(int d1, int d2, double x);
        double fDistributionMean(int d1, int d2);
        double fDistributionVariance(int d1, int d2);
    }

    void printFormulas();
}

#endif // DISTRIBUTIONS_H
