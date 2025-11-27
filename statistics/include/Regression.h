#ifndef REGRESSION_H
#define REGRESSION_H

#include <vector>
#include <string>

namespace Stats {

struct PolynomialRegressionResult {
    std::vector<double> coefficients; // a0, a1, a2... for y = a0 + a1*x + a2*x^2...
    double rSquared;
    double standardError;
};

struct MultipleRegressionResult {
    std::vector<double> coefficients; // b0, b1, b2... for y = b0 + b1*x1 + b2*x2...
    double rSquared;
    double standardError;
};

class Regression {
public:
    // Polynomial Regression
    static PolynomialRegressionResult polynomialRegression(
        const std::vector<double>& x, 
        const std::vector<double>& y, 
        int degree);
        
    // Multiple Linear Regression
    // x is a vector of vectors, where each inner vector is a feature column (x1, x2, etc.)
    static MultipleRegressionResult multipleLinearRegression(
        const std::vector<std::vector<double>>& x, 
        const std::vector<double>& y);
        
    // Prediction helpers
    static double predictPolynomial(double x, const std::vector<double>& coeffs);
    static double predictMultiple(const std::vector<double>& x, const std::vector<double>& coeffs);
};

} // namespace Stats

#endif // REGRESSION_H
