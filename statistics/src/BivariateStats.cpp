#include "BivariateStats.h"
#include "UnivariateStats.h"
#include "Logger.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace Stats {

// Correlation Measures

double BivariateStats::pearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) throw std::invalid_argument("Datasets must have the same size");
    if (x.empty()) throw std::invalid_argument("Cannot compute correlation of empty datasets");
    
    log("Computing Pearson correlation coefficient...");
    
    size_t n = x.size();
    double meanX = UnivariateStats::mean(x);
    double meanY = UnivariateStats::mean(y);
    
    double sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
    
    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - meanX;
        double dy = y[i] - meanY;
        sumXY += dx * dy;
        sumX2 += dx * dx;
        sumY2 += dy * dy;
    }
    
    if (sumX2 == 0 || sumY2 == 0) {
        return 0.0; // No variation in one or both variables
    }
    
    double r = sumXY / std::sqrt(sumX2 * sumY2);
    
    logStep("Pearson r", r);
    logStep("R-squared", r * r);
    
    return r;
}

double BivariateStats::spearmanCorrelation(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) throw std::invalid_argument("Datasets must have the same size");
    if (x.empty()) throw std::invalid_argument("Cannot compute correlation of empty datasets");
    
    log("Computing Spearman rank correlation...");
    
    // Rank the data
    auto ranksX = rankData(x);
    auto ranksY = rankData(y);
    
    // Convert ranks to double vectors
    std::vector<double> rankXDouble(ranksX.begin(), ranksX.end());
    std::vector<double> rankYDouble(ranksY.begin(), ranksY.end());
    
    // Compute Pearson correlation on ranks
    double rho = pearsonCorrelation(rankXDouble, rankYDouble);
    
    logStep("Spearman rho", rho);
    
    return rho;
}

double BivariateStats::kendallTau(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) throw std::invalid_argument("Datasets must have the same size");
    if (x.empty()) throw std::invalid_argument("Cannot compute correlation of empty datasets");
    
    log("Computing Kendall's tau...");
    
    size_t n = x.size();
    int concordant = 0;
    int discordant = 0;
    
    // Count concordant and discordant pairs
    for (size_t i = 0; i < n - 1; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double diffX = x[j] - x[i];
            double diffY = y[j] - y[i];
            
            if ((diffX > 0 && diffY > 0) || (diffX < 0 && diffY < 0)) {
                concordant++;
            } else if ((diffX > 0 && diffY < 0) || (diffX < 0 && diffY > 0)) {
                discordant++;
            }
            // Tied pairs are ignored
        }
    }
    
    double totalPairs = n * (n - 1) / 2.0;
    double tau = (concordant - discordant) / totalPairs;
    
    logStep("Concordant pairs", concordant);
    logStep("Discordant pairs", discordant);
    logStep("Kendall's tau", tau);
    
    return tau;
}

// Covariance

double BivariateStats::covariance(const std::vector<double>& x, const std::vector<double>& y, bool sample) {
    if (x.size() != y.size()) throw std::invalid_argument("Datasets must have the same size");
    if (x.empty()) throw std::invalid_argument("Cannot compute covariance of empty datasets");
    if (sample && x.size() == 1) throw std::invalid_argument("Sample covariance undefined for single value");
    
    log(sample ? "Computing sample covariance..." : "Computing population covariance...");
    
    double meanX = UnivariateStats::mean(x);
    double meanY = UnivariateStats::mean(y);
    
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum += (x[i] - meanX) * (y[i] - meanY);
    }
    
    double divisor = sample ? (x.size() - 1) : x.size();
    double cov = sum / divisor;
    
    logStep("Covariance", cov);
    
    return cov;
}

// Linear Regression

BivariateStats::RegressionResult BivariateStats::linearRegression(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) throw std::invalid_argument("Datasets must have the same size");
    if (x.size() < 2) throw std::invalid_argument("Need at least 2 points for regression");
    
    log("Computing linear regression (least squares)...");
    
    size_t n = x.size();
    double meanX = UnivariateStats::mean(x);
    double meanY = UnivariateStats::mean(y);
    
    double sumXY = 0.0, sumX2 = 0.0;
    
    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - meanX;
        double dy = y[i] - meanY;
        sumXY += dx * dy;
        sumX2 += dx * dx;
    }
    
    if (sumX2 == 0) {
        throw std::invalid_argument("X values have no variation - cannot fit regression line");
    }
    
    RegressionResult result;
    result.slope = sumXY / sumX2;
    result.intercept = meanY - result.slope * meanX;
    result.correlation = pearsonCorrelation(x, y);
    result.rSquared = result.correlation * result.correlation;
    
    logStep("Slope (b1)", result.slope);
    logStep("Intercept (b0)", result.intercept);
    logStep("R-squared", result.rSquared);
    
    // Calculate predictions and residuals
    result.predictions.reserve(n);
    result.residuals.reserve(n);
    double sse = 0.0; // Sum of squared errors
    
    for (size_t i = 0; i < n; ++i) {
        double pred = result.intercept + result.slope * x[i];
        double resid = y[i] - pred;
        result.predictions.push_back(pred);
        result.residuals.push_back(resid);
        sse += resid * resid;
    }
    
    // Standard error of the estimate
    result.standardError = std::sqrt(sse / (n - 2));
    
    // Standard errors of slope and intercept
    double sumX2Centered = sumX2;
    result.slopeStdError = result.standardError / std::sqrt(sumX2Centered);
    
    double sumX2Total = 0.0;
    for (double val : x) {
        sumX2Total += val * val;
    }
    result.interceptStdError = result.standardError * std::sqrt(sumX2Total / (n * sumX2Centered));
    
    logStep("Standard error", result.standardError);
    logStep("Slope std. error", result.slopeStdError);
    
    return result;
}

double BivariateStats::predict(double x, double slope, double intercept) {
    return intercept + slope * x;
}

// Contingency Table Analysis

BivariateStats::ContingencyTable BivariateStats::chiSquareTest(const std::vector<std::vector<int>>& table) {
    if (table.empty() || table[0].empty()) {
        throw std::invalid_argument("Contingency table cannot be empty");
    }
    
    log("Computing chi-square test...");
    
    size_t rows = table.size();
    size_t cols = table[0].size();
    
    // Calculate row and column totals
    std::vector<int> rowTotals(rows, 0);
    std::vector<int> colTotals(cols, 0);
    int grandTotal = 0;
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            rowTotals[i] += table[i][j];
            colTotals[j] += table[i][j];
            grandTotal += table[i][j];
        }
    }
    
    // Calculate expected frequencies and chi-square
    ContingencyTable result;
    result.observed = table;
    result.expected.resize(rows, std::vector<double>(cols));
    result.chiSquare = 0.0;
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            double expected = (static_cast<double>(rowTotals[i]) * colTotals[j]) / grandTotal;
            result.expected[i][j] = expected;
            
            if (expected > 0) {
                double diff = table[i][j] - expected;
                result.chiSquare += (diff * diff) / expected;
            }
        }
    }
    
    result.degreesOfFreedom = (rows - 1) * (cols - 1);
    result.pValue = 1.0 - chiSquareCDF(result.chiSquare, result.degreesOfFreedom);
    
    logStep("Chi-square statistic", result.chiSquare);
    logStep("Degrees of freedom", result.degreesOfFreedom);
    logStep("P-value (approximate)", result.pValue);
    
    return result;
}

// Hypothesis Testing

double BivariateStats::tStatisticForCorrelation(double r, int n) {
    if (n < 3) throw std::invalid_argument("Need at least 3 observations for t-test");
    if (std::abs(r) >= 1.0) return INFINITY;
    
    return r * std::sqrt(n - 2) / std::sqrt(1 - r * r);
}

bool BivariateStats::isCorrelationSignificant(double r, int n, double alpha) {
    (void)alpha; // Unused in simplified approximation
    double t = std::abs(tStatisticForCorrelation(r, n));
    
    // Critical value for two-tailed test (approximate for common cases)
    // For alpha=0.05 and df>30, critical value is approximately 1.96
    // This is a simplified approximation
    double criticalValue = 1.96; // For alpha=0.05, df>30
    
    if (n - 2 <= 30) {
        // Rough approximation for smaller samples
        criticalValue = 2.0 + 0.3 / std::sqrt(n - 2);
    }
    
    return t > criticalValue;
}

void BivariateStats::printRegressionSummary(const RegressionResult& result, const std::string& xName, const std::string& yName) {
    std::cout << "\n========== Regression Analysis ==========" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "\nRegression Equation:" << std::endl;
    std::cout << yName << " = " << result.intercept << " + " << result.slope << " * " << xName << std::endl;
    
    std::cout << "\n--- Goodness of Fit ---" << std::endl;
    std::cout << "R-squared:              " << result.rSquared << std::endl;
    std::cout << "Correlation (r):        " << result.correlation << std::endl;
    std::cout << "Standard Error:         " << result.standardError << std::endl;
    
    std::cout << "\n--- Coefficients ---" << std::endl;
    std::cout << "Intercept:              " << result.intercept 
              << " (SE = " << result.interceptStdError << ")" << std::endl;
    std::cout << "Slope:                  " << result.slope 
              << " (SE = " << result.slopeStdError << ")" << std::endl;
    
    std::cout << "\n--- Residuals ---" << std::endl;
    double meanResidual = UnivariateStats::mean(result.residuals);
    double stdResidual = UnivariateStats::standardDeviation(result.residuals);
    std::cout << "Mean residual:          " << meanResidual << std::endl;
    std::cout << "Std. dev. residuals:    " << stdResidual << std::endl;
    
    std::cout << "========================================\n" << std::endl;
}

// Helper functions

std::vector<size_t> BivariateStats::rankData(const std::vector<double>& data) {
    size_t n = data.size();
    std::vector<std::pair<double, size_t>> indexed;
    indexed.reserve(n);
    
    for (size_t i = 0; i < n; ++i) {
        indexed.push_back({data[i], i});
    }
    
    std::sort(indexed.begin(), indexed.end());
    
    std::vector<size_t> ranks(n);
    
    // Handle ties by assigning average rank
    size_t i = 0;
    while (i < n) {
        size_t j = i;
        while (j < n && indexed[j].first == indexed[i].first) {
            ++j;
        }
        
        double avgRank = (i + j - 1) / 2.0 + 1; // Average rank (1-indexed)
        
        for (size_t k = i; k < j; ++k) {
            ranks[indexed[k].second] = static_cast<size_t>(avgRank);
        }
        
        i = j;
    }
    
    return ranks;
}

double BivariateStats::chiSquareCDF(double x, int df) {
    // Very simplified approximation using Wilson-Hilferty transformation
    // This is just for approximate p-values
    if (df <= 0) return 0.0;
    if (x <= 0) return 0.0;
    
    // For better accuracy, a proper chi-square CDF implementation would be needed
    // This is a placeholder approximation
    double z = std::pow(x / df, 1.0/3.0) - (1.0 - 2.0/(9.0*df));
    z /= std::sqrt(2.0/(9.0*df));
    
    // Approximate normal CDF
    return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
}

} // namespace Stats
