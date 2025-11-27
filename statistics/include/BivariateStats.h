#ifndef BIVARIATE_STATS_H
#define BIVARIATE_STATS_H

#include <vector>
#include <string>
#include <utility>

namespace Stats {

class BivariateStats {
public:
    // Correlation Measures
    static double pearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y);
    static double spearmanCorrelation(const std::vector<double>& x, const std::vector<double>& y);
    static double kendallTau(const std::vector<double>& x, const std::vector<double>& y);
    
    // Covariance
    static double covariance(const std::vector<double>& x, const std::vector<double>& y, bool sample = true);
    
    // Linear Regression
    struct RegressionResult {
        double slope;
        double intercept;
        double rSquared;
        double correlation;
        std::vector<double> residuals;
        std::vector<double> predictions;
        double standardError;
        double slopeStdError;
        double interceptStdError;
    };
    static RegressionResult linearRegression(const std::vector<double>& x, const std::vector<double>& y);
    
    // Prediction
    static double predict(double x, double slope, double intercept);
    
    // Contingency Table Analysis
    struct ContingencyTable {
        std::vector<std::vector<int>> observed;
        std::vector<std::vector<double>> expected;
        double chiSquare;
        int degreesOfFreedom;
        double pValue; // Approximate (requires chi-square distribution table)
    };
    static ContingencyTable chiSquareTest(const std::vector<std::vector<int>>& table);
    
    // Hypothesis Testing
    static double tStatisticForCorrelation(double r, int n);
    static bool isCorrelationSignificant(double r, int n, double alpha = 0.05);
    
    // Utilities
    static void printRegressionSummary(const RegressionResult& result, const std::string& xName = "X", const std::string& yName = "Y");
    
private:
    // Helper functions
    static std::vector<size_t> rankData(const std::vector<double>& data);
    static double chiSquareCDF(double x, int df); // Approximate
};

} // namespace Stats

#endif // BIVARIATE_STATS_H
