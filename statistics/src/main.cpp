#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include "Dataset.h"
#include "UnivariateStats.h"
#include "BivariateStats.h"
#include "DataIO.h"
#include "Logger.h"
#include "DataTransform.h"
#include "HypothesisTesting.h"
#include "FrequencyTable.h"
#include "GroupedDataset.h"
#include "Regression.h"
#include "TimeSeries.h"

using namespace Stats;
namespace fs = std::filesystem;

void clearInput() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

std::string selectDataFile() {
    std::string dataPath = "data";
    if (!fs::exists(dataPath)) {
        std::cout << "Error: 'data' directory not found." << std::endl;
        return "";
    }

    std::vector<std::string> files;
    std::cout << "\nAvailable Data Files:" << std::endl;
    int i = 1;
    for (const auto& entry : fs::directory_iterator(dataPath)) {
        if (entry.path().extension() == ".csv") {
            std::cout << i << ". " << entry.path().filename().string() << std::endl;
            files.push_back(entry.path().string());
            i++;
        }
    }

    if (files.empty()) {
        std::cout << "No CSV files found in 'data' directory." << std::endl;
        return "";
    }

    std::cout << "Select file (1-" << files.size() << "): ";
    int choice;
    if (!(std::cin >> choice) || choice < 1 || choice > static_cast<int>(files.size())) {
        clearInput();
        std::cout << "Invalid selection." << std::endl;
        return "";
    }
    clearInput();

    return files[choice - 1];
}

void printMenu() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "    Statistics Analysis CLI Tool" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "--- Data Management ---" << std::endl;
    std::cout << "1.  Load Data (from data/ folder)" << std::endl;
    std::cout << "2.  View Current Dataset" << std::endl;
    std::cout << "3.  Save Results" << std::endl;
    std::cout << "\n--- Analysis Modules ---" << std::endl;
    std::cout << "4.  Univariate Statistics" << std::endl;
    std::cout << "5.  Grouped Data Analysis" << std::endl;
    std::cout << "6.  Bivariate & Regression" << std::endl;
    std::cout << "7.  Time Series Analysis" << std::endl;
    std::cout << "8.  Hypothesis Testing" << std::endl;
    std::cout << "\n--- Utilities ---" << std::endl;
    std::cout << "9.  Data Transformations" << std::endl;
    std::cout << "10. Settings (Verbose Mode: " << (verboseMode ? "ON" : "OFF") << ")" << std::endl;
    std::cout << "0.  Exit" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Select option: ";
}

int main() {
    Dataset mainDataset("Data");
    Dataset xData("X"), yData("Y");
    bool hasBivariateData = false;
    
    std::cout << "Statistics Analysis Tool - Professional Edition\n";
    std::cout << "Developed in C++17 | High-Performance Statistical Computing\n" << std::endl;
    
    int choice;
    while (true) {
        printMenu();
        
        if (!(std::cin >> choice)) {
            clearInput();
            std::cout << "Invalid input. Please enter a number." << std::endl;
            continue;
        }
        clearInput();
        
        try {
            switch (choice) {
                case 0:
                    std::cout << "Exiting... Thank you for using Statistics Analysis Tool!" << std::endl;
                    return 0;
                    
                case 1: { // Load Data
                    std::string filename = selectDataFile();
                    if (filename.empty()) break;
                    
                    // Auto-detect header and columns
                    std::cout << "\n--- Analyzing file: " << filename << " ---" << std::endl;
                    
                    auto datasets = DataIO::loadMultiColumnCSV(filename, true); // Assume header
                    
                    if (datasets.empty()) {
                        std::cout << "Error: Could not load data from file." << std::endl;
                        break;
                    }
                    
                    int numCols = datasets.size();
                    std::cout << "Detected " << numCols << " column(s)" << std::endl;
                    
                    if (numCols == 1) {
                        // Univariate data
                        mainDataset = datasets[0];
                        hasBivariateData = false;
                        std::cout << "Loaded as: UNIVARIATE" << std::endl;
                        mainDataset.print();
                    } else if (numCols >= 2) {
                        // Bivariate/Multivariate data
                        xData = datasets[0];
                        yData = datasets[1];
                        hasBivariateData = true;
                        std::cout << "Loaded as: BIVARIATE (using first 2 columns)" << std::endl;
                        if (numCols > 2) {
                            std::cout << "Note: " << (numCols - 2) << " additional column(s) ignored." << std::endl;
                        }
                        xData.print();
                        yData.print();
                    }
                    
                    std::cout << "✓ Data loaded successfully!\n" << std::endl;
                    break;
                }

                case 2: { // View Dataset
                    if (hasBivariateData) {
                        std::cout << "\nBivariate Data:" << std::endl;
                        xData.print();
                        yData.print();
                        std::cout << "Sample size: " << xData.size() << " pairs" << std::endl;
                    } else if (!mainDataset.empty()) {
                        mainDataset.print();
                    } else {
                        std::cout << "No data loaded." << std::endl;
                    }
                    break;
                }

                case 3: { // Save Results
                    std::cout << "Enter filename to save (in data/ folder): ";
                    std::string filename;
                    std::getline(std::cin, filename);
                    
                    // Ensure it saves to data/ if not specified
                    if (filename.find("data/") == std::string::npos) {
                        filename = "data/" + filename;
                    }
                    
                    if (hasBivariateData) {
                        DataIO::saveCSV(filename, xData.getData(), yData.getData(), 
                                      xData.getName(), yData.getName());
                    } else {
                        DataIO::saveCSV(filename, mainDataset);
                    }
                    break;
                }
                
                case 4: { // Univariate Statistics
                    if (!hasBivariateData && mainDataset.empty()) {
                        std::cout << "No data loaded. Please load or enter data first." << std::endl;
                        break;
                    }
                    
                    const auto& data = hasBivariateData ? xData.getData() : mainDataset.getData();
                    
                    std::cout << "\n--- Univariate Statistics ---" << std::endl;
                    std::cout << "1. Full Summary" << std::endl;
                    std::cout << "2. Central Tendency (mean, median, mode)" << std::endl;
                    std::cout << "3. Dispersion (variance, std dev, range, IQR)" << std::endl;
                    std::cout << "4. Distribution Shape (skewness, kurtosis)" << std::endl;
                    std::cout << "5. Quantiles and Percentiles" << std::endl;
                    std::cout << "6. Frequency Distribution" << std::endl;
                    std::cout << "Select: ";
                    
                    int subChoice;
                    std::cin >> subChoice;
                    clearInput();
                    
                    if (subChoice == 1) {
                        UnivariateStats::printSummary(data);
                    } else if (subChoice == 2) {
                        std::cout << std::fixed << std::setprecision(6);
                        std::cout << "Mean:             " << UnivariateStats::mean(data) << std::endl;
                        std::cout << "Median:           " << UnivariateStats::median(data) << std::endl;
                        
                        auto modes = UnivariateStats::mode(data);
                        std::cout << "Mode:             ";
                        if (modes.size() == data.size()) {
                            std::cout << "No mode (all unique)" << std::endl;
                        } else {
                            for (size_t i = 0; i < modes.size() && i < 5; ++i) {
                                std::cout << modes[i] << " ";
                            }
                            std::cout << std::endl;
                        }
                        
                        if (data.size() > 0 && *std::min_element(data.begin(), data.end()) > 0) {
                            std::cout << "Geometric Mean:   " << UnivariateStats::geometricMean(data) << std::endl;
                        }
                    } else if (subChoice == 3) {
                        std::cout << std::fixed << std::setprecision(6);
                        std::cout << "Range:            " << UnivariateStats::range(data) << std::endl;
                        std::cout << "IQR:              " << UnivariateStats::interquartileRange(data) << std::endl;
                        std::cout << "Variance:         " << UnivariateStats::variance(data) << std::endl;
                        std::cout << "Std. Deviation:   " << UnivariateStats::standardDeviation(data) << std::endl;
                        std::cout << "MAD:              " << UnivariateStats::meanAbsoluteDeviation(data) << std::endl;
                    } else if (subChoice == 4) {
                        if (data.size() >= 3) {
                            std::cout << std::fixed << std::setprecision(6);
                            std::cout << "Skewness:         " << UnivariateStats::skewness(data) << std::endl;
                            if (data.size() >= 4) {
                                std::cout << "Excess Kurtosis:  " << UnivariateStats::kurtosis(data) << std::endl;
                            }
                        } else {
                            std::cout << "Need at least 3 values for distribution shape analysis." << std::endl;
                        }
                    } else if (subChoice == 5) {
                        auto fns = UnivariateStats::fiveNumberSummary(data);
                        std::cout << std::fixed << std::setprecision(6);
                        std::cout << "Minimum:          " << fns.minimum << std::endl;
                        std::cout << "Q1 (25%):         " << fns.q1 << std::endl;
                        std::cout << "Median (50%):     " << fns.median << std::endl;
                        std::cout << "Q3 (75%):         " << fns.q3 << std::endl;
                        std::cout << "Maximum:          " << fns.maximum << std::endl;
                        
                        std::cout << "\nEnter custom percentile (0-100): ";
                        double p;
                        if (std::cin >> p) {
                            std::cout << p << "th percentile: " 
                                    << UnivariateStats::percentile(data, p) << std::endl;
                        }
                        clearInput();
                    } else if (subChoice == 6) {
                        auto freq = UnivariateStats::frequencyDistribution(data);
                        std::cout << "\nValue\t\tFrequency" << std::endl;
                        std::cout << "----------------------------" << std::endl;
                        
                        int count = 0;
                        for (const auto& pair : freq) {
                            std::cout << pair.first << "\t\t" << pair.second << std::endl;
                            if (++count > 20) {
                                std::cout << "... (" << (freq.size() - 20) << " more)" << std::endl;
                                break;
                            }
                        }
                    }
                    break;
                }
                
                case 5: { // Grouped Data Analysis
                    if (!hasBivariateData && mainDataset.empty()) {
                        std::cout << "No data loaded. Please load or enter data first." << std::endl;
                        break;
                    }
                    
                    const auto& data = hasBivariateData ? xData.getData() : mainDataset.getData();
                    
                    std::cout << "\n--- Grouped Data Analysis ---" << std::endl;
                    std::cout << "1. Generate Frequency Table" << std::endl;
                    std::cout << "2. Calculate Grouped Statistics" << std::endl;
                    std::cout << "Select: ";
                    
                    int subChoice;
                    std::cin >> subChoice;
                    clearInput();
                    
                    if (subChoice == 1) {
                        std::cout << "Enter number of classes (0 for auto): ";
                        int k;
                        std::cin >> k;
                        clearInput();
                        
                        auto table = FrequencyTable::createFromData(data, k);
                        table.print();
                    } else if (subChoice == 2) {
                        std::cout << "Enter number of classes (0 for auto): ";
                        int k;
                        std::cin >> k;
                        clearInput();
                        
                        auto table = FrequencyTable::createFromData(data, k);
                        GroupedDataset grouped(table);
                        
                        std::cout << std::fixed << std::setprecision(6);
                        std::cout << "\n--- Grouped Statistics ---" << std::endl;
                        std::cout << "Mean:             " << grouped.mean() << std::endl;
                        std::cout << "Median:           " << grouped.median() << std::endl;
                        
                        auto modes = grouped.mode();
                        std::cout << "Mode:             ";
                        if (modes.empty()) std::cout << "None";
                        else for (double m : modes) std::cout << m << " ";
                        std::cout << std::endl;
                        
                        std::cout << "Variance:         " << grouped.variance() << std::endl;
                        std::cout << "Std. Deviation:   " << grouped.standardDeviation() << std::endl;
                    }
                    break;
                }

                case 6: { // Bivariate & Regression Analysis
                    if (!hasBivariateData) {
                        std::cout << "No bivariate data loaded. Please load X and Y data first." << std::endl;
                        break;
                    }
                    
                    if (xData.size() != yData.size()) {
                        std::cout << "Error: X and Y must have the same size." << std::endl;
                        break;
                    }
                    
                    std::cout << "\n--- Bivariate & Regression Analysis ---" << std::endl;
                    std::cout << "1. Correlation Analysis" << std::endl;
                    std::cout << "2. Simple Linear Regression" << std::endl;
                    std::cout << "3. Polynomial Regression" << std::endl;
                    std::cout << "4. Multiple Linear Regression (Requires extra features)" << std::endl;
                    std::cout << "Select: ";
                    
                    int subChoice;
                    std::cin >> subChoice;
                    clearInput();
                    
                    const auto& x = xData.getData();
                    const auto& y = yData.getData();
                    
                    if (subChoice == 1) {
                        // ... Existing correlation code ...
                        std::cout << std::fixed << std::setprecision(6);
                        std::cout << "\n--- Correlation Coefficients ---" << std::endl;
                        double pearson = BivariateStats::pearsonCorrelation(x, y);
                        std::cout << "Pearson r:        " << pearson << std::endl;
                        std::cout << "R-squared:        " << pearson * pearson << std::endl;
                        std::cout << "Spearman rho:     " << BivariateStats::spearmanCorrelation(x, y) << std::endl;
                        std::cout << "Kendall tau:      " << BivariateStats::kendallTau(x, y) << std::endl;
                    } else if (subChoice == 2) {
                         // ... Existing linear regression code ...
                        auto result = BivariateStats::linearRegression(x, y);
                        BivariateStats::printRegressionSummary(result, xData.getName(), yData.getName());
                    } else if (subChoice == 3) {
                        std::cout << "Enter polynomial degree: ";
                        int degree;
                        std::cin >> degree;
                        clearInput();
                        
                        try {
                            auto result = Regression::polynomialRegression(x, y, degree);
                            std::cout << "\n--- Polynomial Regression (Degree " << degree << ") ---" << std::endl;
                            std::cout << "R-squared:      " << result.rSquared << std::endl;
                            std::cout << "Std. Error:     " << result.standardError << std::endl;
                            std::cout << "Coefficients:   ";
                            for (size_t i = 0; i < result.coefficients.size(); ++i) {
                                std::cout << "a" << i << "=" << result.coefficients[i] << " ";
                            }
                            std::cout << std::endl;
                        } catch (const std::exception& e) {
                            std::cout << "Error: " << e.what() << std::endl;
                        }
                    } else if (subChoice == 4) {
                        std::cout << "Multiple regression requires loading multiple X columns." << std::endl;
                        std::cout << "Currently only simple X-Y data is loaded." << std::endl;
                        // In a full implementation, we would allow selecting multiple columns from a loaded CSV
                    }
                    break;
                }

                case 7: { // Time Series Analysis
                    if (!hasBivariateData && mainDataset.empty()) {
                        std::cout << "No data loaded." << std::endl;
                        break;
                    }
                    
                    const auto& data = hasBivariateData ? xData.getData() : mainDataset.getData();
                    
                    std::cout << "\n--- Time Series Analysis ---" << std::endl;
                    std::cout << "1. Autocorrelation (ACF)" << std::endl;
                    std::cout << "Select: ";
                    
                    int subChoice;
                    std::cin >> subChoice;
                    clearInput();
                    
                    if (subChoice == 1) {
                        std::cout << "Enter max lag: ";
                        int maxLag;
                        std::cin >> maxLag;
                        clearInput();
                        
                        auto acf = TimeSeries::calculateACF(data, maxLag);
                        std::cout << "\nLag\tAutocorrelation" << std::endl;
                        for (size_t k = 0; k < acf.size(); ++k) {
                            std::cout << k << "\t" << acf[k] << std::endl;
                        }
                    }
                    break;
                }

                case 9: { // Data Transformations
                    if (!hasBivariateData && mainDataset.empty()) {
                        std::cout << "No data loaded. Please load or enter data first." << std::endl;
                        break;
                    }
                    
                    const auto& data = hasBivariateData ? xData.getData() : mainDataset.getData();
                    
                    std::cout << "\n--- Data Transformations ---" << std::endl;
                    std::cout << "1. Z-score Standardization" << std::endl;
                    std::cout << "2. Min-Max Normalization" << std::endl;
                    std::cout << "3. Trimmed Mean" << std::endl;
                    std::cout << "4. Log Transform" << std::endl;
                    std::cout << "5. Square Root Transform" << std::endl;
                    std::cout << "Select: ";
                    
                    int subChoice;
                    std::cin >> subChoice;
                    clearInput();
                    
                    if (subChoice == 1) {
                        auto zScores = DataTransform::zScoreStandardize(data);
                        std::cout << "\nZ-scores (first 10):" << std::endl;
                        for (size_t i = 0; i < zScores.size() && i < 10; ++i) {
                            std::cout << std::fixed << std::setprecision(4) << zScores[i] << " ";
                        }
                        if (zScores.size() > 10) std::cout << "...";
                        std::cout << std::endl;
                    } else if (subChoice == 2) {
                        std::cout << "Enter new min (default 0): ";
                        double newMin = 0.0;
                        if (!(std::cin >> newMin)) {
                            newMin = 0.0;
                            clearInput();
                        }
                        std::cout << "Enter new max (default 1): ";
                        double newMax = 1.0;
                        if (!(std::cin >> newMax)) {
                            newMax = 1.0;
                            clearInput();
                        }
                        clearInput();
                        
                        auto normalized = DataTransform::minMaxNormalize(data, newMin, newMax);
                        std::cout << "\nNormalized values (first 10):" << std::endl;
                        for (size_t i = 0; i < normalized.size() && i < 10; ++i) {
                            std::cout << std::fixed << std::setprecision(4) << normalized[i] << " ";
                        }
                        if (normalized.size() > 10) std::cout << "...";
                        std::cout << std::endl;
                    } else if (subChoice == 3) {
                        std::cout << "Enter trim proportion (0-0.5, default 0.1): ";
                        double trim = 0.1;
                        if (!(std::cin >> trim)) {
                            trim = 0.1;
                        }
                        clearInput();
                        
                        double tMean = DataTransform::trimmedMean(data, trim);
                        std::cout << std::fixed << std::setprecision(6);
                        std::cout << "\nTrimmed Mean " << (trim*100) << "% trim): " << tMean << std::endl;
                        std::cout << "Regular Mean:        " << UnivariateStats::mean(data) << std::endl;
                    } else if (subChoice == 4) {
                        auto transformed = DataTransform::logTransform(data);
                        std::cout << "\nLog-transformed values (first 10):" << std::endl;
                        for (size_t i = 0; i < transformed.size() && i < 10; ++i) {
                            std::cout << std::fixed << std::setprecision(4) << transformed[i] << " ";
                        }
                        if (transformed.size() > 10) std::cout << "...";
                        std::cout << std::endl;
                    } else if (subChoice == 5) {
                        auto transformed = DataTransform::sqrtTransform(data);
                        std::cout << "\nSquare root transformed values (first 10):" << std::endl;
                        for (size_t i = 0; i < transformed.size() && i < 10; ++i) {
                            std::cout << std::fixed << std::setprecision(4) << transformed[i] << " ";
                        }
                        if (transformed.size() > 10) std::cout << "...";
                        std::cout << std::endl;
                    }
                    break;
                }
                
                case 8: { // Hypothesis Testing
                    if (!hasBivariateData && mainDataset.empty()) {
                        std::cout << "No data loaded." << std::endl;
                        break;
                    }
                    
                    std::cout << "\n--- Hypothesis Testing ---" << std::endl;
                    std::cout << "1. One-Sample t-test" << std::endl;
                    std::cout << "2. Two-Sample t-test (independent)" << std::endl;
                    std::cout << "3. Paired t-test" << std::endl;
                    std::cout << "4. Confidence Interval for Mean" << std::endl;
                    std::cout << "Select: ";
                    
                    int subChoice;
                    std::cin >> subChoice;
                    clearInput();
                    
                    if (subChoice == 1) {
                        const auto& data = hasBivariateData ? xData.getData() : mainDataset.getData();
                        
                        std::cout << "Enter hypothesized mean (μ₀): ";
                        double mu0;
                        std::cin >> mu0;
                        clearInput();
                        
                        auto result = HypothesisTesting::oneSampleTTest(data, mu0);
                        HypothesisTesting::printTTestResult(result, "One-Sample t-test");
                    } else if (subChoice == 2) {
                        // For two-sample, we need to load two datasets. 
                        // If we have bivariate data loaded, we can use X and Y.
                        // If not, we might need to ask user to load another file?
                        // Given the constraint "obliterate any other data input method", 
                        // we must rely on loaded data.
                        
                        if (!hasBivariateData) {
                            std::cout << "Error: Two-sample test requires bivariate data (two columns) to be loaded." << std::endl;
                            break;
                        }
                        
                        std::cout << "\nUsing loaded X as Group 1 and Y as Group 2." << std::endl;
                        
                        std::cout << "Assume equal variances? (y/n): ";
                        char eq;
                        std::cin >> eq;
                        clearInput();
                        
                        auto result = HypothesisTesting::twoSampleTTest(
                            xData.getData(), yData.getData(), eq == 'y');
                        HypothesisTesting::printTTestResult(result, "Two-Sample t-test");
                    } else if (subChoice == 3) {
                        if (!hasBivariateData) {
                            std::cout << "Error: Paired test requires bivariate data (two columns) to be loaded." << std::endl;
                            break;
                        }
                        
                        std::cout << "\nUsing loaded X as 'Before' and Y as 'After'." << std::endl;
                        
                        auto result = HypothesisTesting::pairedTTest(
                            xData.getData(), yData.getData());
                        HypothesisTesting::printTTestResult(result, "Paired t-test");
                    } else if (subChoice == 4) {
                        const auto& data = hasBivariateData ? xData.getData() : mainDataset.getData();
                        
                        std::cout << "Enter confidence level (0-1, default 0.95): ";
                        double conf = 0.95;
                        if (!(std::cin >> conf)) {
                            conf = 0.95;
                        }
                        clearInput();
                        
                        auto ci = HypothesisTesting::confidenceInterval(data, conf);
                        std::cout << std::fixed << std::setprecision(6);
                        std::cout << "\n" << (conf*100) << "% Confidence Interval: [" 
                                  << ci.first << ", " << ci.second << "]" << std::endl;
                        std::cout << "Sample Mean: " << UnivariateStats::mean(data) << std::endl;
                    }
                    break;
                }
                
                case 10: { // Toggle Verbose Mode
                    verboseMode = !verboseMode;
                    std::cout << "Verbose mode is now " << (verboseMode ? "ON" : "OFF") << std::endl;
                    std::cout << "Verbose mode shows step-by-step calculation details." << std::endl;
                    break;
                }
                

                
                default:
                    std::cout << "Invalid option. Please try again." << std::endl;
                    break;
            }
        } catch (const std::exception& e) {
            std::cout << "Error: " << e.what() << std::endl;
        }
    }
    
    return 0;
}
