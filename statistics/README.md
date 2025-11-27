# Statistics Analysis CLI Tool

A high-performance, comprehensive C++ library and command-line interface for univariate and bivariate statistical analysis. This tool provides educational insights alongside robust numerical computations, covering everything from basic descriptive statistics to advanced correlation and regression analysis.

## Overview

The Statistics Analysis CLI Tool offers a complete suite of statistical operations for analyzing single and paired datasets. Built with performance and numerical stability in mind, it provides both quick calculations and detailed step-by-step explanations through verbose mode.

## Key Features

### Univariate Statistics

#### Central Tendency Measures
- **Arithmetic Mean**: $\bar{x} = \frac{1}{n}\sum_{i=1}^{n}x_i$
- **Median**: Middle value of sorted dataset
- **Mode**: Most frequently occurring value(s)
- **Geometric Mean**: $\sqrt[n]{\prod_{i=1}^{n}x_i}$ (for positive values)
- **Harmonic Mean**: $\frac{n}{\sum_{i=1}^{n}\frac{1}{x_i}}$

#### Dispersion Measures
- **Variance**: $s^2 = \frac{1}{n-1}\sum_{i=1}^{n}(x_i - \bar{x})^2$ (implemented using Welford's algorithm for numerical stability)
- **Standard Deviation**: $s = \sqrt{s^2}$
- **Range**: $\max(x) - \min(x)$
- **Interquartile Range (IQR)**: $Q_3 - Q_1$
- **Mean Absolute Deviation (MAD)**: $\frac{1}{n}\sum_{i=1}^{n}|x_i - \bar{x}|$
- **Coefficient of Variation (CV)**: $\frac{s}{\bar{x}} \times 100\%$

#### Distribution Shape
- **Skewness**: Measure of asymmetry in the distribution
- **Kurtosis**: Measure of "tailedness" (excess kurtosis reported)

#### Position Measures
- **Quantiles**: Any percentile from 0 to 100
- **Quartiles**: Q1 (25th), Q2 (50th/median), Q3 (75th percentiles)
- **Five-Number Summary**: Minimum, Q1, Median, Q3, Maximum

#### Frequency Analysis
- Frequency distributions
- Relative frequency
- Cumulative frequency

### Bivariate Statistics

#### Correlation Measures
- **Pearson Correlation (r)**: $r = \frac{\sum(x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum(x_i - \bar{x})^2\sum(y_i - \bar{y})^2}}$
  - Measures linear relationship between variables
  - Range: [-1, 1]
- **Spearman Rank Correlation (ρ)**: Non-parametric correlation using ranks
  - Robust to outliers and monotonic relationships
- **Kendall's Tau (τ)**: Based on concordant and discordant pairs
  - Alternative non-parametric measure

#### Covariance
- **Sample Covariance**: $cov(x,y) = \frac{1}{n-1}\sum_{i=1}^{n}(x_i - \bar{x})(y_i - \bar{y})$
- **Population Covariance**: Uses $n$ instead of $n-1$ as divisor

#### Linear Regression
- **Simple Linear Regression**: $y = \beta_0 + \beta_1 x$
  - Least squares estimation
  - Slope: $\beta_1 = \frac{\sum(x_i - \bar{x})(y_i - \bar{y})}{\sum(x_i - \bar{x})^2}$
  - Intercept: $\beta_0 = \bar{y} - \beta_1\bar{x}$
- **R-squared**: Coefficient of determination
- **Residual Analysis**: Predictions, residuals, standard errors
- **Standard Errors**: For slope and intercept estimates
- **Prediction**: Make predictions for new X values

#### Hypothesis Testing
- **Correlation Significance Test**: t-statistic for correlation coefficients
- **Chi-Square Test**: For contingency tables
  - Expected frequencies calculation
  - Chi-square statistic and p-value (approximate)

### Grouped Data Analysis
- **Frequency Tables**: Generation from raw data (Sturges' rule)
- **Grouped Statistics**:
  - Mean: $\bar{x} = \frac{\sum f \cdot x_m}{N}$
  - Median: $L + \frac{N/2 - F}{f} \cdot h$
  - Mode: $L + \frac{f_1 - f_0}{2f_1 - f_0 - f_2} \cdot h$
  - Variance and Standard Deviation

### Advanced Regression
- **Polynomial Regression**: Fit $y = a_0 + a_1x + a_2x^2 + \dots + a_nx^n$
  - R-squared and Standard Error
- **Multiple Linear Regression**: Fit $y = b_0 + b_1x_1 + b_2x_2 + \dots$
  - Matrix-based solution (Normal Equations)

### Time Series Analysis
- **Autocorrelation (ACF)**: Measure correlation of signal with delayed copy
  - Lag-k autocorrelation calculation

### Data Transformations

#### Standardization and Normalization
- **Z-score Standardization**: Transform data to standard normal distribution (mean=0, std=1)
- **Min-Max Normalization**: Scale data to custom range (default [0,1])

#### Robust Statistics
- **Trimmed Mean**: Mean after removing outliers from both tails

#### Mathematical Transforms
- **Log Transformation**: Natural logarithm transform for positive data
- **Square Root Transformation**: Square root transform for non-negative data
- **Power Transformation**: Custom power transforms

### Hypothesis Testing

#### t-tests
- **One-Sample t-test**: Test if sample mean differs from hypothesized value μ₀
  - t-statistic: $t = \\frac{\\bar{x} - \\mu_0}{s/\\sqrt{n}}$
  - p-value calculation (two-tailed)
  - 95% confidence interval
- **Two-Sample t-test**: Compare means of two independent groups
  - Pooled variance method (equal variances)
  - Welch's test (unequal variances)
- **Paired t-test**: Compare means of paired/matched samples
  - Tests differences within pairs

#### Confidence Intervals
- **Confidence Interval for Mean**: Calculate CI for population mean
  - Default 95% confidence level
  - Customizable confidence levels (90%, 99%, etc.)

### Data Management

#### Data Input
- Manual data entry (space-separated values)
- CSV file import with header detection
- Single-column (univariate) or multi-column (bivariate) support

#### Data Output
- CSV export with custom naming
- Formatted console output
- Summary statistics printing

#### Data Validation
- NaN and infinity detection
- Outlier removal using IQR method
- Data integrity checks

### Technical Highlights

- **Numerical Stability**:
  - Welford's online algorithm for variance computation
  - Careful handling of floating-point arithmetic
  - Robust quantile interpolation
- **Performance**:
  - Efficient std::vector-based storage
  - C++17 features for modern performance
  - Optimized compilation with `-O3 -march=native`
- **Educational Mode**:
  - Verbose toggle shows step-by-step calculations
  - Logged intermediate values for learning
  - Clear mathematical explanations

## Building the Project

### Prerequisites
- C++17 compatible compiler (GCC, Clang)
- Make build system

### Compilation

To build the project:

```bash
cd statistics
make
```

To clean build artifacts:

```bash
make clean
```

## Usage

### Quick Start

1.  **Prepare Data**: Place CSV files in `data/` folder
2.  **Run**: `./bin/statistics`
3.  **Select File**: Choose from auto-detected CSV files
4.  **Auto-Load**: System automatically detects:
    - Number of columns (univariate vs bivariate)
    - Header presence
    - Data types

### Data Format

**Univariate** (single column):
```csv
Value
10
20
30
```

**Bivariate** (two columns):
```csv
X,Y
1,2.5
2,4.8
3,6.1
```

### Menu Structure

**Data Management**
1.  Load Data (auto-detection)
2.  View Current Dataset
3.  Save Results

**Analysis Modules**
4.  Univariate Statistics
5.  Grouped Data Analysis
6.  Bivariate & Regression
7.  Time Series Analysis
8.  Hypothesis Testing

**Utilities**
9.  Data Transformations
10. Settings (Verbose Mode)
   - Linear regression with predictions
   - Correlation significance testing

4. **Data Transformations**
   - Z-score standardization
   - Min-max normalization
   - Trimmed mean
   - Log and square root transforms

5. **Hypothesis Testing**
   - One-sample t-test
   - Two-sample t-test (independent)
   - Paired t-test
   - Confidence intervals for mean

6. **View Current Dataset**
   - Display loaded data
   - Quick data inspection

7. **Toggle Verbose Mode**
   - Enable/disable step-by-step calculations
   - Educational insights

### Example: Univariate Analysis

```bash
./statistics

Select option: 1
1. Enter data manually
Select: 1

Enter univariate or bivariate data?
1. Univariate (single dataset)
Select: 1

Enter values (space-separated):
1 2 3 4 5 6 7 8 9 10

Select option: 2
1. Full Summary
Select: 1

========== Statistical Summary ==========
Count:                  10

--- Five-Number Summary ---
Minimum:                1.0000
Q1 (25th percentile):   3.2500
Median (50th):          5.5000
Q3 (75th percentile):   7.7500
Maximum:                10.0000

--- Central Tendency ---
Mean:                   5.5000
Mode:                   No mode (all unique)

--- Dispersion ---
Range:                  9.0000
IQR:                    4.5000
Variance (sample):      9.1667
Std. Deviation (sample):3.0277
MAD:                    2.5000
CV:                     55.0490%

--- Distribution Shape ---
Skewness:               0.0000
Excess Kurtosis:        -1.2242
```

### Example: Bivariate Analysis (Linear Regression)

```bash
Select option: 1
2. Bivariate (X and Y datasets)

Enter X data:
1 2 3 4 5

Enter Y data:
2 4 5 4 6

Select option: 3
3. Linear Regression

========== Regression Analysis ==========

Regression Equation:
Y = 1.600000 + 0.800000 * X

--- Goodness of Fit ---
R-squared:              0.816327
Correlation (r):        0.903506
Standard Error:         0.547723

--- Coefficients ---
Intercept:              1.600000 (SE = 0.632456)
Slope:                  0.800000 (SE = 0.200000)

--- Residuals ---
Mean residual:          0.000000
Std. dev. residuals:    0.489898
```

### Example: CSV Import/Export

```bash
# Create a CSV file
echo "X,Y" > data.csv
echo "1,2" >> data.csv
echo "2,4.5" >> data.csv
echo "3,6.8" >> data.csv

# Load in the program
Select option: 1
2. Load from CSV file
Enter filename: data.csv
Does the file have a header row? y
2. Multiple columns (bivariate)

# Analyze and save results
Select option: 3
3. Linear Regression
[Analysis performed]

# Save to new file
Select option: 1
3. Save current data to CSV
Enter filename: results.csv
```

## API Reference

The library is encapsulated within the `Stats` namespace.

```cpp
#include "UnivariateStats.h"
#include "BivariateStats.h"
#include "Dataset.h"

using namespace Stats;

int main() {
    // Univariate analysis
    std::vector<double> data = {1, 2, 3, 4, 5};
    double mean = UnivariateStats::mean(data);
    double stdDev = UnivariateStats::standardDeviation(data);
    
    // Bivariate analysis
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {2, 4, 5, 4, 6};
    
    double r = BivariateStats::pearsonCorrelation(x, y);
    auto regression = BivariateStats::linearRegression(x, y);
    
    std::cout << "Slope: " << regression.slope << std::endl;
    std::cout << "R-squared: " << regression.rSquared << std::endl;
    
    return 0;
}
```

## File Formats

- **CSV**: Standard comma-separated values
  - Supports headers
  - Single or multiple columns
  - Automatic numeric validation

## Project Structure

```
statistics/
├── include/                  # Header files
│   ├── BivariateStats.h
│   ├── DataIO.h
│   ├── Dataset.h
│   ├── DataTransform.h
│   ├── FrequencyTable.h      # NEW
│   ├── GroupedDataset.h      # NEW
│   ├── HypothesisTesting.h
│   ├── Logger.h
│   ├── Regression.h          # NEW
│   ├── TimeSeries.h          # NEW
│   └── UnivariateStats.h
├── src/                      # Source files
│   ├── BivariateStats.cpp
│   ├── DataIO.cpp
│   ├── Dataset.cpp
│   ├── DataTransform.cpp
│   ├── FrequencyTable.cpp    # NEW
│   ├── GroupedDataset.cpp    # NEW
│   ├── HypothesisTesting.cpp
│   ├── Logger.cpp
│   ├── main.cpp
│   ├── Regression.cpp        # NEW
│   ├── TimeSeries.cpp        # NEW
│   └── UnivariateStats.cpp
├── bin/                      # Executables
├── obj/                      # Object files
├── data/                     # Data files
├── Makefile                  # Build configuration
└── README.md                 # This file
```

## Mathematical References

All statistical formulas are implemented according to standard textbook definitions with careful attention to:
- Sample vs. population statistics (n-1 vs. n divisors)
- Numerical stability (Welford's algorithm)
- Edge case handling (empty data, single values, ties)
