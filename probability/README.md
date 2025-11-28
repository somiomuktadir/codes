# Combinatorics & Probability Calculator

A comprehensive C++ command-line tool for calculating combinatorics, probability, and statistical distributions. Designed for educational purposes and quick calculations.

## Features

### Combinatorics
- **Factorial**: Calculate n!
- **Permutations**: nPr (ordered selection)
- **Combinations**: nCr (unordered selection)
- **Permutations with Repetition**: Arrangements of multisets
- **Circular Permutations**: Arrangements in a circle
- **Derangements**: Permutations with no fixed points (!n)
- **Catalan Numbers**: Sequence of natural numbers in combinatorial math
- **Stirling Numbers (2nd Kind)**: Partitioning a set into k non-empty subsets
- **Bell Numbers**: Number of partitions of a set
- **Fibonacci Numbers**: Sequence where F(n) = F(n-1) + F(n-2)
- **Lucas Numbers**: Similar to Fibonacci but starts with 2, 1
- **Integer Partitions**: Ways to write n as sum of positive integers
- **Multinomial Coefficients**: Permutations of a multiset

### Probability
- **Basic Probability**: Favorable / Total outcomes
- **Complementary Probability**: 1 - P(E)
- **Conditional Probability**: P(A|B)
- **Bayes' Theorem**: Calculate posterior probability
- **Union Probability**: P(A U B)
- **Intersection (Independent)**: P(A n B) for independent events
- **Intersection (Dependent)**: P(A n B) = P(A) × P(B|A) for dependent events
- **Total Probability**: Law of total probability for partitioned events
- **Expected Value**: Mean of a discrete random variable
- **Variance**: Measure of spread for a discrete random variable

### Distributions
**Discrete:**
- **Binomial**: Number of successes in n independent Bernoulli trials (PMF, CDF, Mean, Variance)
- **Poisson**: Number of events in a fixed interval (PMF, CDF, Mean, Variance)
- **Geometric**: Number of trials to get the first success (PMF, CDF, Mean, Variance)
- **Hypergeometric**: Sampling without replacement (PMF, Mean, Variance)
- **Negative Binomial**: Trials to achieve r successes (PMF, Mean, Variance)

**Continuous:**
- **Uniform**: Constant probability density (PDF, CDF, Mean, Variance, StdDev)
- **Normal (Gaussian)**: Bell curve distribution (PDF, CDF, Variance, StdDev)
- **Exponential**: Time between events in a Poisson process (PDF, CDF, Mean, Variance, StdDev)
- **Gamma**: Generalization of Exponential distribution (PDF, Mean, Variance, StdDev)
- **Beta**: Defined on interval [0, 1] (PDF, Mean, Variance, StdDev)
- **Chi-Square**: Sum of squared standard normals (PDF, Mean, Variance, StdDev)
- **Student's t**: Estimating mean of normally distributed population (PDF, Mean, Variance, StdDev)
- **F-Distribution**: Ratio of two chi-square distributions (PDF, Mean, Variance, StdDev)

### Verbose Mode
- **Toggle Verbose Logging**: Enable/disable step-by-step calculation details
- **Formula Display**: Shows mathematical formulas being applied
- **Intermediate Steps**: Displays parameter substitution and calculation progression
- **Educational Tool**: Helps users understand the computation process

## Building the Project

### Prerequisites
- C++17 compatible compiler (GCC, Clang)
- Make build system

### Compilation

To build the project:

```bash
cd combinatorics
make
```

To clean build artifacts:

```bash
make clean
```

## Usage

Run the executable from the `probability` directory:

```bash
./bin/probability
```

### Main Menu Options
1. **Combinatorics** - Access all combinatorial calculations
2. **Probability** - Basic and advanced probability calculations
3. **Distributions** - Discrete and continuous distribution functions
4. **Toggle Verbose Mode** - Enable/disable step-by-step calculation display
0. **Exit** - Close the application

### Using Verbose Mode
When verbose mode is enabled, the calculator will display:
- Mathematical formulas being applied
- Parameter substitution steps
- Intermediate calculation results
- Final answers with context

This feature is particularly helpful for educational purposes and understanding how calculations are performed.

Follow the on-screen menu to select modules and perform calculations.

## Project Structure

```
combinatorics/
├── include/          # Header files
├── src/              # Source files
├── bin/              # Executables
├── obj/              # Object files
├── data/             # Data files (placeholder)
├── Makefile          # Build configuration
└── README.md         # This file
```
