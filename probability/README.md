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
- **Intersection**: P(A n B) for independent events
- **Expected Value**: Mean of a discrete random variable
- **Variance**: Measure of spread for a discrete random variable

### Distributions
**Discrete:**
- **Binomial**: Number of successes in n independent Bernoulli trials
- **Poisson**: Number of events in a fixed interval
- **Geometric**: Number of trials to get the first success
- **Hypergeometric**: Sampling without replacement
- **Negative Binomial**: Trials to achieve r successes

**Continuous:**
- **Uniform**: Constant probability density
- **Normal (Gaussian)**: Bell curve distribution
- **Exponential**: Time between events in a Poisson process
- **Gamma**: Generalization of Exponential distribution
- **Beta**: Defined on interval [0, 1]
- **Chi-Square**: Sum of squared standard normals
- **Student's t**: Estimating mean of normally distributed population
- **F-Distribution**: Ratio of two chi-square distributions

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

Run the executable from the `combinatorics` directory:

```bash
./bin/combinatorics
```

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
