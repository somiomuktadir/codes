
A collection of C programs demonstrating core computer science concepts, from memory management to advanced algorithms.

-   **`01_basics_revision/`**: Foundational programming constructs, string manipulation, and basic algorithmic logic.
-   **`02_pointers_and_memory/`**: Advanced memory management, including pointer arithmetic, function pointers, and dynamic allocation strategies.
-   **`03_data_structures/`**: Implementation of essential data structures such as Linked Lists, Stacks, Queues, Binary Trees, and Hash Tables.
-   **`04_file_io_and_system/`**: System-level programming, including file operations, CSV parsing, and logging mechanisms.
-   **`05_advanced_algorithms/`**: Implementation of sorting algorithms (Merge Sort, Quick Sort), Graph algorithms (Dijkstra), and Dynamic Programming solutions.
-   **`06_data_analysis/`**: Statistical analysis tools and data processing utilities.

### 2. Linear Algebra Library (`lin_alg/`)

A high-performance C++ library designed for advanced linear algebra operations. This project features a custom CLI tool and supports:

-   **Matrix Decompositions**: LU, QR, Cholesky, Eigen, and SVD.
-   **Linear Solvers**: Direct (Gaussian, Cholesky) and Iterative (CG, GMRES) solvers.
-   **Sparse Matrix Support**: Efficient CSR/CSC storage and operations for large-scale data.
-   **Analysis**: PCA, rank computation, and stability analysis.
-   **Optimization**: Integrated support for BLAS and OpenMP for parallel execution.

*For detailed documentation, refer to the [Linear Algebra README](lin_alg/README.md).*

### 3. Statistics Analysis Library (`statistics/`)

A comprehensive C++ library for statistical analysis covering both univariate and bivariate statistics. Features an interactive CLI tool for:

-   **Univariate Statistics**: Mean, median, mode, variance, standard deviation, skewness, kurtosis, quantiles, and frequency distributions.
-   **Bivariate Statistics**: Pearson/Spearman/Kendall correlation, covariance, and linear regression with predictions.
-   **Hypothesis Testing**: One-sample, two-sample, and paired t-tests with confidence intervals and p-values.
-   **Data Transformations**: Z-score standardization, min-max normalization, trimmed mean, and mathematical transforms.
-   **File I/O**: CSV import/export with header support.
-   **Educational Mode**: Verbose toggle for step-by-step calculation insights.

*For detailed documentation, refer to the [Statistics README](statistics/README.md).*

### 4. Database Management (`mysql/`)

A series of advanced SQL scripts designed for the Sakila sample database, demonstrating complex business logic and analytics.

-   **`01_advanced_sakila_queries.sql`**: Complex joins, multi-level aggregations, and reporting queries.
-   **`02_window_functions_analytics.sql`**: Advanced analytics using Window Functions and Common Table Expressions (CTEs).
-   **`03_procedures_and_triggers.sql`**: Stored procedures and triggers for automating business logic and maintaining data integrity.
-   **`sakila_insights.sql`**: Comprehensive analytical reports on database contents.

### 5. Probability & Combinatorics Library (`probability/`)

A comprehensive C++ toolkit for combinatorial mathematics and probability theory with verbose logging capabilities. Features an interactive menu for:

-   **Combinatorics**: Permutations, combinations, number sequences (Fibonacci, Catalan, Lucas), derangements, and partitions.
-   **Probability**: Basic probability, Bayes' theorem, conditional probability, union/intersection (both independent and dependent), total probability, expected value, and variance.
-   **Distributions**: Full statistical support for discrete (Binomial, Poisson, Geometric, Hypergeometric, Negative Binomial) and continuous (Uniform, Normal, Exponential, Gamma, Beta, Chi-Square, Student's t, F-Distribution) distributions with PDF/PMF, CDF, Mean, Variance, and Standard Deviation.
-   **Verbose Mode**: Toggle step-by-step calculation display showing formulas, parameter substitution, and intermediate results for educational purposes.

*For detailed documentation, refer to the [Probability README](probability/README.md).*

## Usage

### Compiling C Programs
To compile and execute a specific C program:

```bash
gcc c/05_advanced_algorithms/merge_sort.c -o merge_sort
./merge_sort
```

### Running the Linear Algebra Tool
Navigate to the `lin_alg` directory and build the project:

```bash
cd lin_alg
make
./linear_algebra
```

### Running the Statistics Tool
Navigate to the `statistics` directory and build the project:

```bash
cd statistics
make
./statistics
```

### Running the Probability Tool
Navigate to the `probability` directory and build the project:

```bash
cd probability
make
./bin/probability
```

### Executing SQL Scripts
Ensure a MySQL instance with the Sakila schema is running. Source the scripts using the MySQL command line client or a GUI tool like Workbench.

```sql
source mysql/01_advanced_sakila_queries.sql;
```

---
*Maintained by Muktadir for academic coursework.*
