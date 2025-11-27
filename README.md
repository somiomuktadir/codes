# University Code Repository

This repository serves as a comprehensive archive of programming assignments, algorithms, and system implementations developed for university-level computer science coursework. It encompasses a variety of domains including system programming in C, high-performance numerical computing in C++, and advanced database management with MySQL.

## Repository Structure

The codebase is organized into three primary modules, each focusing on a specific technical domain.

### 1. System Programming & Algorithms (`c/`)

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

### 3. Database Management (`mysql/`)

A series of advanced SQL scripts designed for the Sakila sample database, demonstrating complex business logic and analytics.

-   **`01_advanced_sakila_queries.sql`**: Complex joins, multi-level aggregations, and reporting queries.
-   **`02_window_functions_analytics.sql`**: Advanced analytics using Window Functions and Common Table Expressions (CTEs).
-   **`03_procedures_and_triggers.sql`**: Stored procedures and triggers for automating business logic and maintaining data integrity.
-   **`sakila_insights.sql`**: Comprehensive analytical reports on database contents.

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

### Executing SQL Scripts
Ensure a MySQL instance with the Sakila schema is running. Source the scripts using the MySQL command line client or a GUI tool like Workbench.

```sql
source mysql/01_advanced_sakila_queries.sql;
```

---
*Maintained by Muktadir for academic coursework.*
