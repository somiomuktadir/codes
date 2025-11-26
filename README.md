# C and SQL Code Repository

This repository contains a collection of C codes, C++ programs, and SQL scripts designed for a university-level course.

## Directory Structure

### C Programming (`c/`)

The C programs are organized into the following categories:

- **`01_basics_revision/`**: Fundamental concepts, loops, arrays, strings, and basic algorithms.
- **`02_pointers_and_memory/`**: Pointer arithmetic, function pointers, dynamic memory allocation, and memory management.
- **`03_data_structures/`**: Implementations of Linked Lists, Stacks, Queues, Trees, and Hash Tables.
- **`04_file_io_and_system/`**: File handling, CSV parsing, logging systems, and system calls.
- **`05_advanced_algorithms/`**: Sorting algorithms (Merge Sort, Quick Sort), Graph algorithms (Dijkstra), and Dynamic Programming.
- **`hello.c`**: Basic Hello World program.
- **`06_data_analysis/`**: Data analysis programs and datasets.

### Linear Algebra Library (`lin_alg/`)

A C++ library for linear algebra operations:

- **`Matrix.cpp` / `Matrix.h`**: Matrix class implementation.
- **`VectorOps.h`**: Vector operations.
- **`LinearSolver.cpp` / `LinearSolver.h`**: Solvers for linear systems.
- **`Decomposer.cpp` / `Decomposer.h`**: Matrix decomposition algorithms (LU, QR, etc.).
- **`Analysis.cpp` / `Analysis.h`**: Data analysis tools.
- **`main.cpp`**: Entry point for testing and usage.

### MySQL Scripts (`mysql/`)

The SQL scripts demonstrate advanced database concepts using the Sakila sample database:

- **`01_advanced_sakila_queries.sql`**: Complex joins, aggregations, and reporting queries.
- **`02_window_functions_analytics.sql`**: Analytical queries using Window Functions and CTEs.
- **`03_procedures_and_triggers.sql`**: Stored procedures and triggers for business logic.
- **`sakila_insights.sql`**: Deep dive insights into the Sakila database.
- **`exercise.sql`**: Practice exercises.
- **`table_data.sql`**: Scripts for populating table data.
- **`another.sql`**: Additional SQL examples.

## Usage

To compile and run a C program:

```bash
gcc c/01_basics_revision/prime_generator.c -o prime_generator
./prime_generator
```

To run SQL scripts, you will need a MySQL database with the Sakila schema installed.
