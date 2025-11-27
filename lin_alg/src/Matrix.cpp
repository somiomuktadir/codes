#include "Matrix.h"
#include "Decomposer.h"
#include <fstream>
#include <sstream>
#include <algorithm>

// Optional BLAS support
#ifdef USE_BLAS
extern "C" {
    void cblas_dgemm(const int Order, const int TransA, const int TransB,
                     const int M, const int N, const int K,
                     const double alpha, const double *A, const int lda,
                     const double *B, const int ldb,
                     const double beta, double *C, const int ldc);
}
#define CblasRowMajor 101
#define CblasNoTrans 111
#endif

namespace LinAlg {

// Constructors
Matrix::Matrix(int r, int c, double initialValue) : rows(r), cols(c) {
    data.resize(r * c, initialValue);
}

Matrix::Matrix(const std::vector<std::vector<double>>& d) {
    if (d.empty()) {
        rows = 0;
        cols = 0;
    } else {
        rows = d.size();
        cols = d[0].size();
        data.resize(rows * cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                data[i * cols + j] = d[i][j];
            }
        }
    }
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) {
    rows = list.size();
    cols = (rows > 0) ? list.begin()->size() : 0;
    data.resize(rows * cols);
    
    int i = 0;
    for (const auto& row : list) {
        int j = 0;
        for (double val : row) {
            data[i * cols + j] = val;
            ++j;
        }
        ++i;
    }
}

// Accessors
int Matrix::getRows() const {
    return rows;
}

int Matrix::getCols() const {
    return cols;
}

double& Matrix::at(int r, int c) {
    if (r < 0 || r >= rows || c < 0 || c >= cols) {
        throw std::out_of_range("Matrix index out of bounds");
    }
    return data[r * cols + c];
}

const double& Matrix::at(int r, int c) const {
    if (r < 0 || r >= rows || c < 0 || c >= cols) {
        throw std::out_of_range("Matrix index out of bounds");
    }
    return data[r * cols + c];
}

double& Matrix::operator()(int r, int c) {
    return data[r * cols + c];
}

const double& Matrix::operator()(int r, int c) const {
    return data[r * cols + c];
}

// Comparison operators
bool Matrix::operator==(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) return false;
    for (int i = 0; i < rows * cols; ++i) {
        if (std::abs(data[i] - other.data[i]) > epsilon()) return false;
    }
    return true;
}

bool Matrix::operator!=(const Matrix& other) const {
    return !(*this == other);
}

// Arithmetic operations
Matrix Matrix::operator+(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; ++i) {
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction");
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; ++i) {
        result.data[i] = data[i] - other.data[i];
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols != other.rows) {
        throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
    }
    Matrix result(rows, other.cols);
    
#ifdef USE_BLAS
    // Use BLAS cblas_dgemm for matrix multiplication: C = alpha*A*B + beta*C
    // result = 1.0 * this * other + 0.0 * result
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                rows, other.cols, cols,
                1.0, data.data(), cols,
                other.data.data(), other.cols,
                0.0, result.data.data(), result.cols);
#else
    // Native implementation with optional OpenMP parallelization
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(rows > 100 && other.cols > 100)
    #endif
    for (int i = 0; i < rows; ++i) {
        for (int k = 0; k < cols; ++k) {
            double aik = data[i * cols + k];
            for (int j = 0; j < other.cols; ++j) {
                result.data[i * other.cols + j] += aik * other.data[k * other.cols + j];
            }
        }
    }
#endif
    
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; ++i) {
        result.data[i] = data[i] * scalar;
    }
    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(cols, rows);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[j * rows + i] = data[i * cols + j];
        }
    }
    return result;
}

// Element-wise operations
Matrix Matrix::hadamard(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for Hadamard product");
    }
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; ++i) {
        result.data[i] = data[i] * other.data[i];
    }
    return result;
}

Matrix Matrix::applyFunction(std::function<double(double)> func) const {
    Matrix result(rows, cols);
    for (int i = 0; i < rows * cols; ++i) {
        result.data[i] = func(data[i]);
    }
    return result;
}

// Matrix manipulation
void Matrix::resize(int newRows, int newCols, double fillValue) {
    std::vector<double> newData(newRows * newCols, fillValue);
    int minRows = std::min(rows, newRows);
    int minCols = std::min(cols, newCols);
    
    for (int i = 0; i < minRows; ++i) {
        for (int j = 0; j < minCols; ++j) {
            newData[i * newCols + j] = data[i * cols + j];
        }
    }
    
    data = std::move(newData);
    rows = newRows;
    cols = newCols;
}

Matrix Matrix::submatrix(int startRow, int startCol, int numRows, int numCols) const {
    if (startRow < 0 || startCol < 0 || 
        startRow + numRows > rows || startCol + numCols > cols) {
        throw std::out_of_range("Submatrix indices out of bounds");
    }
    
    Matrix result(numRows, numCols);
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            result.data[i * numCols + j] = data[(startRow + i) * cols + (startCol + j)];
        }
    }
    return result;
}

Matrix Matrix::hstack(const Matrix& A, const Matrix& B) {
    if (A.rows != B.rows) {
        throw std::invalid_argument("Matrices must have same number of rows for hstack");
    }
    
    Matrix result(A.rows, A.cols + B.cols);
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            result.data[i * result.cols + j] = A.data[i * A.cols + j];
        }
        for (int j = 0; j < B.cols; ++j) {
            result.data[i * result.cols + A.cols + j] = B.data[i * B.cols + j];
        }
    }
    return result;
}

Matrix Matrix::vstack(const Matrix& A, const Matrix& B) {
    if (A.cols != B.cols) {
        throw std::invalid_argument("Matrices must have same number of columns for vstack");
    }
    
    Matrix result(A.rows + B.rows, A.cols);
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.cols; ++j) {
            result.data[i * result.cols + j] = A.data[i * A.cols + j];
        }
    }
    for (int i = 0; i < B.rows; ++i) {
        for (int j = 0; j < B.cols; ++j) {
            result.data[(A.rows + i) * result.cols + j] = B.data[i * B.cols + j];
        }
    }
    return result;
}

// File I/O
void Matrix::saveCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    file << std::setprecision(15);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file << data[i * cols + j];
            if (j < cols - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
}

Matrix Matrix::loadCSV(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for reading: " + filename);
    }
    
    std::vector<std::vector<double>> tempData;
    std::string line;
    
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;
        
        while (std::getline(ss, cell, ',')) {
            row.push_back(std::stod(cell));
        }
        
        if (!row.empty()) {
            tempData.push_back(row);
        }
    }
    
    file.close();
    
    if (tempData.empty()) {
        throw std::runtime_error("Empty CSV file or invalid format");
    }
    
    return Matrix(tempData);
}

// Utilities
Matrix Matrix::identity(int n) {
    Matrix result(n, n);
    for (int i = 0; i < n; ++i) {
        result.data[i * n + i] = 1.0;
    }
    return result;
}

void Matrix::print() const {
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < rows; ++i) {
        std::cout << "[ ";
        for (int j = 0; j < cols; ++j) {
            std::cout << data[i * cols + j] << " ";
        }
        std::cout << "]" << std::endl;
    }
}

Matrix operator*(double scalar, const Matrix& mat) {
    return mat * scalar;
}

double Matrix::trace() const {
    if (rows != cols) throw std::invalid_argument("Matrix must be square for trace");
    double sum = 0;
    for (int i = 0; i < rows; ++i) sum += data[i * cols + i];
    return sum;
}

int Matrix::rank() const {
    // Gaussian elimination to row echelon form
    Matrix M = *this;
    int pivotRow = 0;
    for (int col = 0; col < cols && pivotRow < rows; ++col) {
        int sel = pivotRow;
        for (int i = pivotRow + 1; i < rows; ++i) {
            if (std::abs(M.data[i * cols + col]) > std::abs(M.data[sel * cols + col])) {
                sel = i;
            }
        }
        
        if (std::abs(M.data[sel * cols + col]) < epsilon()) continue;
        
        // Swap rows
        for (int j = 0; j < cols; ++j) {
            std::swap(M.data[sel * cols + j], M.data[pivotRow * cols + j]);
        }
        
        for (int i = pivotRow + 1; i < rows; ++i) {
            double factor = M.data[i * cols + col] / M.data[pivotRow * cols + col];
            for (int j = col; j < cols; ++j) {
                M.data[i * cols + j] -= factor * M.data[pivotRow * cols + j];
            }
        }
        pivotRow++;
    }
    return pivotRow;
}

Matrix Matrix::rref() const {
    Matrix M = *this;
    int lead = 0;
    int rowCount = rows;
    int colCount = cols;

    for (int r = 0; r < rowCount; ++r) {
        if (colCount <= lead) break;
        int i = r;
        while (M(i, lead) == 0) {
            i++;
            if (rowCount == i) {
                i = r;
                lead++;
                if (colCount == lead) return M;
            }
        }

        for (int k = 0; k < colCount; ++k) {
            std::swap(M(i, k), M(r, k));
        }

        double val = M(r, lead);
        for (int k = 0; k < colCount; ++k) {
            M(r, k) /= val;
        }

        for (int i = 0; i < rowCount; ++i) {
            if (i != r) {
                double val = M(i, lead);
                for (int k = 0; k < colCount; ++k) {
                    M(i, k) -= val * M(r, k);
                }
            }
        }
        lead++;
    }
    return M;
}

double Matrix::conditionNumber() const {
    // Condition number is ratio of max singular value to min singular value
    // We can use SVD for this.
    try {
        auto [U, S, V] = Decomposer::SVD(*this);
        double maxSv = 0;
        double minSv = 1e300; // Large number
        
        int minDim = std::min(rows, cols);
        for(int i=0; i<minDim; ++i) {
            double s = S(i,i);
            if(s > maxSv) maxSv = s;
            if (s < minSv && s > epsilon()) minSv = s; // Avoid zero singular values for singular matrices
        }
        
        if (minSv > 1e299) return 1.0/0.0; // Infinity (singular)
        if (minSv < epsilon()) return 1.0/0.0;
        
        return maxSv / minSv;
    } catch (...) {
        return -1.0; // Error
    }
}

} // namespace LinAlg
