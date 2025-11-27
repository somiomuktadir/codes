#ifndef MATRIX_MARKET_H
#define MATRIX_MARKET_H

#include "Matrix.h"
#include "SparseMatrix.h"
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace LinAlg {

// Matrix Market (.mtx) file format support
// Compliant with NIST Matrix Market specification
class MatrixMarket {
public:
    enum class Format { COORDINATE, ARRAY };
    enum class Field { REAL, INTEGER, COMPLEX, PATTERN };
    enum class Symmetry { GENERAL, SYMMETRIC, SKEW_SYMMETRIC, HERMITIAN };
    
    // Load matrix from Matrix Market file
    // Returns dense matrix (can convert to sparse if needed)
    static Matrix loadDense(const std::string& filename);
    
    // Load sparse matrix from Matrix Market file (coordinate format only)
    static SparseMatrixCSR loadSparse(const std::string& filename);
    
    // Save dense matrix to Matrix Market file (array format)
    static void saveDense(const std::string& filename, const Matrix& mat);
    
    // Save sparse matrix to Matrix Market file (coordinate format)
    static void saveSparse(const std::string& filename, const SparseMatrixCSR& mat);
    
    // Query file format without loading
    static Format getFormat(const std::string& filename);
    
private:
    // Parse header line
    static void parseHeader(const std::string& line, Format& format, 
                           Field& field, Symmetry& symmetry);
    
    // Skip comment lines (starting with %)
    static void skipComments(std::ifstream& file);
    
    // Write header
    static void writeHeader(std::ofstream& file, Format format, 
                           Field field, Symmetry symmetry);
};

} // namespace LinAlg

#endif // MATRIX_MARKET_H
