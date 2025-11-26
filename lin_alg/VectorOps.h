#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>

namespace VectorOps {

    inline double dot(const std::vector<double>& v1, const std::vector<double>& v2) {
        if (v1.size() != v2.size()) {
            throw std::invalid_argument("Vector dimensions must match for dot product");
        }
        double result = 0.0;
        for (size_t i = 0; i < v1.size(); ++i) {
            result += v1[i] * v2[i];
        }
        return result;
    }

    inline double norm(const std::vector<double>& v) {
        return std::sqrt(dot(v, v));
    }

    inline std::vector<double> normalize(const std::vector<double>& v) {
        double n = norm(v);
        if (n == 0) return v;
        std::vector<double> result = v;
        for (double& val : result) {
            val /= n;
        }
        return result;
    }

    inline std::vector<double> cross(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != 3 || b.size() != 3) {
            throw std::invalid_argument("Cross product only defined for 3D vectors");
        }
        return {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        };
    }
    
    inline void print(const std::vector<double>& v) {
        std::cout << "( ";
        for (double val : v) {
            std::cout << val << " ";
        }
        std::cout << ")" << std::endl;
    }
}

#endif // VECTOR_OPS_H
