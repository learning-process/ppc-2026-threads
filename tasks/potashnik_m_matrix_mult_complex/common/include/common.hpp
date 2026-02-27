#pragma once

#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace potashnik_m_matrix_mult_complex {

struct Complex {
  double real;
  double imaginary;

  explicit Complex(double real = 0.0, double imaginary = 0.0) {
    this->real = real;
    this->imaginary = imaginary;
  }

  bool operator==(const Complex &comp) const {
    return comp.real == real && comp.imaginary == imaginary;
  }

  bool operator!=(const Complex &comp) const {
    return !(comp == *this);
  }

  Complex operator*(const Complex &comp) const {
    double tmp_real = (real * comp.real) - (imaginary * comp.imaginary);
    double tmp_imaginary = (real * comp.imaginary) + (imaginary * comp.real);
    return Complex{tmp_real, tmp_imaginary};
  }

  Complex operator+(const Complex &comp) const {
    double tmp_real = real + comp.real;
    double tmp_imaginary = imaginary + comp.imaginary;
    return Complex{tmp_real, tmp_imaginary};
  }

  Complex &operator+=(const Complex &comp) {
    real += comp.real;
    imaginary += comp.imaginary;
    return *this;
  }
};

struct CCSMatrix {
  std::vector<Complex> val;
  std::vector<size_t> row_ind;
  std::vector<size_t> col_ptr;
  size_t height = 0;
  size_t width = 0;

  CCSMatrix() = default;

  explicit CCSMatrix(const std::vector<std::vector<Complex>> &matr) {
    height = matr.size();
    if (!matr.empty()) {
      width = matr[0].size();
    } else {
      width = 0;
    }

    if (height == 0 || width == 0) {
      return;
    }

    for (size_t col = 0; col < width; ++col) {
      for (size_t row = 0; row < height; ++row) {
        if (matr[row][col] != Complex(0.0)) {
          val.push_back(matr[row][col]);
          row_ind.push_back(row);
          col_ptr.push_back(col);
        }
      }
    }
  }

  [[nodiscard]] size_t Count() const {
    return val.size();
  }

  bool compare(const CCSMatrix &matr) {
    if (height != matr.height) {
      return false;
    }
    if (width != matr.height) {
      return false;
    }

    if (val.size() != matr.val.size()) {
      return false;
    }
    if (row_ind.size() != matr.row_ind.size()) {
      return false;
    }
    if (col_ptr.size() != matr.col_ptr.size()) {
      return false;
    }

    for (size_t i = 0; i < val.size(); i++) {
      if (val[i] != matr.val[i]) {
        return false;
      }
    }
    for (size_t i = 0; i < row_ind.size(); i++) {
      if (row_ind[i] != matr.row_ind[i]) {
        return false;
      }
    }
    for (size_t i = 0; i < col_ptr.size(); i++) {
      if (col_ptr[i] != matr.col_ptr[i]) {
        return false;
      }
    }

    return true;
  }
};

using InType = std::tuple<CCSMatrix, CCSMatrix>;
using OutType = CCSMatrix;
using TestType = std::tuple<size_t, size_t, size_t>;  // n, m, k. Matrix_1: n*m, Matrix_2: m*k
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace potashnik_m_matrix_mult_complex
