#pragma once

#include <cstddef>  // for size_t
#include <map>
#include <stdexcept>  // for std::invalid_argument
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

struct Complex {
  double re;  // real
  double im;  // imaginary

  explicit Complex(double re = 0.0, double im = 0.0) : re(re), im(im) {}

  bool operator==(const Complex &other) const {
    return other.re == re && other.im == im;
  }

  bool operator!=(const Complex &other) const {
    return !(other == *this);
  }

  Complex operator*(const Complex &other) const {
    return Complex{(re * other.re) - (im * other.im), (re * other.im) + (im * other.re)};
  }

  Complex operator+(const Complex &other) const {
    return Complex{re + other.re, im + other.im};
  }

  Complex &operator+=(const Complex &other) {
    re += other.re;
    im += other.im;
    return *this;
  }
};

struct SparseMatrix {
  std::vector<Complex> val;
  std::vector<size_t> row_ind;
  std::vector<size_t> col_ind;
  size_t height = 0;
  size_t width = 0;

  SparseMatrix() = default;

  explicit SparseMatrix(const std::vector<std::vector<Complex>> &matr)
      : height(matr.size()), width(matr.empty() ? 0 : matr[0].size()) {
    if (height == 0 || width == 0) {
      return;
    }

    for (size_t col = 0; col < width; ++col) {
      for (size_t row = 0; row < height; ++row) {
        if (matr[row][col] != Complex(0.0)) {
          val.push_back(matr[row][col]);
          row_ind.push_back(row);
          col_ind.push_back(col);
        }
      }
    }
  }

  [[nodiscard]] size_t Count() const {
    return val.size();
  }

  SparseMatrix operator*(const SparseMatrix &matr_b) const {
    if (width != matr_b.height) {
      throw std::invalid_argument("Incompatible matrix dimensions for multiplication");
    }

    std::map<std::pair<size_t, size_t>, Complex> mp;  // <row, col> -> val

    for (size_t i = 0; i < Count(); ++i) {
      size_t row_a = row_ind[i];
      size_t col_a = col_ind[i];
      Complex val_a = val[i];

      for (size_t j = 0; j < matr_b.Count(); ++j) {
        size_t row_b = matr_b.row_ind[j];
        size_t col_b = matr_b.col_ind[j];
        Complex val_b = matr_b.val[j];

        if (col_a == row_b) {
          mp[{row_a, col_b}] += val_a * val_b;
        }
      }
    }

    SparseMatrix res;
    res.width = matr_b.width;
    res.height = height;
    for (const auto &[key, value] : mp) {
      res.val.push_back(value);
      res.row_ind.push_back(key.first);
      res.col_ind.push_back(key.second);
    }

    return res;
  }
};

using InType = std::tuple<SparseMatrix, SparseMatrix>;
using OutType = SparseMatrix;
using TestType = std::tuple<size_t, size_t, size_t>;  // n, m, k. Matrix_1: n*m, Matrix_2: m*k
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zavyalov_a_compl_sparse_matr_mult
