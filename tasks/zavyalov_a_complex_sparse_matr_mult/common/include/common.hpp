#pragma once

#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

struct Complex {
  double re;  // real
  double im;  // imaginary
  Complex(double _re = 0.0, double _im = 0.0) {
    re = _re;
    im = _im;
  }
  bool operator==(Complex other) {
    return other.re == re && other.im == im;
  }
  bool operator!=(Complex other) {
    return !(other == *this);
  }
  Complex operator*(Complex other) {
    return Complex(re * other.re - im * other.im, re * other.im + im * other.re);
  }
  Complex operator+(Complex other) {
    return Complex(re + other.re, im + other.im);
  }
  Complex& operator+=(Complex other) {
    re += other.re;
    im += other.im;
    return *this;
  }
};

struct Sparse_matrix {
  std::vector<Complex> val;
  std::vector<size_t> row_ind;
  std::vector<size_t> col_ind;
  size_t height = 0;
  size_t width = 0;
  Sparse_matrix() {}
  Sparse_matrix(std::vector<std::vector<Complex>> matr) {
    height = matr.size();
    width = matr[0].size();
    for (size_t col = 0; col < width; col++) {
      for (size_t row = 0; row < height; row++) {
        if (matr[row][col] != Complex(0.0)) {
          val.push_back(matr[row][col]);
          row_ind.push_back(row);
          col_ind.push_back(col);
        }
      }
    }
  }
  size_t count() {
    return val.size();
  }
  Sparse_matrix operator*(Sparse_matrix& matr_b) {
    if (width != matr_b.height) {
      throw "Incompatible matrix dimensions for multiplication";
    }

    std::map<std::pair<size_t, size_t>, Complex> mp;  // <col_b, row_a> -> val_c

    for (size_t i = 0; i < count(); i++) {
      size_t row_a = row_ind[i];
      size_t col_a = col_ind[i];
      Complex val_a = val[i];

      for (size_t j = 0; j < matr_b.count(); j++) {
        size_t row_b = matr_b.row_ind[j];
        size_t col_b = matr_b.col_ind[j];
        Complex val_b = matr_b.val[j];

        if (col_a == row_b) {
          mp[{col_b, row_a}] += val_a * val_b;
        }
      }
    }

    Sparse_matrix res;
    res.width = matr_b.width;
    res.height = height;
    for (const auto& p : mp) {
      res.val.push_back(p.second);
      res.col_ind.push_back(p.first.first);
      res.row_ind.push_back(p.first.second);
    }

    return res;
  }
};

using InType = std::tuple<Sparse_matrix, Sparse_matrix>;
using OutType = Sparse_matrix;
using TestType = std::tuple<size_t, size_t, size_t>;  // n, m, k. Matrix_1: n*m, Matrix_2: m*k
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace zavyalov_a_compl_sparse_matr_mult
