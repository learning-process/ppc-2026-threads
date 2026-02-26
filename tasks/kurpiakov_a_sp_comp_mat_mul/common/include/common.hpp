#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

template <typename T>
class Complex {
 public:
  T _re;
  T _im;

  Complex() : _re(T(0)), _im(T(0)) {}
  Complex(T r, T i) : _re(r), _im(i) {}

  Complex operator+(const Complex& other) const {
    return {_re + other._re, _im + other._im};
  }

  Complex operator-(const Complex& other) const {
    return {_re - other._re, _im - other._im};
  }

  Complex operator*(const Complex& other) const {
    return {_re * other._re - _im * other._im, _re * other._im + _im * other._re};
  }

  Complex& operator+=(const Complex& other) {
    _re += other._re;
    _im += other._im;
    return *this;
  }

  bool operator==(const Complex& other) const {
    constexpr double kEps = 1e-9;
    return std::abs(_re - other._re) < kEps && std::abs(_im - other._im) < kEps;
  }

  bool operator!=(const Complex& other) const {
    return !(*this == other);
  }
};

template <typename T>
class CSRMatrix {
 public:
  int _rows;
  int _cols;
  std::vector<Complex<T>> _values;
  std::vector<int> _col_indices;
  std::vector<int> _row_ptr;

  CSRMatrix() : _rows(0), _cols(0), _row_ptr(1, 0) {}

  CSRMatrix(int r, int c) : _rows(r), _cols(c), _row_ptr(r + 1, 0) {}

  CSRMatrix(int r, int c, std::vector<Complex<T>> vals, std::vector<int> col_idx, std::vector<int> rp)
      : _rows(r), _cols(c), _values(std::move(vals)), _col_indices(std::move(col_idx)), _row_ptr(std::move(rp)) {}

  bool operator==(const CSRMatrix& other) const {
    if (_rows != other._rows || _cols != other._cols) {
      return false;
    }
    if (_row_ptr != other._row_ptr || _col_indices != other._col_indices) {
      return false;
    }
    if (_values.size() != other._values.size()) {
      return false;
    }
    for (size_t i = 0; i < _values.size(); ++i) {
      if (_values[i] != other._values[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const CSRMatrix& other) const {
    return !(*this == other);
  }

  CSRMatrix Multiply(const CSRMatrix& other) const {
    if (_cols != other._rows) {
      return {};
    }

    CSRMatrix result(_rows, other._cols);

    std::vector<Complex<T>> row_acc(other._cols);
    std::vector<bool> row_used(other._cols, false);

    for (int i = 0; i < _rows; ++i) {
      std::vector<int> used_cols;
      used_cols.reserve(other._cols);

      for (int ja = _row_ptr[i]; ja < _row_ptr[i + 1]; ++ja) {
        int ka = _col_indices[ja];
        const Complex<T>& a_val = _values[ja];

        for (int jb = other._row_ptr[ka]; jb < other._row_ptr[ka + 1]; ++jb) {
          int cb = other._col_indices[jb];
          const Complex<T>& b_val = other._values[jb];

          if (!row_used[cb]) {
            row_used[cb] = true;
            row_acc[cb] = Complex<T>();
            used_cols.push_back(cb);
          }
          row_acc[cb] += a_val * b_val;
        }
      }

      std::sort(used_cols.begin(), used_cols.end());

      for (int c : used_cols) {
        result._values.push_back(row_acc[c]);
        result._col_indices.push_back(c);
        row_used[c] = false;
      }
      result._row_ptr[i + 1] = static_cast<int>(result._values.size());
    }

    return result;
  }

  std::vector<Complex<T>> ToDense() const {
    std::vector<Complex<T>> dense(_rows * _cols);
    for (int i = 0; i < _rows; ++i) {
      for (int j = _row_ptr[i]; j < _row_ptr[i + 1]; ++j) {
        dense[i * _cols + _col_indices[j]] = _values[j];
      }
    }
    return dense;
  }
};

using ComplexD = Complex<double>;
using SparseMatrix = CSRMatrix<double>;
using InType = std::pair<SparseMatrix, SparseMatrix>;
using OutType = SparseMatrix;
using TestType = std::tuple<InType, std::string, OutType>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kurpiakov_a_sp_comp_mat_mul
