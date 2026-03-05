#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

template <typename T>
class Complex {
 public:
  T re;
  T im;

  Complex() : re(T(0)), im(T(0)) {}
  Complex(T r, T i) : re(r), im(i) {}

  Complex operator+(const Complex &other) const {
    return {re + other.re, im + other.im};
  }

  Complex operator-(const Complex &other) const {
    return {re - other.re, im - other.im};
  }

  Complex operator*(const Complex &other) const {
    return {(re * other.re) - (im * other.im), (re * other.im) + (im * other.re)};
  }

  Complex &operator+=(const Complex &other) {
    re += other.re;
    im += other.im;
    return *this;
  }

  bool operator==(const Complex &other) const {
    constexpr double kEps = 1e-9;
    return std::abs(re - other.re) < kEps && std::abs(im - other.im) < kEps;
  }

  bool operator!=(const Complex &other) const {
    return !(*this == other);
  }
};

template <typename T>
class CSRMatrix {
 public:
  int rows;
  int cols;
  std::vector<Complex<T>> values;
  std::vector<int> col_indices;
  std::vector<int> row_ptr;

  CSRMatrix() : rows(0), cols(0), row_ptr(1, 0) {}

  CSRMatrix(int r, int c) : rows(r), cols(c), row_ptr(r + 1, 0) {}

  CSRMatrix(int r, int c, std::vector<Complex<T>> vals, std::vector<int> col_idx, std::vector<int> rp)
      : rows(r), cols(c), values(std::move(vals)), col_indices(std::move(col_idx)), row_ptr(std::move(rp)) {}

  bool operator==(const CSRMatrix &other) const {
    if (rows != other.rows || cols != other.cols) {
      return false;
    }
    if (row_ptr != other.row_ptr || col_indices != other.col_indices) {
      return false;
    }
    if (values.size() != other.values.size()) {
      return false;
    }
    for (size_t i = 0; i < values.size(); ++i) {
      if (values[i] != other.values[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const CSRMatrix &other) const {
    return !(*this == other);
  }

  [[nodiscard]] CSRMatrix Multiply(const CSRMatrix &other) const {
    if (cols != other.rows) {
      return {};
    }

    CSRMatrix result(rows, other.cols);

    std::vector<Complex<T>> row_acc(other.cols);
    std::vector<bool> row_used(other.cols, false);

    for (int i = 0; i < rows; ++i) {
      std::vector<int> used_cols;
      used_cols.reserve(other.cols);

      for (int ja = row_ptr[i]; ja < row_ptr[i + 1]; ++ja) {
        int ka = col_indices[ja];
        const Complex<T> &a_val = values[ja];

        for (int jb = other.row_ptr[ka]; jb < other.row_ptr[ka + 1]; ++jb) {
          int cb = other.col_indices[jb];
          const Complex<T> &b_val = other.values[jb];

          if (!row_used[cb]) {
            row_used[cb] = true;
            row_acc[cb] = Complex<T>();
            used_cols.push_back(cb);
          }
          row_acc[cb] += a_val * b_val;
        }
      }

      std::ranges::sort(used_cols);

      for (int c : used_cols) {
        result.values.push_back(row_acc[c]);
        result.col_indices.push_back(c);
        row_used[c] = false;
      }
      result.row_ptr[i + 1] = static_cast<int>(result.values.size());
    }

    return result;
  }

 private:
  static void ProcessRow(int i, const CSRMatrix &self, const CSRMatrix &other, std::vector<T> &acc_re,
                          std::vector<T> &acc_im, std::vector<bool> &local_used,
                          std::vector<std::vector<Complex<T>>> &row_values,
                          std::vector<std::vector<int>> &row_col_indices) {
    std::vector<int> used_cols;
    used_cols.reserve(other.cols);

    for (int ja = self.row_ptr[i]; ja < self.row_ptr[i + 1]; ++ja) {
      int ka = self.col_indices[ja];
      T a_re = self.values[ja].re;
      T a_im = self.values[ja].im;
      int jb_start = other.row_ptr[ka];
      int jb_end = other.row_ptr[ka + 1];

      for (int jb = jb_start; jb < jb_end; ++jb) {
        int cb = other.col_indices[jb];
        if (!local_used[cb]) {
          local_used[cb] = true;
          acc_re[cb] = T(0);
          acc_im[cb] = T(0);
          used_cols.push_back(cb);
        }
        acc_re[cb] += (a_re * other.values[jb].re) - (a_im * other.values[jb].im);
        acc_im[cb] += (a_re * other.values[jb].im) + (a_im * other.values[jb].re);
      }
    }

    std::ranges::sort(used_cols);

    row_values[i].reserve(used_cols.size());
    row_col_indices[i].reserve(used_cols.size());

    for (int c : used_cols) {
      row_values[i].emplace_back(acc_re[c], acc_im[c]);
      row_col_indices[i].push_back(c);
      local_used[c] = false;
    }
  }

 public:
  [[nodiscard]] CSRMatrix OMPMultiply(const CSRMatrix &other) const {
    if (cols != other.rows) {
      return {};
    }

    CSRMatrix result(rows, other.cols);
    std::vector<std::vector<Complex<T>>> row_values(rows);
    std::vector<std::vector<int>> row_col_indices(rows);
    int nrows = rows;
    int ncols = other.cols;

#pragma omp parallel default(none) shared(other, row_values, row_col_indices, nrows, ncols)
    {
      std::vector<T> acc_re(ncols, T(0));
      std::vector<T> acc_im(ncols, T(0));
      std::vector<bool> local_used(ncols, false);

#pragma omp for schedule(dynamic)
      for (int i = 0; i < nrows; ++i) {
        ProcessRow(i, *this, other, acc_re, acc_im, local_used, row_values, row_col_indices);
      }
    }

    int total_nnz = 0;
    for (int i = 0; i < rows; ++i) {
      total_nnz += static_cast<int>(row_values[i].size());
    }

    result.values.reserve(total_nnz);
    result.col_indices.reserve(total_nnz);

    for (int i = 0; i < rows; ++i) {
      result.values.insert(result.values.end(), row_values[i].begin(), row_values[i].end());
      result.col_indices.insert(result.col_indices.end(), row_col_indices[i].begin(), row_col_indices[i].end());
      result.row_ptr[i + 1] = static_cast<int>(result.values.size());
    }

    return result;
  }

  [[nodiscard]] std::vector<Complex<T>> ToDense() const {
    std::vector<Complex<T>> dense(rows * cols);
    for (int i = 0; i < rows; ++i) {
      for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
        dense[(i * cols) + col_indices[j]] = values[j];
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
