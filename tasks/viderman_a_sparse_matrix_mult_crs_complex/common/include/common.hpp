#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

using Complex = std::complex<double>;

struct CRSMatrix {
  int rows = 0;
  int cols = 0;
  std::vector<int> row_ptr;
  std::vector<int> col_indices;
  std::vector<Complex> values;

  CRSMatrix() = default;

  CRSMatrix(int r, int c) : rows(r), cols(c) {
    if (r >= 0) {
      row_ptr.assign(r + 1, 0);
    }
  }

  bool IsValid() const {
    if (rows < 0 || cols < 0) {
      return false;
    }

    if (static_cast<int>(row_ptr.size()) != rows + 1) {
      return false;
    }

    if (row_ptr[0] != 0) {
      return false;
    }

    if (col_indices.size() != values.size()) {
      return false;
    }

    const int nnz = static_cast<int>(values.size());

    if (row_ptr.back() != nnz) {
      return false;
    }

    for (int i = 0; i < rows; ++i) {
      if (row_ptr[i] > row_ptr[i + 1]) {
        return false;
      }

      if (row_ptr[i] < 0 || row_ptr[i] > nnz) {
        return false;
      }
    }

    for (int i = 0; i < rows; ++i) {
      for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
        int col = col_indices[j];

        if (col < 0 || col >= cols) {
          return false;
        }

        if (j + 1 < row_ptr[i + 1]) {
          if (col_indices[j] >= col_indices[j + 1]) {
            return false;
          }
        }
      }
    }

    return true;
  }

  [[nodiscard]] std::size_t NonZeros() const {
    return values.size();
  }
};

using InType = std::tuple<CRSMatrix, CRSMatrix>;
using OutType = CRSMatrix;
using BaseTask = ppc::task::Task<InType, OutType>;

constexpr double kEpsilon = 1e-14;

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
