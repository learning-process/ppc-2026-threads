#pragma once

#include <complex>
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
    row_ptr.resize(r + 1, 0);
  }

  bool IsValid() const {
    if (rows < 0 || cols < 0) {
      return false;
    }
    if (row_ptr.size() != static_cast<size_t>(rows + 1)) {
      return false;
    }
    if (col_indices.size() != values.size()) {
      return false;
    }

    for (int i = 0; i < rows; ++i) {
      if (row_ptr[i] > row_ptr[i + 1]) {
        return false;
      }
    }

    for (size_t i = 0; i < col_indices.size(); ++i) {
      if (col_indices[i] < 0 || col_indices[i] >= cols) {
        return false;
      }
    }

    for (int i = 0; i < rows; ++i) {
      for (int j = row_ptr[i]; j < row_ptr[i + 1] - 1; ++j) {
        if (col_indices[j] >= col_indices[j + 1]) {
          return false;
        }
      }
    }

    return true;
  }

  size_t NonZeros() const {
    return values.size();
  }
};

using InType = std::tuple<CRSMatrix, CRSMatrix>;
using OutType = CRSMatrix;
using TestType = std::tuple<int, int, int, int, int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

constexpr double EPSILON = 1e-14;

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
