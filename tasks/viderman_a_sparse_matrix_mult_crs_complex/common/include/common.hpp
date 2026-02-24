#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <ranges>
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

  [[nodiscard]] bool is_valid() const {
    return has_valid_shape() && has_valid_row_ptr_size() && has_matching_value_sizes() && has_monotonic_row_ptr() &&
           has_valid_col_indices() && has_sorted_rows();
  }

  [[nodiscard]] size_t non_zeros() const {
    return values.size();
  }

 private:
  [[nodiscard]] bool has_valid_shape() const {
    return rows >= 0 && cols >= 0;
  }

  [[nodiscard]] bool has_valid_row_ptr_size() const {
    const auto expected = static_cast<std::size_t>(rows) + 1;
    return row_ptr.size() == expected;
  }

  [[nodiscard]] bool has_matching_value_sizes() const {
    return col_indices.size() == values.size();
  }

  [[nodiscard]] bool has_monotonic_row_ptr() const {
    for (int i = 0; i < rows; ++i) {
      if (row_ptr[i] > row_ptr[i + 1]) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] bool has_valid_col_indices() const {
    return std::ranges::all_of(col_indices, [this](int col) { return col >= 0 && col < cols; });
  }

  [[nodiscard]] bool has_sorted_rows() const {
    for (int i = 0; i < rows; ++i) {
      for (int j = row_ptr[i]; j < row_ptr[i + 1] - 1; ++j) {
        if (col_indices[j] >= col_indices[j + 1]) {
          return false;
        }
      }
    }
    return true;
  }
};

using InType = std::tuple<CRSMatrix, CRSMatrix>;
using OutType = CRSMatrix;
using TestType = std::tuple<int, int, int, int, int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

constexpr double kEpsilon = 1e-14;

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
