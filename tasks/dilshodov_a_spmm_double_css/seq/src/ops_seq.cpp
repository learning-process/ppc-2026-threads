#include "dilshodov_a_spmm_double_css/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <map>
#include <utility>

#include "dilshodov_a_spmm_double_css/common/include/common.hpp"

namespace dilshodov_a_spmm_double_css {

namespace {
constexpr double kEps = 1e-10;

bool HasValidDimensions(const SparseMatrixCCS &m) {
  return m.rows_count > 0 && m.cols_count > 0;
}

bool HasValidContainers(const SparseMatrixCCS &m) {
  if (m.col_ptrs.size() != static_cast<size_t>(m.cols_count) + 1) {
    return false;
  }
  if (m.row_indices.size() != m.values.size()) {
    return false;
  }
  if (m.col_ptrs.empty() || m.col_ptrs.front() != 0) {
    return false;
  }
  if (std::cmp_not_equal(m.col_ptrs.back(), m.values.size())) {
    return false;
  }

  return true;
}

bool HasValidColumnOrdering(const SparseMatrixCCS &m) {
  for (int j = 0; j < m.cols_count; ++j) {
    if (m.col_ptrs[j] > m.col_ptrs[j + 1]) {
      return false;
    }
    int prev_row = -1;
    for (int idx = m.col_ptrs[j]; idx < m.col_ptrs[j + 1]; ++idx) {
      const int row = m.row_indices[idx];
      if (row < 0 || row >= m.rows_count) {
        return false;
      }
      if (row <= prev_row) {
        return false;
      }
      prev_row = row;
    }
  }

  return true;
}

bool IsValidCCS(const SparseMatrixCCS &m) {
  return HasValidDimensions(m) && HasValidContainers(m) && HasValidColumnOrdering(m);
}

}  // namespace

DilshodovASpmmDoubleCssSeq::DilshodovASpmmDoubleCssSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DilshodovASpmmDoubleCssSeq::ValidationImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  return IsValidCCS(matrix_a) && IsValidCCS(matrix_b) && matrix_a.cols_count == matrix_b.rows_count;
}

bool DilshodovASpmmDoubleCssSeq::PreProcessingImpl() {
  GetOutput() = SparseMatrixCCS{};
  return true;
}

bool DilshodovASpmmDoubleCssSeq::RunImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  matrix_c.rows_count = matrix_a.rows_count;
  matrix_c.cols_count = matrix_b.cols_count;
  matrix_c.col_ptrs.assign(static_cast<size_t>(matrix_c.cols_count) + 1, 0);
  matrix_c.row_indices.clear();
  matrix_c.values.clear();

  for (int col_b = 0; col_b < matrix_b.cols_count; ++col_b) {
    std::map<int, double> accumulator;

    for (int idx_b = matrix_b.col_ptrs[col_b]; idx_b < matrix_b.col_ptrs[col_b + 1]; ++idx_b) {
      const int pivot_row = matrix_b.row_indices[idx_b];
      const double pivot_value = matrix_b.values[idx_b];

      for (int idx_a = matrix_a.col_ptrs[pivot_row]; idx_a < matrix_a.col_ptrs[pivot_row + 1]; ++idx_a) {
        const int result_row = matrix_a.row_indices[idx_a];
        accumulator[result_row] += matrix_a.values[idx_a] * pivot_value;
      }
    }

    for (const auto &[row, value] : accumulator) {
      if (std::abs(value) > kEps) {
        matrix_c.row_indices.push_back(row);
        matrix_c.values.push_back(value);
      }
    }

    matrix_c.col_ptrs[col_b + 1] = static_cast<int>(matrix_c.values.size());
  }

  matrix_c.non_zeros = static_cast<int>(matrix_c.values.size());
  return true;
}

bool DilshodovASpmmDoubleCssSeq::PostProcessingImpl() {
  GetOutput().non_zeros = static_cast<int>(GetOutput().values.size());
  return true;
}

}  // namespace dilshodov_a_spmm_double_css
