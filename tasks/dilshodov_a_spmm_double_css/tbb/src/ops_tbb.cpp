#include "dilshodov_a_spmm_double_css/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

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
  if (m.col_ptrs.back() < 0) {
    return false;
  }
  if (static_cast<size_t>(m.col_ptrs.back()) != m.values.size()) {
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

struct ScratchData {
  std::vector<double> acc;
  std::vector<int> marker;
  std::vector<int> used_rows;
};

void AccumulateColumnProduct(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs, int rhs_col, ScratchData &scratch,
                             std::vector<std::pair<int, double>> &col) {
  scratch.used_rows.clear();
  col.clear();

  for (int idx_b = rhs.col_ptrs[rhs_col]; idx_b < rhs.col_ptrs[rhs_col + 1]; ++idx_b) {
    const int k = rhs.row_indices[idx_b];
    const double rhs_value = rhs.values[idx_b];

    for (int idx_a = lhs.col_ptrs[k]; idx_a < lhs.col_ptrs[k + 1]; ++idx_a) {
      const int row = lhs.row_indices[idx_a];
      const double value = lhs.values[idx_a] * rhs_value;

      if (scratch.marker[row] != rhs_col) {
        scratch.marker[row] = rhs_col;
        scratch.acc[row] = value;
        scratch.used_rows.push_back(row);
      } else {
        scratch.acc[row] += value;
      }
    }
  }

  col.reserve(scratch.used_rows.size());
  for (int row : scratch.used_rows) {
    const double value = scratch.acc[row];
    if (std::abs(value) > kEps) {
      col.emplace_back(row, value);
    }
  }

  std::ranges::sort(col, {}, &std::pair<int, double>::first);
}

void BuildOutputFromColumns(const std::vector<std::vector<std::pair<int, double>>> &column_results,
                            SparseMatrixCCS &out) {
  int offset = 0;
  for (int col = 0; col < out.cols_count; ++col) {
    for (const auto &[row, value] : column_results[col]) {
      out.row_indices.push_back(row);
      out.values.push_back(value);
      ++offset;
    }
    out.col_ptrs[col + 1] = offset;
  }
}

}  // namespace

DilshodovASpmmDoubleCssTbb::DilshodovASpmmDoubleCssTbb(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DilshodovASpmmDoubleCssTbb::ValidationImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  return IsValidCCS(lhs) && IsValidCCS(rhs) && lhs.cols_count == rhs.rows_count;
}

bool DilshodovASpmmDoubleCssTbb::PreProcessingImpl() {
  GetOutput() = SparseMatrixCCS{};
  return true;
}

bool DilshodovASpmmDoubleCssTbb::RunImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  auto &out = GetOutput();

  out.rows_count = lhs.rows_count;
  out.cols_count = rhs.cols_count;
  out.col_ptrs.assign(static_cast<size_t>(out.cols_count) + 1, 0);
  out.row_indices.clear();
  out.values.clear();

  std::vector<std::vector<std::pair<int, double>>> column_results(rhs.cols_count);

  ScratchData exemplar;
  exemplar.acc.assign(lhs.rows_count, 0.0);
  exemplar.marker.assign(lhs.rows_count, -1);

  tbb::enumerable_thread_specific<ScratchData> tls(exemplar);

  tbb::parallel_for(tbb::blocked_range<int>(0, rhs.cols_count), [&](const tbb::blocked_range<int> &range) {
    auto &scratch = tls.local();
    for (int rhs_col = range.begin(); rhs_col < range.end(); ++rhs_col) {
      AccumulateColumnProduct(lhs, rhs, rhs_col, scratch, column_results[rhs_col]);
    }
  });

  BuildOutputFromColumns(column_results, out);

  out.non_zeros = static_cast<int>(out.values.size());
  return true;
}

bool DilshodovASpmmDoubleCssTbb::PostProcessingImpl() {
  GetOutput().non_zeros = static_cast<int>(GetOutput().values.size());
  return true;
}

}  // namespace dilshodov_a_spmm_double_css
