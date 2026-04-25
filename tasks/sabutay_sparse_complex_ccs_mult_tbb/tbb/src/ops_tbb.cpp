#include "../include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "../../common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"

namespace sabutay_sparse_complex_ccs_mult_tbb {

namespace {

constexpr double kEps = 1e-14;

void AccumulateColumnContributions(const CCS &a, const CCS &b, int col, std::vector<int> &rows,
                                   std::vector<int> &marker, std::vector<std::complex<double>> &acc) {
  for (int k = b.col_ptr[col]; k < b.col_ptr[col + 1]; ++k) {
    const std::complex<double> b_value = b.values[k];
    const int b_row = b.row_ind[k];

    for (int a_pos = a.col_ptr[b_row]; a_pos < a.col_ptr[b_row + 1]; ++a_pos) {
      const int a_row = a.row_ind[a_pos];
      acc[a_row] += b_value * a.values[a_pos];
      if (marker[a_row] == -1) {
        rows.push_back(a_row);
        marker[a_row] = 1;
      }
    }
  }
}

void FlushAccumulatedValues(const std::vector<int> &rows, std::vector<int> &marker,
                            std::vector<std::complex<double>> &acc, std::vector<int> &out_rows,
                            std::vector<std::complex<double>> &out_values) {
  out_rows.reserve(rows.size());
  out_values.reserve(rows.size());

  for (const int row : rows) {
    if (std::abs(acc[row]) > kEps) {
      out_values.push_back(acc[row]);
      out_rows.push_back(row);
    }
    acc[row] = std::complex<double>(0.0, 0.0);
    marker[row] = -1;
  }
}

bool IsValidCCS(const CCS &matrix) {
  if (matrix.m < 0 || matrix.n < 0) {
    return false;
  }
  if (matrix.col_ptr.size() != static_cast<std::size_t>(matrix.n) + 1) {
    return false;
  }
  if (matrix.row_ind.size() != matrix.values.size()) {
    return false;
  }
  if (matrix.col_ptr.empty() || matrix.col_ptr[0] != 0) {
    return false;
  }
  if (static_cast<std::size_t>(matrix.col_ptr.back()) != matrix.row_ind.size()) {
    return false;
  }
  for (int j = 0; j < matrix.n; ++j) {
    if (matrix.col_ptr[j] > matrix.col_ptr[j + 1]) {
      return false;
    }
  }
  return std::ranges::all_of(matrix.row_ind, [&](int row) { return row >= 0 && row < matrix.m; });
}

}  // namespace

SabutayASparseComplexCcsMultTBB::SabutayASparseComplexCcsMultTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void SabutayASparseComplexCcsMultTBB::SpMM(const CCS &a, const CCS &b, CCS &c) {
  c.m = a.m;
  c.n = b.n;
  c.col_ptr.assign(b.n + 1, 0);
  c.row_ind.clear();
  c.values.clear();

  std::vector<std::vector<int>> local_row_ind(b.n);
  std::vector<std::vector<std::complex<double>>> local_values(b.n);
  std::vector<int> local_sizes(b.n, 0);

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, b.n), [&](const oneapi::tbb::blocked_range<int> &range) {
    std::vector<int> rows;
    std::vector<int> marker(a.m, -1);
    std::vector<std::complex<double>> acc(a.m);

    for (int j = range.begin(); j < range.end(); ++j) {
      rows.clear();
      local_row_ind[j].clear();
      local_values[j].clear();

      AccumulateColumnContributions(a, b, j, rows, marker, acc);
      FlushAccumulatedValues(rows, marker, acc, local_row_ind[j], local_values[j]);

      local_sizes[j] = static_cast<int>(local_values[j].size());
    }
  });

  for (int j = 0; j < b.n; ++j) {
    c.col_ptr[j + 1] = c.col_ptr[j] + local_sizes[j];
  }

  c.row_ind.reserve(c.col_ptr.back());
  c.values.reserve(c.col_ptr.back());

  for (int j = 0; j < b.n; ++j) {
    c.row_ind.insert(c.row_ind.end(), local_row_ind[j].begin(), local_row_ind[j].end());
    c.values.insert(c.values.end(), local_values[j].begin(), local_values[j].end());
  }
}

bool SabutayASparseComplexCcsMultTBB::ValidationImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  return IsValidCCS(a) && IsValidCCS(b) && a.n == b.m;
}

bool SabutayASparseComplexCcsMultTBB::PreProcessingImpl() {
  return true;
}

bool SabutayASparseComplexCcsMultTBB::RunImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  CCS &c = GetOutput();

  SpMM(a, b, c);

  return true;
}

bool SabutayASparseComplexCcsMultTBB::PostProcessingImpl() {
  return true;
}

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
